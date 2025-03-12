using Easy.Common.Extensions;
using FlashLFQ.Interfaces;
using MzLibUtil;
using NetSerializer;
using Readers;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibUtil;
using FlashLFQ.PeakIndexing;
using Chemistry;

namespace FlashLFQ
{
    public class TimsTofIndexingEngine : IIndexingEngine
    {
        private List<IndexedTimsTofPeak>[] _indexedPeaks;
        private readonly Serializer _serializer;
        private const int BinsPerDalton = 100;
        private readonly int _maxThreads;
        public SpectraFileInfo FileInfo { get; init; }
        public Ms1ScanInfo[] Ms1ScanInfoArray { get; private set; }
        /// <summary>
        /// This int represents the number of tims scans that are averaged to create a single IMS scan
        /// Every spectrum from the same frame will have a different timsIndex that is at least
        /// ImsResolution scans apart from the previous one and no more than 2x ImsResolution apart from the previous one
        /// </summary>
        public int ImsResolution { get; private set; }

        /// <summary>
        /// Used to convert the tofIndices stored in the .d file to m/z values
        /// </summary>
        public double[] MzLookupArray { get; set; }
        /// <summary>
        /// Used to convert scan number to 1/K0 values
        /// </summary>
        public double[] OneOverK0LookupArray { get; set; }

        public TimsTofIndexingEngine(SpectraFileInfo fileInfo, int maxThreads)
        {
            var messageTypes = new List<Type>
            {
                typeof(List<IndexedTimsTofPeak>[]),
                typeof(List<IndexedTimsTofPeak>),
                typeof(IndexedTimsTofPeak)
            };
            _serializer = new Serializer(messageTypes);
            _maxThreads = maxThreads;
            FileInfo = fileInfo;
        }

        public bool IndexPeaks(SpectraFileInfo fileInfo, bool silent)
        {
            TimsTofFileReader file = new TimsTofFileReader(fileInfo.FullFilePathWithExtension);
            file.InitiateDynamicConnection();
            MzLookupArray = file.GetMzLookupTable();

            Ms1ScanInfoArray = new Ms1ScanInfo[file.NumberOfMs1Frames];
            var frameArray = file.GetMs1InfoFrameByFrame(out int scansPerSpectra, maxThreads: _maxThreads);
            ImsResolution = scansPerSpectra;
            _indexedPeaks = new List<IndexedTimsTofPeak>[(int)Math.Ceiling(file.ScanWindow.Maximum) * BinsPerDalton + 1];
            file.CloseDynamicConnection();

            for (int i = 0; i < frameArray.Length; i++)
            {
                var ms1Scan = frameArray[i];

                for (int j = 0; j < ms1Scan.TimsScanIdxMs1SpectraList.Count; j++)
                {
                    var spectrum = ms1Scan.TimsScanIdxMs1SpectraList[j].Spectrum;
                    int timsIndex = ms1Scan.TimsScanIdxMs1SpectraList[j].ScanIdx;

                    // for every mz peak, create an IonMobilityPeak and assign it to the appropriate TraceableTimsTofPeak
                    for (int spectrumIdx = 0; spectrumIdx < spectrum.Size; spectrumIdx++)
                    {
                        if (spectrum.YArray[spectrumIdx] < 150) continue;
                        double peakMz = MzLookupArray[spectrum.XArray[spectrumIdx]];
                        int roundedMz = (int)Math.Round(peakMz * BinsPerDalton, 0);
                        _indexedPeaks[roundedMz] ??= new List<IndexedTimsTofPeak>(frameArray.Length / 100);
                        _indexedPeaks[roundedMz].Add(new IndexedTimsTofPeak(spectrum.XArray[spectrumIdx], timsIndex, spectrum.YArray[spectrumIdx], i));
                    }
                }

                Ms1ScanInfoArray[i] = new Ms1ScanInfo((int)ms1Scan.FrameId, i, ms1Scan.RetentionTime);
                frameArray[i] = null;
            }

            if (_indexedPeaks == null || _indexedPeaks.Length == 0)
            {
                if (!silent)
                {
                    Console.WriteLine("FlashLFQ Error: The file " + fileInfo.FilenameWithoutExtension + " contained no MS1 peaks!");
                }

                return false;
            }

            //throw new MzLibException("Done indexing");

            return true;
            
        }



        public IndexedIonMobilityPeak GetIndexedPeak(double theorMass, int zeroBasedScanIndex, Tolerance tolerance, int chargeState, int? timsIndex = null)
        {
            IndexedMassSpectralPeak bestPeak = null;
            double floorMz = tolerance.GetMinimumValue(theorMass).ToMz(chargeState);
            double ceilingMz = tolerance.GetMaximumValue(theorMass).ToMz(chargeState);

            int floorBin = (int)Math.Floor(floorMz * BinsPerDalton);
            int ceilingBin = (int)Math.Ceiling(ceilingMz * BinsPerDalton);

            //We can have multiple peaks, each corresponding to different IMS scans in the same frame
            Dictionary<int, List<IndexedTimsTofPeak>> peaksInFrame = new();
            for (int j = floorBin; j <= ceilingBin; j++)
            {
                if (j < _indexedPeaks.Length && _indexedPeaks[j] != null)
                {
                    List<IndexedTimsTofPeak> bin = _indexedPeaks[j];
                    if (bin == null || bin.Count == 0) return null;
                  
                    int index = BinarySearchForIndexedPeak(bin, zeroBasedScanIndex);

                    for (int i = index; i < bin.Count; i++)
                    {
                        IndexedTimsTofPeak peak = bin[i];

                        if (peak.ZeroBasedMs1FrameIndex > zeroBasedScanIndex)
                            break;

                        double expMass = MzLookupArray[peak.TofIndex].ToMass(chargeState);
                        if (tolerance.Within(expMass, theorMass) && peak.ZeroBasedMs1FrameIndex == zeroBasedScanIndex)
                        {
                            if (peaksInFrame.TryGetValue(peak.TimsIndex, out var peaksAtIms))
                                peaksAtIms.Add(peak);
                            else
                                peaksInFrame[peak.TimsIndex] = new List<IndexedTimsTofPeak> { peak };
                        }
                    }
                }
            }

            return MergeTimsTofPeaks(peaksInFrame, theorMass.ToMz(chargeState));
        }

        private IndexedIonMobilityPeak MergeTimsTofPeaks(Dictionary<int, List<IndexedTimsTofPeak>> peaksInFrame, double targetMz, int? targetTimsScanIndex = null)
        {
            if (!peaksInFrame.Any()) return null;

            IndexedTimsTofPeak apex = default(IndexedTimsTofPeak);
            List<IndexedTimsTofPeak> peaksByTimsScanIndex = new();
            foreach (var peakList in peaksInFrame.OrderBy(kvp => kvp.Key).Select(kvp => kvp.Value))
            {
                if(peakList.Count > 1)
                {
                    var bestPeak = peakList.MinBy(p => Math.Abs(targetMz - MzLookupArray[p.TofIndex]));
                    peaksByTimsScanIndex.Add(bestPeak);
                    if(bestPeak.Intensity > apex.Intensity) 
                        apex = bestPeak;
                }
                else
                {
                    peaksByTimsScanIndex.Add(peakList[0]);
                    if(peakList[0].Intensity > apex.Intensity)
                          apex = peakList[0];
                }
            }

            int apexIndex = peaksByTimsScanIndex.IndexOf(apex);

            int leftIndex = apexIndex - 1;
            int previousTimsIndex = apex.TimsIndex;
            while (leftIndex >= 0)
            {
                if (previousTimsIndex - peaksByTimsScanIndex[leftIndex].TimsIndex < 1.5 * ImsResolution)
                {
                    previousTimsIndex = peaksByTimsScanIndex[leftIndex].TimsIndex;
                    leftIndex--;
                }
                else
                {
                    leftIndex++;
                    break;
                }
            }
            if(leftIndex < 0) leftIndex = 0; // if we are at the beginning of the list, we don't want to go out of bounds (leftIndex--)

            int rightIndex = apexIndex + 1;
            previousTimsIndex = apex.TimsIndex;
            while (rightIndex < peaksByTimsScanIndex.Count)
            {
                if (peaksByTimsScanIndex[rightIndex].TimsIndex - previousTimsIndex < 1.5 * ImsResolution)
                {
                    previousTimsIndex = peaksByTimsScanIndex[rightIndex].TimsIndex;
                    rightIndex++;
                }
                else
                {
                    rightIndex--;
                    break;
                }
            }
            if(rightIndex >= peaksByTimsScanIndex.Count) rightIndex = peaksByTimsScanIndex.Count - 1; // if we are at the end of the list, we don't want to go out of bounds (rightIndex++)


            return new IndexedIonMobilityPeak(
                MzLookupArray[apex.TofIndex], 
                peaksByTimsScanIndex[leftIndex..(rightIndex+1)].Sum(p => p.Intensity),
                apex.ZeroBasedMs1FrameIndex,
                Ms1ScanInfoArray[apex.ZeroBasedMs1FrameIndex].RetentionTime,
                ionMobilityValue: apex.TimsIndex);
        }

        private int BinarySearchForIndexedPeak(List<IndexedTimsTofPeak> bin, int zeroBasedScanIndex)
        {
            int min = 0;
            int max = bin.Count - 1;

            while (min <= max)
            {
                int mid = (min + max) / 2;
                if (bin[mid].ZeroBasedMs1FrameIndex == zeroBasedScanIndex)
                {
                    while (mid - 1 >=0 && bin[mid - 1].ZeroBasedMs1FrameIndex == zeroBasedScanIndex)
                    {
                        mid--;
                    }
                    return mid;
                }

                if (bin[mid].ZeroBasedMs1FrameIndex < zeroBasedScanIndex)
                {
                    min = mid + 1;
                }
                else
                {
                    max = mid - 1;
                }
            }

            return min;
        }

        public void ClearIndex()
        {
            _indexedPeaks = null;
            GC.Collect();
        }

        public void SerializeIndex(SpectraFileInfo file)
        {
            string dir = Path.GetDirectoryName(file.FullFilePathWithExtension);
            string indexPath = Path.Combine(dir, file.FilenameWithoutExtension + ".ind");

            using (var indexFile = File.Create(indexPath))
            {
                _serializer.Serialize(indexFile, _indexedPeaks);
            }
        }

        public void DeserializeIndex(SpectraFileInfo file)
        {
            string dir = Path.GetDirectoryName(file.FullFilePathWithExtension);
            string indexPath = Path.Combine(dir, file.FilenameWithoutExtension + ".ind");

            using (var indexFile = File.OpenRead(indexPath))
            {
                _indexedPeaks = (List<IndexedTimsTofPeak>[])_serializer.Deserialize(indexFile);
            }

            File.Delete(indexPath);
        }

        public IIndexedPeak GetIndexedPeak(double theorMass, int zeroBasedScanIndex, Tolerance tolerance, int chargeState)
        {
            return GetIndexedPeak(theorMass, zeroBasedScanIndex, tolerance, chargeState, timsIndex: null);
        }

    }
}
