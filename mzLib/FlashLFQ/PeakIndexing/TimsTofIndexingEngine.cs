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
using MassSpectrometry;
using System.CodeDom;

namespace FlashLFQ
{
    public class TimsTofIndexingEngine : IIndexingEngine
    {
        private List<IndexedTimsTofPeak>[] _indexedPeaks;
        private readonly Serializer _serializer;
        private const int BinsPerDalton = 100;
        private readonly int _maxThreads;
        public SpectraFileInfo FileInfo { get; init; }
        public Ms1ScanInfo[] ScanInfoArray { get; private set; }
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

        public TimsTofIndexingEngine(SpectraFileInfo fileInfo) : this(fileInfo, -1, 24) { }

        public TimsTofIndexingEngine(SpectraFileInfo fileInfo, int maxThreads, int spectraPerFrame = 8)
        {
            var messageTypes = new List<Type>
            {
                typeof(List<IndexedTimsTofPeak>[]),
                typeof(List<IndexedTimsTofPeak>),
                typeof(IndexedTimsTofPeak)
            };
            _serializer = new Serializer(messageTypes);
            if(maxThreads <= 0) maxThreads = Environment.ProcessorCount - 1;
            _maxThreads = maxThreads;
            FileInfo = fileInfo;
            SpectraPerFrame = spectraPerFrame;
        }

        public int SpectraPerFrame { get; init; }

        public bool IndexPeaks(bool silent)
        {
            TimsTofFileReader file = new TimsTofFileReader(FileInfo.FullFilePathWithExtension);
            file.InitiateDynamicConnection(ppmToleranceForCentroiding: 5);
            MzLookupArray = file.GetMzLookupTable();

            ScanInfoArray = new Ms1ScanInfo[file.NumberOfMs1Frames];
            var frameArray = file.GetMs1InfoFrameByFrame(SpectraPerFrame, out int scansPerSpectra, maxThreads: _maxThreads);
            ImsResolution = scansPerSpectra;
            _indexedPeaks = new List<IndexedTimsTofPeak>[(int)Math.Ceiling(file.ScanWindow.Maximum) * BinsPerDalton + 1];
            file.CloseDynamicConnection();

            //int noiseBaseline = GetApproximateNoiseLevel(frameArray);

            for (int i = 0; i < frameArray.Length; i++)
            {
                var ms1Scan = frameArray[i];

                for (int j = 0; j < ms1Scan.ComponentSpectra.Count; j++)
                {
                    var spectrum = ms1Scan.ComponentSpectra[j];
                    int timsIndex = ms1Scan.ComponentSpectra[j].ZeroIndexedIonMobilityScanNumber;

                    // for every mz peak, create an IonMobilityPeak and assign it to the appropriate TraceableTimsTofPeak
                    for (int spectrumIdx = 0; spectrumIdx < spectrum.Size; spectrumIdx++)
                    {
                        //if (spectrum.YArray[spectrumIdx] < noiseBaseline) continue;
                        double peakMz = MzLookupArray[spectrum.XArray[spectrumIdx]];
                        int roundedMz = (int)Math.Floor(peakMz * BinsPerDalton);
                        _indexedPeaks[roundedMz] ??= new List<IndexedTimsTofPeak>(frameArray.Length / 100);
                        _indexedPeaks[roundedMz].Add(new IndexedTimsTofPeak(spectrum.XArray[spectrumIdx], timsIndex, spectrum.YArray[spectrumIdx], i));
                    }
                }

                ScanInfoArray[i] = new Ms1ScanInfo((int)ms1Scan.FrameId, i, ms1Scan.RetentionTime);
                frameArray[i] = null;
            }

            if (_indexedPeaks == null || _indexedPeaks.Length == 0)
            {
                if (!silent)
                {
                    Console.WriteLine("FlashLFQ Error: The file " + FileInfo.FilenameWithoutExtension + " contained no MS1 peaks!");
                }

                return false;
            }

            //throw new MzLibException("Done indexing");

            return true;
            
        }

        public static int GetApproximateNoiseLevel(TimsDataScan[] frameArray)
        {
            //Grab 20 frames, create an intensity histogram, and find the most frequent intensity level. That's the noise floor
            int maxNoiseIntensity = 20000;
            int[] intensityFrequencyArray = new int[maxNoiseIntensity]; //implicit assumption that the noise floor is below 300
            for(int i = 0; i < frameArray.Length; i += frameArray.Length/20)
            {
                foreach(var scanSpectrumTuple in frameArray[i].ComponentSpectra)
                {
                    var intensitySpectrum = scanSpectrumTuple.YArray;
                    for(int j = 0; j < intensitySpectrum.Length; j++)
                    {
                        if (intensitySpectrum[j] < maxNoiseIntensity)
                            intensityFrequencyArray[(int)intensitySpectrum[j]]++;
                    }
                }
            }
            var frequencySortedList = intensityFrequencyArray.Where(x => x > 0).OrderByDescending(x => x).ToList();
            int point1PercentileFrequency = frequencySortedList[frequencySortedList.Count / 1000];
            return intensityFrequencyArray.IndexOf(point1PercentileFrequency);
        }

        public IndexedIonMobilityPeak GetIndexedPeak(double mz, int zeroBasedScanIndex, PpmTolerance tolerance,  int? timsIndex = null)
        {
            IndexedMassSpectralPeak bestPeak = null;
            double floorMz = tolerance.GetMinimumValue(mz);
            double ceilingMz = tolerance.GetMaximumValue(mz);

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

                        double expMz = MzLookupArray[peak.TofIndex];
                        if (tolerance.Within(expMz, mz) && peak.ZeroBasedMs1FrameIndex == zeroBasedScanIndex)
                        {
                            if (peaksInFrame.TryGetValue(peak.TimsIndex, out var peaksAtIms))
                                peaksAtIms.Add(peak);
                            else
                                peaksInFrame[peak.TimsIndex] = new List<IndexedTimsTofPeak> { peak };
                        }
                    }
                }
            }

            return MergeTimsTofPeaks(peaksInFrame, mz, timsIndex);
        }

        private PpmTolerance mergeTolerance = new PpmTolerance(15);

        private IndexedIonMobilityPeak MergeTimsTofPeaks(Dictionary<int, List<IndexedTimsTofPeak>> peaksInFrame, double targetMz, int? targetTimsScanIndex = null)
        {
            if (!peaksInFrame.Any()) return null;

            List<IndexedTimsTofPeak> peaksByTimsScanIndex = new();
            foreach (var peakList in peaksInFrame.OrderBy(kvp => kvp.Key).Select(kvp => kvp.Value))
            {
                if (peakList.Count > 1)
                {
                    
                    var bestPeak = peakList.MinBy(p => Math.Abs(targetMz - MzLookupArray[p.TofIndex]));
                    peaksByTimsScanIndex.Add(bestPeak);
                    if (peakList.Count >= 2)
                    {
                        peakList.Remove(bestPeak);
                        uint minTofIndex = bestPeak.TofIndex - 2;
                        uint maxTofIndex = bestPeak.TofIndex + 2;
                        var otherPeaks = peakList.Where(p => p.TofIndex >= minTofIndex && p.TofIndex <= maxTofIndex).ToList();
                        if(otherPeaks.IsNotNullOrEmpty())
                            peaksByTimsScanIndex.AddRange(otherPeaks);
                    }
                }
                else
                {
                    peaksByTimsScanIndex.Add(peakList[0]);
                }
            }

            return MergeTimsTofPeaksHelper(peaksByTimsScanIndex, targetTimsScanIndex);
        }

        private IndexedIonMobilityPeak MergeTimsTofPeaksHelper(List<IndexedTimsTofPeak> peaksByTimsScanIndex, int? targetTimsScanIndex)
        {
            if(peaksByTimsScanIndex.Count < 1) return null;

            IndexedTimsTofPeak peakCenter = default(IndexedTimsTofPeak);
            if(targetTimsScanIndex == null)
            {
                peakCenter = peaksByTimsScanIndex.MaxBy(p => p.Intensity);
            }
            else
            {
                peakCenter = peaksByTimsScanIndex.MinBy(p => Math.Abs(p.TimsIndex - targetTimsScanIndex.Value));
                if (Math.Abs(peakCenter.TimsIndex - (int)targetTimsScanIndex) > 2.5 * ImsResolution)
                    return null;
            }

            int centerIndex = peaksByTimsScanIndex.IndexOf(peakCenter);
            int leftIndex = centerIndex - 1;
            int previousTimsIndex = peakCenter.TimsIndex;
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
            if (leftIndex < 0) leftIndex = 0; // if we are at the beginning of the list, we don't want to go out of bounds (leftIndex--)

            int rightIndex = centerIndex + 1;
            previousTimsIndex = peakCenter.TimsIndex;
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
            if (rightIndex >= peaksByTimsScanIndex.Count) rightIndex = peaksByTimsScanIndex.Count - 1; // if we are at the end of the list, we don't want to go out of bounds (rightIndex++)

            var apex = peaksByTimsScanIndex[leftIndex..(rightIndex + 1)].MaxBy(p => p.Intensity);
            // TODO: can probably calculate the summed intensity and find the apex in the loops above

            return new IndexedIonMobilityPeak(
                GetWeightedAverageMz(peaksByTimsScanIndex[leftIndex..(rightIndex + 1)], out int summedIntensity),
                summedIntensity,
                peakCenter.ZeroBasedMs1FrameIndex,
                ScanInfoArray[peakCenter.ZeroBasedMs1FrameIndex].RetentionTime,
                ionMobilityValues: peaksByTimsScanIndex[leftIndex..(rightIndex + 1)].Select(p => p.TimsIndex).ToHashSet(),
                apexIonMobilityValue: apex.TimsIndex);
        }

        public double GetWeightedAverageMz(List<IndexedTimsTofPeak> timsTofPeaks, out int summedIntensity)
        {
            summedIntensity = timsTofPeaks.Sum(p => p.Intensity);
            double averagedMz = 0;
            foreach(var peak in timsTofPeaks)
            {
                double weight = (double)peak.Intensity / (double)summedIntensity;
                averagedMz += weight * MzLookupArray[peak.TofIndex];
            }
            return averagedMz;
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

        public void SerializeIndex()
        {
            string dir = Path.GetDirectoryName(FileInfo.FullFilePathWithExtension);
            string indexPath = Path.Combine(dir, FileInfo.FilenameWithoutExtension + ".ind");

            using (var indexFile = File.Create(indexPath))
            {
                _serializer.Serialize(indexFile, _indexedPeaks);
            }
        }

        public void DeserializeIndex()
        {
            string dir = Path.GetDirectoryName(FileInfo.FullFilePathWithExtension);
            string indexPath = Path.Combine(dir, FileInfo.FilenameWithoutExtension + ".ind");

            using (var indexFile = File.OpenRead(indexPath))
            {
                _indexedPeaks = (List<IndexedTimsTofPeak>[])_serializer.Deserialize(indexFile);
            }

            File.Delete(indexPath);
        }

        public IIndexedMzPeak GetIndexedPeak(double mz, int zeroBasedScanIndex, PpmTolerance tolerance)
        {
            return GetIndexedPeak(mz, zeroBasedScanIndex, tolerance, timsIndex: null);
        }

    }
}
