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

            //ConcurrentDictionary<int, Dictionary<int, List<TraceableTimsTofPeak>>> binToFramePeakDict = new();

            
            var frameArray = file.GetMs1InfoFrameByFrame(out int scansPerSpectra, maxThreads: _maxThreads);
            ImsResolution = scansPerSpectra;
            _indexedPeaks = new List<IndexedTimsTofPeak>[(int)Math.Ceiling(file.ScanWindow.Maximum) * BinsPerDalton + 1];

            // Read in a chunk of frames
            // Iterate through individual scans, creating ims peaks
            // Group imsPeaks int TraceableTimsTofPeaks
            // Process TraceableTimsTofPeaks into IndexedTimsTofPeaks
            // Sort the IndexedTimsTofPeaks into the appropriate bins of either _indexedPeaks or an _indexedPeaks subarray

            for (int i = 0; i < frameArray.Length; i++)
            {
                var ms1Scan = frameArray[i];
                //Dictionary<int, List<TraceableTimsTofPeak>> binPeakDictionary = new();

                for (int j = 0; j < ms1Scan.TimsScanIdxMs1SpectraList.Count; j++)
                {
                    var spectrum = ms1Scan.TimsScanIdxMs1SpectraList[j].Spectrum;
                    int timsIndex = ms1Scan.TimsScanIdxMs1SpectraList[j].ScanIdx;

                    // for every mz peak, create an IonMobilityPeak and assign it to the appropriate TraceableTimsTofPeak
                    for (int spectrumIdx = 0; spectrumIdx < spectrum.Size; spectrumIdx++)
                    {
                        double peakMz = MzLookupArray[spectrum.XArray[spectrumIdx]];
                        int roundedMz = (int)Math.Round(peakMz * BinsPerDalton, 0);
                        _indexedPeaks[roundedMz] ??= new List<IndexedTimsTofPeak>(frameArray.Length / 100);
                        _indexedPeaks[roundedMz].Add(new IndexedTimsTofPeak(spectrum.XArray[spectrumIdx], timsIndex, spectrum.YArray[spectrumIdx], i));


                        //if (binPeakDictionary.TryGetValue(roundedMz, out var framePeaks))
                        //{
                        //    var matchingPeak = framePeaks
                        //        .MinBy(p => Math.Abs(peakMz - MzLookupArray[p.TofIndex])); // This could probably be optimized
                        //    if (tolerance.Within(matchingPeak.TofIndex, peakMz)) matchingPeak.AddIonMobilityPeak(ionMobilityPeak);
                        //    else framePeaks.Add(new TraceableTimsTofPeak(i, ms1Scan.RetentionTime, ionMobilityPeak));
                        //}
                        //else
                        //{
                        //    binPeakDictionary[roundedMz] = new List<TraceableTimsTofPeak> { new TraceableTimsTofPeak(i, ms1Scan.RetentionTime, ionMobilityPeak) };
                        //}
                    }
                }

                //foreach(var kvp in binPeakDictionary)
                //{
                //    _indexedPeaks[kvp.Key] ??= new List<IndexedTimsTofPeak>();
                //    foreach(var traceablePeak in kvp.Value)
                //    {
                //        var ittPeaks = traceablePeak.GetIndexedPeaks();
                //        if (ittPeaks.IsNotNullOrEmpty())
                //            _indexedPeaks[kvp.Key].AddRange(ittPeaks);
                //    }
                    
                //}

                Ms1ScanInfoArray[i] = new Ms1ScanInfo((int)ms1Scan.FrameId, i, ms1Scan.RetentionTime);
                frameArray[i] = null;
            }

            //_indexedPeaks = new List<IndexedTimsTofPeak>[(int)Math.Ceiling(file.ScanWindow.Maximum) * BinsPerDalton + 1];
           

            if (_indexedPeaks == null || _indexedPeaks.Length == 0)
            {
                if (!silent)
                {
                    Console.WriteLine("FlashLFQ Error: The file " + fileInfo.FilenameWithoutExtension + " contained no MS1 peaks!");
                }

                return false;
            }

            throw new MzLibException("Done indexing");

            return true;
        }



        public IndexedMassSpectralPeak GetIndexedPeak(double theorMass, int zeroBasedScanIndex, Tolerance tolerance, int chargeState, int? timsIndex = null)
        {
            IndexedMassSpectralPeak bestPeak = null;
            double floorMz = tolerance.GetMinimumValue(theorMass).ToMz(chargeState);
            double ceilingMz = tolerance.GetMaximumValue(theorMass).ToMz(chargeState);

            int floorBin = (int)Math.Floor(floorMz);
            int ceilingBin = (int)Math.Ceiling(ceilingMz);


            //List<IndexedTimsTofPeak> peaksInFrame = new(); 

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

            return bestPeak;
        }

        private 

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
            Ms1ScanInfoArray = null;
        }

        public void SerializeIndex(SpectraFileInfo file)
        {
            using (var stream = File.OpenWrite(file.FilenameWithoutExtension + ".ind"))
            {
                _serializer.Serialize(stream, _indexedPeaks);
            }
        }

        public void DeserializeIndex(SpectraFileInfo file)
        {
            using (var stream = File.OpenRead(file.FilenameWithoutExtension + ".ind"))
            {
                _indexedPeaks = (List<IndexedTimsTofPeak>[])_serializer.Deserialize(stream);
            }
        }

        public IndexedMassSpectralPeak GetIndexedPeak(double theorMass, int zeroBasedScanIndex, Tolerance tolerance, int chargeState)
        {
            // Implementation of GetIndexedPeak method
            throw new NotImplementedException();
        }

    }
}
