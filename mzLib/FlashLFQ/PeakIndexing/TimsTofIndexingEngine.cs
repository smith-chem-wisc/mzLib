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
using FlashLFQ.PeakIndexing; 

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
            PpmTolerance tolerance = new PpmTolerance(20);

            ConcurrentDictionary<int, Dictionary<int, List<TraceableTimsTofPeak>>> binToFramePeakDict = new();

            var frameArray = file.GetMs1InfoFrameByFrame(maxThreads: _maxThreads);
            _indexedPeaks = new List<IndexedTimsTofPeak>[(int)Math.Ceiling(file.ScanWindow.Maximum) * BinsPerDalton + 1];

            // Read in a chunk of frames
            // Iterate through individual scans, creating ims peaks
            // Group imsPeaks int TraceableTimsTofPeaks
            // Process TraceableTimsTofPeaks into IndexedTimsTofPeaks
            // Sort the IndexedTimsTofPeaks into the appropriate bins of either _indexedPeaks or an _indexedPeaks subarray

            for (int i = 0; i < frameArray.Length; i++)
            {
                var ms1Scan = frameArray[i];
                Dictionary<int, List<TraceableTimsTofPeak>> binPeakDictionary = new();

                for (int j = 0; j < ms1Scan.TimsScanIdxMs1SpectraList.Count; j++)
                {
                    var spectrum = ms1Scan.TimsScanIdxMs1SpectraList[j].Spectrum;
                    int timsIndex = ms1Scan.TimsScanIdxMs1SpectraList[j].ScanIdx;

                    // for every mz peak, create an IonMobilityPeak and assign it to the appropriate TraceableTimsTofPeak
                    for (int spectrumIdx = 0; spectrumIdx < spectrum.Size; spectrumIdx++)
                    {
                        var ionMobilityPeak = new IonMobilityPeak(spectrum.XArray[spectrumIdx], spectrum.YArray[spectrumIdx], timsIndex);
                        double peakMz = MzLookupArray[ionMobilityPeak.TofIndex];
                        int roundedMz = (int)Math.Round(peakMz * BinsPerDalton, 0);

                        if (binPeakDictionary.TryGetValue(roundedMz, out var framePeaks))
                        {
                            var matchingPeak = framePeaks
                                .MinBy(p => Math.Abs(peakMz - MzLookupArray[p.TofIndex])); // This could probably be optimized
                            if (tolerance.Within(matchingPeak.TofIndex, peakMz)) matchingPeak.AddIonMobilityPeak(ionMobilityPeak);
                            else framePeaks.Add(new TraceableTimsTofPeak(i, ms1Scan.RetentionTime, ionMobilityPeak));
                        }
                        else
                        {
                            binPeakDictionary[roundedMz] = new List<TraceableTimsTofPeak> { new TraceableTimsTofPeak(i, ms1Scan.RetentionTime, ionMobilityPeak) };
                        }
                    }
                }

                foreach(var kvp in binPeakDictionary)
                {
                    _indexedPeaks[kvp.Key] ??= new List<IndexedTimsTofPeak>();
                    foreach(var traceablePeak in kvp.Value)
                    {
                        var ittPeaks = traceablePeak.GetIndexedPeaks();
                        if (ittPeaks.IsNotNullOrEmpty())
                            _indexedPeaks[kvp.Key].AddRange(ittPeaks);
                    }
                    
                }

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

            return true;
        }



        public IndexedMassSpectralPeak GetIndexedPeak(double theorMass, int zeroBasedScanIndex, Tolerance tolerance, int chargeState, int? timsIndex = null)
        {
            // Here, we'll need to convert between tofIndices for the stored peaks and mz values
            throw new NotImplementedException();
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
