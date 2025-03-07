using Chemistry;
using Readers;
using MassSpectrometry;
using MzLibUtil;
using NetSerializer;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using FlashLFQ.PeakIndexing;
using Easy.Common.Extensions;
using System.Collections.Concurrent;
using System.Threading.Tasks;

namespace FlashLFQ
{
    public class PeakIndexingEngine
    {
        private List<IndexedMassSpectralPeak>[] _indexedPeaks;
        private readonly Serializer _serializer;
        private const int BinsPerDalton = 100;
        private readonly int _maxThreads;

        public PeakIndexingEngine(int maxThreads)
        {
            var messageTypes = new List<Type>
            {
                typeof(List<IndexedMassSpectralPeak>[]), 
                typeof(List<IndexedMassSpectralPeak>),
                typeof(IndexedMassSpectralPeak),
                typeof(List<IndexedTimsTofPeak>[]),
                typeof(List<IndexedTimsTofPeak>),
                typeof(IndexedTimsTofPeak)
            };
            _serializer = new Serializer(messageTypes);
            _maxThreads = maxThreads;
        }

        public bool IndexMassSpectralPeaks(SpectraFileInfo fileInfo, bool silent, Dictionary<SpectraFileInfo, Ms1ScanInfo[]> _ms1Scans)
        {
            if (!silent)
            {
                Console.WriteLine("Reading spectra file");
            }

            MsDataScan[] msDataScans = null;

            // read spectra file
            string fileName = fileInfo.FullFilePathWithExtension;
            var reader = MsDataFileReader.GetDataFile(fileName);
            if (reader is TimsTofFileReader)
                return IndexTimsTofPeaks((TimsTofFileReader)reader, fileInfo, silent, _ms1Scans);
            reader.LoadAllStaticData();
            // retrieve only the ms1s. 
            msDataScans = reader.GetMS1Scans().Where(i => i.MsnOrder == 1)
                .Select(i => i)
                .OrderBy(i => i.OneBasedScanNumber)
                .ToArray(); 
            
            if (!msDataScans.Any(p => p != null))
            {
                _indexedPeaks = Array.Empty<List<IndexedMassSpectralPeak>>();
                return false;
            }

            _indexedPeaks = new List<IndexedMassSpectralPeak>[(int)Math.Ceiling(msDataScans.Where(p => p != null
                && p.MassSpectrum.LastX != null).Max(p => p.MassSpectrum.LastX.Value) * BinsPerDalton) + 1];

            int scanIndex = 0;
            List<Ms1ScanInfo> scanInfo = new List<Ms1ScanInfo>();

            for (int i = 0; i < msDataScans.Length; i++)
            {
                if (msDataScans[i] == null)
                {
                    continue;
                }

                scanInfo.Add(new Ms1ScanInfo(msDataScans[i].OneBasedScanNumber, scanIndex, msDataScans[i].RetentionTime));

                for (int j = 0; j < msDataScans[i].MassSpectrum.XArray.Length; j++)
                {
                    int roundedMz = (int)Math.Round(msDataScans[i].MassSpectrum.XArray[j] * BinsPerDalton, 0);
                    if (_indexedPeaks[roundedMz] == null)
                    {
                        _indexedPeaks[roundedMz] = new List<IndexedMassSpectralPeak>();
                    }

                    _indexedPeaks[roundedMz].Add(new IndexedMassSpectralPeak(msDataScans[i].MassSpectrum.XArray[j],
                        msDataScans[i].MassSpectrum.YArray[j], scanIndex, msDataScans[i].RetentionTime));
                }

                scanIndex++;
            }

            _ms1Scans.Add(fileInfo, scanInfo.ToArray());

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

        public bool IndexTimsTofPeaks(TimsTofFileReader file, SpectraFileInfo fileInfo, bool silent, Dictionary<SpectraFileInfo, Ms1ScanInfo[]> _ms1Scans)
        {
            file.InitiateDynamicConnection();

            Ms1ScanInfo[] scanInfoArray = new Ms1ScanInfo[file.NumberOfMs1Frames];
            PpmTolerance tolerance = new PpmTolerance(20);
            //Dictionary<int, List<TraceableTimsTofPeak>> roundedMzObservedPeakDict = new();

            //int zeroBasedMs1FrameIndex = 0;
            ConcurrentDictionary<int, Dictionary<int, List<IndexedMassSpectralPeak>>> frameObservedPeaksDict = new();
            // foreach frame, build a collection of TraceableTimsTofPeaks that will be added to _indexedPeaks
            //foreach (TimsDataScan ms1Scan in file.GetMs1InfoScanByScan())

            //
            var frameCollection = new BlockingCollection<TimsDataScan>(new ConcurrentQueue<TimsDataScan>());
            file.GetMs1InfoFrameByFrame(frameCollection, maxThreads: _maxThreads);

            //file.GetMs1InfoFrameByFrame(maxThreads: _maxThreads);


            Task[] indexingTasks = new Task[_maxThreads];
            //var enumerator = scanEnumerable.GetEnumerator();
            for (int thread = 0; thread < _maxThreads; thread++)
            {
                var indexTask = Task.Factory.StartNew(() =>
                    {
                        while (frameCollection.TryTake(out var ms1Scan))
                        //while(enumerator.)
                        {
                            int zeroBasedMs1FrameIndex = file.Ms1FrameIds.IndexOf(ms1Scan.FrameId);
                            Dictionary<int, List<TraceableTimsTofPeak>> roundedMzObservedPeakDict = new();
                            for (int scanIdx = 0; scanIdx < ms1Scan.TimsScanIdxMs1SpectraList.Count; scanIdx++)
                            {
                                var spectrum = ms1Scan.TimsScanIdxMs1SpectraList[scanIdx].Spectrum;

                                // for every mz peak, create an IonMobilityPeak and assign it to the appropriate TraceableTimsTofPeak
                                for (int spectrumIdx = 0; spectrumIdx < spectrum.Size; spectrumIdx++)
                                {
                                    var ionMobilityPeak = new IonMobilityPeak(spectrum.XArray[spectrumIdx], (int)spectrum.YArray[spectrumIdx], scanIdx + 1);
                                    int roundedMz = (int)Math.Round(ionMobilityPeak.Mz * BinsPerDalton, 0);

                                    if (roundedMzObservedPeakDict.TryGetValue(roundedMz, out var traceablePeaks))
                                    {
                                        var matchingPeak = traceablePeaks
                                            .MinBy(p => Math.Abs(spectrum.XArray[spectrumIdx] - p.Mz));
                                        if (tolerance.Within(matchingPeak.Mz, spectrum.XArray[spectrumIdx])) matchingPeak.AddIonMobilityPeak(ionMobilityPeak);
                                        else traceablePeaks.Add(new TraceableTimsTofPeak(zeroBasedMs1FrameIndex, ms1Scan.RetentionTime, ionMobilityPeak));
                                    }
                                    else
                                    {
                                        roundedMzObservedPeakDict[roundedMz] = new List<TraceableTimsTofPeak> {
                                new TraceableTimsTofPeak(zeroBasedMs1FrameIndex, ms1Scan.RetentionTime, ionMobilityPeak) };
                                    }
                                }
                            }

                            Dictionary<int, List<IndexedMassSpectralPeak>> peakDict = new();
                            foreach (var kvp in roundedMzObservedPeakDict)
                            {
                                List<IndexedMassSpectralPeak> imps = new();
                                foreach (var traceablePeak in kvp.Value)
                                {
                                    var peaksFromTraceable = traceablePeak.GetIndexedPeaks();
                                    if (peaksFromTraceable.IsNotNullOrEmpty())
                                    {
                                        imps.AddRange(peaksFromTraceable);
                                    }
                                }
                                peakDict[kvp.Key] = imps;
                            }

                            frameObservedPeaksDict.TryAdd(
                                (int)ms1Scan.FrameId,
                                peakDict);

                            scanInfoArray[zeroBasedMs1FrameIndex] = new Ms1ScanInfo((int)ms1Scan.FrameId, zeroBasedMs1FrameIndex, ms1Scan.RetentionTime);
                        }
                    });
                indexingTasks[thread] = indexTask;
            }

            Task.WaitAll(indexingTasks);
            _indexedPeaks = new List<IndexedMassSpectralPeak>[(int)Math.Ceiling(file.ScanWindow.Maximum) * BinsPerDalton + 1];
            foreach (var frameDictPair in frameObservedPeaksDict.OrderBy(kvp => kvp.Key))
            {
                foreach (var kvp in frameDictPair.Value)
                {

                    if (kvp.Value.IsNotNullOrEmpty())
                    {
                        if (_indexedPeaks[kvp.Key] == null)
                            _indexedPeaks[kvp.Key] = new List<IndexedMassSpectralPeak>();
                        else
                        {
                            _indexedPeaks[kvp.Key].AddRange(kvp.Value);
                        }
                    }
                    
                }
            }

            _ms1Scans.Add(fileInfo, scanInfoArray);

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

        public void ClearIndex()
        {
            if (_indexedPeaks != null)
            {
                for (int i = 0; i < _indexedPeaks.Length; i++)
                {
                    if (_indexedPeaks[i] == null)
                    {
                        continue;
                    }

                    _indexedPeaks[i].Clear();
                    _indexedPeaks[i].TrimExcess();
                    _indexedPeaks[i] = null;
                }
            }

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
                _indexedPeaks = (List<IndexedMassSpectralPeak>[])_serializer.Deserialize(indexFile);
            }

            File.Delete(indexPath);
        }

        public IndexedMassSpectralPeak GetIndexedPeak(double theorMass, int zeroBasedScanIndex, Tolerance tolerance, int chargeState, int? timsIndex = null)
        {
            IndexedMassSpectralPeak bestPeak = null;
            int ceilingMz = (int)Math.Ceiling(tolerance.GetMaximumValue(theorMass).ToMz(chargeState) * BinsPerDalton);
            int floorMz = (int)Math.Floor(tolerance.GetMinimumValue(theorMass).ToMz(chargeState) * BinsPerDalton);

            for (int j = floorMz; j <= ceilingMz; j++)
            {
                if (j < _indexedPeaks.Length && _indexedPeaks[j] != null)
                {
                    List<IndexedMassSpectralPeak> bin = _indexedPeaks[j];
                    if (bin == null || bin.Count == 0) return null;
                    int index = BinarySearchForIndexedPeak(bin, zeroBasedScanIndex);

                    for (int i = index; i < bin.Count; i++)
                    {
                        IndexedMassSpectralPeak peak = bin[i];

                        if (peak.ZeroBasedMs1ScanIndex > zeroBasedScanIndex)
                        {
                            break;
                        }

                        double expMass = peak.Mz.ToMass(chargeState);

                        if (tolerance.Within(expMass, theorMass) && peak.ZeroBasedMs1ScanIndex == zeroBasedScanIndex)
                        {
                            if (bestPeak == null) bestPeak = peak;
                            else if (timsIndex != null && peak is IndexedTimsTofPeak && bestPeak is IndexedTimsTofPeak)
                            { 
                                // TODO: Decide on how to choose when considering ion mobility and mass accuracy
                            }
                            else if(Math.Abs(expMass - theorMass) < Math.Abs(bestPeak.Mz.ToMass(chargeState) - theorMass))
                            {
                                bestPeak = peak;
                            }
                        }
                    }
                }
            }

            return bestPeak;
        }

        private int BinarySearchForIndexedPeak(List<IndexedMassSpectralPeak> indexedPeaks, int zeroBasedScanIndex)
        {
            int m = 0;
            int l = 0;
            int r = indexedPeaks.Count - 1;

            while (l <= r)
            {
                m = l + ((r - l) / 2);

                if (r - l < 2)
                {
                    break;
                }
                if (indexedPeaks[m].ZeroBasedMs1ScanIndex < zeroBasedScanIndex)
                {
                    l = m + 1;
                }
                else
                {
                    r = m - 1;
                }
            }

            for (int i = m; i >= 0; i--)
            {
                if (indexedPeaks[i].ZeroBasedMs1ScanIndex < zeroBasedScanIndex)
                {
                    break;
                }

                m--;
            }

            if (m < 0)
            {
                m = 0;
            }

            return m;
        }
    }
}