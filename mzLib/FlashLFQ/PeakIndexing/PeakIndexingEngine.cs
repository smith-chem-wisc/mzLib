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
using FlashLFQ.Interfaces;

namespace FlashLFQ
{
    public class PeakIndexingEngine : IIndexingEngine
    {
        private List<IndexedMassSpectralPeak>[] _indexedPeaks;
        private readonly Serializer _serializer;
        private const int BinsPerDalton = 100;
        private readonly int _maxThreads;
        public Ms1ScanInfo[] Ms1ScanInfoArray { get; private set; }
        public SpectraFileInfo spectraFileInfo { get; private set; }

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

        public bool IndexPeaks(SpectraFileInfo fileInfo, bool silent)
        {
            if (!silent)
            {
                Console.WriteLine("Reading spectra file");
            }

            MsDataScan[] msDataScans = null;

            // read spectra file
            string fileName = fileInfo.FullFilePathWithExtension;
            var reader = MsDataFileReader.GetDataFile(fileName);
            //if (reader is TimsTofFileReader)
            //    return IndexTimsTofPeaks((TimsTofFileReader)reader, fileInfo, silent, _ms1Scans);
            reader.LoadAllStaticData();
            // retrieve only the ms1s. 
            msDataScans = reader.GetMS1Scans()
                .Where(i => i != null && i.MsnOrder == 1)
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
            Ms1ScanInfoArray = new Ms1ScanInfo[msDataScans.Length];
            List<Ms1ScanInfo> scanInfo = new List<Ms1ScanInfo>();

            for (int i = 0; i < msDataScans.Length; i++)
            {

                Ms1ScanInfoArray[i] = new Ms1ScanInfo(msDataScans[i].OneBasedScanNumber, scanIndex, msDataScans[i].RetentionTime);
                //scanInfo.Add(new Ms1ScanInfo(msDataScans[i].OneBasedScanNumber, scanIndex, msDataScans[i].RetentionTime));

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

            //_ms1Scans.Add(fileInfo, scanInfo.ToArray());

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

        //public bool IndexTimsTofPeaks(TimsTofFileReader file, SpectraFileInfo fileInfo, bool silent, Dictionary<SpectraFileInfo, Ms1ScanInfo[]> _ms1Scans)
        //{
        //    file.InitiateDynamicConnection();

        //    Ms1ScanInfo[] scanInfoArray = new Ms1ScanInfo[file.NumberOfMs1Frames];
        //    PpmTolerance tolerance = new PpmTolerance(20);

        //    ConcurrentDictionary<int, Dictionary<int, List<TraceableTimsTofPeak>>> binToFramePeakDict = new();

        //    var frameArray = file.GetMs1InfoFrameByFrame(maxThreads: _maxThreads);

        //    for (int i = 0; i < frameArray.Length; i++)
        //    {
        //        var ms1Scan = frameArray[i];
        //        //Dictionary<int, List<TraceableTimsTofPeak>> roundedMzObservedPeakDict = new();
        //        for (int scanIdx = 0; scanIdx < ms1Scan.TimsScanIdxMs1SpectraList.Count; scanIdx++)
        //        {
        //            var spectrum = ms1Scan.TimsScanIdxMs1SpectraList[scanIdx].Spectrum;

        //            // for every mz peak, create an IonMobilityPeak and assign it to the appropriate TraceableTimsTofPeak
        //            for (int spectrumIdx = 0; spectrumIdx < spectrum.Size; spectrumIdx++)
        //            {
        //                var ionMobilityPeak = new IonMobilityPeak(spectrum.TofIndices[spectrumIdx], (int)spectrum.YArray[spectrumIdx], scanIdx + 1);
        //                int roundedMz = (int)Math.Round(ionMobilityPeak.Mz * BinsPerDalton, 0);

        //                if (binToFramePeakDict.TryGetValue(roundedMz, out var framePeaks))
        //                {
        //                    if (framePeaks.TryGetValue(i, out var traceablePeaks))
        //                    {
        //                        var matchingPeak = traceablePeaks
        //                            .MinBy(p => Math.Abs(spectrum.TofIndices[spectrumIdx] - p.Mz));
        //                        if (tolerance.Within(matchingPeak.Mz, spectrum.TofIndices[spectrumIdx])) matchingPeak.AddIonMobilityPeak(ionMobilityPeak);
        //                        else traceablePeaks.Add(new TraceableTimsTofPeak(i, ms1Scan.RetentionTime, ionMobilityPeak));
        //                    }
        //                    else
        //                    {
        //                        framePeaks[i] = new List<TraceableTimsTofPeak>() { new TraceableTimsTofPeak(i, ms1Scan.RetentionTime, ionMobilityPeak) };
        //                    }
        //                }
        //                else
        //                {
        //                    binToFramePeakDict[roundedMz] = new Dictionary<int, List<TraceableTimsTofPeak>>() {
        //                        { i, new List<TraceableTimsTofPeak>() { new TraceableTimsTofPeak(i, ms1Scan.RetentionTime, ionMobilityPeak) } } };
        //                }
        //            }
        //        }

        //        scanInfoArray[i] = new Ms1ScanInfo((int)ms1Scan.FrameId, i, ms1Scan.RetentionTime);
        //        frameArray[i] = null;
        //    }

        //    _indexedPeaks = new List<IndexedMassSpectralPeak>[(int)Math.Ceiling(file.ScanWindow.Maximum) * BinsPerDalton + 1];
        //    foreach (var binToFramePeakDictKvp in binToFramePeakDict.OrderBy(kvp => kvp.Key))
        //    {
        //        _indexedPeaks[binToFramePeakDictKvp.Key] = binToFramePeakDictKvp.Value
        //            .OrderBy(kvp => kvp.Key)
        //            .SelectMany(kvp => kvp.Value)
        //            .Where(x => x != null)
        //            .AsParallel()
        //            .SelectMany(p => p.GetIndexedPeaks())
        //            .Where(p => p != null)
        //            .ToList();    
        //    }

        //    _ms1Scans.Add(fileInfo, scanInfoArray);

        //    if (_indexedPeaks == null || _indexedPeaks.Length == 0)
        //    {
        //        if (!silent)
        //        {
        //            Console.WriteLine("FlashLFQ Error: The file " + fileInfo.FilenameWithoutExtension + " contained no MS1 peaks!");
        //        }

        //        return false;
        //    }

        //    return true;
        //}

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
                _indexedPeaks = (List<IndexedMassSpectralPeak>[])_serializer.Deserialize(indexFile);
            }

            File.Delete(indexPath);
        }

        public IIndexedPeak GetIndexedPeak(double theorMass, int zeroBasedScanIndex, Tolerance tolerance, int chargeState)
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
                            if (bestPeak == null) 
                                bestPeak = peak;
                            
                            else if(Math.Abs(expMass - theorMass) < Math.Abs(bestPeak.Mz.ToMass(chargeState) - theorMass))
                                bestPeak = peak;
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