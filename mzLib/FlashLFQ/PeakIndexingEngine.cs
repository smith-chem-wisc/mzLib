using Chemistry;
using Readers;
using MassSpectrometry;
using MzLibUtil;
using NetSerializer;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MzLibUtil.SparseMatrix;

namespace FlashLFQ
{
    public class PeakIndexingEngine
    {
        private List<IndexedMassSpectralPeak>[] _indexedPeaks;
        private readonly Serializer _serializer;
        private const int BinsPerDalton = 100;

        public PeakIndexingEngine()
        {
            var messageTypes = new List<Type>
            {
                typeof(List<IndexedMassSpectralPeak>[]), 
                typeof(List<IndexedMassSpectralPeak>),
                typeof(IndexedMassSpectralPeak),
                typeof(List<IndexedTimsTofPeak>[]),
                typeof(List<IndexedTimsTofPeak>),
                typeof(IndexedTimsTofPeak),
                typeof(List<IonMobilityPeak>),
                typeof(IonMobilityPeak)
            };
            _serializer = new Serializer(messageTypes);
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
                _indexedPeaks = new List<IndexedMassSpectralPeak>[0];
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
            // create the indexed peaks array
            _indexedPeaks = new List<IndexedMassSpectralPeak>[(int)Math.Ceiling(file.ScanWindow.Maximum) * BinsPerDalton + 1];
            file.InitiateDynamicConnection();

            Ms1ScanInfo[] scanInfoArray = new Ms1ScanInfo[file.NumberOfMs1Frames];

            // set the default list length to 1/25th of the number of frames
            int defaultListLength = file.NumberOfMs1Frames / 25;

            // Populate the _indexedPeaks array with the peaks from the TimsTofFileReader
            int zeroBasedMs1FrameIndex = 0;
            HashSet<int> observedRoundedMzs = new();
            // foreach frame...
            foreach (TimsDataScan ms1Scan in file.GetMs1InfoScanByScan())
            {
                observedRoundedMzs.Clear();
                // for each scan in the frame...
                for (int scanIdx = 0; scanIdx < file.NumberOfScansPerFrame; scanIdx++)
                {
                    var spectrum = ms1Scan.Ms1SpectraIndexedByZeroBasedScanNumber[scanIdx];
                    if (spectrum == null) continue; // If there are no peaks in the spectrum, continue to the next scan
                    // for each peak in the spectrum
                    for (int spectrumIdx = 0; spectrumIdx < spectrum.Size; spectrumIdx++)
                    {
                        int roundedMz = (int)Math.Round(spectrum.XArray[spectrumIdx] * BinsPerDalton, 0);
                        observedRoundedMzs.Add(roundedMz);
                        // If the list of IndexedMassSpectralPeaks doesn't exist for the given mz, create it and add a new TimsIndexedMassSpectralPeak
                        if (_indexedPeaks[roundedMz] == null)
                        {
                            _indexedPeaks[roundedMz] = new List<IndexedMassSpectralPeak>(defaultListLength);
                            _indexedPeaks[roundedMz].Add(new IndexedTimsTofPeak(spectrum.XArray[spectrumIdx], zeroBasedMs1FrameIndex, ms1Scan.RetentionTime,
                                new IonMobilityPeak(scanIdx + 1, spectrum.YArray[spectrumIdx])));
                        }

                        // Otherwise, check if the list already contains a peak for the given scan index. If it does, add a new IonMobilityPeak.
                        else if (_indexedPeaks[roundedMz][^1].ZeroBasedMs1ScanIndex == zeroBasedMs1FrameIndex)
                            ((IndexedTimsTofPeak)_indexedPeaks[roundedMz][^1]).AddIonMobilityPeak(new IonMobilityPeak(scanIdx + 1, spectrum.YArray[spectrumIdx]));

                        // If it doesn't, create and add a new TimsIndexedMassSpectralPeak.
                        else
                            _indexedPeaks[roundedMz].Add(new IndexedTimsTofPeak(spectrum.XArray[spectrumIdx], zeroBasedMs1FrameIndex, ms1Scan.RetentionTime,
                                new IonMobilityPeak(scanIdx + 1, spectrum.YArray[spectrumIdx])));
                    }
                }
                foreach (int roundedMz in observedRoundedMzs)
                    ((IndexedTimsTofPeak)_indexedPeaks[roundedMz][^1]).IonMobilityPeaks.TrimExcess();
                scanInfoArray[zeroBasedMs1FrameIndex] = new Ms1ScanInfo((int)ms1Scan.FrameId, zeroBasedMs1FrameIndex, ms1Scan.RetentionTime);
                zeroBasedMs1FrameIndex++;
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

        public IndexedMassSpectralPeak GetIndexedPeak(double theorMass, int zeroBasedScanIndex, Tolerance tolerance, int chargeState)
        {
            IndexedMassSpectralPeak bestPeak = null;
            int ceilingMz = (int)Math.Ceiling(tolerance.GetMaximumValue(theorMass).ToMz(chargeState) * BinsPerDalton);
            int floorMz = (int)Math.Floor(tolerance.GetMinimumValue(theorMass).ToMz(chargeState) * BinsPerDalton);

            for (int j = floorMz; j <= ceilingMz; j++)
            {
                if (j < _indexedPeaks.Length && _indexedPeaks[j] != null)
                {
                    List<IndexedMassSpectralPeak> bin = _indexedPeaks[j];
                    int index = BinarySearchForIndexedPeak(bin, zeroBasedScanIndex);

                    for (int i = index; i < bin.Count; i++)
                    {
                        IndexedMassSpectralPeak peak = bin[i];

                        if (peak.ZeroBasedMs1ScanIndex > zeroBasedScanIndex)
                        {
                            break;
                        }

                        double expMass = peak.Mz.ToMass(chargeState);

                        if (tolerance.Within(expMass, theorMass) && peak.ZeroBasedMs1ScanIndex == zeroBasedScanIndex
                            && (bestPeak == null || Math.Abs(expMass - theorMass) < Math.Abs(bestPeak.Mz.ToMass(chargeState) - theorMass)))
                        {
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