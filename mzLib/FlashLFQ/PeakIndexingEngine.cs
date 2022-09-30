using Chemistry;
using Readers;
using MassSpectrometry;
using MzLibUtil;
using NetSerializer;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

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
                typeof(List<IndexedMassSpectralPeak>[]), typeof(List<IndexedMassSpectralPeak>),
                typeof(IndexedMassSpectralPeak)
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
            var reader = ReaderCreator.CreateReader(fileName); 
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