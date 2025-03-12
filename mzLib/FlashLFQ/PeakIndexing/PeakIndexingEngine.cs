using Chemistry;
using Readers;
using MassSpectrometry;
using MzLibUtil;
using NetSerializer;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using FlashLFQ.Interfaces;

namespace FlashLFQ
{
    public class PeakIndexingEngine : IIndexingEngine
    {
        private List<IndexedMassSpectralPeak>[] _indexedPeaks;
        private readonly Serializer _serializer;
        private const int BinsPerDalton = 100;
        public Ms1ScanInfo[] Ms1ScanInfoArray { get; private set; }
        public SpectraFileInfo SpectraFile { get; init; }

        public PeakIndexingEngine(SpectraFileInfo file)
        {
            var messageTypes = new List<Type>
            {
                typeof(List<IndexedMassSpectralPeak>[]), typeof(List<IndexedMassSpectralPeak>),
                typeof(IndexedMassSpectralPeak)
            };
            _serializer = new Serializer(messageTypes);
            SpectraFile = file;
        }

        public bool IndexPeaks(bool silent)
        {
            if (!silent)
            {
                Console.WriteLine("Reading spectra file");
            }

            // read spectra file
            string fileName = SpectraFile.FullFilePathWithExtension;
            var reader = MsDataFileReader.GetDataFile(fileName); 
            reader.LoadAllStaticData();
            // retrieve only the ms1s. 
            MsDataScan[] msDataScans = reader.GetMS1Scans()
                .Where(i => i != null && i.MsnOrder == 1)
                .OrderBy(i => i.OneBasedScanNumber)
                .ToArray(); 
            
            if (msDataScans.All(p => p == null))
            {
                _indexedPeaks = null;
                return false;
            }

            _indexedPeaks = new List<IndexedMassSpectralPeak>[(int)Math.Ceiling(msDataScans.Where(p => p != null
                && p.MassSpectrum.LastX != null).Max(p => p.MassSpectrum.LastX.Value) * BinsPerDalton) + 1];
            Ms1ScanInfoArray = new Ms1ScanInfo[msDataScans.Length];

            List<Ms1ScanInfo> scanInfo = new List<Ms1ScanInfo>();

            for (int scanIndex = 0; scanIndex < msDataScans.Length; scanIndex++)
            {
                Ms1ScanInfoArray[scanIndex] = new Ms1ScanInfo(msDataScans[scanIndex].OneBasedScanNumber, scanIndex, msDataScans[scanIndex].RetentionTime);

                for (int j = 0; j < msDataScans[scanIndex].MassSpectrum.XArray.Length; j++)
                {
                    int roundedMz = (int)Math.Round(msDataScans[scanIndex].MassSpectrum.XArray[j] * BinsPerDalton, 0);
                    _indexedPeaks[roundedMz] ??= new List<IndexedMassSpectralPeak>();
                    _indexedPeaks[roundedMz].Add(
                        new IndexedMassSpectralPeak(
                            msDataScans[scanIndex].MassSpectrum.XArray[j],
                            msDataScans[scanIndex].MassSpectrum.YArray[j], 
                            scanIndex, 
                            msDataScans[scanIndex].RetentionTime));
                }

            }

            if (_indexedPeaks == null || _indexedPeaks.Length == 0)
            {
                if (!silent)
                {
                    Console.WriteLine("FlashLFQ Error: The file " + SpectraFile.FilenameWithoutExtension + " contained no MS1 peaks!");
                }

                return false;
            }

            return true;
        }

        public void ClearIndex()
        {
            _indexedPeaks = null;
            GC.Collect();
        }

        public void SerializeIndex()
        {
            string dir = Path.GetDirectoryName(SpectraFile.FullFilePathWithExtension);
            string indexPath = Path.Combine(dir, SpectraFile.FilenameWithoutExtension + ".ind");

            using (var indexFile = File.Create(indexPath))
            {
                _serializer.Serialize(indexFile, _indexedPeaks);
            }
            ClearIndex();
        }

        public void DeserializeIndex()
        {
            string dir = Path.GetDirectoryName(SpectraFile.FullFilePathWithExtension);
            string indexPath = Path.Combine(dir, SpectraFile.FilenameWithoutExtension + ".ind");

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