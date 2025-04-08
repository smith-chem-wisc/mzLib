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
    public class PeakIndexingEngine : IndexingEngine<IndexedMassSpectralPeak>, IFlashLfqIndexingEngine
    {
        private readonly Serializer _serializer;
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

        public bool IndexPeaks()
        {
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
                IndexedPeaks = null;
                return false;
            }

            return IndexPeaks(msDataScans);
        }

        public void ClearIndex()
        {
            IndexedPeaks = null;
            GC.Collect();
        }

        public void SerializeIndex()
        {
            string dir = Path.GetDirectoryName(SpectraFile.FullFilePathWithExtension);
            string indexPath = Path.Combine(dir, SpectraFile.FilenameWithoutExtension + ".ind");

            using (var indexFile = File.Create(indexPath))
            {
                _serializer.Serialize(indexFile, IndexedPeaks);
            }
            ClearIndex();
        }

        public void DeserializeIndex()
        {
            string dir = Path.GetDirectoryName(SpectraFile.FullFilePathWithExtension);
            string indexPath = Path.Combine(dir, SpectraFile.FilenameWithoutExtension + ".ind");

            using (var indexFile = File.OpenRead(indexPath))
            {
                IndexedPeaks = (List<IndexedMassSpectralPeak>[])_serializer.Deserialize(indexFile);
            }

            File.Delete(indexPath);
        }
    }
}