﻿using FlashLFQ.Interfaces;
using MassSpectrometry;
using MzLibUtil;
using NetSerializer;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;

namespace FlashLFQ
{
    public class PeakIndexingEngine : IndexingEngine<IndexedMassSpectralPeak>, IFlashLfqIndexingEngine
    {
        private readonly Serializer _serializer;
        public SpectraFileInfo SpectraFile { get; private set; }
        internal PeakIndexingEngine()
        {
            var messageTypes = new List<Type>
            {
                typeof(List<IndexedMassSpectralPeak>[]), typeof(List<IndexedMassSpectralPeak>),
                typeof(IndexedMassSpectralPeak)
            };
            _serializer = new Serializer(messageTypes);
        }

        /// <summary>
        /// This factory method returns an IndexingEngine instance where the peaks in all MS1 scans have been indexed. 
        /// This method ignores MS2 scans when indexing
        /// </summary>
        public static PeakIndexingEngine? InitializeIndexingEngine(SpectraFileInfo file)
        {
            // read spectra file
            string fileName = file.FullFilePathWithExtension;
            var reader = MsDataFileReader.GetDataFile(fileName);
            reader.LoadAllStaticData();

            var peakIndexingEngine = InitializeIndexingEngine(reader);
            if(peakIndexingEngine != null) peakIndexingEngine.SpectraFile = file;
            return peakIndexingEngine;
        }

        /// <summary>
        /// This factory method returns an IndexingEngine instance where the peaks in all MS1 scans have been indexed. 
        /// This method ignores MS2 scans when indexing
        /// </summary>
        public static PeakIndexingEngine? InitializeIndexingEngine(MsDataFile dataFile)
        {
            var scanArray = dataFile.GetMS1Scans()
                .Where(i => i != null && i.MsnOrder == 1)
                .OrderBy(i => i.OneBasedScanNumber)
                .ToArray();
            return InitializeIndexingEngine(scanArray);
        }

        /// <summary>
        /// Read in all spectral peaks from the scanArray, index the peaks based on mass and retention time, 
        /// and store them in a jagged array of Lists containing all peaks within a particular mass range
        /// </summary>
        /// <param name="scanArray">An array of raw data scans</param>
        public static PeakIndexingEngine? InitializeIndexingEngine(MsDataScan[] scanArray)
        {
            PeakIndexingEngine newEngine = new();
            if (newEngine.IndexPeaks(scanArray))
                return newEngine;
            return null;
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

        /// <summary>  
        /// Prune the index engine to remove any unnecessary data or entries for the better memory usage.  
        /// </summary>  
        public void PruneIndex(List<float> targetMz)
        {
            PpmTolerance ppmTolerance = new PpmTolerance(10); // Default tolerance, can be adjusted as needed
            if (IndexedPeaks == null || targetMz == null || !targetMz.Any())
                return;
            List<IndexedMassSpectralPeak>[] indexedPeaks = new List<IndexedMassSpectralPeak>[IndexedPeaks.Length];
            var maxIndex = IndexedPeaks.Length - 1;

            foreach (var mz in targetMz)
            {
                int ceilingMz = (int)Math.Ceiling(ppmTolerance.GetMaximumValue(mz) * BinsPerDalton);
                int floorMz = (int)Math.Floor(ppmTolerance.GetMinimumValue(mz) * BinsPerDalton);
                if (ceilingMz > maxIndex || floorMz > maxIndex)
                {
                    // If the mz is out of bounds, skip it
                    continue;
                }

                for (int i = floorMz; i <= ceilingMz; i++)
                {
                    if (indexedPeaks[i] == null)
                    {
                        indexedPeaks[i] = IndexedPeaks[i];
                    }

                }
            }
            IndexedPeaks = indexedPeaks;
        }
    }
}