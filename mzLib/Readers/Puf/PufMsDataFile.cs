using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Readers.Puf
{
    /// <summary>
    /// MsDataFile implementation for a directory of PUF files.
    /// Each PUF experiment is mapped to an MsDataScan.
    /// </summary>
    public class PufMsDataFile : MsDataFile
    {
        private readonly List<string> _pufFilePaths = new();
        public readonly List<PufResultFile> PufFiles = new();
        protected MsDataScan[] IndexedScans { get; set; }

        public PufMsDataFile(string directoryPath) : base(directoryPath)
        {
            if (!Directory.Exists(directoryPath))
                throw new DirectoryNotFoundException($"PUF directory not found: {directoryPath}");

            _pufFilePaths = Directory.GetFiles(directoryPath, "*.puf").OrderBy(f => f).ToList();
        }

        public override MsDataFile LoadAllStaticData(FilteringParams filteringParams = null, int maxThreads = 1)
        {
            foreach (var file in _pufFilePaths)
            {
                var pufFile = new PufResultFile(file);
                pufFile.LoadResults();
                PufFiles.Add(pufFile);
            }


            SourceFile = GetSourceFile();

            // ensures that if a scan (OneBasedScanNumber) does not exist,
            // the final scans array will contain a null value  
            // this unique case is due to the nature of loading MGF files
            var orderedScans = PufFiles.SelectMany(p => p.Scans)
                .OrderBy(x => x.OneBasedScanNumber)
                .ToArray();

            var indexedScans = new MsDataScan[orderedScans[^1].OneBasedScanNumber];
            foreach (var scan in orderedScans)
                indexedScans[scan.OneBasedScanNumber - 1] = scan;

            IndexedScans = indexedScans;
            Scans = orderedScans;
            return this;
        }

        public override SourceFile GetSourceFile()
        {
            return new SourceFile("no nativeID format", "PUF directory", null, null, FilePath, "PUFdir");
        }
        public override MsDataScan GetOneBasedScan(int scanNumber)
        {
            var scan = IndexedScans[scanNumber - 1];

            // If attempt to get scan is erroneously by index and not scan number return that.
            // This is done in mzml writing. 
            if (scan == null && scanNumber < Scans.Length)
                scan = Scans[scanNumber - 1];

            return scan;
        }

        public override MsDataScan GetOneBasedScanFromDynamicConnection(int oneBasedScanNumber, IFilteringParams filterParams = null)
        {
            // For simplicity, just return from Scans (no dynamic connection for PUF)
            return GetOneBasedScan(oneBasedScanNumber);
        }

        public override void CloseDynamicConnection()
        {
            // No dynamic connection for PUF
        }

        public override void InitiateDynamicConnection()
        {
            // No dynamic connection for PUF
        }
    }
}
