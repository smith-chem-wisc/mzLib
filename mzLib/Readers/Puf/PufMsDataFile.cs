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

        public PufMsDataFile(string directoryPath) : base(directoryPath)
        {
            if (!Directory.Exists(directoryPath))
                throw new DirectoryNotFoundException($"PUF directory not found: {directoryPath}");

            SourceFile = GetSourceFile();
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

            Scans = PufFiles.SelectMany(p => p.Scans)
                .OrderBy(p => p.OneBasedScanNumber)
                .ToArray();

            return this;
        }

        public override SourceFile GetSourceFile()
        {
            return new SourceFile("no nativeID format", "PUF directory", null, null, FilePath, "PUFdir");
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
