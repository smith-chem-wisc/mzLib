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
            var scans = new List<MsDataScan>();
            foreach (var file in _pufFilePaths)
            {
                var pufFile = new PufResultFile(file);
                pufFile.LoadResults();
                PufFiles.Add(pufFile);

                foreach (var exp in pufFile.DataSet.Experiments)
                {
                    // Map PufMsMsExperiment to MsDataScan
                    var inst = exp.InstrumentData;
                    if (inst == null)
                        continue;

                    // Use the first intact as precursor, if available
                    var intact = inst.Intacts.FirstOrDefault();
                    double precursorMz = intact?.MzMonoisotopic ?? 0;
                    double intensity = intact?.Intensity ?? 0;
                    double retentionTime = 0; // PUF does not store RT, set to 0 or parse from comment if available

                    int scanNumber = int.Parse(exp.Id);

                    var dissociationType = Enum.TryParse<DissociationType>(inst.FragmentationMethod, true, out var dissType)
                        ? dissType
                        : DissociationType.Unknown;

                    var spectrum = (NeutralMassSpectrum)inst;
                    var scan = new MsDataScan(
                        spectrum,
                        scanNumber,
                        2,
                        true,
                        Polarity.Unknown,
                        retentionTime,
                        new MzLibUtil.MzRange(spectrum.FirstX!.Value, spectrum.LastX!.Value),
                        null,
                        MZAnalyzerType.Unknown,
                        spectrum.SumOfAllY,
                        null, 
                        null,
                        "",
                        precursorMz,
                        1,
                        intensity,
                        precursorMz,
                        1,
                        dissociationType, null, precursorMz, null, exp.Comment);

                    scans.Add(scan);
                }
            }

            Scans = scans.OrderBy(p => p.OneBasedScanNumber).ToArray();
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
