using System;
using System.IO;
using System.Linq;
using MassSpectrometry;

namespace Readers
{
    public sealed class MzmlImspExportService : IMzmlImspExportService
    {
        private readonly IImspExportService imspExportService;

        public MzmlImspExportService()
            : this(new ImspExportService())
        {
        }

        public MzmlImspExportService(IImspExportService imspExportService)
        {
            this.imspExportService = imspExportService ?? throw new ArgumentNullException(nameof(imspExportService));
        }

        public string ConvertToImspFile(string inputPath, string? outputPath = null,
            int binsPerDalton = ImspExportService.DefaultBinsPerDalton,
            double intensityThreshold = ImspExportService.DefaultIntensityThreshold)
        {
            ValidateMzmlPath(inputPath);
            outputPath = string.IsNullOrWhiteSpace(outputPath)
                ? Path.ChangeExtension(inputPath, ".imsp")
                : outputPath;

            var ms1Scans = LoadMs1Scans(inputPath);
            imspExportService.WriteFile(ms1Scans, outputPath, binsPerDalton, intensityThreshold);
            return outputPath;
        }

        public byte[] ConvertToImspBytes(string inputPath, int binsPerDalton = ImspExportService.DefaultBinsPerDalton,
            double intensityThreshold = ImspExportService.DefaultIntensityThreshold)
        {
            ValidateMzmlPath(inputPath);
            var ms1Scans = LoadMs1Scans(inputPath);
            return imspExportService.WriteBytes(ms1Scans, binsPerDalton, intensityThreshold);
        }

        private static void ValidateMzmlPath(string inputPath)
        {
            if (string.IsNullOrWhiteSpace(inputPath))
            {
                throw new ArgumentException("Input path is required.", nameof(inputPath));
            }

            if (!string.Equals(Path.GetExtension(inputPath), ".mzML", StringComparison.OrdinalIgnoreCase))
            {
                throw new NotSupportedException("Only mzML input files are supported for IMSP export.");
            }
        }

        private static MsDataScan[] LoadMs1Scans(string inputPath)
        {
            var reader = MsDataFileReader.GetDataFile(inputPath);
            reader.LoadAllStaticData();
            return reader.GetMS1Scans()
                .OrderBy(s => s.OneBasedScanNumber)
                .ToArray();
        }
    }
}
