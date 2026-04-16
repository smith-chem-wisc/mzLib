namespace Readers
{
    public interface IMzmlImspExportService
    {
        string ConvertToImspFile(string inputPath, string? outputPath = null, int binsPerDalton = MassSpectrometry.ImspExportService.DefaultBinsPerDalton,
            double intensityThreshold = MassSpectrometry.ImspExportService.DefaultIntensityThreshold);

        byte[] ConvertToImspBytes(string inputPath, int binsPerDalton = MassSpectrometry.ImspExportService.DefaultBinsPerDalton,
            double intensityThreshold = MassSpectrometry.ImspExportService.DefaultIntensityThreshold);
    }
}
