using System.Collections.Generic;
using System.IO;

namespace MassSpectrometry
{
    public interface IImspExportService
    {
        int WriteFile(IEnumerable<MsDataScan> ms1Scans, string outputPath, int binsPerDalton = ImspExportService.DefaultBinsPerDalton,
            double intensityThreshold = ImspExportService.DefaultIntensityThreshold);

        byte[] WriteBytes(IEnumerable<MsDataScan> ms1Scans, int binsPerDalton = ImspExportService.DefaultBinsPerDalton,
            double intensityThreshold = ImspExportService.DefaultIntensityThreshold);

        int WriteStream(IEnumerable<MsDataScan> ms1Scans, Stream outputStream, int binsPerDalton = ImspExportService.DefaultBinsPerDalton,
            double intensityThreshold = ImspExportService.DefaultIntensityThreshold);
    }
}
