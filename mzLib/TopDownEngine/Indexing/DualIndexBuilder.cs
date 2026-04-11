using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Linq;

namespace TopDownEngine.Indexing;

public class DualIndexBuilder
{
    public (PeakIndexingEngine fine, ThickIndexView thick) Build(MsDataFile file)
    {
        if (file == null)
        {
            throw new ArgumentNullException(nameof(file));
        }

        MsDataScan[] scans = file.GetMS1Scans()
            .Where(scan => scan != null && scan.MsnOrder == 1)
            .OrderBy(scan => scan.OneBasedScanNumber)
            .ToArray();

        PeakIndexingEngine fine = PeakIndexingEngine.InitializeIndexingEngine(scans)
            ?? throw new MzLibException("Error: Could not initialize fine peak index from input file.");

        ThickIndexView thick = new(fine);
        return (fine, thick);
    }
}
