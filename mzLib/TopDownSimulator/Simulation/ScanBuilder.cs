using System;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;
using Readers;

namespace TopDownSimulator.Simulation;

/// <summary>
/// Converts a rasterized simulation grid into mzLib scan/file objects.
/// </summary>
public sealed class ScanBuilder
{
    public MsDataScan[] BuildMs1Scans(
        RasterizedScanGrid grid,
        bool isCentroid = false,
        Polarity polarity = Polarity.Positive,
        MZAnalyzerType analyzer = MZAnalyzerType.Orbitrap,
        string scanFilter = "synthetic")
    {
        int nScans = grid.ScanTimes.Length;
        int nMz = grid.MzGrid.Length;
        var scans = new MsDataScan[nScans];
        var mzRange = new MzRange(grid.MzGrid[0], grid.MzGrid[^1]);

        for (int s = 0; s < nScans; s++)
        {
            var intensities = new double[nMz];
            double tic = 0;
            for (int b = 0; b < nMz; b++)
            {
                double intensity = grid.Intensities[s, b];
                intensities[b] = intensity;
                tic += intensity;
            }

            scans[s] = new MsDataScan(
                massSpectrum: new MzSpectrum((double[])grid.MzGrid.Clone(), intensities, false),
                oneBasedScanNumber: s + 1,
                msnOrder: 1,
                isCentroid: isCentroid,
                polarity: polarity,
                retentionTime: grid.ScanTimes[s],
                scanWindowRange: mzRange,
                scanFilter: scanFilter,
                mzAnalyzer: analyzer,
                totalIonCurrent: tic,
                injectionTime: 1.0,
                noiseData: null,
                nativeId: $"scan={s + 1}");
        }

        return scans;
    }

    public GenericMsDataFile BuildMsDataFile(MsDataScan[] scans, SourceFile? sourceFile = null)
    {
        sourceFile ??= new SourceFile(
            nativeIdFormat: "scan number only nativeID format",
            massSpectrometerFileFormat: "mzML format",
            checkSum: null,
            fileChecksumType: null,
            id: "TopDownSimulator");

        return new GenericMsDataFile(scans, sourceFile);
    }
}
