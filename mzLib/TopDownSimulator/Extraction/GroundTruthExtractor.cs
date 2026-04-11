using System.Linq;
using MassSpectrometry;
using MzLibUtil;
using TopDownSimulator.Model;

namespace TopDownSimulator.Extraction;

/// <summary>
/// Pulls a <see cref="ProteoformGroundTruth"/> tensor from a pre-built
/// <see cref="PeakIndexingEngine"/> for a given proteoform ID.
///
/// For each MS1 scan within (rtCenter ± rtHalfWidth) and each charge in
/// [minCharge, maxCharge], the extractor queries the index at every theoretical
/// isotopologue m/z using a ppm tolerance and records the matched peak intensity.
/// </summary>
public sealed class GroundTruthExtractor
{
    private readonly IndexingEngine<IndexedMassSpectralPeak> _index;
    private readonly double _ppmTolerance;

    public GroundTruthExtractor(IndexingEngine<IndexedMassSpectralPeak> index, double ppmTolerance = 10.0)
    {
        _index = index;
        _ppmTolerance = ppmTolerance;
    }

    public ProteoformGroundTruth Extract(
        double monoisotopicMass,
        double rtCenter,
        double rtHalfWidth,
        int minCharge,
        int maxCharge)
    {
        var scanInfos = _index.ScanInfoArray!
            .Where(s => s != null && s.MsnOrder == 1
                        && s.RetentionTime >= rtCenter - rtHalfWidth
                        && s.RetentionTime <= rtCenter + rtHalfWidth)
            .OrderBy(s => s.ZeroBasedScanIndex)
            .ToArray();

        int nScans = scanInfos.Length;
        int nCharges = maxCharge - minCharge + 1;

        var kernel = new IsotopeEnvelopeKernel(monoisotopicMass);
        int nIsotopologues = kernel.IsotopologueCount;

        var centroidMzs = new double[nCharges][];
        for (int z = minCharge; z <= maxCharge; z++)
            centroidMzs[z - minCharge] = kernel.CentroidMzs(z);

        var intensities = new double[nCharges][][];
        var chargeXics = new double[nCharges][];
        var tolerance = new PpmTolerance(_ppmTolerance);

        for (int c = 0; c < nCharges; c++)
        {
            intensities[c] = new double[nIsotopologues][];
            for (int i = 0; i < nIsotopologues; i++)
                intensities[c][i] = new double[nScans];
            chargeXics[c] = new double[nScans];
        }

        int[] scanIdxArray = new int[nScans];
        double[] scanTimes = new double[nScans];
        for (int s = 0; s < nScans; s++)
        {
            scanIdxArray[s] = scanInfos[s].ZeroBasedScanIndex;
            scanTimes[s] = scanInfos[s].RetentionTime;
        }

        for (int c = 0; c < nCharges; c++)
        {
            for (int i = 0; i < nIsotopologues; i++)
            {
                double mz = centroidMzs[c][i];
                for (int s = 0; s < nScans; s++)
                {
                    var peak = _index.GetIndexedPeak(mz, scanIdxArray[s], tolerance);
                    if (peak == null) continue;
                    double intensity = peak.Intensity;
                    intensities[c][i][s] = intensity;
                    chargeXics[c][s] += intensity;
                }
            }
        }

        return new ProteoformGroundTruth
        {
            MonoisotopicMass = monoisotopicMass,
            RetentionTimeCenter = rtCenter,
            MinCharge = minCharge,
            MaxCharge = maxCharge,
            ZeroBasedScanIndices = scanIdxArray,
            ScanTimes = scanTimes,
            CentroidMzs = centroidMzs,
            IsotopologueIntensities = intensities,
            ChargeXics = chargeXics,
        };
    }
}
