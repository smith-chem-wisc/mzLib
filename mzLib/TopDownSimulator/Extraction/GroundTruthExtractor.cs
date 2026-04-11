using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;
using TopDownSimulator.Model;

namespace TopDownSimulator.Extraction;

/// <summary>
/// Pulls a <see cref="ProteoformGroundTruth"/> tensor from a pre-built peak
/// <see cref="IndexingEngine{T}"/> for a given proteoform.
///
/// For each MS1 scan within (rtCenter ± rtHalfWidth) and each charge in
/// [minCharge, maxCharge], the extractor collects every peak that falls within
/// ±<c>mzWindowHalfWidth</c> of each theoretical isotopologue centroid, and also
/// records the single closest peak's intensity in the scalar
/// <see cref="ProteoformGroundTruth.IsotopologueIntensities"/> tensor.
/// </summary>
public sealed class GroundTruthExtractor
{
    private readonly IndexingEngine<IndexedMassSpectralPeak> _index;
    private readonly Dictionary<int, MzSpectrum> _spectraByScanIndex;
    private readonly double _ppmTolerance;
    private readonly double _mzWindowHalfWidth;

    public GroundTruthExtractor(
        IndexingEngine<IndexedMassSpectralPeak> index,
        MsDataScan[] scans,
        double ppmTolerance = 10.0,
        double mzWindowHalfWidth = 0.05)
    {
        _index = index;
        _ppmTolerance = ppmTolerance;
        _mzWindowHalfWidth = mzWindowHalfWidth;

        _spectraByScanIndex = new Dictionary<int, MzSpectrum>(scans.Length);
        for (int i = 0; i < scans.Length; i++)
        {
            if (scans[i] != null)
                _spectraByScanIndex[i] = scans[i].MassSpectrum;
        }
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
        var peakWindows = new PeakSample[nCharges][][][];
        var chargeXics = new double[nCharges][];

        for (int c = 0; c < nCharges; c++)
        {
            intensities[c] = new double[nIsotopologues][];
            peakWindows[c] = new PeakSample[nIsotopologues][][];
            for (int i = 0; i < nIsotopologues; i++)
            {
                intensities[c][i] = new double[nScans];
                peakWindows[c][i] = new PeakSample[nScans][];
            }
            chargeXics[c] = new double[nScans];
        }

        int[] scanIdxArray = new int[nScans];
        double[] scanTimes = new double[nScans];
        for (int s = 0; s < nScans; s++)
        {
            scanIdxArray[s] = scanInfos[s].ZeroBasedScanIndex;
            scanTimes[s] = scanInfos[s].RetentionTime;
        }

        var windowTolerance = new AbsoluteTolerance(_mzWindowHalfWidth);
        var ppmTolerance = new PpmTolerance(_ppmTolerance);

        for (int s = 0; s < nScans; s++)
        {
            int scanIdx = scanIdxArray[s];
            if (!_spectraByScanIndex.TryGetValue(scanIdx, out var spectrum))
                continue;

            for (int c = 0; c < nCharges; c++)
            {
                for (int i = 0; i < nIsotopologues; i++)
                {
                    double centroid = centroidMzs[c][i];

                    // 1) Collect every peak within the absolute m/z window.
                    var indices = spectrum.GetPeakIndicesWithinTolerance(centroid, windowTolerance);
                    PeakSample[] samples;
                    if (indices.Count == 0)
                    {
                        samples = System.Array.Empty<PeakSample>();
                    }
                    else
                    {
                        samples = new PeakSample[indices.Count];
                        for (int k = 0; k < indices.Count; k++)
                        {
                            int idx = indices[k];
                            samples[k] = new PeakSample(spectrum.XArray[idx], spectrum.YArray[idx]);
                        }
                    }
                    peakWindows[c][i][s] = samples;

                    // 2) Scalar intensity = peak closest to centroid and within ppm tolerance.
                    double bestIntensity = 0;
                    double bestDelta = double.MaxValue;
                    for (int k = 0; k < samples.Length; k++)
                    {
                        if (!ppmTolerance.Within(samples[k].Mz, centroid)) continue;
                        double delta = System.Math.Abs(samples[k].Mz - centroid);
                        if (delta < bestDelta)
                        {
                            bestDelta = delta;
                            bestIntensity = samples[k].Intensity;
                        }
                    }
                    intensities[c][i][s] = bestIntensity;
                    chargeXics[c][s] += bestIntensity;
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
            IsotopologuePeakWindows = peakWindows,
            ChargeXics = chargeXics,
            MzWindowHalfWidth = _mzWindowHalfWidth,
        };
    }
}
