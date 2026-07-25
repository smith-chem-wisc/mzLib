using System;
using System.Collections.Generic;
using MassSpectrometry;
using MzLibUtil;
using TopDownSimulator.Model;

namespace TopDownSimulator.Extraction;

/// <summary>
/// Pulls a <see cref="ProteoformGroundTruth"/> tensor straight from the MS1 scans of a run.
///
/// For each MS1 scan within (rtCenter ± rtHalfWidth) and each charge in
/// [minCharge, maxCharge], the extractor records the intensity of the peak closest to each
/// theoretical isotopologue centroid, and collects the full set of peaks around each centroid at
/// the apex (charge, scan) cell.
/// </summary>
/// <remarks>
/// This deliberately does not use <c>IndexingEngine</c>. All lookups are binary searches inside a
/// single scan's <see cref="MzSpectrum"/>, so a binned peak index buys nothing here and building
/// one over a full run costs seconds and hundreds of megabytes.
/// </remarks>
public sealed class GroundTruthExtractor
{
    private readonly int[] _ms1ScanIndices;
    private readonly double[] _ms1RetentionTimes;
    private readonly MzSpectrum[] _ms1Spectra;
    private readonly double _ppmTolerance;
    private readonly double _mzWindowHalfWidth;

    public GroundTruthExtractor(
        MsDataScan[] scans,
        double ppmTolerance = 10.0,
        double mzWindowHalfWidth = 0.05)
    {
        if (scans is null)
            throw new ArgumentNullException(nameof(scans));

        _ppmTolerance = ppmTolerance;
        _mzWindowHalfWidth = mzWindowHalfWidth;

        var indices = new List<int>(scans.Length);
        var retentionTimes = new List<double>(scans.Length);
        var spectra = new List<MzSpectrum>(scans.Length);

        // Scan index is the position in this array, matching what the peak index used to report as
        // ZeroBasedScanIndex. Ascending iteration keeps the output ordered by scan index.
        for (int i = 0; i < scans.Length; i++)
        {
            if (scans[i] is null || scans[i].MsnOrder != 1)
                continue;

            indices.Add(i);
            retentionTimes.Add(scans[i].RetentionTime);
            spectra.Add(scans[i].MassSpectrum);
        }

        _ms1ScanIndices = indices.ToArray();
        _ms1RetentionTimes = retentionTimes.ToArray();
        _ms1Spectra = spectra.ToArray();
    }

    public ProteoformGroundTruth Extract(
        double monoisotopicMass,
        double rtCenter,
        double rtHalfWidth,
        int minCharge,
        int maxCharge)
    {
        double rtLow = rtCenter - rtHalfWidth;
        double rtHigh = rtCenter + rtHalfWidth;

        int nScans = 0;
        for (int i = 0; i < _ms1RetentionTimes.Length; i++)
        {
            if (_ms1RetentionTimes[i] >= rtLow && _ms1RetentionTimes[i] <= rtHigh)
                nScans++;
        }

        int nCharges = maxCharge - minCharge + 1;

        var kernel = new IsotopeEnvelopeKernel(monoisotopicMass);
        int nIsotopologues = kernel.IsotopologueCount;

        var centroidMzs = new double[nCharges][];
        for (int z = minCharge; z <= maxCharge; z++)
            centroidMzs[z - minCharge] = kernel.CentroidMzs(z);

        var scanIdxArray = new int[nScans];
        var scanTimes = new double[nScans];
        var windowSpectra = new MzSpectrum[nScans];
        for (int i = 0, s = 0; i < _ms1RetentionTimes.Length; i++)
        {
            if (_ms1RetentionTimes[i] < rtLow || _ms1RetentionTimes[i] > rtHigh)
                continue;

            scanIdxArray[s] = _ms1ScanIndices[i];
            scanTimes[s] = _ms1RetentionTimes[i];
            windowSpectra[s] = _ms1Spectra[i];
            s++;
        }

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
                for (int s = 0; s < nScans; s++)
                    peakWindows[c][i][s] = Array.Empty<PeakSample>();
            }
            chargeXics[c] = new double[nScans];
        }

        var windowTolerance = new AbsoluteTolerance(_mzWindowHalfWidth);
        var ppmTolerance = new PpmTolerance(_ppmTolerance);

        // Pass 1: scalar intensities and charge XICs over the whole window. Peak samples are read
        // but not retained — materializing all nScans x nCharges x nIsotopologues of them is the
        // dominant allocation in this class and only one (charge, scan) column is ever consumed.
        for (int s = 0; s < nScans; s++)
        {
            var spectrum = windowSpectra[s];
            if (spectrum is null)
                continue;

            for (int c = 0; c < nCharges; c++)
            {
                for (int i = 0; i < nIsotopologues; i++)
                {
                    double centroid = centroidMzs[c][i];
                    var indices = spectrum.GetPeakIndicesWithinTolerance(centroid, windowTolerance);

                    double bestIntensity = 0;
                    double bestDelta = double.MaxValue;
                    for (int k = 0; k < indices.Count; k++)
                    {
                        int idx = indices[k];
                        double mz = spectrum.XArray[idx];
                        if (!ppmTolerance.Within(mz, centroid)) continue;

                        double delta = Math.Abs(mz - centroid);
                        if (delta < bestDelta)
                        {
                            bestDelta = delta;
                            bestIntensity = spectrum.YArray[idx];
                        }
                    }

                    intensities[c][i][s] = bestIntensity;
                    chargeXics[c][s] += bestIntensity;
                }
            }
        }

        // Pass 2: collect the peak windows for the apex cell only. The selection matches
        // EnvelopeWidthFitter's, which is the sole consumer.
        int apexCharge = 0, apexScan = 0;
        double apexValue = double.NegativeInfinity;
        for (int c = 0; c < nCharges; c++)
        {
            for (int s = 0; s < nScans; s++)
            {
                double v = chargeXics[c][s];
                if (v > apexValue)
                {
                    apexValue = v;
                    apexCharge = c;
                    apexScan = s;
                }
            }
        }

        if (apexValue > 0 && windowSpectra[apexScan] is not null)
        {
            var spectrum = windowSpectra[apexScan];
            for (int i = 0; i < nIsotopologues; i++)
            {
                double centroid = centroidMzs[apexCharge][i];
                var indices = spectrum.GetPeakIndicesWithinTolerance(centroid, windowTolerance);
                if (indices.Count == 0)
                    continue;

                var samples = new PeakSample[indices.Count];
                for (int k = 0; k < indices.Count; k++)
                {
                    int idx = indices[k];
                    samples[k] = new PeakSample(spectrum.XArray[idx], spectrum.YArray[idx]);
                }

                peakWindows[apexCharge][i][apexScan] = samples;
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
            ApexChargeIndex = apexValue > 0 ? apexCharge : -1,
            ApexScanIndex = apexValue > 0 ? apexScan : -1,
        };
    }
}
