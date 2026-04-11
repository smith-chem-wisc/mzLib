using System;
using TopDownSimulator.Extraction;

namespace TopDownSimulator.Fitting;

public enum WidthFitMode
{
    /// <summary>Width was recovered from a multi-point profile peak shape.</summary>
    Profile,
    /// <summary>Input was centroided (≤1 peak per isotopologue window); returned the fallback.</summary>
    CentroidedFallback,
}

public sealed record FittedEnvelopeWidth(double SigmaMz, WidthFitMode Mode, int PeaksUsed);

/// <summary>
/// Estimates the per-peak m/z Gaussian width σ_m for a proteoform by looking at the
/// strongest charge's apex scan. When the input is profile data, σ_m is recovered
/// as the weighted-pooled second central moment of the peak samples around each
/// theoretical centroid. When the input is centroided (only one peak falls in each
/// window), a fallback σ_m is returned.
/// </summary>
public sealed class EnvelopeWidthFitter
{
    private readonly double _fallbackSigmaMz;

    public EnvelopeWidthFitter(double fallbackSigmaMz = 0.01)
    {
        if (fallbackSigmaMz <= 0) throw new ArgumentOutOfRangeException(nameof(fallbackSigmaMz));
        _fallbackSigmaMz = fallbackSigmaMz;
    }

    public FittedEnvelopeWidth Fit(ProteoformGroundTruth truth)
    {
        // 1. Pick the (charge, scan) cell with the highest XIC value.
        int bestC = 0, bestS = 0;
        double bestApex = double.NegativeInfinity;
        for (int c = 0; c < truth.ChargeCount; c++)
        {
            for (int s = 0; s < truth.ScanCount; s++)
            {
                double v = truth.ChargeXics[c][s];
                if (v > bestApex)
                {
                    bestApex = v;
                    bestC = c;
                    bestS = s;
                }
            }
        }

        if (bestApex <= 0)
            return new FittedEnvelopeWidth(_fallbackSigmaMz, WidthFitMode.CentroidedFallback, 0);

        int nIso = truth.CentroidMzs[bestC].Length;

        // 2. Classify as profile vs centroided: max peak count across isotopologues at apex.
        int maxPeaksPerWindow = 0;
        for (int i = 0; i < nIso; i++)
        {
            var w = truth.IsotopologuePeakWindows[bestC][i][bestS];
            if (w.Length > maxPeaksPerWindow) maxPeaksPerWindow = w.Length;
        }
        if (maxPeaksPerWindow <= 1)
            return new FittedEnvelopeWidth(_fallbackSigmaMz, WidthFitMode.CentroidedFallback, 0);

        // 3. Per-isotopologue weighted second central moment, pooled by total intensity.
        double totalWeightedVar = 0;
        double totalWeight = 0;
        int peaksUsed = 0;

        for (int i = 0; i < nIso; i++)
        {
            double centroid = truth.CentroidMzs[bestC][i];
            double leftBoundary = i == 0
                ? double.NegativeInfinity
                : 0.5 * (truth.CentroidMzs[bestC][i - 1] + centroid);
            double rightBoundary = i == nIso - 1
                ? double.PositiveInfinity
                : 0.5 * (centroid + truth.CentroidMzs[bestC][i + 1]);

            var samples = truth.IsotopologuePeakWindows[bestC][i][bestS];
            double sw = 0, swx = 0, swxx = 0;
            int n = 0;
            for (int k = 0; k < samples.Length; k++)
            {
                var p = samples[k];
                if (p.Intensity <= 0) continue;
                if (p.Mz < leftBoundary || p.Mz > rightBoundary) continue;
                sw += p.Intensity;
                swx += p.Intensity * p.Mz;
                swxx += p.Intensity * p.Mz * p.Mz;
                n++;
            }
            if (n < 2 || sw <= 0) continue;

            double mean = swx / sw;
            double variance = swxx / sw - mean * mean;
            if (variance <= 0) continue;

            totalWeightedVar += sw * variance;
            totalWeight += sw;
            peaksUsed += n;
        }

        if (totalWeight <= 0 || peaksUsed < 2)
            return new FittedEnvelopeWidth(_fallbackSigmaMz, WidthFitMode.CentroidedFallback, peaksUsed);

        double sigma = Math.Sqrt(totalWeightedVar / totalWeight);
        return sigma > 0
            ? new FittedEnvelopeWidth(sigma, WidthFitMode.Profile, peaksUsed)
            : new FittedEnvelopeWidth(_fallbackSigmaMz, WidthFitMode.CentroidedFallback, peaksUsed);
    }
}
