using System;
using TopDownSimulator.Extraction;
using TopDownSimulator.Model;

namespace TopDownSimulator.Fitting;

public sealed record FittedAbundance(double Abundance, double Residual, int SamplesUsed);

/// <summary>
/// Given fixed σ_m, an RT profile g(t), and a charge distribution f(z), fits the scalar
/// abundance A that minimizes the squared residual between observed centroid intensities
/// and the forward-model prediction A·g(t_s)·f(z_c)·w_i/(σ_m·√2π), where w_i is the
/// isotopologue weight from the proteoform's averagine kernel.
/// </summary>
public sealed class AbundanceFitter
{
    public FittedAbundance Fit(
        ProteoformGroundTruth truth,
        double sigmaMz,
        EmgProfile rtProfile,
        IChargeStateDistribution chargeDistribution)
    {
        if (sigmaMz <= 0) throw new ArgumentOutOfRangeException(nameof(sigmaMz));

        var kernel = new IsotopeEnvelopeKernel(truth.MonoisotopicMass);
        int nIso = kernel.IsotopologueCount;
        int nC = truth.ChargeCount;
        int nS = truth.ScanCount;

        double peakScale = 1.0 / (sigmaMz * Math.Sqrt(2 * Math.PI));

        double sumOp = 0;  // Σ obs * pred
        double sumPp = 0;  // Σ pred²
        int used = 0;

        for (int s = 0; s < nS; s++)
        {
            double gt = rtProfile.Evaluate(truth.ScanTimes[s]);
            for (int c = 0; c < nC; c++)
            {
                int z = truth.MinCharge + c;
                double fz = chargeDistribution.Evaluate(z);
                for (int i = 0; i < nIso; i++)
                {
                    double obs = truth.IsotopologueIntensities[c][i][s];
                    double pred = gt * fz * kernel.Intensity(i) * peakScale;
                    sumOp += obs * pred;
                    sumPp += pred * pred;
                    if (pred > 0) used++;
                }
            }
        }

        if (sumPp <= 0)
            throw new InvalidOperationException("Cannot fit abundance: model predicts zero everywhere.");

        double A = sumOp / sumPp;

        // Compute residual sum-of-squares under the fitted A.
        double rss = 0;
        for (int s = 0; s < nS; s++)
        {
            double gt = rtProfile.Evaluate(truth.ScanTimes[s]);
            for (int c = 0; c < nC; c++)
            {
                int z = truth.MinCharge + c;
                double fz = chargeDistribution.Evaluate(z);
                for (int i = 0; i < nIso; i++)
                {
                    double obs = truth.IsotopologueIntensities[c][i][s];
                    double pred = A * gt * fz * kernel.Intensity(i) * peakScale;
                    double d = obs - pred;
                    rss += d * d;
                }
            }
        }

        return new FittedAbundance(A, rss, used);
    }
}
