#nullable enable
using System;
using TopDownSimulator.Extraction;
using TopDownSimulator.Model;

namespace TopDownSimulator.Fitting;

public sealed record FittedAbundance(double Abundance, double Residual, int SamplesUsed);

/// <summary>
/// Given a fixed peak-width model, an RT profile g(t), and a charge distribution f(z), fits the
/// scalar abundance A that minimizes the squared residual between observed centroid intensities
/// and the forward-model prediction A·g(t_s)·f(z_c)·φ(mz_i; M, z_c).
/// </summary>
/// <remarks>
/// The envelope term is the kernel evaluated at the isotopologue's own m/z, not its normalized
/// weight over a lone Gaussian's peak height. The two agree only while neighbouring isotopologues
/// are fully resolved; once σ_m approaches the 1.00235/z spacing — which is exactly what happens at
/// high m/z and high charge — the rendered peak carries its neighbours' tails as well, and fitting
/// against the isolated-peak form would under-estimate abundance there. Using the kernel means an
/// abundance fitted under a width model reproduces the peak height rendered under that same model.
/// </remarks>
public sealed class AbundanceFitter
{
    public FittedAbundance Fit(
        ProteoformGroundTruth truth,
        IPeakWidthModel widthModel,
        EmgProfile rtProfile,
        IChargeStateDistribution chargeDistribution)
    {
        if (widthModel is null) throw new ArgumentNullException(nameof(widthModel));

        var kernel = new IsotopeEnvelopeKernel(truth.MonoisotopicMass);
        int nIso = kernel.IsotopologueCount;
        int nC = truth.ChargeCount;
        int nS = truth.ScanCount;

        // The envelope shape does not depend on scan time, so it is evaluated once per
        // (charge, isotopologue) rather than once per sample.
        var shape = new double[nC][];
        for (int c = 0; c < nC; c++)
        {
            int z = truth.MinCharge + c;
            shape[c] = new double[nIso];
            for (int i = 0; i < nIso; i++)
                shape[c][i] = kernel.Evaluate(truth.CentroidMzs[c][i], z, widthModel);
        }

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
                    double pred = gt * fz * shape[c][i];
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
                    double pred = A * gt * fz * shape[c][i];
                    double d = obs - pred;
                    rss += d * d;
                }
            }
        }

        return new FittedAbundance(A, rss, used);
    }

    public FittedAbundance Fit(
        ProteoformGroundTruth truth,
        double sigmaMz,
        EmgProfile rtProfile,
        IChargeStateDistribution chargeDistribution) =>
        Fit(truth, new ConstantPeakWidth(sigmaMz), rtProfile, chargeDistribution);
}
