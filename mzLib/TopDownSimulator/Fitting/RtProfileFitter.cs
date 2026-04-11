using System;
using TopDownSimulator.Extraction;
using TopDownSimulator.Model;

namespace TopDownSimulator.Fitting;

public sealed record FittedRtProfile(EmgProfile Profile, double TotalIntensity, int ScansUsed);

/// <summary>
/// Fits an EMG retention-time profile (μ, σ, τ) to the summed-charge XIC of a proteoform
/// using a closed-form method-of-moments estimator.
///
/// For an exponentially-modified Gaussian:
///   mean     = μ + τ
///   variance = σ² + τ²
///   m₃       = 2·τ³
///
/// so given sample central moments (m₁, m₂, m₃) of the weighted XIC we recover
/// τ = (m₃/2)^(1/3), σ² = m₂ − τ², μ = m₁ − τ. When m₃ ≤ 0 the data is symmetric or
/// left-skewed and we fall back to τ=0 (pure Gaussian).
/// </summary>
public sealed class RtProfileFitter
{
    public FittedRtProfile Fit(ProteoformGroundTruth truth)
    {
        int nScans = truth.ScanCount;
        int nCharges = truth.ChargeCount;
        if (nScans == 0)
            throw new InvalidOperationException("Cannot fit RT profile: ground truth has no scans.");

        // Sum across charges to build a single per-scan XIC.
        double[] totals = new double[nScans];
        double sumI = 0;
        int nonZero = 0;
        for (int s = 0; s < nScans; s++)
        {
            double v = 0;
            for (int c = 0; c < nCharges; c++)
                v += truth.ChargeXics[c][s];
            totals[s] = v;
            if (v > 0)
            {
                sumI += v;
                nonZero++;
            }
        }
        if (sumI <= 0 || nonZero < 2)
            throw new InvalidOperationException("Cannot fit RT profile: XIC is empty.");

        double[] times = truth.ScanTimes;

        // First moment (mean RT).
        double m1 = 0;
        for (int s = 0; s < nScans; s++)
            m1 += totals[s] * times[s];
        m1 /= sumI;

        // Second and third central moments.
        double m2 = 0, m3 = 0;
        for (int s = 0; s < nScans; s++)
        {
            double d = times[s] - m1;
            m2 += totals[s] * d * d;
            m3 += totals[s] * d * d * d;
        }
        m2 /= sumI;
        m3 /= sumI;

        double tau, sigma, mu;
        if (m3 <= 0)
        {
            // Symmetric or left-skewed → fall back to pure Gaussian.
            tau = 0;
            sigma = Math.Sqrt(Math.Max(m2, 0));
            mu = m1;
        }
        else
        {
            tau = Math.Cbrt(m3 / 2.0);
            double sigmaSq = m2 - tau * tau;
            if (sigmaSq <= 0)
            {
                // τ overshot; clamp so σ stays small but positive.
                tau = Math.Sqrt(Math.Max(m2, 0));
                sigma = 1e-6;
            }
            else
            {
                sigma = Math.Sqrt(sigmaSq);
            }
            mu = m1 - tau;
        }

        return new FittedRtProfile(
            new EmgProfile(Mu: mu, Sigma: sigma, Tau: tau),
            TotalIntensity: sumI,
            ScansUsed: nonZero);
    }
}
