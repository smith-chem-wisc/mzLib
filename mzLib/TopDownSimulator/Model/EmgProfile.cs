using System;

namespace TopDownSimulator.Model;

/// <summary>
/// Exponentially-modified Gaussian RT elution profile g_p(t).
/// Parameters: μ (center), σ (Gaussian width), τ (exponential tail).
/// PDF is normalized so ∫ g(t) dt = 1; callers scale by abundance A_p.
/// </summary>
public sealed record EmgProfile(double Mu, double Sigma, double Tau)
{
    public double Evaluate(double t)
    {
        if (Tau <= 0)
        {
            double z0 = (t - Mu) / Sigma;
            return Math.Exp(-0.5 * z0 * z0) / (Sigma * Math.Sqrt(2 * Math.PI));
        }

        double arg = (Mu - t) / Tau + Sigma * Sigma / (2 * Tau * Tau);
        double z = (t - Mu - Sigma * Sigma / Tau) / (Sigma * Math.Sqrt(2));
        return Math.Exp(arg) / (2 * Tau) * Erfc(-z);
    }

    private static double Erfc(double x) => 1.0 - Erf(x);

    // Abramowitz & Stegun 7.1.26 — accurate to ~1.5e-7, adequate for RT modeling.
    private static double Erf(double x)
    {
        double sign = Math.Sign(x);
        x = Math.Abs(x);
        const double a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
        const double a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;
        double t = 1.0 / (1.0 + p * x);
        double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.Exp(-x * x);
        return sign * y;
    }
}
