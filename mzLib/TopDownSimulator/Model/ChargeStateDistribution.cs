using System;

namespace TopDownSimulator.Model;

/// <summary>
/// f_p(z) — the relative abundance of charge state z for a single proteoform.
/// Implementations are normalized so Σ_z f(z) = 1 over a caller-supplied charge range.
/// </summary>
public interface IChargeStateDistribution
{
    double Evaluate(int z);
}

/// <summary>
/// Gaussian-on-charge: f(z) ∝ exp(-(z - μ_z)² / (2 σ_z²)).
/// Default Phase 1 charge-state distribution; upgrade to skew-normal in Phase 2
/// if residuals are systematically asymmetric.
/// </summary>
public sealed record GaussianChargeDistribution(double MuZ, double SigmaZ) : IChargeStateDistribution
{
    public double Evaluate(int z)
    {
        double d = (z - MuZ) / SigmaZ;
        return Math.Exp(-0.5 * d * d) / (SigmaZ * Math.Sqrt(2 * Math.PI));
    }
}
