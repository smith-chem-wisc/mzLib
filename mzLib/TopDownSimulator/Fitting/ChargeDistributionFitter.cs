using System;
using TopDownSimulator.Extraction;
using TopDownSimulator.Model;

namespace TopDownSimulator.Fitting;

public sealed record FittedChargeDistribution(
    GaussianChargeDistribution Distribution,
    double[] ApexIntensitiesByCharge,
    int ChargesUsed);

/// <summary>
/// Fits a Gaussian-on-charge distribution (μ_z, σ_z) to per-charge XIC apex intensities
/// using a weighted method-of-moments estimator. When fewer than two charges carry
/// signal, σ_z falls back to a user-supplied default (1.0 by default).
/// </summary>
public sealed class ChargeDistributionFitter
{
    private readonly double _fallbackSigmaZ;

    public ChargeDistributionFitter(double fallbackSigmaZ = 1.0)
    {
        if (fallbackSigmaZ <= 0) throw new ArgumentOutOfRangeException(nameof(fallbackSigmaZ));
        _fallbackSigmaZ = fallbackSigmaZ;
    }

    public FittedChargeDistribution Fit(ProteoformGroundTruth truth)
    {
        int nCharges = truth.ChargeCount;
        int nScans = truth.ScanCount;

        // Per-charge apex intensity from its XIC.
        double[] apex = new double[nCharges];
        int chargesUsed = 0;
        for (int c = 0; c < nCharges; c++)
        {
            double best = 0;
            for (int s = 0; s < nScans; s++)
            {
                if (truth.ChargeXics[c][s] > best) best = truth.ChargeXics[c][s];
            }
            apex[c] = best;
            if (best > 0) chargesUsed++;
        }

        if (chargesUsed == 0)
            throw new InvalidOperationException("Cannot fit charge distribution: all charge XICs are zero.");

        double sumW = 0, sumWz = 0;
        for (int c = 0; c < nCharges; c++)
        {
            int z = truth.MinCharge + c;
            sumW += apex[c];
            sumWz += apex[c] * z;
        }
        double muZ = sumWz / sumW;

        double sigmaZ;
        if (chargesUsed < 2)
        {
            sigmaZ = _fallbackSigmaZ;
        }
        else
        {
            double sumWdd = 0;
            for (int c = 0; c < nCharges; c++)
            {
                int z = truth.MinCharge + c;
                double d = z - muZ;
                sumWdd += apex[c] * d * d;
            }
            double var = sumWdd / sumW;
            sigmaZ = var > 0 ? Math.Sqrt(var) : _fallbackSigmaZ;
        }

        return new FittedChargeDistribution(
            new GaussianChargeDistribution(muZ, sigmaZ),
            apex,
            chargesUsed);
    }
}
