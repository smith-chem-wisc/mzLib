#nullable enable
using System;

namespace TopDownSimulator.Model;

/// <summary>
/// σ_m as a function of m/z — the Gaussian width of a single simulated peak.
/// </summary>
/// <remarks>
/// Peak width relative to isotopologue spacing (1.00235 / z) is what decides whether an envelope
/// is resolved or merged, so it is the difficulty knob of any feature-finding benchmark built on
/// simulated data. A single scalar σ makes every charge state at every m/z equally easy, which is
/// not what an instrument does.
/// </remarks>
public interface IPeakWidthModel
{
    /// <summary>σ_m of a peak centred at <paramref name="mz"/>.</summary>
    double SigmaAt(double mz);

    /// <summary>
    /// An upper bound on σ over [<paramref name="minMz"/>, <paramref name="maxMz"/>]. Callers use
    /// this to size evaluation windows and grid padding conservatively; under-estimating it
    /// silently truncates envelope tails.
    /// </summary>
    double MaxSigmaOver(double minMz, double maxMz);
}

/// <summary>
/// σ_m = constant, independent of m/z. The historical behaviour of the simulator, kept so existing
/// fixtures stay meaningful and so a run can be made deliberately uniform.
/// </summary>
public sealed record ConstantPeakWidth(double SigmaMz) : IPeakWidthModel
{
    public double SigmaMz { get; init; } = Validate(SigmaMz);

    public double SigmaAt(double mz) => SigmaMz;

    public double MaxSigmaOver(double minMz, double maxMz) => SigmaMz;

    public override string ToString() => $"Constant(sigma={SigmaMz:R})";

    private static double Validate(double sigmaMz) =>
        sigmaMz > 0 && double.IsFinite(sigmaMz)
            ? sigmaMz
            : throw new ArgumentOutOfRangeException(nameof(sigmaMz), sigmaMz, "σ_m must be finite and positive.");
}

/// <summary>
/// σ_m = K · (m/z)^1.5, the Orbitrap width law.
/// </summary>
/// <remarks>
/// Resolving power falls as R ∝ (m/z)^−1/2 and R ≡ (m/z)/FWHM, so FWHM ∝ (m/z)^1.5 and
/// σ = FWHM / (2√(2 ln 2)) = K · (m/z)^1.5. One free parameter, fitted from data by
/// <see cref="Fitting.PeakWidthModelFitter"/>.
/// </remarks>
public sealed record OrbitrapPeakWidth(double K) : IPeakWidthModel
{
    /// <summary>The exponent of the width law. Fixed rather than fitted; see PeakWidthModelFitter.</summary>
    public const double Exponent = 1.5;

    /// <summary>FWHM = σ · 2√(2 ln 2).</summary>
    public const double FwhmPerSigma = 2.354820045030949;

    public double K { get; init; } = Validate(K);

    public double SigmaAt(double mz) => mz <= 0 ? 0 : K * Math.Pow(mz, Exponent);

    // σ increases monotonically with m/z, so the bound is attained at the top of the span.
    public double MaxSigmaOver(double minMz, double maxMz) => SigmaAt(Math.Max(minMz, maxMz));

    /// <summary>The model that passes through a single (m/z, σ) observation.</summary>
    public static OrbitrapPeakWidth FromSigmaAt(double mz, double sigmaMz)
    {
        if (mz <= 0 || !double.IsFinite(mz))
            throw new ArgumentOutOfRangeException(nameof(mz), mz, "m/z must be finite and positive.");
        if (sigmaMz <= 0 || !double.IsFinite(sigmaMz))
            throw new ArgumentOutOfRangeException(nameof(sigmaMz), sigmaMz, "σ_m must be finite and positive.");

        return new OrbitrapPeakWidth(sigmaMz / Math.Pow(mz, Exponent));
    }

    public override string ToString() => $"Orbitrap(k={K:R})";

    private static double Validate(double k) =>
        k > 0 && double.IsFinite(k)
            ? k
            : throw new ArgumentOutOfRangeException(nameof(k), k, "k must be finite and positive.");
}
