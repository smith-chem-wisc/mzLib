#nullable enable
namespace TopDownSimulator.Model;

/// <summary>
/// Parameter vector θ_p for a single proteoform in the forward model
/// I(s, b) = Σ_p A_p · g_p(t_s) · Σ_z f_p(z) · φ(b; M_p, z, σ_m) + ε.
/// </summary>
public sealed record ProteoformModel(
    double MonoisotopicMass,
    double Abundance,
    EmgProfile RtProfile,
    IChargeStateDistribution ChargeDistribution,
    string? Identifier = null);
