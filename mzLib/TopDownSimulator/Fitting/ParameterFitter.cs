#nullable enable
using TopDownSimulator.Extraction;
using TopDownSimulator.Model;

namespace TopDownSimulator.Fitting;

public sealed record FittedProteoform(
    ProteoformModel Model,
    double SigmaMz,
    WidthFitMode WidthMode,
    int WidthPeaksUsed,
    double Residual);

/// <summary>
/// Top-level Phase 1 fitter. Runs the per-factor fitters in order —
/// envelope width → RT profile → charge distribution → abundance —
/// and packages the result into a <see cref="ProteoformModel"/>.
/// </summary>
public sealed class ParameterFitter
{
    private readonly EnvelopeWidthFitter _widthFitter;
    private readonly RtProfileFitter _rtFitter;
    private readonly ChargeDistributionFitter _chargeFitter;
    private readonly AbundanceFitter _abundanceFitter;

    public ParameterFitter(
        EnvelopeWidthFitter? widthFitter = null,
        RtProfileFitter? rtFitter = null,
        ChargeDistributionFitter? chargeFitter = null,
        AbundanceFitter? abundanceFitter = null)
    {
        _widthFitter = widthFitter ?? new EnvelopeWidthFitter();
        _rtFitter = rtFitter ?? new RtProfileFitter();
        _chargeFitter = chargeFitter ?? new ChargeDistributionFitter();
        _abundanceFitter = abundanceFitter ?? new AbundanceFitter();
    }

    public FittedProteoform Fit(ProteoformGroundTruth truth, string? identifier = null)
    {
        var width = _widthFitter.Fit(truth);
        var rt = _rtFitter.Fit(truth);
        var charge = _chargeFitter.Fit(truth);
        var abundance = _abundanceFitter.Fit(truth, width.SigmaMz, rt.Profile, charge.Distribution);

        var model = new ProteoformModel(
            MonoisotopicMass: truth.MonoisotopicMass,
            Abundance: abundance.Abundance,
            RtProfile: rt.Profile,
            ChargeDistribution: charge.Distribution,
            Identifier: identifier);

        return new FittedProteoform(
            model,
            SigmaMz: width.SigmaMz,
            WidthMode: width.Mode,
            WidthPeaksUsed: width.PeaksUsed,
            Residual: abundance.Residual);
    }
}
