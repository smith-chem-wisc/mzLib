#nullable enable
using System;
using System.Collections.Generic;
using TopDownSimulator.Extraction;
using TopDownSimulator.Model;

namespace TopDownSimulator.Fitting;

public sealed record FittedProteoform(
    ProteoformModel Model,
    double SigmaMz,
    WidthFitMode WidthMode,
    int WidthPeaksUsed,
    double Residual,
    IReadOnlyList<PeakWidthMeasurement>? WidthMeasurements = null,
    int WidthWindowsTooCoarselySampled = 0)
{
    /// <summary>
    /// The (m/z, σ) observations this record contributed. Pooled across every record they are the
    /// input to <see cref="PeakWidthModelFitter"/>; <see cref="SigmaMz"/> alone cannot support a
    /// width law.
    /// </summary>
    public IReadOnlyList<PeakWidthMeasurement> WidthMeasurements { get; init; } =
        WidthMeasurements ?? Array.Empty<PeakWidthMeasurement>();
}

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

    /// <summary>
    /// Fits one proteoform. When <paramref name="widthModel"/> is supplied the abundance is fitted
    /// under it, so the fitted number is consistent with what a simulation under that same model
    /// will render; otherwise the abundance is fitted under this record's own measured σ.
    /// </summary>
    public FittedProteoform Fit(
        ProteoformGroundTruth truth,
        string? identifier = null,
        IPeakWidthModel? widthModel = null)
    {
        var width = _widthFitter.Fit(truth);
        var rt = _rtFitter.Fit(truth);
        var charge = _chargeFitter.Fit(truth);
        var abundance = _abundanceFitter.Fit(
            truth,
            widthModel ?? new ConstantPeakWidth(width.SigmaMz),
            rt.Profile,
            charge.Distribution);

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
            Residual: abundance.Residual,
            WidthMeasurements: width.Measurements,
            WidthWindowsTooCoarselySampled: width.WindowsTooCoarselySampled);
    }
}
