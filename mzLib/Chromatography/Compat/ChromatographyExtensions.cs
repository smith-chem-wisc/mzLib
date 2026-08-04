using Omics;

// This file intentionally declares a namespace that does not match its assembly. It preserves the
// pre-Chromatography type locations for downstream consumers (notably MetaMorpheus) that still
// reference Proteomics.RetentionTimePrediction.SSRCalc3. It lives in the Chromatography project
// rather than Proteomics so that Proteomics no longer has to reference Chromatography — that edge
// was what dragged TorchSharp and libtorch into the dependency graph of every project built on
// Proteomics, which is to say almost all of them.
namespace Proteomics.RetentionTimePrediction;

/// <summary>
/// Backwards-compatibility shims for the type layout that predated the Chromatography namespace.
/// </summary>
/// <remarks>
/// Prefer <see cref="Chromatography.RetentionTimePrediction.SSRCalc.SSRCalc3"/> and its
/// <c>ScoreSequence(string)</c> overload directly. These shims are scheduled for removal in the
/// next major version.
/// </remarks>
[Obsolete("Use Chromatography.RetentionTimePrediction.SSRCalc.SSRCalc3.ScoreSequence(string) directly. This shim will be removed in the next major version.")]
public static class ChromatographyExtensions
{
    /// <summary>
    /// Scores a biopolymer's base sequence with SSRCalc3.
    /// </summary>
    /// <remarks>
    /// The parameter was widened from <c>PeptideWithSetModifications</c> to
    /// <see cref="IBioPolymerWithSetMods"/> so this shim depends only on Omics. Existing call sites
    /// that pass a <c>PeptideWithSetModifications</c> continue to compile unchanged.
    /// </remarks>
    public static double ScoreSequence(this Chromatography.RetentionTimePrediction.SSRCalc.SSRCalc3 predictor,
        IBioPolymerWithSetMods peptide) => predictor.ScoreSequence(peptide.BaseSequence);
}

/// <inheritdoc cref="ChromatographyExtensions"/>
[Obsolete("Use Chromatography.RetentionTimePrediction.SSRCalc.SSRCalc3 directly. This shim will be removed in the next major version.")]
public class SSRCalc3 : Chromatography.RetentionTimePrediction.SSRCalc.SSRCalc3
{
    public SSRCalc3(string name, Column column) : base(name, column) { }
}
