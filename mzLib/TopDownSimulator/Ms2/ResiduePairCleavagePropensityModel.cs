#nullable enable
using System;

namespace TopDownSimulator.Ms2;

/// <summary>
/// <see cref="IBondCleavageModel"/> backed by the published amino-acid-pair propensity table.
/// </summary>
/// <remarks>
/// <para>
/// The table values are already fractions in [0, 1], so the default configuration is the identity:
/// the probability handed to the cascade is exactly the published propensity. Rescaling is
/// available but deliberately explicit — <see cref="PropensityScale"/> defaults to 1 and
/// <see cref="MaximumCleavageProbability"/> defaults to 1, so neither does anything unless you ask.
/// </para>
/// <para>
/// <see cref="MissingPairPropensity"/> covers the two ways the table can come up empty: an
/// unrecognized residue letter (selenocysteine, ambiguity codes, 'X') and a pair whose propensity
/// was never measured. Every pair in the shipped table is defined, so only the former can occur
/// unless a custom table is supplied.
/// </para>
/// </remarks>
public sealed class ResiduePairCleavagePropensityModel : IBondCleavageModel
{
    public ResiduePairCleavagePropensityModel(
        ResiduePairCleavagePropensityTable? table = null,
        double propensityScale = 1.0,
        double missingPairPropensity = 0.0,
        double maximumCleavageProbability = 1.0)
    {
        if (propensityScale < 0)
            throw new ArgumentOutOfRangeException(nameof(propensityScale), "Scale must be non-negative.");
        if (missingPairPropensity < 0 || missingPairPropensity > 1)
            throw new ArgumentOutOfRangeException(nameof(missingPairPropensity), "Must be in [0, 1].");
        if (maximumCleavageProbability < 0 || maximumCleavageProbability > 1)
            throw new ArgumentOutOfRangeException(nameof(maximumCleavageProbability), "Must be in [0, 1].");

        Table = table ?? ResiduePairCleavagePropensityTable.Default;
        PropensityScale = propensityScale;
        MissingPairPropensity = missingPairPropensity;
        MaximumCleavageProbability = maximumCleavageProbability;
    }

    public ResiduePairCleavagePropensityTable Table { get; }

    /// <summary>
    /// Multiplier applied to every looked-up propensity before clamping. 1.0 (the default) means
    /// the published fractions are used as-is; raising it makes the precursor fragment further up
    /// the backbone, lowering it leaves more surviving precursor.
    /// </summary>
    public double PropensityScale { get; }

    /// <summary>Probability used when the table has no measured value for a bond. Default 0.</summary>
    public double MissingPairPropensity { get; }

    /// <summary>Upper clamp applied after scaling. Default 1.0, i.e. no clamping in practice.</summary>
    public double MaximumCleavageProbability { get; }

    public double GetCleavageProbability(string baseSequence, int zeroBasedBondIndex)
    {
        if (string.IsNullOrEmpty(baseSequence))
            throw new ArgumentException("A base sequence is required.", nameof(baseSequence));
        if (zeroBasedBondIndex < 0 || zeroBasedBondIndex > baseSequence.Length - 2)
            throw new ArgumentOutOfRangeException(nameof(zeroBasedBondIndex),
                $"Bond index must be in [0, {baseSequence.Length - 2}] for a sequence of length {baseSequence.Length}.");

        char nSide = baseSequence[zeroBasedBondIndex];
        char cSide = baseSequence[zeroBasedBondIndex + 1];

        double propensity = Table.TryGetPropensity(nSide, cSide, out double looked)
            ? looked
            : MissingPairPropensity;

        return Math.Clamp(propensity * PropensityScale, 0.0, MaximumCleavageProbability);
    }
}
