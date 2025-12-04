using System;

namespace Omics
{
    /// <summary>
    /// Hypothesis for a spectral match: a single candidate assignment for a spectrum
    /// (for example a specific biopolymer sequence with localized set modifications).
    /// Implementations represent a scored candidate and must provide value-equality
    /// semantics via <see cref="IEquatable{T}"/>.
    /// </summary>
    public interface ISpectralMatchHypothesis : IEquatable<ISpectralMatchHypothesis>
    {
        /// <summary>
        /// Numeric score for the hypothesis. Higher values typically indicate better confidence.
        /// Callers should document whether higher-is-better semantics apply.
        /// </summary>
        double Score { get; }

        /// <summary>
        /// The specific biopolymer (peptide/proteoform/oligo) with localized set modifications
        /// that this hypothesis assigns to the spectrum.
        /// </summary>
        IBioPolymerWithSetMods SpecificBioPolymer { get; }

        /// <summary>
        /// True when this hypothesis derives from a decoy sequence (used for FDR estimation).
        /// </summary>
        bool IsDecoy { get; }

        /// <summary>
        /// Full sequence string for this hypothesis (including modification annotations).
        /// </summary>
        string FullSequence { get; }

        /// <summary>
        /// Optional per-notch q-value or other localized confidence metric. Null when not available.
        /// </summary>
        double? QValueNotch { get; }

        /// <summary>
        /// Determines whether this hypothesis is equal to <paramref name="other"/>.
        /// Implementations should compare identity-relevant fields (sequence, mods, location, etc.)
        /// so that logically equivalent hypotheses are treated as equal.
        /// </summary>
        /// <param name="other">Other hypothesis to compare.</param>
        /// <returns>true if equivalent; otherwise false.</returns>
        bool Equals(ISpectralMatchHypothesis? other);

    }
}