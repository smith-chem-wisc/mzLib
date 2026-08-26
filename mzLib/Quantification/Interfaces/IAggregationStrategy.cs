using System;

namespace Quantification.Interfaces
{
    /// <summary>
    /// Combines several intensity values into one.
    ///
    /// This is deliberately separate from <see cref="ICollapseStrategy"/>, which decides *which*
    /// samples are combined. Keeping the two apart means each strategy type still describes exactly
    /// one behaviour, so reading <see cref="QuantificationParameters"/> is enough to know how the data
    /// was processed — the same reason <see cref="INormalizationStrategy"/> and
    /// <see cref="IRollUpStrategy"/> have one implementation per behaviour rather than one
    /// parameterized implementation.
    /// </summary>
    public interface IAggregationStrategy
    {
        string Name { get; }

        /// <summary>
        /// Reduces <paramref name="values"/> to a single value. An empty span aggregates to 0.
        /// </summary>
        /// <remarks>
        /// Implementations must not mutate <paramref name="values"/>; callers reuse the buffer.
        /// Zeros are included: the caller decides whether a zero means "not observed" and should be
        /// excluded before calling.
        /// </remarks>
        double Aggregate(ReadOnlySpan<double> values);
    }
}
