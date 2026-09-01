using System;

namespace Quantification.Interfaces
{
    /// <summary>
    /// Collapse strategies decide *which* samples are combined, reducing the number of columns in a
    /// QuantMatrix. An example is combining the fractions of a sample into a single column.
    /// This is a column-wise operation.
    /// Incoming matrix has dimensions (p x s) where s = number of samples and p = number of peptide/proteins
    /// Outgoing matrix has dimensions (p x S) where S = number of collapsed samples, S &lt;= s
    ///
    /// *How* the combined values are reduced to one is the separate concern of
    /// <see cref="IAggregationStrategy"/>, so that each strategy type describes exactly one behaviour
    /// and <see cref="QuantificationParameters"/> fully describes how the data was processed.
    /// </summary>
    public interface ICollapseStrategy
    {
        string Name { get; }

        /// <summary>
        /// Whether <see cref="CollapseSamples{T}"/> reads the aggregation strategy it is handed.
        /// </summary>
        /// <remarks>
        /// A strategy that groups nothing never reduces two values to one, so it has no use for an
        /// aggregator -- and requiring the caller to supply one anyway would mean naming a behaviour
        /// that does not happen. <see cref="QuantificationEngine"/> reads this to decide whether a null
        /// <see cref="QuantificationParameters.CollapseAggregationStrategy"/> is a configuration error
        /// or simply unused.
        ///
        /// Defaults to true, which is the safe answer: an implementation that does aggregate gets the
        /// check without opting in, and only a genuine no-op has to say so.
        /// </remarks>
        bool RequiresAggregation => true;

        /// <summary>
        /// Combines the columns this strategy groups together, reducing each group to a single column
        /// with <paramref name="aggregation"/>.
        /// </summary>
        /// <param name="quantMatrix">The matrix whose columns are to be collapsed.</param>
        /// <param name="aggregation">How the values within a group are reduced to one.</param>
        QuantMatrix<T> CollapseSamples<T>(QuantMatrix<T> quantMatrix, IAggregationStrategy aggregation)
            where T : IEquatable<T>;
    }
}
