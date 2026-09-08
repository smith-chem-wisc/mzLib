using System;
using System.Buffers;
using System.Collections.Generic;
using Quantification.Interfaces;

namespace Quantification.Strategies
{
    /// <summary>
    /// Rolls lower-level rows up to higher-level rows by combining each column with an
    /// <see cref="IAggregationStrategy"/>.
    ///
    /// Which rows are combined is not a strategy decision here -- the caller's map already says so --
    /// so the only behaviour a roll-up carries is how the values are combined, and that is exactly
    /// what IAggregationStrategy names. Splitting it out is the same separation that
    /// <see cref="ICollapseStrategy"/> and IAggregationStrategy make for sample collapsing, and it
    /// keeps <see cref="QuantificationParameters"/> readable as a full description of the processing:
    /// <see cref="SumRollUp"/> still means sum and nothing else.
    ///
    /// Abstract, because there is no default way to combine intensities. A roll-up has to be named by
    /// a concrete subclass -- <see cref="SumRollUp"/>, <see cref="MedianRollUp"/> -- so that reading
    /// <see cref="QuantificationParameters"/> says what was actually done to compute the intensities.
    /// Adding one is a two-line subclass naming its aggregation.
    ///
    /// Zeros are excluded before aggregating -- exactly the zero sentinel, not everything that is not
    /// positive. The matrix uses 0 for "not observed", so a row never measured in a sample should not
    /// pull a median down or a mean toward zero; a negative value, by contrast, is data, and excluding
    /// it would silently change a sum. Since adding zero changes nothing, <see cref="SumRollUp"/>
    /// produces exactly what it always did, for any input.
    ///
    /// <see cref="MedianRollUp"/> keeps its results for non-negative input, which is what intensities
    /// are. It differs from its previous implementation only where a value is negative: that used
    /// <c>v &gt; 0</c> and so treated a negative as missing, which conflated "not measured" with
    /// "measured, and less than zero".
    ///
    /// Note that the exclusion is done HERE, by the roll-up, and not by the aggregation itself --
    /// <see cref="IAggregationStrategy.Aggregate"/> is documented as including zeros and leaving that
    /// choice to its caller. So the same aggregation passed to
    /// <see cref="QuantificationParameters.CollapseAggregationStrategy"/> will see zeros, and behave
    /// differently there. An aggregation whose meaning depends on zeros being excluded is therefore
    /// not recommended as a collapse aggregation strategy.
    /// </summary>
    public abstract class AggregatingRollUp : IRollUpStrategy
    {
        private static readonly List<int> EmptyIndices = new List<int>();

        private readonly IAggregationStrategy _aggregation;

        protected AggregatingRollUp(IAggregationStrategy aggregation)
        {
            _aggregation = aggregation ?? throw new ArgumentNullException(nameof(aggregation));
        }

        public string Name => _aggregation.Name + " Roll-Up";

        public QuantMatrix<THigh> RollUp<TLow, THigh>(QuantMatrix<TLow> matrix, Dictionary<THigh, List<int>> map)
            where TLow : IEquatable<TLow>
            where THigh : IEquatable<THigh>
        {
            var rolledUpMatrix = new QuantMatrix<THigh>(map.Keys, matrix.ColumnKeys, matrix.ExperimentalDesign);

            var pool = ArrayPool<double>.Shared;
            double[] aggregated = pool.Rent(matrix.ColumnCount);

            // Sized by the largest group, not by the row count. The map is caller-supplied, so a group
            // may list an index more than once, and a list longer than the matrix has rows would run
            // off the end of a row-count-sized buffer part way through a roll-up.
            int largestGroup = 0;
            foreach (var kvp in map)
            {
                if (kvp.Value != null && kvp.Value.Count > largestGroup)
                {
                    largestGroup = kvp.Value.Count;
                }
            }

            double[] observed = pool.Rent(Math.Max(1, largestGroup));

            try
            {
                foreach (var kvp in map)
                {
                    List<int> lowIndices = kvp.Value ?? EmptyIndices;

                    for (int col = 0; col < matrix.ColumnCount; col++)
                    {
                        // Gather the observed values for this column. A rented buffer is not zeroed,
                        // so the count -- not the buffer length -- bounds what the aggregator sees.
                        int observedCount = 0;
                        foreach (int lowIndex in lowIndices)
                        {
                            double value = matrix.Matrix[lowIndex, col];

                            // Exactly the zero sentinel, not everything non-positive. A negative value
                            // is data, not a marker for missing, and dropping it would silently change
                            // a sum.
                            if (value != 0)
                            {
                                observed[observedCount++] = value;
                            }
                        }

                        aggregated[col] = _aggregation.Aggregate(
                            new ReadOnlySpan<double>(observed, 0, observedCount));
                    }

                    rolledUpMatrix.SetRow(kvp.Key, aggregated);
                }
            }
            finally
            {
                pool.Return(aggregated, clearArray: true);
                pool.Return(observed, clearArray: true);
            }

            return rolledUpMatrix;
        }
    }
}
