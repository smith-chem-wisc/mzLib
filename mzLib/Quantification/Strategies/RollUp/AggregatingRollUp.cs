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
    /// Zeros are excluded before aggregating. The matrix uses 0 for "not observed", so a row that was
    /// never measured in a sample should not pull a median down or a mean toward zero. Summing is
    /// unaffected by the exclusion, which is why <see cref="SumRollUp"/> and <see cref="MedianRollUp"/>
    /// both keep their previous results.
    /// </summary>
    public class AggregatingRollUp : IRollUpStrategy
    {
        private readonly IAggregationStrategy _aggregation;

        public AggregatingRollUp(IAggregationStrategy aggregation)
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
            double[] observed = pool.Rent(Math.Max(1, matrix.RowCount));

            try
            {
                foreach (var kvp in map)
                {
                    List<int> lowIndices = kvp.Value;

                    for (int col = 0; col < matrix.ColumnCount; col++)
                    {
                        // Gather the observed values for this column. A rented buffer is not zeroed,
                        // so the count -- not the buffer length -- bounds what the aggregator sees.
                        int observedCount = 0;
                        foreach (int lowIndex in lowIndices)
                        {
                            double value = matrix.Matrix[lowIndex, col];
                            if (value > 0)
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
