using System;
using System.Buffers;
using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;
using Quantification.Interfaces;

namespace Quantification.Strategies
{
    /// <summary>
    /// How values are combined when several samples collapse into one.
    /// </summary>
    public enum AggregationType
    {
        Median,
        Average,
        Sum,
        Max
    }

    /// <summary>
    /// Which sample dimension to collapse.
    /// </summary>
    public enum CollapseDimension
    {
        /// <summary>
        /// Collapse fractions within the same Condition + BiologicalReplicate + TechnicalReplicate.
        /// </summary>
        Fraction,

        /// <summary>
        /// Collapse technical replicates within the same Condition + BiologicalReplicate + Fraction.
        /// </summary>
        TechnicalReplicate,

        /// <summary>
        /// Collapse biological replicates within the same Condition + TechnicalReplicate + Fraction.
        /// Less common, but useful for condition-level summaries.
        /// </summary>
        BiologicalReplicate
    }

    /// <summary>
    /// Collapses samples along one named dimension, grouping by the dimensions that remain and
    /// aggregating each group with the chosen <see cref="AggregationType"/>.
    ///
    /// This is the general form of <see cref="SumCollapse"/> and <see cref="MeanCollapse"/>, which
    /// group by Condition + BiologicalReplicate only and therefore collapse every other dimension at
    /// once. Applying this strategy twice expresses the usual bottom-up recipe — sum fractions, then
    /// take the median across technical replicates:
    ///
    /// <code>
    /// var summed = new SampleCollapseStrategy(CollapseDimension.Fraction, AggregationType.Sum)
    ///     .CollapseSamples(matrix);
    /// var final = new SampleCollapseStrategy(CollapseDimension.TechnicalReplicate, AggregationType.Median)
    ///     .CollapseSamples(summed);
    /// </code>
    ///
    /// Note that samples agreeing on every remaining dimension are combined even if they are distinct
    /// channels of the same file — two reference channels sharing a condition and biological replicate
    /// collapse together, which is the same behaviour <see cref="SumCollapse"/> has.
    /// </summary>
    public class SampleCollapseStrategy : ICollapseStrategy
    {
        public string Name => $"Collapse_{Dimension}_{Aggregation}";

        public AggregationType Aggregation { get; }
        public CollapseDimension Dimension { get; }

        private static ArrayPool<double> ArrayPool => ArrayPool<double>.Shared;

        public SampleCollapseStrategy(
            CollapseDimension dimension = CollapseDimension.Fraction,
            AggregationType aggregation = AggregationType.Median)
        {
            Dimension = dimension;
            Aggregation = aggregation;
        }

        public QuantMatrix<T> CollapseSamples<T>(QuantMatrix<T> quantMatrix) where T : IEquatable<T>
        {
            var columnGroups = GroupColumnsByDimension(quantMatrix.ColumnKeys);

            var collapsedColumns = columnGroups
                .Select(g => CreateCollapsedSampleInfo(g[0].Column, Dimension))
                .ToList();

            var collapsedMatrix = new QuantMatrix<T>(
                quantMatrix.RowKeys,
                collapsedColumns,
                quantMatrix.ExperimentalDesign);

            // One scratch buffer for the collapsed row, one for the values feeding a single group.
            int widestGroup = columnGroups.Count == 0 ? 0 : columnGroups.Max(g => g.Count);
            double[] collapsedValues = ArrayPool.Rent(Math.Max(columnGroups.Count, 1));
            double[] groupValues = ArrayPool.Rent(Math.Max(widestGroup, 1));

            try
            {
                for (int rowIdx = 0; rowIdx < quantMatrix.RowCount; rowIdx++)
                {
                    for (int groupIdx = 0; groupIdx < columnGroups.Count; groupIdx++)
                    {
                        var group = columnGroups[groupIdx];
                        for (int i = 0; i < group.Count; i++)
                            groupValues[i] = quantMatrix.Matrix[rowIdx, group[i].Index];

                        collapsedValues[groupIdx] = Aggregate(groupValues, group.Count);
                    }

                    collapsedMatrix.SetRow(quantMatrix.RowKeys[rowIdx], collapsedValues);
                }
            }
            finally
            {
                ArrayPool.Return(collapsedValues, clearArray: true);
                ArrayPool.Return(groupValues, clearArray: true);
            }

            return collapsedMatrix;
        }

        /// <summary>
        /// Groups the columns by the dimensions that are NOT being collapsed, preserving the original
        /// column order both between groups and within each group.
        /// </summary>
        private List<List<(ISampleInfo Column, int Index)>> GroupColumnsByDimension(
            IReadOnlyList<ISampleInfo> columns)
        {
            var order = new List<(string, int, int)>();
            var groups = new Dictionary<(string, int, int), List<(ISampleInfo, int)>>();

            for (int index = 0; index < columns.Count; index++)
            {
                var column = columns[index];
                var key = Dimension switch
                {
                    CollapseDimension.Fraction =>
                        (column.Condition ?? string.Empty, column.BiologicalReplicate, column.TechnicalReplicate),
                    CollapseDimension.TechnicalReplicate =>
                        (column.Condition ?? string.Empty, column.BiologicalReplicate, column.Fraction),
                    CollapseDimension.BiologicalReplicate =>
                        (column.Condition ?? string.Empty, column.TechnicalReplicate, column.Fraction),
                    _ => throw new ArgumentOutOfRangeException(nameof(Dimension))
                };

                if (!groups.TryGetValue(key, out var members))
                {
                    members = new List<(ISampleInfo, int)>();
                    groups[key] = members;
                    order.Add(key);
                }
                members.Add((column, index));
            }

            return order.Select(k => groups[k]).ToList();
        }

        /// <summary>
        /// Aggregates the first <paramref name="count"/> entries of <paramref name="values"/>.
        /// The buffer may be longer than <paramref name="count"/> because it comes from an array pool.
        /// </summary>
        private double Aggregate(double[] values, int count)
        {
            if (count == 0) return 0.0;

            switch (Aggregation)
            {
                case AggregationType.Sum:
                {
                    double sum = 0;
                    for (int i = 0; i < count; i++) sum += values[i];
                    return sum;
                }
                case AggregationType.Average:
                {
                    double sum = 0;
                    for (int i = 0; i < count; i++) sum += values[i];
                    return sum / count;
                }
                case AggregationType.Max:
                {
                    double max = values[0];
                    for (int i = 1; i < count; i++) if (values[i] > max) max = values[i];
                    return max;
                }
                case AggregationType.Median:
                    return Median(values, count);
                default:
                    throw new ArgumentOutOfRangeException(nameof(Aggregation));
            }
        }

        private static double Median(double[] values, int count)
        {
            var sorted = new double[count];
            Array.Copy(values, sorted, count);
            Array.Sort(sorted);
            int mid = count / 2;
            return count % 2 == 0 ? (sorted[mid - 1] + sorted[mid]) / 2.0 : sorted[mid];
        }

        /// <summary>
        /// Builds the <see cref="ISampleInfo"/> that represents a collapsed group, zeroing the
        /// dimension that was collapsed and keeping the rest.
        /// </summary>
        /// <remarks>
        /// An <see cref="IsobaricQuantSampleInfo"/> collapses to another
        /// <see cref="IsobaricQuantSampleInfo"/> so that the channel label, plex, reporter m/z and
        /// reference-channel flag survive. Collapsing an isobaric column into a plain
        /// <see cref="SpectraFileInfo"/> would silently disable
        /// <see cref="ReferenceChannelNormalization"/> for everything downstream.
        /// </remarks>
        internal static ISampleInfo CreateCollapsedSampleInfo(ISampleInfo source, CollapseDimension dimension)
        {
            (int biorep, int techrep, int fraction) = dimension switch
            {
                CollapseDimension.Fraction => (source.BiologicalReplicate, source.TechnicalReplicate, 0),
                CollapseDimension.TechnicalReplicate => (source.BiologicalReplicate, 0, source.Fraction),
                CollapseDimension.BiologicalReplicate => (0, source.TechnicalReplicate, source.Fraction),
                _ => throw new ArgumentOutOfRangeException(nameof(dimension))
            };

            string label = $"{source.Condition}_Bio{biorep}_Frac{fraction}_Tech{techrep}_Collapsed";

            if (source is IsobaricQuantSampleInfo isobaric)
            {
                return new IsobaricQuantSampleInfo(
                    fullFilePathWithExtension: label,
                    condition: isobaric.Condition,
                    biologicalReplicate: biorep,
                    technicalReplicate: techrep,
                    fraction: fraction,
                    plexId: isobaric.PlexId,
                    channelLabel: isobaric.ChannelLabel,
                    reporterIonMz: isobaric.ReporterIonMz,
                    isReferenceChannel: isobaric.IsReferenceChannel);
            }

            return new SpectraFileInfo(
                fullFilePathWithExtension: label,
                condition: source.Condition,
                biorep: biorep,
                techrep: techrep,
                fraction: fraction);
        }
    }
}
