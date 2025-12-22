using System.Collections.Immutable;
using MassSpectrometry;
using MathNet.Numerics.Statistics;
using Quantification.Interfaces;

namespace Quantification
{
    /// <summary>
    /// Defines how values should be aggregated when collapsing samples.
    /// </summary>
    public enum AggregationType
    {
        Median,
        Average,
        Sum
    }

    /// <summary>
    /// Defines which sample dimension to collapse.
    /// </summary>
    public enum CollapseDimension
    {
        /// <summary>
        /// Collapse fractions within the same Condition + BiologicalReplicate + TechnicalReplicate.
        /// </summary>
        Fraction,

        /// <summary>
        /// Collapse technical replicates within the same Condition + BiologicalReplicate.
        /// </summary>
        TechnicalReplicate,

        /// <summary>
        /// Collapse biological replicates within the same Condition.
        /// (Less common, but useful for generating condition-level summaries)
        /// </summary>
        BiologicalReplicate
    }

    /// <summary>
    /// Collapses samples along a specified dimension (Fraction, TechnicalReplicate, or BiologicalReplicate).
    /// Groups samples by the remaining dimensions, then aggregates intensities using the specified method.
    /// </summary>
    public class SampleCollapseStrategy : ICollapseStrategy
    {
        public string Name => $"Collapse_{Dimension}_{Aggregation}";
        public AggregationType Aggregation { get; }
        public CollapseDimension Dimension { get; }

        public SampleCollapseStrategy(
            CollapseDimension dimension = CollapseDimension.Fraction,
            AggregationType aggregation = AggregationType.Median)
        {
            Dimension = dimension;
            Aggregation = aggregation;
        }

        public QuantMatrix<T> CollapseSamples<T>(QuantMatrix<T> quantMatrix) where T : IEquatable<T>
        {
            // Group columns based on which dimension we're collapsing
            var columnGroups = GroupColumnsByDimension(quantMatrix.ColumnKeys);

            // Create new collapsed column keys (one per group)
            var collapsedColumns = columnGroups
                .Select(g => CreateCollapsedSampleInfo(g.First().Column, Dimension))
                .ToList();

            // Create new matrix with collapsed columns
            var collapsedMatrix = new QuantMatrix<T>(
                quantMatrix.RowKeys,
                collapsedColumns,
                quantMatrix.ExperimentalDesign);

            // For each row, aggregate values across the collapsed dimension
            for (int rowIdx = 0; rowIdx < quantMatrix.RowKeys.Count; rowIdx++)
            {
                var collapsedValues = new double[columnGroups.Count];

                for (int groupIdx = 0; groupIdx < columnGroups.Count; groupIdx++)
                {
                    var group = columnGroups[groupIdx];
                    var values = group.Select(x => quantMatrix.Matrix[rowIdx, x.Index]).ToList();
                    collapsedValues[groupIdx] = Aggregate(values);
                }

                collapsedMatrix.SetRow(quantMatrix.RowKeys[rowIdx], collapsedValues);
            }

            return collapsedMatrix;
        }

        private List<IGrouping<object, (ISampleInfo Column, int Index)>> GroupColumnsByDimension(
            ImmutableList<ISampleInfo> columns)
        {
            var indexed = columns.Select((col, index) => (Column: col, Index: index));

            return Dimension switch
            {
                // Collapse fractions: group by Condition + BioRep + TechRep
                CollapseDimension.Fraction => indexed
                    .GroupBy(x => (object)new
                    {
                        x.Column.Condition,
                        x.Column.BiologicalReplicate,
                        x.Column.TechnicalReplicate
                    })
                    .ToList(),

                // Collapse tech reps: group by Condition + BioRep + Fraction
                CollapseDimension.TechnicalReplicate => indexed
                    .GroupBy(x => (object)new
                    {
                        x.Column.Condition,
                        x.Column.BiologicalReplicate,
                        x.Column.Fraction
                    })
                    .ToList(),

                // Collapse bio reps: group by Condition + TechRep + Fraction
                CollapseDimension.BiologicalReplicate => indexed
                    .GroupBy(x => (object)new
                    {
                        x.Column.Condition,
                        x.Column.TechnicalReplicate,
                        x.Column.Fraction
                    })
                    .ToList(),

                _ => throw new ArgumentOutOfRangeException(nameof(Dimension))
            };
        }

        private double Aggregate(List<double> values)
        {
            if (values.Count == 0) return 0.0;

            return Aggregation switch
            {
                AggregationType.Median => values.Median(),
                AggregationType.Average => values.Average(),
                AggregationType.Sum => values.Sum(),
                _ => throw new ArgumentOutOfRangeException(nameof(Aggregation))
            };
        }

        private static ISampleInfo CreateCollapsedSampleInfo(ISampleInfo source, CollapseDimension dimension)
        {
            // Set the collapsed dimension to 0 to indicate it has been collapsed
            return dimension switch
            {
                CollapseDimension.Fraction => new SpectraFileInfo(
                    fullFilePathWithExtension: $"{source.Condition}_Bio{source.BiologicalReplicate}_Tech{source.TechnicalReplicate}_Collapsed",
                    condition: source.Condition,
                    biorep: source.BiologicalReplicate,
                    techrep: source.TechnicalReplicate,
                    fraction: 0),

                CollapseDimension.TechnicalReplicate => new SpectraFileInfo(
                    fullFilePathWithExtension: $"{source.Condition}_Bio{source.BiologicalReplicate}_Frac{source.Fraction}_Collapsed",
                    condition: source.Condition,
                    biorep: source.BiologicalReplicate,
                    techrep: 0,
                    fraction: source.Fraction),

                CollapseDimension.BiologicalReplicate => new SpectraFileInfo(
                    fullFilePathWithExtension: $"{source.Condition}_Tech{source.TechnicalReplicate}_Frac{source.Fraction}_Collapsed",
                    condition: source.Condition,
                    biorep: 0,
                    techrep: source.TechnicalReplicate,
                    fraction: source.Fraction),

                _ => throw new ArgumentOutOfRangeException(nameof(dimension))
            };
        }
    }
}