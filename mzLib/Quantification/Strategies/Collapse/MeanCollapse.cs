using System.Collections.Immutable;
using System.Linq;
using Quantification.Interfaces;

namespace Quantification.Strategies
{
    /// <summary>
    /// Collapses technical replicates and fractions by averaging (rather than summing) intensities.
    /// Samples are grouped by <c>Condition + '_' + BiologicalReplicate</c>, matching the grouping
    /// used by <see cref="SumCollapse"/>. The collapsed value is the arithmetic mean across all
    /// members of each group, using the total group count as the denominator (not the count of
    /// non-zero values). This treats zero as "observed zero" rather than "missing".
    ///
    /// Compared to <see cref="SumCollapse"/>, the mean keeps absolute intensities on the same scale
    /// as individual measurements, making them easier to interpret.
    /// </summary>
    public class MeanCollapse : ICollapseStrategy
    {
        public string Name => "Mean Collapse";

        public QuantMatrix<T> CollapseSamples<T>(QuantMatrix<T> quantMatrix) where T : IEquatable<T>
        {
            var groupings = quantMatrix.ColumnKeys.GroupBy(s => s.Condition + '_' + s.BiologicalReplicate);
            var collapsedMatrix = new QuantMatrix<T>(
                quantMatrix.RowKeys,
                groupings.Select(g => g.First()).ToImmutableList(),
                quantMatrix.ExperimentalDesign);

            foreach (var group in groupings)
            {
                int groupSize = group.Count();
                double[] summedValues = new double[quantMatrix.RowKeys.Count];

                foreach (var sample in group)
                {
                    int columnIndex = quantMatrix.ColumnKeys.IndexOf(sample);
                    for (int rowIndex = 0; rowIndex < quantMatrix.RowKeys.Count; rowIndex++)
                        summedValues[rowIndex] += quantMatrix.Matrix[rowIndex, columnIndex];
                }

                // Divide by group size to produce the mean.
                double[] meanValues = summedValues.Select(v => v / groupSize).ToArray();
                collapsedMatrix.SetColumn(group.First(), meanValues);
            }

            return collapsedMatrix;
        }
    }
}
