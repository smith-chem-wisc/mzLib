using MassSpectrometry;
using Quantification.Interfaces;
using System.Collections.Immutable;

namespace Quantification.Strategies
{
    /// <summary>
    /// This strategy will collapse by summing intensities. All technical replicates and fractions for each Condition_BiologicalReplicate will be summed together.
    /// </summary>
    public class SumCollapse : ICollapseStrategy
    {
        public string Name => "Sum Collapse";

        public QuantMatrix<T> CollapseSamples<T>(QuantMatrix<T> quantMatrix) where T : IEquatable<T>
        {
            // Pre-compute column index map to avoid O(n) IndexOf calls
            var colIndexMap = new Dictionary<ISampleInfo, int>(quantMatrix.ColumnCount);
            for (int i = 0; i < quantMatrix.ColumnKeys.Count; i++)
                colIndexMap[quantMatrix.ColumnKeys[i]] = i;

            var groupings = quantMatrix.ColumnKeys.GroupBy(s => s.Condition + '_' + s.BiologicalReplicate);
            var collapsedMatrix = new QuantMatrix<T>(quantMatrix.RowKeys, groupings.Select(g => g.First()).ToImmutableList(), quantMatrix.ExperimentalDesign);

            foreach (var group in groupings)
            {
                double[] summedValues = new double[quantMatrix.RowKeys.Count];
                foreach (var sample in group)
                {
                    int columnIndex = colIndexMap[sample];
                    for (int rowIndex = 0; rowIndex < quantMatrix.RowKeys.Count; rowIndex++)
                    {
                        summedValues[rowIndex] += quantMatrix.Matrix[rowIndex, columnIndex];
                    }
                }
                collapsedMatrix.SetColumn(group.First(), summedValues);
            }

            return collapsedMatrix;
        }
    }
}
