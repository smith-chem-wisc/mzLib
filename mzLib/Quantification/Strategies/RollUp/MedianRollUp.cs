using Quantification.Interfaces;
using System.Buffers;
using System.Collections.Generic;
using System.Linq;

namespace Quantification.Strategies
{
    /// <summary>
    /// Rolls up lower-level entities (e.g., PSMs) to higher-level entities (e.g., peptides or proteins)
    /// by taking the median of non-zero values in each column. Zero values are treated as missing
    /// data and excluded from the median calculation. If all values for a column are zero, the
    /// rolled-up value is zero.
    ///
    /// Compared to <see cref="SumRollUp"/>, median roll-up is more robust to outliers and is
    /// not biased toward entities with more identified lower-level features.
    /// </summary>
    public class MedianRollUp : IRollUpStrategy
    {
        public string Name => "Median Roll-Up";
        public ArrayPool<double> ArrayPool => ArrayPool<double>.Shared;

        public QuantMatrix<THigh> RollUp<TLow, THigh>(QuantMatrix<TLow> matrix, Dictionary<THigh, List<int>> map)
            where TLow : IEquatable<TLow>
            where THigh : IEquatable<THigh>
        {
            var rolledUpMatrix = new QuantMatrix<THigh>(map.Keys, matrix.ColumnKeys, matrix.ExperimentalDesign);

            foreach (var kvp in map)
            {
                THigh highKey = kvp.Key;
                List<int> lowIndices = kvp.Value;

                double[] medianValues = ArrayPool.Rent(matrix.ColumnCount);

                for (int sampleIndex = 0; sampleIndex < matrix.ColumnCount; sampleIndex++)
                {
                    var nonZeroValues = new List<double>();
                    foreach (int i in lowIndices)
                    {
                        double v = matrix.Matrix[i, sampleIndex];
                        if (v > 0) nonZeroValues.Add(v);
                    }

                    medianValues[sampleIndex] = nonZeroValues.Count > 0
                        ? Median(nonZeroValues)
                        : 0.0;
                }

                rolledUpMatrix.SetRow(highKey, medianValues);
                ArrayPool.Return(medianValues);
            }

            return rolledUpMatrix;
        }

        private static double Median(List<double> values)
        {
            var arr = values.ToArray();
            Array.Sort(arr);
            int mid = arr.Length / 2;
            return arr.Length % 2 == 0
                ? (arr[mid - 1] + arr[mid]) / 2.0
                : arr[mid];
        }
    }
}
