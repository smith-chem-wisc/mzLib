using System;
using System.Collections.Generic;
using System.Linq;
using Quantification.Interfaces;

namespace Quantification.Strategies
{
    /// <summary>
    /// Normalizes intensities by equalizing the median log2 intensity across all columns.
    /// Each column is shifted in log2 space so that its median equals the global median of all
    /// column medians. Zero values (missing data) are preserved and excluded from median computation.
    /// </summary>
    public class GlobalMedianNormalization : INormalizationStrategy
    {
        public string Name => "Global Median Normalization";

        public QuantMatrix<T> NormalizeIntensities<T>(QuantMatrix<T> quantMatrix) where T : IEquatable<T>
        {
            int rows = quantMatrix.RowKeys.Count;
            int cols = quantMatrix.ColumnKeys.Count;

            // Step 1: Compute per-column log2 medians, ignoring zeros and negatives.
            // Skip columns with fewer than 2 positive values (no meaningful shift possible).
            var columnLog2Medians = new double[cols];
            var log2Vals = new List<double>(rows);
            for (int col = 0; col < cols; col++)
            {
                for (int row = 0; row < rows; row++)
                {
                    double val = quantMatrix.Matrix[row, col];
                    if (val > 0)
                        log2Vals.Add(Math.Log2(val));
                }
                columnLog2Medians[col] = log2Vals.Count > 1 ? Median(log2Vals) : double.NaN;
                log2Vals.Clear();
            }

            // Step 2: Compute the global median across all valid column medians.
            var validMedians = columnLog2Medians.Where(m => !double.IsNaN(m)).ToList();
            if (validMedians.Count == 0)
                return quantMatrix; // Nothing to normalize; return as-is.

            double globalMedian = Median(validMedians);

            // Step 3: Create new matrix and apply per-column shifts.
            var result = new QuantMatrix<T>(quantMatrix.RowKeys, quantMatrix.ColumnKeys, quantMatrix.ExperimentalDesign);

            for (int col = 0; col < cols; col++)
            {
                if (double.IsNaN(columnLog2Medians[col]))
                {
                    // Column has 0 or 1 positive values — copy as-is, no normalization.
                    for (int row = 0; row < rows; row++)
                        result.Matrix[row, col] = quantMatrix.Matrix[row, col];
                    continue;
                }

                // Multiplicative factor corresponding to the log2 shift.
                double factor = Math.Pow(2, globalMedian - columnLog2Medians[col]);

                for (int row = 0; row < rows; row++)
                {
                    double val = quantMatrix.Matrix[row, col];
                    // Zeros remain zero; positive values are scaled.
                    result.Matrix[row, col] = val > 0 ? val * factor : 0.0;
                }
            }

            return result;
        }

        private static double Median(List<double> values)
        {
            var sorted = values.OrderBy(x => x).ToList();
            int mid = sorted.Count / 2;
            return sorted.Count % 2 == 0
                ? (sorted[mid - 1] + sorted[mid]) / 2.0
                : sorted[mid];
        }
    }
}
