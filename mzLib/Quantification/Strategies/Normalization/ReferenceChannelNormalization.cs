using System.Collections.Generic;
using System.Linq;
using MassSpectrometry.ExperimentalDesign;
using Quantification.Interfaces;

namespace Quantification.Strategies
{
    /// <summary>
    /// Normalizes intensities by dividing each channel by the mean of the reference channels
    /// within the same file. Reference channels are identified via
    /// <see cref="IsobaricQuantSampleInfo.IsReferenceChannel"/>. The result converts
    /// absolute intensities into fold-changes relative to the pooled reference.
    ///
    /// After normalization:
    /// - Reference channels are set to 1.0.
    /// - Non-reference channels carry the ratio to the reference mean.
    /// - If both reference channels are zero for a row/file, all channels for that row/file
    ///   are set to zero (missing data; cannot normalize without a reference).
    ///
    /// Columns that are not <see cref="IsobaricQuantSampleInfo"/> are left unchanged.
    /// </summary>
    public class ReferenceChannelNormalization : INormalizationStrategy
    {
        public string Name => "Reference Channel Normalization";

        public QuantMatrix<T> NormalizeIntensities<T>(QuantMatrix<T> quantMatrix) where T : IEquatable<T>
        {
            int rows = quantMatrix.RowKeys.Count;

            // Create the result matrix and pre-fill with the input values.
            var result = new QuantMatrix<T>(quantMatrix.RowKeys, quantMatrix.ColumnKeys, quantMatrix.ExperimentalDesign);
            for (int row = 0; row < rows; row++)
                for (int col = 0; col < quantMatrix.ColumnKeys.Count; col++)
                    result.Matrix[row, col] = quantMatrix.Matrix[row, col];

            // Group isobaric columns by file (FullFilePathWithExtension).
            // Non-isobaric columns are left as-is (already copied above).
            var fileGroups = quantMatrix.ColumnKeys
                .Select((s, i) => (s, i))
                .Where(t => t.s is IsobaricQuantSampleInfo)
                .GroupBy(t => t.s.FullFilePathWithExtension);

            foreach (var fileGroup in fileGroups)
            {
                var fileColumns = fileGroup.ToList();

                var refColIndices = fileColumns
                    .Where(t => ((IsobaricQuantSampleInfo)t.s).IsReferenceChannel)
                    .Select(t => t.i)
                    .ToList();

                var nonRefColIndices = fileColumns
                    .Where(t => !((IsobaricQuantSampleInfo)t.s).IsReferenceChannel)
                    .Select(t => t.i)
                    .ToList();

                // No reference channels in this file — leave values unchanged.
                if (refColIndices.Count == 0)
                    continue;

                for (int row = 0; row < rows; row++)
                {
                    // Reference mean: average of non-zero reference channel values for this row.
                    var refValues = refColIndices
                        .Select(c => quantMatrix.Matrix[row, c])
                        .Where(v => v > 0)
                        .ToList();

                    double refMean = refValues.Count > 0 ? refValues.Average() : 0.0;

                    if (refMean == 0.0)
                    {
                        // No reference signal — zero out all channels for this row/file.
                        foreach (var (_, colIdx) in fileColumns)
                            result.Matrix[row, colIdx] = 0.0;
                        continue;
                    }

                    // Reference channels → 1.0 by definition.
                    foreach (int refCol in refColIndices)
                        result.Matrix[row, refCol] = 1.0;

                    // Non-reference channels → ratio to reference mean.
                    foreach (int nonRefCol in nonRefColIndices)
                        result.Matrix[row, nonRefCol] = quantMatrix.Matrix[row, nonRefCol] / refMean;
                }
            }

            return result;
        }
    }
}
