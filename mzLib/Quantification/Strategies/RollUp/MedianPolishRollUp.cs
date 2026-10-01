using System;
using System.Collections.Generic;
using System.Linq;
using FlashLFQ;
using Quantification.Interfaces;

namespace Quantification.Strategies
{
    /// <summary>
    /// Rolls peptides up to proteins with FlashLFQ's median polish protein quantification
    /// (<see cref="FlashLfqResults.MedianPolishProteinIntensities"/>). Each group's rows are log2-transformed
    /// and decomposed into an overall effect, a per-row effect (roughly each peptide's ionization efficiency)
    /// and a per-column effect (the change between samples). The protein intensity in a column is
    /// 2^(overall + column effect) times the number of rows in the group.
    ///
    /// Unlike the <see cref="AggregatingRollUp"/> strategies, columns are not independent here: every
    /// column of a group is fitted together, which is what lets median polish separate a peptide that
    /// ionizes poorly from a protein that is less abundant.
    ///
    /// How this compares with running median polish inside FlashLFQ:
    ///  - FlashLFQ collapses fractions (most intense fraction) and technical replicates (mean) before
    ///    polishing. Here each column of the incoming matrix is one sample; collapsing is left to
    ///    <see cref="QuantificationParameters.CollapseStrategy"/>, which runs before this roll-up.
    ///  - A value that is not positive is not observed, since it has no logarithm. The matrix uses 0 for
    ///    that; a negative value is treated the same way.
    ///  - FlashLFQ reports NaN for a sample where the protein was observed but is unquantifiable. The
    ///    matrix has one marker for "no value", 0, and downstream strategies (e.g.
    ///    <see cref="GlobalMedianNormalization"/>) treat anything else as a measurement, so such a sample
    ///    rolls up to 0 here.
    ///  - A row with no observed value in any column is left out, as FlashLFQ leaves out a peptide it never
    ///    quantified. Kept, it would add nothing to the fit but still count toward the number of rows the
    ///    result scales by.
    ///  - A row index listed more than once in a group counts once. It is still one peptide, and a
    ///    duplicate would both weight the fit toward it and inflate the row count the result scales by.
    /// </summary>
    public class MedianPolishRollUp : IRollUpStrategy
    {
        public string Name => "Median Polish Roll-Up";

        public QuantMatrix<THigh> RollUp<TLow, THigh>(QuantMatrix<TLow> matrix, Dictionary<THigh, List<int>> map)
            where TLow : IEquatable<TLow>
            where THigh : IEquatable<THigh>
        {
            var rolledUpMatrix = new QuantMatrix<THigh>(map.Keys, matrix.ColumnKeys, matrix.ExperimentalDesign);

            foreach (var kvp in map)
            {
                if (kvp.Value == null || kvp.Value.Count == 0)
                {
                    // Nothing to polish; the row stays at the 0 it was created with.
                    continue;
                }

                double[][] rows = kvp.Value
                    .Distinct()
                    .Select(lowIndex => matrix.GetRow(lowIndex))
                    .Where(row => row.Any(value => value > 0))
                    .ToArray();

                if (rows.Length == 0)
                {
                    continue;
                }

                double[] proteinIntensities = FlashLfqResults.MedianPolishProteinIntensities(rows);

                for (int col = 0; col < proteinIntensities.Length; col++)
                {
                    if (double.IsNaN(proteinIntensities[col]))
                    {
                        proteinIntensities[col] = 0;
                    }
                }

                rolledUpMatrix.SetRow(kvp.Key, proteinIntensities);
            }

            return rolledUpMatrix;
        }
    }
}
