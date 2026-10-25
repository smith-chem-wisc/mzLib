using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Statistics;
using Quantification.Interfaces;

namespace Quantification.Strategies
{
    /// <summary>
    /// Rolls peptides up to proteins with median polish, the protein quantification FlashLFQ uses
    /// (FlashLfqResults.CalculateProteinResultsMedianPolish calls <see cref="QuantifyGroup"/>). Each group's
    /// rows are log2-transformed and decomposed into an overall effect, a per-row effect (roughly each
    /// peptide's ionization efficiency) and a per-column effect (the change between samples). The protein
    /// intensity in a column is 2^(overall + column effect) times the number of rows in the group.
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

                double[] proteinIntensities = QuantifyGroup(rows);

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

        /// <summary>
        /// Quantifies one protein from the intensities of its peptides using the median polish algorithm.
        /// This is the per-protein step of FlashLFQ's median polish protein quantification, on a plain
        /// peptide-by-sample table, so FlashLFQ and this roll-up share one implementation.
        /// </summary>
        /// <param name="peptideIntensities">One row per peptide, one column per sample, all rows the same length.
        /// Values are un-logged intensities; a value that is not positive means the peptide was not observed.</param>
        /// <returns>One protein intensity per sample: 0 where no peptide was observed, NaN where peptides were
        /// observed but the protein is unquantifiable in that sample, and the median polish estimate otherwise.</returns>
        public static double[] QuantifyGroup(double[][] peptideIntensities)
        {
            int numPeptides = peptideIntensities.Length;
            int numSamples = numPeptides == 0 ? 0 : peptideIntensities[0].Length;
            double[] proteinIntensities = new double[numSamples];

            if (numPeptides == 0)
            {
                return proteinIntensities;
            }

            // set up peptide intensity table
            // top row is the column effects, left column is the row effects
            // the other cells are log2-transformed peptide intensity measurements
            // if a value is missing, it will be filled with NaN
            double[][] peptideIntensityMatrix = new double[numPeptides + 1][];
            peptideIntensityMatrix[0] = new double[numSamples + 1];
            bool[] sampleIsObserved = new bool[numSamples];

            for (int p = 0; p < numPeptides; p++)
            {
                peptideIntensityMatrix[p + 1] = new double[numSamples + 1];

                for (int s = 0; s < numSamples; s++)
                {
                    double sampleIntensity = peptideIntensities[p][s];

                    if (sampleIntensity > 0)
                    {
                        peptideIntensityMatrix[p + 1][s + 1] = Math.Log(sampleIntensity, 2);
                        sampleIsObserved[s] = true;
                    }
                    else
                    {
                        peptideIntensityMatrix[p + 1][s + 1] = double.NaN;
                    }
                }
            }

            // if there are any peptides that have only one measurement, mark them as NaN
            // unless we have ONLY peptides with one measurement
            var peptidesWithMoreThanOneMmt = peptideIntensityMatrix.Skip(1).Count(row => row.Skip(1).Count(cell => !double.IsNaN(cell)) > 1);
            if (peptidesWithMoreThanOneMmt > 0)
            {
                for (int i = 1; i < peptideIntensityMatrix.Length; i++)
                {
                    int validValueCount = peptideIntensityMatrix[i].Count(p => !double.IsNaN(p) && p != 0);

                    if (validValueCount < 2 && numSamples >= 2)
                    {
                        for (int j = 1; j < peptideIntensityMatrix[0].Length; j++)
                        {
                            peptideIntensityMatrix[i][j] = double.NaN;
                        }
                    }
                }
            }

            // do median polish protein quantification
            // row effects in a protein can be considered ~ relative ionization efficiency
            // column effects are differences between conditions
            MedianPolish(peptideIntensityMatrix);

            double overallEffect = peptideIntensityMatrix[0][0];
            double[] columnEffects = peptideIntensityMatrix[0].Skip(1).ToArray();
            double referenceProteinIntensity = Math.Pow(2, overallEffect) * numPeptides;

            // check for unquantifiable proteins; these are proteins w/ quantified peptides, but
            // the protein is still not quantifiable because there are not peptides to compare across runs.
            // the column effect can be 0 in some cases. sometimes it's a valid value and sometimes it's not.
            int possibleUnquantifiableSampleCount = 0;
            for (int s = 0; s < numSamples; s++)
            {
                if (sampleIsObserved[s] && columnEffects[s] == 0)
                {
                    possibleUnquantifiableSampleCount++;
                }
            }

            for (int s = 0; s < numSamples; s++)
            {
                if (!sampleIsObserved[s])
                {
                    continue;
                }

                if (possibleUnquantifiableSampleCount > 1 && columnEffects[s] == 0)
                {
                    proteinIntensities[s] = double.NaN;
                }
                else
                {
                    // this step un-logs the protein "intensity". in reality this value is more like a fold-change
                    // than an intensity, but unlike a fold-change it's not relative to a particular sample.
                    // by multiplying this value by the reference protein intensity calculated earlier, then we get
                    // a protein intensity value
                    proteinIntensities[s] = Math.Pow(2, columnEffects[s]) * referenceProteinIntensity;
                }
            }

            return proteinIntensities;
        }

        /// <summary>
        /// Fits the median polish model to <paramref name="table"/> in place. The top row accumulates the
        /// column effects and the left column the row effects, with the overall effect at [0][0]; the
        /// other cells are log-scale measurements (NaN for missing) and are left holding the residuals.
        /// </summary>
        public static void MedianPolish(double[][] table, int maxIterations = 10, double improvementCutoff = 0.0001)
        {
            // technically, this is weighted mean polish and not median polish.
            // but it should give similar results while being more robust to issues
            // arising from missing values.
            // the weights are inverse square difference to median.

            // subtract overall effect
            List<double> allValues = table.SelectMany(p => p.Where(p => !double.IsNaN(p) && p != 0)).ToList();

            if (allValues.Any())
            {
                double overallEffect = allValues.Median();
                table[0][0] += overallEffect;

                for (int r = 1; r < table.Length; r++)
                {
                    for (int c = 1; c < table[0].Length; c++)
                    {
                        table[r][c] -= overallEffect;
                    }
                }
            }

            double sumAbsoluteResiduals = double.MaxValue;

            for (int i = 0; i < maxIterations; i++)
            {
                // subtract row effects
                for (int r = 0; r < table.Length; r++)
                {
                    List<double> rowValues = table[r].Skip(1).Where(p => !double.IsNaN(p)).ToList();

                    if (rowValues.Any())
                    {
                        double rowMedian = rowValues.Median();
                        double[] weights = rowValues.Select(p => 1.0 / Math.Max(0.0001, Math.Pow(p - rowMedian, 2))).ToArray();
                        double rowEffect = rowValues.Sum(p => p * weights[rowValues.IndexOf(p)]) / weights.Sum();
                        table[r][0] += rowEffect;

                        for (int c = 1; c < table[0].Length; c++)
                        {
                            table[r][c] -= rowEffect;
                        }
                    }
                }

                // subtract column effects
                for (int c = 0; c < table[0].Length; c++)
                {
                    List<double> colValues = table.Skip(1).Select(p => p[c]).Where(p => !double.IsNaN(p)).ToList();

                    if (colValues.Any())
                    {
                        double colMedian = colValues.Median();
                        double[] weights = colValues.Select(p => 1.0 / Math.Max(0.0001, Math.Pow(p - colMedian, 2))).ToArray();
                        double colEffect = colValues.Sum(p => p * weights[colValues.IndexOf(p)]) / weights.Sum();
                        table[0][c] += colEffect;

                        for (int r = 1; r < table.Length; r++)
                        {
                            table[r][c] -= colEffect;
                        }
                    }
                }

                // calculate sum of absolute residuals and end the algorithm if it is not improving
                double iterationSumAbsoluteResiduals = table.Skip(1).SelectMany(p => p.Skip(1)).Where(p => !double.IsNaN(p)).Sum(p => Math.Abs(p));

                if (Math.Abs((iterationSumAbsoluteResiduals - sumAbsoluteResiduals) / sumAbsoluteResiduals) < improvementCutoff)
                {
                    break;
                }

                sumAbsoluteResiduals = iterationSumAbsoluteResiduals;
            }
        }
    }
}
