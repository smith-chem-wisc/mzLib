using Easy.Common.Extensions;
using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace FlashLFQ
{
    /// <summary>
    /// Median polish protein quantification strategy: uses the median polish algorithm
    /// to calculate protein quantities in each biological replicate
    /// See https://mgimond.github.io/ES218/Week11a.html for an example of the median polish algorithm.
    /// </summary>
    public class MedianPolishProteinQuantificationStrategy<TPeptide, TProteinGroup, TSpectraFile> 
        : IProteinQuantificationStrategy<TPeptide, TProteinGroup, TSpectraFile>
        where TPeptide : IQuantifiablePeptide
        where TProteinGroup : IQuantifiableProteinGroup
        where TSpectraFile : IQuantifiableSpectraFile
    {
        private readonly int _maxIterations;
        private readonly double _improvementCutoff;

        public MedianPolishProteinQuantificationStrategy(
            int maxIterations = 10, 
            double improvementCutoff = 0.0001)
        {
            _maxIterations = maxIterations;
            _improvementCutoff = improvementCutoff;
        }

        public void QuantifyProteins(
            Dictionary<TProteinGroup, List<TPeptide>> proteinGroupToPeptides,
            List<TSpectraFile> spectraFiles,
            bool useSharedPeptides)
        {
            // Quantify each protein
            foreach (var kvp in proteinGroupToPeptides)
            {
                TProteinGroup proteinGroup = kvp.Key;
                List<TPeptide> peptidesForThisProtein = kvp.Value;

                // Set up peptide intensity table
                int numSamples = spectraFiles
                    .Select(p => p.Condition + p.BiologicalReplicate)
                    .Distinct()
                    .Count();

                double[][] peptideIntensityMatrix = new double[peptidesForThisProtein.Count + 1][];
                for (int i = 0; i < peptideIntensityMatrix.Length; i++)
                {
                    peptideIntensityMatrix[i] = new double[numSamples + 1];
                }

                // Populate matrix with log2-transformed peptide intensities
                int sampleN = 0;
                foreach (var group in spectraFiles.GroupBy(p => p.Condition).OrderBy(p => p.Key))
                {
                    foreach (var sample in group.GroupBy(p => p.BiologicalReplicate).OrderBy(p => p.Key))
                    {
                        foreach (TPeptide peptide in peptidesForThisProtein)
                        {
                            double sampleIntensity = CalculateSampleIntensity(peptide, sample.ToList());

                            if (sampleIntensity == 0)
                            {
                                sampleIntensity = double.NaN;
                            }
                            else
                            {
                                sampleIntensity = Math.Log(sampleIntensity, 2);
                            }

                            peptideIntensityMatrix[peptidesForThisProtein.IndexOf(peptide) + 1][sampleN + 1] = sampleIntensity;
                        }

                        sampleN++;
                    }
                }

                // Filter out peptides with only one measurement (unless all peptides have one measurement)
                FilterSingleMeasurementPeptides(peptideIntensityMatrix, numSamples);

                // Do median polish protein quantification
                MedianPolish(peptideIntensityMatrix, _maxIterations, _improvementCutoff);

                // Extract results and set protein intensities
                double overallEffect = peptideIntensityMatrix[0][0];
                double[] columnEffects = peptideIntensityMatrix[0].Skip(1).ToArray();
                double referenceProteinIntensity = Math.Pow(2, overallEffect) * peptidesForThisProtein.Count;

                // Check for unquantifiable proteins
                List<string> possibleUnquantifiableSample = IdentifyUnquantifiableProteins(
                    spectraFiles, peptidesForThisProtein, columnEffects);

                // Set the sample protein intensities
                SetProteinIntensities(
                    proteinGroup, 
                    spectraFiles, 
                    peptidesForThisProtein, 
                    columnEffects, 
                    referenceProteinIntensity, 
                    possibleUnquantifiableSample);
            }
        }

        private double CalculateSampleIntensity(TPeptide peptide, List<TSpectraFile> sampleFiles)
        {
            double sampleIntensity = 0;
            double highestFractionIntensity = 0;

            // The fraction with the highest intensity is used as the sample intensity
            foreach (var fraction in sampleFiles.GroupBy(p => p.Fraction))
            {
                double fractionIntensity = 0;
                int replicatesWithValidValues = 0;

                foreach (TSpectraFile replicate in fraction.OrderBy(p => p.TechnicalReplicate))
                {
                    double replicateIntensity = peptide.GetIntensity(replicate);

                    if (replicateIntensity > 0)
                    {
                        fractionIntensity += replicateIntensity;
                        replicatesWithValidValues++;
                    }
                }

                if (replicatesWithValidValues > 0)
                {
                    fractionIntensity /= replicatesWithValidValues;
                }

                if (fractionIntensity > highestFractionIntensity)
                {
                    highestFractionIntensity = fractionIntensity;
                    sampleIntensity = highestFractionIntensity;
                }
            }

            return sampleIntensity;
        }

        private void FilterSingleMeasurementPeptides(double[][] peptideIntensityMatrix, int numSamples)
        {
            var peptidesWithMoreThanOneMmt = peptideIntensityMatrix
                .Skip(1)
                .Count(row => row.Skip(1).Count(cell => !double.IsNaN(cell)) > 1);

            if (peptidesWithMoreThanOneMmt > 0)
            {
                for (int i = 1; i < peptideIntensityMatrix.Length; i++)
                {
                    int validValueCount = peptideIntensityMatrix[i]
                        .Count(p => !double.IsNaN(p) && p != 0);

                    if (validValueCount < 2 && numSamples >= 2)
                    {
                        for (int j = 1; j < peptideIntensityMatrix[0].Length; j++)
                        {
                            peptideIntensityMatrix[i][j] = double.NaN;
                        }
                    }
                }
            }
        }

        private List<string> IdentifyUnquantifiableProteins(
            List<TSpectraFile> spectraFiles,
            List<TPeptide> peptidesForThisProtein,
            double[] columnEffects)
        {
            List<string> possibleUnquantifiableSample = new List<string>();
            int sampleN = 0;

            foreach (var group in spectraFiles.GroupBy(p => p.Condition).OrderBy(p => p.Key))
            {
                foreach (var sample in group.GroupBy(p => p.BiologicalReplicate).OrderBy(p => p.Key))
                {
                    bool isMissingValue = true;

                    foreach (TSpectraFile spectraFile in sample)
                    {
                        if (peptidesForThisProtein.Any(p => p.GetIntensity(spectraFile) != 0))
                        {
                            isMissingValue = false;
                            break;
                        }
                    }

                    if (!isMissingValue && columnEffects[sampleN] == 0)
                    {
                        possibleUnquantifiableSample.Add(group.Key + "_" + sample.Key);
                    }

                    sampleN++;
                }
            }

            return possibleUnquantifiableSample;
        }

        private void SetProteinIntensities(
            TProteinGroup proteinGroup,
            List<TSpectraFile> spectraFiles,
            List<TPeptide> peptidesForThisProtein,
            double[] columnEffects,
            double referenceProteinIntensity,
            List<string> possibleUnquantifiableSample)
        {
            int sampleN = 0;

            foreach (var group in spectraFiles.GroupBy(p => p.Condition).OrderBy(p => p.Key))
            {
                foreach (var sample in group.GroupBy(p => p.BiologicalReplicate).OrderBy(p => p.Key))
                {
                    double columnEffect = columnEffects[sampleN];
                    double sampleProteinIntensity = Math.Pow(2, columnEffect) * referenceProteinIntensity;

                    // Check if this is a valid value
                    bool isMissingValue = true;

                    foreach (TSpectraFile spectraFile in sample)
                    {
                        if (peptidesForThisProtein.Any(p => p.GetIntensity(spectraFile) != 0))
                        {
                            isMissingValue = false;
                            break;
                        }
                    }

                    if (!isMissingValue)
                    {
                        if (possibleUnquantifiableSample.Count > 1 && 
                            possibleUnquantifiableSample.Contains(group.Key + "_" + sample.Key))
                        {
                            proteinGroup.SetIntensity(sample.First(), double.NaN);
                        }
                        else
                        {
                            proteinGroup.SetIntensity(sample.First(), sampleProteinIntensity);
                        }
                    }

                    sampleN++;
                }
            }
        }

        /// <summary>
        /// Performs the median polish algorithm on a data table
        /// Technically, this is weighted mean polish and not median polish, but it gives similar results 
        /// while being more robust to issues arising from missing values.
        /// </summary>
        private void MedianPolish(double[][] table, int maxIterations, double improvementCutoff)
        {
            // Subtract overall effect
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
                // Subtract row effects
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

                // Subtract column effects
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

                // Calculate sum of absolute residuals and end if not improving
                double iterationSumAbsoluteResiduals = table.Skip(1)
                    .SelectMany(p => p.Skip(1))
                    .Where(p => !double.IsNaN(p))
                    .Sum(p => Math.Abs(p));

                if (Math.Abs((iterationSumAbsoluteResiduals - sumAbsoluteResiduals) / sumAbsoluteResiduals) < improvementCutoff)
                {
                    break;
                }

                sumAbsoluteResiduals = iterationSumAbsoluteResiduals;
            }
        }
    }
}
