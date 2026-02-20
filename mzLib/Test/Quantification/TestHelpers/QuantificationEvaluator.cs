using System;
using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;
using MassSpectrometry.ExperimentalDesign;
using Omics.BioPolymerGroup;
using Quantification;

namespace Test.Quantification.TestHelpers;

/// <summary>
/// Helper utilities for evaluating quantification results against known ground truth.
/// Computes per-condition mean intensities and fold changes from a protein matrix.
/// </summary>
public static class QuantificationEvaluator
{
    /// <summary>
    /// Calculates the mean intensity for each condition across all columns in the protein matrix
    /// that belong to that condition (reference channels are excluded).
    /// </summary>
    /// <param name="proteinMatrix">The protein-level quantity matrix.</param>
    /// <param name="protein">The protein group to look up.</param>
    /// <returns>
    /// Dictionary mapping condition name → mean intensity (NaN if no non-zero values exist
    /// for that condition).
    /// </returns>
    public static Dictionary<string, double> GetMeanIntensityByCondition(
        QuantMatrix<IBioPolymerGroup> proteinMatrix,
        IBioPolymerGroup protein)
    {
        var row = proteinMatrix.GetRow(protein);
        var columnKeys = proteinMatrix.ColumnKeys;

        // Group column indices by condition, skipping reference channels
        var conditionValues = new Dictionary<string, List<double>>();

        for (int i = 0; i < columnKeys.Count; i++)
        {
            if (columnKeys[i] is not IsobaricQuantSampleInfo isoInfo)
                continue;

            if (isoInfo.IsReferenceChannel)
                continue;

            string cond = isoInfo.Condition;
            if (!conditionValues.ContainsKey(cond))
                conditionValues[cond] = new List<double>();

            if (row[i] > 0)
                conditionValues[cond].Add(row[i]);
        }

        return conditionValues.ToDictionary(
            kvp => kvp.Key,
            kvp => kvp.Value.Count > 0 ? kvp.Value.Average() : double.NaN);
    }

    /// <summary>
    /// Calculates the fold change between two conditions: numerator / denominator.
    /// Returns NaN if either condition is missing or has value 0 / NaN.
    /// </summary>
    public static double CalculateFoldChange(
        Dictionary<string, double> meanIntensities,
        string numeratorCondition,
        string denominatorCondition)
    {
        if (!meanIntensities.TryGetValue(numeratorCondition, out double num) || double.IsNaN(num) || num == 0)
            return double.NaN;
        if (!meanIntensities.TryGetValue(denominatorCondition, out double den) || double.IsNaN(den) || den == 0)
            return double.NaN;
        return num / den;
    }

    /// <summary>
    /// For each protein in the matrix, computes the fold change between two conditions,
    /// returning a list of (accession, foldChange) pairs where fold change is finite.
    /// </summary>
    public static List<(string accession, double foldChange)> ComputeFoldChanges(
        QuantMatrix<IBioPolymerGroup> proteinMatrix,
        string numeratorCondition,
        string denominatorCondition)
    {
        var results = new List<(string, double)>();
        foreach (var protein in proteinMatrix.RowKeys)
        {
            var meansByCondition = GetMeanIntensityByCondition(proteinMatrix, protein);
            double fc = CalculateFoldChange(meansByCondition, numeratorCondition, denominatorCondition);
            if (!double.IsNaN(fc) && !double.IsInfinity(fc))
                results.Add((protein.BioPolymerGroupName, fc));
        }
        return results;
    }

    /// <summary>
    /// Returns the median of a sequence of values.
    /// </summary>
    public static double Median(IEnumerable<double> values)
    {
        var sorted = values.OrderBy(x => x).ToList();
        if (sorted.Count == 0) return double.NaN;
        int mid = sorted.Count / 2;
        return sorted.Count % 2 == 0
            ? (sorted[mid - 1] + sorted[mid]) / 2.0
            : sorted[mid];
    }
}
