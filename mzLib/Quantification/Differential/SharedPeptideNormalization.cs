using System;
using System.Collections.Generic;
using System.Linq;

namespace Quantification.Differential;

/// <summary>A measured per-sample normalization (<c>DEF-DIFF-NORM</c>).</summary>
/// <param name="Setting">
/// <see cref="SharedPeptideNormalization.SharedPeptideMedian"/> or
/// <see cref="SharedPeptideNormalization.SharedPeptideMedianHalf"/>: which reference set the shifts were measured on.
/// </param>
/// <param name="ReferenceSetSize">How many peptides were in the reference set.</param>
/// <param name="MinimumReferenceSetSize">Below this many peptides in every sample, the half-samples set was used (GR-25).</param>
/// <param name="PerSampleShift">Each sample's shift in log2 units, subtracted from its values; NaN when undefined.</param>
/// <param name="Warnings">Which samples have no shift, and why.</param>
public sealed record NormalizationResult(string Setting, int ReferenceSetSize, int MinimumReferenceSetSize,
    IReadOnlyDictionary<string, double> PerSampleShift, IReadOnlyList<string> Warnings);

/// <summary>
/// The shared-peptide median normalization of GR-2 / GR-12 (<c>DEF-DIFF-NORM</c>): each sample is shifted by its median
/// log2 difference from the across-sample reference, over the peptides with a value in every sample.
/// </summary>
/// <remarks>
/// <para>
/// <b>Reference.</b> A peptide's reference is the mean of its log2 values over the samples where it has one (the log of
/// the geometric mean, as in DESeq's median-of-ratios, Anders and Huber 2010). A sample's shift is the median, over the
/// reference peptides with a value in it, of its value minus the reference. The median of an even count is the mean of
/// the middle two.
/// </para>
/// <para>
/// <b>Reference set.</b> Peptides of decoy, contaminant or entrapment groups never enter it (<c>DEF-DIFF-ROWSET</c>). It
/// is the remaining peptides with a value in every sample; when fewer than
/// <see cref="DefaultMinimumPeptides"/> (GR-25), it is those with a value in at least half the samples and in at least 2
/// (<see cref="SharedPeptideMedianHalf"/>).
/// </para>
/// <para>
/// <b>Scope.</b> The shifts are measured over every sample of the table given, so measure one stratum at a time
/// (<see cref="ObservationTable.ForSamples"/>; GR-5). They are measured whatever the setting, because the metadata always
/// reports them; a PTM enrichment records <see cref="None"/> and is simply not shifted (GR-12).
/// </para>
/// </remarks>
public static class SharedPeptideNormalization
{
    /// <summary><c>DEF-DIFF-NORM</c>: shifted on the peptides with a value in every sample.</summary>
    public const string SharedPeptideMedian = "shared_peptide_median";

    /// <summary><c>DEF-DIFF-NORM</c>: shifted on the peptides with a value in at least half the samples.</summary>
    public const string SharedPeptideMedianHalf = "shared_peptide_median_half";

    /// <summary><c>DEF-DIFF-NORM</c>: not shifted.</summary>
    public const string None = "none";

    /// <summary>GR-25: the every-sample set must hold at least this many peptides.</summary>
    public const int DefaultMinimumPeptides = 100;

    /// <summary>Measures each sample's shift over every sample of <paramref name="table"/>.</summary>
    /// <param name="table">One stratum's table.</param>
    /// <param name="minimumPeptides">The every-sample set's floor (GR-25); at least 1.</param>
    public static NormalizationResult Measure(ObservationTable table, int minimumPeptides = DefaultMinimumPeptides)
    {
        ArgumentNullException.ThrowIfNull(table);
        ArgumentOutOfRangeException.ThrowIfLessThan(minimumPeptides, 1);

        int samples = table.Samples.Count;
        var eligible = Enumerable.Range(0, table.Peptides.Count)
            .Where(p => IsReferenceEligible(table, table.Peptides[p]))
            .Select(p => (Peptide: p, Present: Enumerable.Range(0, samples).Count(s => !double.IsNaN(table.Log2Intensity(p, s)))))
            .ToList();

        string setting = SharedPeptideMedian;
        var reference = eligible.Where(e => e.Present == samples && samples > 0).Select(e => e.Peptide).ToList();
        if (reference.Count < minimumPeptides)
        {
            setting = SharedPeptideMedianHalf;
            int needed = Math.Max(2, (samples + 1) / 2);
            reference = eligible.Where(e => e.Present >= needed).Select(e => e.Peptide).ToList();
        }

        // Each sample's differences from the peptide's mean over the samples that have it.
        var differences = Enumerable.Range(0, samples).Select(_ => new List<double>()).ToArray();
        foreach (int p in reference)
        {
            double sum = 0;
            int count = 0;
            for (int s = 0; s < samples; s++)
            {
                double value = table.Log2Intensity(p, s);
                if (double.IsNaN(value)) continue;
                sum += value;
                count++;
            }

            double mean = sum / count;
            for (int s = 0; s < samples; s++)
            {
                double value = table.Log2Intensity(p, s);
                if (!double.IsNaN(value))
                    differences[s].Add(value - mean);
            }
        }

        var shifts = new Dictionary<string, double>(StringComparer.Ordinal);
        var unshifted = new List<string>();
        for (int s = 0; s < samples; s++)
        {
            shifts[table.Samples[s]] = Median(differences[s]);
            if (differences[s].Count == 0)
                unshifted.Add(table.Samples[s]);
        }

        var warnings = new List<string>();
        if (unshifted.Count > 0)
            warnings.Add($"No reference peptide ({reference.Count}, {setting}) has a value in sample(s) " +
                $"{string.Join(", ", unshifted)}, so their shift is undefined and they are not shifted.");

        return new NormalizationResult(setting, reference.Count, minimumPeptides, shifts, warnings);
    }

    /// <summary>
    /// The table with each sample's shift subtracted from its values. A sample whose shift is NaN is left as it is (the
    /// result's warnings say why).
    /// </summary>
    /// <exception cref="ArgumentException">The result has no shift for one of the table's samples.</exception>
    public static ObservationTable Apply(ObservationTable table, NormalizationResult result)
    {
        ArgumentNullException.ThrowIfNull(table);
        ArgumentNullException.ThrowIfNull(result);

        var shifts = new double[table.Samples.Count];
        for (int s = 0; s < shifts.Length; s++)
        {
            if (!result.PerSampleShift.TryGetValue(table.Samples[s], out double shift))
                throw new ArgumentException($"The normalization has no shift for sample '{table.Samples[s]}'.", nameof(result));
            shifts[s] = double.IsNaN(shift) ? 0 : shift;
        }

        return table.WithLog2((p, s) => table.Log2Intensity(p, s) - shifts[s]);
    }

    /// <summary>Whether a peptide may be in the reference set: none of its groups is a decoy, contaminant or entrapment group.</summary>
    public static bool IsReferenceEligible(ObservationTable table, ObservedPeptide peptide)
    {
        ArgumentNullException.ThrowIfNull(table);
        ArgumentNullException.ThrowIfNull(peptide);

        foreach (string name in peptide.ProteinGroups)
            if (table.ProteinGroups.TryGetValue(name, out var group) && (group.IsDecoy || group.IsContaminant || group.IsEntrapment))
                return false;
        return true;
    }

    /// <summary>The median; for an even count the mean of the middle two; NaN when empty.</summary>
    private static double Median(List<double> values)
    {
        if (values.Count == 0)
            return double.NaN;

        var sorted = values.ToArray();
        Array.Sort(sorted);
        int middle = sorted.Length / 2;
        return sorted.Length % 2 == 1 ? sorted[middle] : (sorted[middle - 1] + sorted[middle]) / 2;
    }
}
