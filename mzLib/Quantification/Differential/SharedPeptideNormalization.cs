using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Statistics;

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
/// (<see cref="SharedPeptideMedianHalf"/>). The half set has the same floor: below it no sample is shifted, the setting
/// is <see cref="None"/>, and a warning says why. "The samples" here are those with a value on at least one such
/// peptide: a lost run has no shift, and it does not count towards either set.
/// </para>
/// <para>
/// <b>The half set is settled.</b> Its peptides are missing from different samples, so one pass measures each sample
/// against means over different samples and pulls the shifts together. The pass is repeated on the shifted values until
/// the shifts stop moving, and the shifts are then centred to average 0. If they have not settled within 1,000 passes,
/// the first pass's shifts are used and a warning says so. On the every-sample set one pass already leaves every sample
/// with the same median difference, so it is not repeated.
/// </para>
/// <para>
/// <b>Samples the reference does not link.</b> The shifts are fixed only up to one constant per group of samples that
/// reference peptides link (two samples are linked when a chain of peptides, each with a value in both of a pair, joins
/// them). When the half set falls into two or more such groups, each group is settled and centred on its own, and a
/// warning says the groups have no common scale: nothing in the data says how far apart they are.
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

    /// <summary>The half set repeats its pass until no shift moves by more than this, in log2 units.</summary>
    private const double SettledWithin = 1e-10;

    /// <summary>The most passes the half set makes before it stops and says so.</summary>
    private const int MaximumPasses = 1000;

    /// <summary>Measures each sample's shift over every sample of <paramref name="table"/>.</summary>
    /// <param name="table">One stratum's table.</param>
    /// <param name="minimumPeptides">The floor on either reference set (GR-25); at least 1.</param>
    public static NormalizationResult Measure(ObservationTable table, int minimumPeptides = DefaultMinimumPeptides)
        => Measure(table, minimumPeptides, MaximumPasses);

    /// <param name="maximumPasses">The half set's limit on passes, lowered by a test to reach its warning.</param>
    /// <param name="minimumHalfSetPeptides">
    /// The half set's floor when not <paramref name="minimumPeptides"/>, lowered by a test to measure a few peptides.
    /// </param>
    internal static NormalizationResult Measure(ObservationTable table, int minimumPeptides, int maximumPasses,
        int? minimumHalfSetPeptides = null)
    {
        ArgumentNullException.ThrowIfNull(table);
        ArgumentOutOfRangeException.ThrowIfLessThan(minimumPeptides, 1);

        int samples = table.Samples.Count;
        var eligible = Enumerable.Range(0, table.Peptides.Count)
            .Where(p => IsReferenceEligible(table, table.Peptides[p]))
            .Select(p => (Peptide: p, Present: Enumerable.Range(0, samples).Count(s => !double.IsNaN(table.Log2Intensity(p, s)))))
            .ToList();

        // A sample with no value on any eligible peptide (a lost run) can have no shift. It does not count towards "every
        // sample" or "half the samples": otherwise one lost run takes the whole stratum off the every-sample set.
        int withValues = Enumerable.Range(0, samples)
            .Count(s => eligible.Any(e => !double.IsNaN(table.Log2Intensity(e.Peptide, s))));

        string setting = SharedPeptideMedian;
        var reference = eligible.Where(e => e.Present == withValues && withValues > 0).Select(e => e.Peptide).ToList();
        var warnings = new List<string>();
        if (reference.Count < minimumPeptides)
        {
            int inEverySample = reference.Count;
            setting = SharedPeptideMedianHalf;
            int needed = Math.Max(2, (withValues + 1) / 2);
            reference = eligible.Where(e => e.Present >= needed).Select(e => e.Peptide).ToList();
            int halfSetFloor = minimumHalfSetPeptides ?? minimumPeptides;
            if (reference.Count < halfSetFloor)
            {
                // Fewer than the floor is not trusted on either set (review round 4: one peptide shifted four samples).
                warnings.Add($"{inEverySample} reference peptide(s) have a value in every sample and {reference.Count} in " +
                    $"at least half, under the floor of {halfSetFloor}, so no sample is shifted ({None}).");
                var none = table.Samples.ToDictionary(s => s, _ => double.NaN, StringComparer.Ordinal);
                return new NormalizationResult(None, reference.Count, minimumPeptides, none, warnings);
            }
        }

        double[] shift = MedianDifferences(table, reference, new double[samples]);
        if (setting == SharedPeptideMedianHalf)
        {
            var groups = LinkedGroups(table, reference, samples);
            if (groups.Count > 1)
                warnings.Add($"The {SharedPeptideMedianHalf} reference peptides link the samples in {groups.Count} groups " +
                    $"that share none ({string.Join("; ", groups.Select(g => string.Join(", ", g.Select(i => table.Samples[i]))))}). " +
                    "Each group is settled and centred on its own, so the groups have no common scale.");
            Settle(table, reference, shift, groups, maximumPasses, warnings);
        }

        var shifts = new Dictionary<string, double>(StringComparer.Ordinal);
        var unshifted = new List<string>();
        for (int s = 0; s < samples; s++)
        {
            shifts[table.Samples[s]] = shift[s];
            if (double.IsNaN(shift[s]))
                unshifted.Add(table.Samples[s]);
        }

        if (unshifted.Count > 0)
            warnings.Insert(0, $"No reference peptide ({reference.Count}, {setting}) has a value in sample(s) " +
                $"{string.Join(", ", unshifted)}, so their shift is undefined and they are not shifted.");

        return new NormalizationResult(setting, reference.Count, minimumPeptides, shifts, warnings);
    }

    /// <summary>
    /// One pass: each sample's median, over the reference peptides with a value in it, of its value minus its shift minus
    /// the peptide's mean over the samples that have it, each also less its shift. NaN for a sample with no reference value.
    /// </summary>
    private static double[] MedianDifferences(ObservationTable table, List<int> reference, double[] shift)
    {
        int samples = shift.Length;
        double Shifted(int p, int s) => table.Log2Intensity(p, s) - (double.IsNaN(shift[s]) ? 0 : shift[s]);

        var differences = Enumerable.Range(0, samples).Select(_ => new List<double>()).ToArray();
        foreach (int p in reference)
        {
            double sum = 0;
            int count = 0;
            for (int s = 0; s < samples; s++)
            {
                double value = Shifted(p, s);
                if (double.IsNaN(value)) continue;
                sum += value;
                count++;
            }

            double mean = sum / count;
            for (int s = 0; s < samples; s++)
            {
                double value = Shifted(p, s);
                if (!double.IsNaN(value))
                    differences[s].Add(value - mean);
            }
        }

        // MathNet's median: the mean of the middle two for an even count, NaN when empty.
        return Array.ConvertAll(differences, d => d.Median());
    }

    /// <summary>
    /// The half set's shifts, settled. Its peptides are missing from different samples, so one pass measures each sample
    /// against peptide means taken over different samples, and the shifts come out pulled together. Repeating the pass on
    /// the shifted values, and adding up the steps, until no step moves a shift by more than <see cref="SettledWithin"/>
    /// removes that. The steps after the first are centred within each linked group, because the shifts are fixed only up
    /// to one constant per group; centring across groups instead turns that constant into a step that never shrinks
    /// (review round 4). If the steps have not settled within <paramref name="maximumPasses"/>, the first pass is used. The
    /// result is centred to average 0 within each group.
    /// </summary>
    private static void Settle(ObservationTable table, List<int> reference, double[] shift, List<int[]> groups,
        int maximumPasses, List<string> warnings)
    {
        double[] firstPass = (double[])shift.Clone();
        bool settled = false;
        for (int pass = 2; pass <= maximumPasses; pass++)
        {
            double[] step = MedianDifferences(table, reference, shift);
            Centre(step, groups);
            settled = step.All(d => double.IsNaN(d) || Math.Abs(d) <= SettledWithin);
            if (settled)
                break;
            for (int s = 0; s < shift.Length; s++)
                if (!double.IsNaN(step[s]))
                    shift[s] += step[s];
        }

        if (!settled)
        {
            Array.Copy(firstPass, shift, shift.Length);
            warnings.Add($"The {SharedPeptideMedianHalf} shifts did not settle within {maximumPasses} passes, so the first " +
                "pass's shifts are used.");
        }
        Centre(shift, groups);
    }

    /// <summary>
    /// The groups of samples the reference peptides link, each in sample order and the groups in order of their first
    /// sample. A sample with no reference value is in none.
    /// </summary>
    private static List<int[]> LinkedGroups(ObservationTable table, List<int> reference, int samples)
    {
        var parent = Enumerable.Range(0, samples).ToArray();
        int Root(int s)
        {
            while (parent[s] != s)
                s = parent[s] = parent[parent[s]];
            return s;
        }

        var hasValue = new bool[samples];
        foreach (int p in reference)
        {
            int first = -1;
            for (int s = 0; s < samples; s++)
            {
                if (double.IsNaN(table.Log2Intensity(p, s))) continue;
                hasValue[s] = true;
                if (first < 0)
                    first = s;
                else
                    parent[Root(s)] = Root(first);
            }
        }
        return Enumerable.Range(0, samples).Where(s => hasValue[s]).GroupBy(Root).Select(g => g.ToArray()).ToList();
    }

    /// <summary>Within each group, subtracts the mean of the group's defined values from each; NaN stays NaN.</summary>
    private static void Centre(double[] values, List<int[]> groups)
    {
        foreach (int[] group in groups)
        {
            var defined = group.Where(s => !double.IsNaN(values[s])).ToList();
            if (defined.Count == 0)
                continue;
            double mean = defined.Average(s => values[s]);
            foreach (int s in defined)
                values[s] -= mean;
        }
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
}
