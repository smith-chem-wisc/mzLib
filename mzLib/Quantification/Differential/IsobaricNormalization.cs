using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Statistics;

namespace Quantification.Differential;

/// <summary>One channel's shift in one plex run (<c>DEF-DIFF-NORM</c>, isobaric).</summary>
/// <param name="Plex">The plex.</param>
/// <param name="TechnicalReplicate">The plex's technical replicate.</param>
/// <param name="SampleId">The channel's sample.</param>
/// <param name="Values">How many values the channel's median was taken over.</param>
/// <param name="Shift">The shift in log2 units, subtracted from the channel's values in that plex run; NaN when undefined.</param>
public sealed record IsobaricChannelShift(int Plex, int TechnicalReplicate, string SampleId, int Values, double Shift);

/// <summary>A measured isobaric normalization (<c>DEF-DIFF-NORM</c>).</summary>
/// <param name="Setting"><see cref="IsobaricNormalization.GlobalChannelMedian"/>, or <see cref="IsobaricNormalization.None"/> when nothing could be measured.</param>
/// <param name="Shifts">One per sample channel of every plex run with rows, ordered by plex, technical replicate and sample.</param>
/// <param name="Warnings">Which channels have no shift, and why.</param>
public sealed record IsobaricNormalizationResult(string Setting, IReadOnlyList<IsobaricChannelShift> Shifts,
    IReadOnlyList<string> Warnings);

/// <summary>
/// The global channel-median normalization of an isobaric dataset (REQ-QP-15 PLAN 6f rule 3; MSstatsTMT 2.20's
/// <c>global_norm = TRUE</c>): each channel's median log2 intensity in each plex run is moved to the median of those
/// medians, over every plex run.
/// </summary>
/// <remarks>
/// <para>
/// <b>Medians.</b> A channel's median in a plex run is taken over that run's rows with a value in the channel. Rows of
/// decoy, contaminant or entrapment groups never enter it (<c>DEF-DIFF-ROWSET</c>, as
/// <see cref="SharedPeptideNormalization"/>). The median of an even count is the mean of the middle two.
/// </para>
/// <para>
/// <b>Target.</b> The median of every channel median that exists. Reference and unused channels are not in the table, so
/// they do not enter it; MSstatsTMT's target includes them. The target is one constant for every value, so it changes no
/// comparison.
/// </para>
/// <para>
/// <b>Reference channels.</b> GR-20's reference-channel shift is not applied here: it moves an ion's values in one plex
/// run by one amount in every channel, which a model with a term per ion per plex run absorbs (PLAN 6h).
/// </para>
/// <para>
/// <b>Never throws on data.</b> A channel with no value in a plex run has a NaN shift and is left as it is; a warning says
/// so. With no channel median at all, the setting is <see cref="None"/>.
/// </para>
/// </remarks>
public static class IsobaricNormalization
{
    /// <summary><c>DEF-DIFF-NORM</c>: each channel's median in each plex run moved to the median of those medians.</summary>
    public const string GlobalChannelMedian = "global_channel_median";

    /// <summary><c>DEF-DIFF-NORM</c>: not shifted.</summary>
    public const string None = "none";

    /// <summary>Measures every channel's shift in every plex run of <paramref name="table"/>.</summary>
    public static IsobaricNormalizationResult Measure(IsobaricObservationTable table)
    {
        ArgumentNullException.ThrowIfNull(table);

        var medians = new List<(int Plex, int TechnicalReplicate, int Sample, int Values, double Median)>();
        foreach (var run in Enumerable.Range(0, table.Rows.Count)
                     .GroupBy(r => (table.Rows[r].Plex, table.Rows[r].TechnicalReplicate))
                     .OrderBy(g => g.Key.Plex).ThenBy(g => g.Key.TechnicalReplicate))
        {
            var eligible = run.Where(r => IsEligible(table, table.Rows[r].Peptide)).ToList();
            foreach (int s in table.SamplesOfPlex(run.Key.Plex))
            {
                var values = eligible.Select(r => table.Log2Intensity(r, s)).Where(v => !double.IsNaN(v)).ToList();
                // MathNet's median: the mean of the middle two for an even count, NaN when empty.
                medians.Add((run.Key.Plex, run.Key.TechnicalReplicate, s, values.Count, values.Median()));
            }
        }

        var defined = medians.Where(m => !double.IsNaN(m.Median)).Select(m => m.Median).ToList();
        double target = defined.Median();
        var shifts = medians.Select(m => new IsobaricChannelShift(m.Plex, m.TechnicalReplicate, table.Samples[m.Sample],
            m.Values, m.Median - target)).ToList();

        var warnings = new List<string>();
        var undefined = shifts.Where(s => double.IsNaN(s.Shift)).ToList();
        if (defined.Count == 0 && shifts.Count > 0)
            warnings.Add("No channel has a value outside decoy, contaminant and entrapment groups, so nothing is shifted.");
        else if (undefined.Count > 0)
            warnings.Add($"{undefined.Count} channel(s) have no value in a plex run outside decoy, contaminant and " +
                "entrapment groups, so they are not shifted there: " +
                string.Join(", ", undefined.Take(5).Select(s => $"'{s.SampleId}' in plex {s.Plex}, technical replicate {s.TechnicalReplicate}")) +
                (undefined.Count > 5 ? ", ..." : "") + ".");

        return new IsobaricNormalizationResult(defined.Count > 0 ? GlobalChannelMedian : None, shifts, warnings);
    }

    /// <summary>
    /// The table with each channel's shift subtracted from its values in each plex run. A NaN shift leaves the channel as
    /// it is (the result's warnings say why).
    /// </summary>
    /// <exception cref="ArgumentException">The result has no shift for a channel of a plex run the table has.</exception>
    public static IsobaricObservationTable Apply(IsobaricObservationTable table, IsobaricNormalizationResult result)
    {
        ArgumentNullException.ThrowIfNull(table);
        ArgumentNullException.ThrowIfNull(result);

        var sampleIndex = table.Samples.Select((id, i) => (id, i)).ToDictionary(t => t.id, t => t.i, StringComparer.Ordinal);
        var shift = new Dictionary<(int, int, int), double>();
        foreach (var s in result.Shifts)
            if (sampleIndex.TryGetValue(s.SampleId, out int i))
                shift[(s.Plex, s.TechnicalReplicate, i)] = double.IsNaN(s.Shift) ? 0 : s.Shift;

        foreach (var row in table.Rows.Select(r => (r.Plex, r.TechnicalReplicate)).Distinct())
            foreach (int s in table.SamplesOfPlex(row.Plex))
                if (!shift.ContainsKey((row.Plex, row.TechnicalReplicate, s)))
                    throw new ArgumentException(
                        $"The normalization has no shift for sample '{table.Samples[s]}' in plex {row.Plex}, technical " +
                        $"replicate {row.TechnicalReplicate}.", nameof(result));

        return table.WithLog2((r, s) =>
            table.Log2Intensity(r, s) - shift[(table.Rows[r].Plex, table.Rows[r].TechnicalReplicate, s)]);
    }

    /// <summary>Whether a row may enter a median: none of its groups is a decoy, contaminant or entrapment group.</summary>
    private static bool IsEligible(IsobaricObservationTable table, ObservedPeptide peptide)
    {
        foreach (string name in peptide.ProteinGroups)
            if (table.ProteinGroups.TryGetValue(name, out var group) && (group.IsDecoy || group.IsContaminant || group.IsEntrapment))
                return false;
        return true;
    }
}
