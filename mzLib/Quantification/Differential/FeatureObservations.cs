using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace Quantification.Differential;

/// <summary>
/// One value a per-protein model fits (STAT1 M5b): its log2 value; its nuisance level (a peptide, or a TMT peptide ion in
/// its plex run); the analysis sample it belongs to (the design row, and the sample counted by the support rules); and
/// its random-effect labels.
/// </summary>
/// <param name="Log2">The normalized log2 value; finite.</param>
/// <param name="Nuisance">The nuisance level: the peptide's full sequence (label-free), or ion and plex run (isobaric).</param>
/// <param name="SampleId">The <see cref="AnalysisSample.SampleId"/> whose design row the value is fitted on.</param>
/// <param name="RandomSample">The inner random intercept's level.</param>
/// <param name="RandomIndividual">The outer random intercept's level; null when the model has only one.</param>
/// <param name="IsMbr">The value includes a match-between-runs peak.</param>
internal readonly record struct FeatureValue(double Log2, string Nuisance, string SampleId, string RandomSample,
    string? RandomIndividual, bool IsMbr);

/// <summary>
/// One protein group's values in one stratum (STAT1 M5b rulings 6, 13): only its usable nuisance levels (a peptide or ion
/// unique to the group, GR-23, with a value in at least the minimum number of the stratum's samples), and how many unique
/// levels there were before that rule, for the decisions record.
/// </summary>
internal sealed class FeatureObservations
{
    internal FeatureObservations(string group, IReadOnlyList<FeatureValue> values, int uniqueLevels)
    {
        Group = group;
        Values = values;
        UniqueLevels = uniqueLevels;
        NuisanceLevels = values.Select(v => v.Nuisance).Distinct(StringComparer.Ordinal).Order(StringComparer.Ordinal).ToArray();
        SampleIds = values.Select(v => v.SampleId).Distinct(StringComparer.Ordinal).Order(StringComparer.Ordinal).ToArray();
    }

    /// <summary>The protein group's name.</summary>
    public string Group { get; }

    /// <summary>The values of the usable nuisance levels, in a fixed order (level, then sample, ordinal).</summary>
    public IReadOnlyList<FeatureValue> Values { get; }

    /// <summary>The usable nuisance levels, ordinal; the first is the reference of the treatment coding.</summary>
    public IReadOnlyList<string> NuisanceLevels { get; }

    /// <summary>Nuisance levels unique to the group, usable or not.</summary>
    public int UniqueLevels { get; }

    /// <summary>The analysis samples with at least one usable value, ordinal.</summary>
    public IReadOnlyList<string> SampleIds { get; }

    /// <summary>How many of <paramref name="sampleIds"/> have at least one usable value.</summary>
    public int SamplesWithValue(IEnumerable<string> sampleIds)
    {
        var mine = SampleIds.ToHashSet(StringComparer.Ordinal);
        return sampleIds.Count(mine.Contains);
    }
}

/// <summary>Builds <see cref="FeatureObservations"/> from either observation table, and the row set (DEF-DIFF-ROWSET).</summary>
internal static class FeatureObservationBuilder
{
    /// <summary>The protein q-value at or below which a target group enters the row set (DEF-PROTSET-1PCT).</summary>
    public const double RowSetQValue = 0.01;

    /// <summary>Why a group named by the table's peptides gets no row (ruling 12), or null when it gets one.</summary>
    public static string? RowSetExclusion(string group, IReadOnlyDictionary<string, ProteinGroupInfo> groups)
    {
        if (!groups.TryGetValue(group, out var info)) return "undescribed";
        if (info.IsDecoy) return "decoy";
        if (info.IsEntrapment) return "entrapment";
        if (info.IsContaminant) return "contaminant";
        if (!(info.QValue <= RowSetQValue)) return "above_1pct_q";
        return null;
    }

    /// <summary>
    /// Label-free: each peptide unique to a group, valued in at least <paramref name="minimumSamples"/> of
    /// <paramref name="samples"/>, gives one value per sample with a value. The random intercepts are the sample, within
    /// the individual when the samples carry one (GR-6).
    /// </summary>
    /// <param name="table">The stratum's normalized table.</param>
    /// <param name="samples">The stratum's analysis samples, by id.</param>
    /// <param name="minimumSamples">The usable-peptide minimum (ruling 6: 2).</param>
    public static IReadOnlyDictionary<string, FeatureObservations> FromTable(ObservationTable table,
        IReadOnlyDictionary<string, AnalysisSample> samples, int minimumSamples)
    {
        var columns = Enumerable.Range(0, table.Samples.Count).Where(s => samples.ContainsKey(table.Samples[s])).ToArray();
        var values = new Dictionary<string, List<FeatureValue>>(StringComparer.Ordinal);
        var unique = new Dictionary<string, int>(StringComparer.Ordinal);
        for (int p = 0; p < table.Peptides.Count; p++)
        {
            var peptide = table.Peptides[p];
            string? group = peptide.UniqueProteinGroup;
            if (group is null) continue;
            unique[group] = unique.GetValueOrDefault(group) + 1;
            var valued = columns.Where(s => !double.IsNaN(table.Log2Intensity(p, s))).ToArray();
            if (valued.Length < minimumSamples) continue;
            if (!values.TryGetValue(group, out var list)) values[group] = list = new List<FeatureValue>();
            foreach (int s in valued)
            {
                string id = table.Samples[s];
                list.Add(new FeatureValue(table.Log2Intensity(p, s), peptide.FullSequence, id, id, samples[id].Individual,
                    table.State(p, s) == ObservationState.MbrTransferred));
            }
        }
        return Collect(values, unique);
    }

    /// <summary>
    /// Isobaric: each ion in its plex run (an <see cref="IsobaricIonRow"/>) whose peptide is unique to a group, valued in at
    /// least <paramref name="minimumSamples"/> of the stratum's channel samples, gives one value per channel with a value.
    /// Random intercepts (ruling 13, model A): the channel's sample, then the channel in each plex run (repeat injections);
    /// when the samples name a real individual (one animal in several plexes), the individual, then the channel's sample.
    /// </summary>
    public static IReadOnlyDictionary<string, FeatureObservations> FromIsobaric(IsobaricObservationTable table,
        IReadOnlyDictionary<string, AnalysisSample> samples, int minimumSamples)
    {
        var values = new Dictionary<string, List<FeatureValue>>(StringComparer.Ordinal);
        var unique = new Dictionary<string, int>(StringComparer.Ordinal);
        for (int r = 0; r < table.Rows.Count; r++)
        {
            var row = table.Rows[r];
            string? group = row.Peptide.UniqueProteinGroup;
            if (group is null) continue;
            unique[group] = unique.GetValueOrDefault(group) + 1;
            var valued = table.SamplesOfPlex(row.Plex)
                .Where(s => samples.ContainsKey(table.Samples[s]) && !double.IsNaN(table.Log2Intensity(r, s))).ToArray();
            if (valued.Length < minimumSamples) continue;
            string run = string.Create(CultureInfo.InvariantCulture, $"plex{row.Plex}.run{row.TechnicalReplicate}");
            string level = string.Create(CultureInfo.InvariantCulture, $"{row.Peptide.FullSequence}/{row.Charge}@{run}");
            if (!values.TryGetValue(group, out var list)) values[group] = list = new List<FeatureValue>();
            foreach (int s in valued)
            {
                string id = table.Samples[s];
                string? animal = samples[id].Individual;
                list.Add(animal is null
                    ? new FeatureValue(table.Log2Intensity(r, s), level, id, id + "@" + run, id, false)
                    : new FeatureValue(table.Log2Intensity(r, s), level, id, id, animal, false));
            }
        }
        return Collect(values, unique);
    }

    private static IReadOnlyDictionary<string, FeatureObservations> Collect(Dictionary<string, List<FeatureValue>> values,
        Dictionary<string, int> unique)
    {
        var result = new Dictionary<string, FeatureObservations>(StringComparer.Ordinal);
        foreach (var (group, count) in unique)
        {
            var list = values.TryGetValue(group, out var v) ? v : new List<FeatureValue>();
            var ordered = list.OrderBy(x => x.Nuisance, StringComparer.Ordinal).ThenBy(x => x.SampleId, StringComparer.Ordinal)
                .ThenBy(x => x.RandomSample, StringComparer.Ordinal).ToArray();
            result[group] = new FeatureObservations(group, ordered, count);
        }
        return result;
    }
}
