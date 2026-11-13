using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Statistics;

namespace Quantification.Differential;

/// <summary>
/// The global shift of a between-levels contrast (GR-27; STAT1 M5b ruling 9; <see cref="DifferentialGlobalShift"/>): how
/// far the whole numerator side sits from the denominator side before normalization, so a reader can see what the
/// normalization took out. Only peptides no decoy, contaminant or entrapment group holds count, as for the normalization's
/// reference set (ruling 9a); a sample with no value on any of them, a lost run, is skipped (ruling 9b, GR-2). No shift is
/// measured for a slope.
/// </summary>
internal static class GlobalShift
{
    /// <summary>
    /// Label-free: over the peptides valued in every kept sample of both sides, the median of (the mean log2 value on the
    /// numerator side minus the mean on the denominator side). No such peptide: NaN and 0.
    /// </summary>
    /// <param name="table">The stratum's table BEFORE normalization.</param>
    public static DifferentialGlobalShift LabelFree(ObservationTable table, IEnumerable<string> numerator,
        IEnumerable<string> denominator, string stratum, string? basis)
    {
        ArgumentNullException.ThrowIfNull(table);
        var column = table.Samples.Select((s, i) => (s, i)).ToDictionary(t => t.s, t => t.i, StringComparer.Ordinal);
        var eligible = Enumerable.Range(0, table.Peptides.Count)
            .Where(p => SharedPeptideNormalization.IsReferenceEligible(table, table.Peptides[p])).ToList();
        List<int> Kept(IEnumerable<string> side) => side.Where(column.ContainsKey).Select(s => column[s])
            .Where(s => eligible.Any(p => !double.IsNaN(table.Log2Intensity(p, s)))).ToList();
        var num = Kept(numerator);
        var den = Kept(denominator);

        var differences = new List<double>();
        if (num.Count > 0 && den.Count > 0)
            foreach (int p in eligible)
            {
                if (num.Any(s => double.IsNaN(table.Log2Intensity(p, s))) || den.Any(s => double.IsNaN(table.Log2Intensity(p, s))))
                    continue;
                differences.Add(num.Average(s => table.Log2Intensity(p, s)) - den.Average(s => table.Log2Intensity(p, s)));
            }
        return Result(differences, stratum, basis);
    }

    /// <summary>
    /// Isobaric (ruling 9c): per ion in its plex run, the mean log2 value over that run's numerator channels minus the mean
    /// over its denominator channels, for ions valued in every one of them; the median over all ions of every plex run that
    /// holds both sides, reported with that ion count. A channel with no value on any eligible ion of a plex run is skipped
    /// in that run. No such ion: NaN and 0.
    /// </summary>
    /// <param name="table">The table BEFORE normalization.</param>
    public static DifferentialGlobalShift Isobaric(IsobaricObservationTable table, IEnumerable<string> numerator,
        IEnumerable<string> denominator, string stratum)
    {
        ArgumentNullException.ThrowIfNull(table);
        var numeratorIds = new HashSet<string>(numerator, StringComparer.Ordinal);
        var denominatorIds = new HashSet<string>(denominator, StringComparer.Ordinal);
        var differences = new List<double>();
        foreach (var run in Enumerable.Range(0, table.Rows.Count)
                     .Where(r => IsobaricNormalization.IsEligible(table, table.Rows[r].Peptide))
                     .GroupBy(r => (table.Rows[r].Plex, table.Rows[r].TechnicalReplicate)))
        {
            var ions = run.ToList();
            var channels = table.SamplesOfPlex(run.Key.Plex)
                .Where(s => ions.Any(r => !double.IsNaN(table.Log2Intensity(r, s)))).ToList();
            var num = channels.Where(s => numeratorIds.Contains(table.Samples[s])).ToList();
            var den = channels.Where(s => denominatorIds.Contains(table.Samples[s])).ToList();
            if (num.Count == 0 || den.Count == 0) continue;
            foreach (int r in ions)
            {
                if (num.Any(s => double.IsNaN(table.Log2Intensity(r, s))) || den.Any(s => double.IsNaN(table.Log2Intensity(r, s))))
                    continue;
                differences.Add(num.Average(s => table.Log2Intensity(r, s)) - den.Average(s => table.Log2Intensity(r, s)));
            }
        }
        return Result(differences, stratum, null);
    }

    private static DifferentialGlobalShift Result(List<double> differences, string stratum, string? basis) =>
        new(stratum, basis, differences.Count == 0 ? double.NaN : differences.Median(), differences.Count);
}
