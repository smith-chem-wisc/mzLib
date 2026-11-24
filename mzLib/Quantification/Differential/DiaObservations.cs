using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using MassSpectrometry;

namespace Quantification.Differential;

/// <summary>
/// One DIA precursor's quantity in one run, as a search engine reported it: the input to <see cref="DiaObservations.Build"/>.
/// It names no engine's format; an adapter maps an engine's precursor table (MetaMorpheus's <c>AllDiaPrecursors.tsv</c>,
/// DIA-NN's report) to it.
/// </summary>
/// <param name="FileName">
/// The spectra file's name without its extension, as the experimental design names the run. A name that ends in
/// <c>-calib</c> also matches the design file without it, as in <see cref="ObservationTable"/>.
/// </param>
/// <param name="FullSequence">The engine's full sequence.</param>
/// <param name="BaseSequence">The unmodified sequence.</param>
/// <param name="Charge">The precursor charge; at least 1.</param>
/// <param name="ProteinGroups">
/// Every protein group the precursor's peptide maps to, named as the producer's protein-group table names them (GR-23), in
/// any order. Never the library's accessions: those are not groups.
/// </param>
/// <param name="Quantity">
/// The quantity when <paramref name="Kind"/> is <see cref="PrecursorQuantityKind.Value"/>, finite and above 0; NaN for the
/// other kinds, never 0.
/// </param>
/// <param name="Kind">What the quantity is (<c>DEF-DIA-KIND</c>).</param>
/// <param name="RunQValue">The precursor's q-value within this run, from 0 to 1.</param>
/// <param name="GlobalQValue">The precursor's q-value over every run of the search, from 0 to 1.</param>
/// <param name="IsDecoy">The precursor is a decoy.</param>
/// <param name="IsEntrapment">The precursor is from an entrapment sequence.</param>
/// <param name="IsContaminant">The precursor is a contaminant's: every protein it maps to is in a contaminant database.</param>
public sealed record DiaPrecursorRun(string FileName, string FullSequence, string BaseSequence, int Charge,
    IReadOnlyList<string> ProteinGroups, double Quantity, PrecursorQuantityKind Kind, double RunQValue, double GlobalQValue,
    bool IsDecoy = false, bool IsEntrapment = false, bool IsContaminant = false);

/// <summary>
/// Builds the differential-statistics input (<see cref="ObservationTable"/>) from DIA precursors (QuantProject dia 067,
/// 072, 076, 078): one row per precursor (full sequence and charge), each its own level in the model, never summed over
/// charges. Build the request with <c>QuantStyle = "dia"</c>.
/// </summary>
/// <remarks>
/// <para>
/// <b>States</b> (<c>DEF-DIA-KIND</c>). A precursor that passes both q-values in a run maps by its kind:
/// <see cref="PrecursorQuantityKind.Value"/> is <see cref="ObservationState.Quantified"/>;
/// <see cref="PrecursorQuantityKind.BelowFloor"/> and <see cref="PrecursorQuantityKind.NotLooked"/> are
/// <see cref="ObservationState.IdentifiedNotQuantified"/>, with no value. The table has no below-floor state of its own
/// yet, so the two kinds are one state here; the precursor table keeps them apart. A precursor that fails either q-value
/// in a run, or has no record there, is <see cref="ObservationState.NotDetected"/> in that run. No cell is ever 0.
/// </para>
/// <para>
/// <b>Which precursors enter.</b> Decoy, entrapment and contaminant precursors never enter, whatever groups they name: a
/// group's own flags keep such groups out of the analysis, but a precursor's group names are the producer's, and a decoy
/// can carry its owner's name. A target precursor that passes the q-values in no run has no row. The protein groups'
/// q-values and flags come from the producer's protein-group table, as for label-free.
/// </para>
/// <para>
/// <b>Runs.</b> Technical replicates of a sample are averaged on the log scale, as for every table. Fractions are
/// refused: two runs of one sample and technical replicate. A design run where no target precursor passes is a run where
/// nothing was reported. When no run of its sample reports anything, the sample is kept with no values
/// (<see cref="ObservationTable.SamplesWithoutValues"/>). When another run of its sample does, it is refused, because the
/// replicate average would leave it out with no record. The table is built under <see cref="QuantBasis.MsmsOnly"/>: DIA
/// has no match-between-runs values.
/// </para>
/// <para>
/// <b>Refusals</b> throw <see cref="ArgumentException"/> naming the precursor and run: a precursor twice in one run, runs
/// of one precursor that disagree on its base sequence or groups, a value that is not finite and above 0, a kind without a
/// value that carries a number, a kind the adapter does not map, a q-value that is not from 0 to 1, a run the design does
/// not have, a run lost inside its sample (naming the sample), and what <c>ObservationTable.Build</c> refuses.
/// </para>
/// </remarks>
public static class DiaObservations
{
    /// <summary>The run and global q-value at or below which a precursor's run enters (<c>DEF-DIFF-ROWSET</c> for <c>dia</c>).</summary>
    public const double DefaultQValueThreshold = 0.01;

    /// <summary>The suffix MetaMorpheus gives a calibrated file's name.</summary>
    private const string CalibratedSuffix = "-calib";

    /// <summary>The table from a search's precursors.</summary>
    /// <param name="runs">The experimental design's runs; each record's file matches one by name.</param>
    /// <param name="precursors">Every precursor in every run the search reported, decoys and failing q-values included or not.</param>
    /// <param name="proteinGroups">The producer's protein groups, with their q-values and flags.</param>
    /// <param name="qValueThreshold">A run enters when its run and global q-values are both at or below this; above 0, at most 1.</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="qValueThreshold"/> is not above 0 and at most 1.</exception>
    /// <exception cref="ArgumentException">See the remarks on the class.</exception>
    public static ObservationTable Build(IEnumerable<ObservationRun> runs, IEnumerable<DiaPrecursorRun> precursors,
        IEnumerable<ProteinGroupInfo> proteinGroups, double qValueThreshold = DefaultQValueThreshold)
    {
        ArgumentNullException.ThrowIfNull(runs);
        ArgumentNullException.ThrowIfNull(precursors);
        ArgumentNullException.ThrowIfNull(proteinGroups);
        if (!(qValueThreshold > 0 && qValueThreshold <= 1))
            throw new ArgumentOutOfRangeException(nameof(qValueThreshold), qValueThreshold,
                "A q-value threshold is above 0 and at most 1.");

        var runList = runs.ToList();
        if (runList.Contains(null!))
            throw new ArgumentException("A run is null.", nameof(runs));
        RefuseFractions(runList);

        var labels = new SortedSet<string>(StringComparer.Ordinal);
        var targets = new List<DiaPrecursorRun>();
        foreach (DiaPrecursorRun precursor in precursors)
        {
            if (precursor is null)
                throw new ArgumentException("A precursor is null.", nameof(precursors));
            labels.Add(precursor.FileName);
            if (precursor.IsDecoy || precursor.IsEntrapment || precursor.IsContaminant)
                continue;
            CheckQValues(precursor);
            targets.Add(precursor);
        }

        // A design run that no record names was searched and reported nothing; whether that is allowed depends on its sample.
        foreach (ObservationRun run in runList)
            if (!Names(labels, run))
                labels.Add(run.FileName);

        var reporting = new HashSet<string>(StringComparer.Ordinal);
        var table = ObservationTable.Build(runList, labels, ToRows(targets, qValueThreshold, reporting), proteinGroups,
            QuantBasis.MsmsOnly);
        RefuseARunLostInsideItsSample(runList, reporting);
        return table;
    }

    /// <summary>
    /// One row per precursor (full sequence and charge), with its value in each run where it passes both q-values. This is
    /// the one place precursors become rows, so a different charge rule changes only this method.
    /// </summary>
    /// <param name="reporting">Filled with the file of every run where at least one precursor passes.</param>
    private static List<PeptideRunValues> ToRows(List<DiaPrecursorRun> targets, double qValueThreshold, HashSet<string> reporting)
    {
        var rows = new List<PeptideRunValues>();
        var ordered = targets.OrderBy(p => p.FullSequence, StringComparer.Ordinal).ThenBy(p => p.Charge)
            .ThenBy(p => p.FileName, StringComparer.Ordinal).ToList();

        for (int start = 0; start < ordered.Count;)
        {
            DiaPrecursorRun first = ordered[start];
            string[] groups = SortedGroups(first);
            var byRun = new Dictionary<string, RunValue>(StringComparer.Ordinal);

            int end = start;
            for (; end < ordered.Count && ordered[end].Charge == first.Charge
                   && string.Equals(ordered[end].FullSequence, first.FullSequence, StringComparison.Ordinal); end++)
            {
                DiaPrecursorRun precursor = ordered[end];
                if (end > start && string.Equals(precursor.FileName, ordered[end - 1].FileName, StringComparison.Ordinal))
                    throw new ArgumentException($"Precursor '{precursor.FullSequence}' at charge {precursor.Charge} is listed " +
                        $"twice in run '{precursor.FileName}'.", "precursors");

                string[] theseGroups = SortedGroups(precursor);
                if (!string.Equals(precursor.BaseSequence, first.BaseSequence, StringComparison.Ordinal)
                    || !theseGroups.SequenceEqual(groups, StringComparer.Ordinal))
                    throw new ArgumentException($"Precursor '{first.FullSequence}' at charge {first.Charge} has base sequence " +
                        $"'{first.BaseSequence}' and groups '{string.Join(';', groups)}' in run '{first.FileName}' but " +
                        $"'{precursor.BaseSequence}' and '{string.Join(';', theseGroups)}' in run '{precursor.FileName}'.",
                        "precursors");

                RunValue value = ToRunValue(precursor);
                if (precursor.RunQValue <= qValueThreshold && precursor.GlobalQValue <= qValueThreshold)
                {
                    byRun[precursor.FileName] = value;
                    reporting.Add(precursor.FileName);
                }
            }

            if (byRun.Count > 0)
                rows.Add(new PeptideRunValues(first.FullSequence, first.BaseSequence, groups, byRun, first.Charge));
            start = end;
        }
        return rows;
    }

    /// <summary>One run's value from the precursor's kind (the mapping in the class remarks).</summary>
    private static RunValue ToRunValue(DiaPrecursorRun precursor)
    {
        switch (precursor.Kind)
        {
            case PrecursorQuantityKind.Value:
                if (double.IsFinite(precursor.Quantity) && precursor.Quantity > 0)
                    return new RunValue(precursor.Quantity, ObservationState.Quantified);
                break;
            case PrecursorQuantityKind.BelowFloor:
            case PrecursorQuantityKind.NotLooked:
                if (double.IsNaN(precursor.Quantity))
                    return new RunValue(double.NaN, ObservationState.IdentifiedNotQuantified);
                break;
            default:
                throw new ArgumentException($"Precursor '{precursor.FullSequence}' at charge {precursor.Charge} in run " +
                    $"'{precursor.FileName}' has quantity kind {precursor.Kind}, which this adapter does not map.", "precursors");
        }

        throw new ArgumentException($"Precursor '{precursor.FullSequence}' at charge {precursor.Charge} in run " +
            $"'{precursor.FileName}' has quantity {Format(precursor.Quantity)} with kind {precursor.Kind}; a Value is finite " +
            "and above 0, and the other kinds carry NaN, never a number.", "precursors");
    }

    private static void CheckQValues(DiaPrecursorRun precursor)
    {
        if (precursor.RunQValue is >= 0 and <= 1 && precursor.GlobalQValue is >= 0 and <= 1)
            return;
        throw new ArgumentException($"Precursor '{precursor.FullSequence}' at charge {precursor.Charge} in run " +
            $"'{precursor.FileName}' has run q-value {Format(precursor.RunQValue)} and global q-value " +
            $"{Format(precursor.GlobalQValue)}; a q-value is from 0 to 1.", "precursors");
    }

    /// <summary>
    /// Refuses a design run where no precursor passes when another run of its sample has one: the replicate average would
    /// leave it out with no record (QuantProject 086). A sample where no run reports anything is kept, without values.
    /// </summary>
    private static void RefuseARunLostInsideItsSample(List<ObservationRun> runs, HashSet<string> reporting)
    {
        foreach (var sample in runs.GroupBy(r => r.SampleId, StringComparer.Ordinal)
                     .OrderBy(g => g.Key, StringComparer.Ordinal))
        {
            var lost = sample.Where(r => !Names(reporting, r)).Select(r => r.FileName).Order(StringComparer.Ordinal).ToList();
            if (lost.Count == 0 || lost.Count == sample.Count())
                continue;

            var others = sample.Where(r => Names(reporting, r)).Select(r => $"'{r.FileName}'").Order(StringComparer.Ordinal);
            throw new ArgumentException($"Run '{lost[0]}' of sample '{sample.Key}' reports no precursor that passes, while " +
                $"the sample's run(s) {string.Join(", ", others)} do, so averaging the sample's replicates would leave it " +
                "out with no record. Remove it from the experimental design to analyse the sample without it.", "runs");
        }
    }

    /// <summary>Whether <paramref name="files"/> names <paramref name="run"/>: by its name, or its name plus <c>-calib</c>.</summary>
    private static bool Names(ISet<string> files, ObservationRun run) =>
        files.Contains(run.FileName) || files.Contains(run.FileName + CalibratedSuffix);

    /// <summary>Refuses fractions: two runs of one sample and technical replicate at different fractions.</summary>
    private static void RefuseFractions(List<ObservationRun> runs)
    {
        var split = runs.GroupBy(r => (r.SampleId, r.TechnicalReplicate))
            .FirstOrDefault(g => g.Select(r => r.Fraction).Distinct().Count() > 1);
        if (split is null)
            return;

        var files = split.Select(r => r.FileName).Order(StringComparer.Ordinal).Select(f => $"'{f}'");
        throw new ArgumentException($"Runs {string.Join(" and ", files)} are different fractions of sample " +
            $"'{split.Key.SampleId}', technical replicate {split.Key.TechnicalReplicate}; DIA fractions are not supported, " +
            "so give each run its own sample or technical replicate.", nameof(runs));
    }

    private static string[] SortedGroups(DiaPrecursorRun precursor) =>
        (precursor.ProteinGroups ?? Array.Empty<string>()).Distinct(StringComparer.Ordinal).Order(StringComparer.Ordinal).ToArray();

    private static string Format(double value) => value.ToString(CultureInfo.InvariantCulture);
}
