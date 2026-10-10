using System;
using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;
using Omics.BioPolymerGroup;

namespace Quantification.Differential;

/// <summary>Which match-between-runs values enter an analysis (ST-15; <c>DEF-DIFF-COLUMNS</c> <c>quant_basis</c>).</summary>
public enum QuantBasis
{
    /// <summary><c>msms_only</c>: no match-between-runs value enters.</summary>
    MsmsOnly,

    /// <summary><c>mbr_kept</c>: the match-between-runs values the producer kept (<c>DEF-MBR-KEPT</c>) enter.</summary>
    MbrKept,
}

/// <summary>The machine names of <see cref="QuantBasis"/>, as <c>DEF-DIFF-COLUMNS</c> writes them.</summary>
public static class QuantBasisNames
{
    /// <summary><c>msms_only</c> or <c>mbr_kept</c>.</summary>
    public static string ToMachineName(this QuantBasis basis) => basis switch
    {
        QuantBasis.MsmsOnly => "msms_only",
        QuantBasis.MbrKept => "mbr_kept",
        _ => throw new ArgumentOutOfRangeException(nameof(basis), basis, null),
    };
}

/// <summary>
/// What a peptide's value in one run, or in one sample, is (STATS-FRAMEWORK section 3.1; <c>DEF-PEP-DT</c>). A cell
/// without a value still says why: states are never collapsed. Declared from most to least informative; a sample with
/// no value takes the most informative state among its runs.
/// </summary>
public enum ObservationState
{
    /// <summary>A value from MS/MS-identified peaks only.</summary>
    Quantified,

    /// <summary>
    /// The value includes a match-between-runs peak. Under <see cref="QuantBasis.MsmsOnly"/> such a value is excluded, so
    /// a cell can carry this state with no value.
    /// </summary>
    MbrTransferred,

    /// <summary>
    /// A peak was found but withheld as ambiguous: shared with another peptide, or the sample's strongest fraction was.
    /// No value.
    /// </summary>
    AmbiguousPeak,

    /// <summary>Identified by MS/MS, but no peak could be quantified. No value.</summary>
    IdentifiedNotQuantified,

    /// <summary>Neither identified by MS/MS nor transferred. No value.</summary>
    NotDetected,
}

/// <summary>What a producer wrote for one peptide in one run.</summary>
/// <param name="Intensity">
/// The intensity on the linear scale. Read only when <paramref name="State"/> is <see cref="ObservationState.Quantified"/>
/// or <see cref="ObservationState.MbrTransferred"/>, and then it must be finite and greater than 0.
/// </param>
/// <param name="State">Where the value came from, or why there is none.</param>
public readonly record struct RunValue(double Intensity, ObservationState State);

/// <summary>
/// One spectra file and its place in the experimental design: the biological sample it belongs to, its fraction and its
/// technical replicate. Grouping comes only from the design; a file name is never parsed for it.
/// </summary>
/// <param name="FileName">
/// The file's name without its extension, as FlashLFQ labels its columns (<c>Intensity_&lt;FileName&gt;</c>). A column
/// labelled with this name plus <c>-calib</c> also matches it, because MetaMorpheus quantifies the calibrated files.
/// </param>
/// <param name="SampleId">The biological sample: the key the design layer gives it (<c>AnalysisSample.SampleId</c>).</param>
/// <param name="Fraction">The fraction; fractions of one technical replicate are summed before the log.</param>
/// <param name="TechnicalReplicate">The technical replicate; replicates of one sample are averaged on the log scale.</param>
public sealed record ObservationRun(string FileName, string SampleId, int Fraction = 0, int TechnicalReplicate = 0)
{
    /// <summary>
    /// The runs of a label-free experimental design: each file's sample is <c>{condition}_{biological replicate + 1}</c>,
    /// the label MetaMorpheus's quantified tables use (<see cref="SampleGroupLabels.ForSample"/>).
    /// </summary>
    public static IReadOnlyList<ObservationRun> FromSpectraFiles(IEnumerable<SpectraFileInfo> files)
    {
        ArgumentNullException.ThrowIfNull(files);
        return files.Select(f => new ObservationRun(f.FilenameWithoutExtension, SampleGroupLabels.ForSample(f),
            f.Fraction, f.TechnicalReplicate)).ToList();
    }
}

/// <summary>
/// A MetaMorpheus protein group as its protein-group table describes it (GR-23): its name, q-value, and whether it is a
/// decoy, contaminant or entrapment group. FlashLFQ's own protein groups carry none of these.
/// </summary>
/// <param name="Name">
/// The group's name exactly as MetaMorpheus writes it: the member accessions joined by <c>|</c>. An opaque key, never
/// re-sorted, because the writer sorts with the current culture.
/// </param>
/// <param name="QValue">The group's <c>Protein QValue</c>.</param>
/// <param name="IsDecoy">The group is a decoy.</param>
/// <param name="IsContaminant">The group is a contaminant.</param>
/// <param name="IsEntrapment">The group comes from an entrapment database.</param>
/// <param name="Genes">The <c>Gene</c> cell as written; null when the table has none.</param>
public sealed record ProteinGroupInfo(string Name, double QValue, bool IsDecoy, bool IsContaminant, bool IsEntrapment = false,
    string? Genes = null)
{
    /// <summary>Every member accession: the name split on <c>|</c>, ordinal-sorted (<c>DEF-DIFF-COLUMNS</c> <c>feature_accessions</c>).</summary>
    public IReadOnlyList<string> Accessions =>
        Name.Split('|', StringSplitOptions.RemoveEmptyEntries).Order(StringComparer.Ordinal).ToList();
}

/// <summary>A peptide (one modified sequence) and the MetaMorpheus protein groups it maps to.</summary>
/// <param name="FullSequence">The engine's full sequence; the peptide's key.</param>
/// <param name="BaseSequence">The unmodified sequence.</param>
/// <param name="ProteinGroups">Every group the peptide's <c>Protein Groups</c> cell names: distinct, ordinal-sorted.</param>
public sealed record ObservedPeptide(string FullSequence, string BaseSequence, IReadOnlyList<string> ProteinGroups)
{
    /// <summary>The placeholder group MetaMorpheus gives a quantified peptide that is in no surviving group.</summary>
    public const string UndefinedProteinGroup = "UNDEFINED";

    /// <summary>
    /// The group this peptide is unique to, or null (GR-23): unique means the cell names exactly one group and it is not
    /// <see cref="UndefinedProteinGroup"/>. Sharing counts against every group MetaMorpheus wrote, whatever its q-value.
    /// </summary>
    public string? UniqueProteinGroup =>
        ProteinGroups.Count == 1 && !string.Equals(ProteinGroups[0], UndefinedProteinGroup, StringComparison.Ordinal)
            ? ProteinGroups[0]
            : null;
}

/// <summary>One peptide's values in each run, as a producer wrote them.</summary>
/// <param name="FullSequence">The engine's full sequence.</param>
/// <param name="BaseSequence">The unmodified sequence.</param>
/// <param name="ProteinGroups">The group names the producer lists for the peptide, in any order.</param>
/// <param name="ByRun">The value per run, keyed by the producer's run label; a run with no entry is <see cref="ObservationState.NotDetected"/>.</param>
public sealed record PeptideRunValues(string FullSequence, string BaseSequence, IReadOnlyList<string> ProteinGroups,
    IReadOnlyDictionary<string, RunValue> ByRun);

/// <summary>One cell of an <see cref="ObservationTable"/>: a peptide in a sample.</summary>
/// <param name="Peptide">The peptide.</param>
/// <param name="SampleId">The sample.</param>
/// <param name="Log2Intensity">log2 of the sample's intensity, or NaN when there is no value. Never 0 for "missing".</param>
/// <param name="State">Where the value came from, or why there is none.</param>
public readonly record struct Observation(ObservedPeptide Peptide, string SampleId, double Log2Intensity, ObservationState State);

/// <summary>
/// The input to differential statistics for one quantification basis (STAT1 milestone M5; STATS-FRAMEWORK section 3.1):
/// every peptide in every biological sample, as log2 intensity or NaN, with a state.
/// </summary>
/// <remarks>
/// <para>
/// <b>Samples, not runs.</b> Fractions of one technical replicate are summed on the linear scale, then the replicates of
/// a sample are averaged on the log2 scale. A sample whose runs gave no value at all is kept, and listed in
/// <see cref="SamplesWithoutValues"/>.
/// </para>
/// <para>
/// <b>Bases.</b> Under <see cref="QuantBasis.MsmsOnly"/> a match-between-runs value never enters a sum. Under
/// <see cref="QuantBasis.MbrKept"/> it does, and the cell's state is <see cref="ObservationState.MbrTransferred"/>.
/// </para>
/// <para>
/// <b>Protein groups</b> are MetaMorpheus's parsimony groups (GR-23), named as the producer names them; their q-values
/// and flags come from <see cref="ProteinGroupInfo"/>.
/// </para>
/// <para>
/// <b>Refusals.</b> Input that cannot be read consistently throws <see cref="ArgumentException"/> naming the run or
/// peptide: two runs with one file name, a producer run that matches no design run or the reverse, a run label claimed
/// twice, a peptide listed twice, a value that is not finite and positive, a group name containing <c>;</c>.
/// </para>
/// <para>
/// <b>Determinism.</b> Peptides are ordered by full sequence and samples by id (ordinal); sums run in fraction order, so
/// the table does not depend on input order.
/// </para>
/// </remarks>
public sealed class ObservationTable
{
    /// <summary>The suffix MetaMorpheus gives a calibrated file's name.</summary>
    private const string CalibratedSuffix = "-calib";

    // [peptide][sample], aligned with Peptides and Samples.
    private readonly double[][] _log2;
    private readonly ObservationState[][] _state;

    private ObservationTable(QuantBasis basis, IReadOnlyList<ObservationRun> runs, IReadOnlyList<string> samples,
        IReadOnlyList<ObservedPeptide> peptides, IReadOnlyDictionary<string, ProteinGroupInfo> proteinGroups,
        double[][] log2, ObservationState[][] state, IReadOnlyList<string> warnings)
    {
        Basis = basis;
        Runs = runs;
        Samples = samples;
        Peptides = peptides;
        ProteinGroups = proteinGroups;
        _log2 = log2;
        _state = state;
        Warnings = warnings;
        SamplesWithoutValues = Enumerable.Range(0, samples.Count)
            .Where(s => log2.All(row => double.IsNaN(row[s])))
            .Select(s => samples[s])
            .ToList();
    }

    /// <summary>The basis this table was built under.</summary>
    public QuantBasis Basis { get; }

    /// <summary>The runs, ordered by file name (ordinal).</summary>
    public IReadOnlyList<ObservationRun> Runs { get; }

    /// <summary>The biological samples, ordered by id (ordinal); every sample of the design, with or without values.</summary>
    public IReadOnlyList<string> Samples { get; }

    /// <summary>The peptides, ordered by full sequence (ordinal).</summary>
    public IReadOnlyList<ObservedPeptide> Peptides { get; }

    /// <summary>The protein groups the producer described, by name.</summary>
    public IReadOnlyDictionary<string, ProteinGroupInfo> ProteinGroups { get; }

    /// <summary>The samples in which no peptide has a value under this basis.</summary>
    public IReadOnlyList<string> SamplesWithoutValues { get; }

    /// <summary>Things a reader should know that did not stop the build (e.g. groups the protein table does not describe).</summary>
    public IReadOnlyList<string> Warnings { get; }

    /// <summary>The log2 intensity of peptide <paramref name="peptide"/> in sample <paramref name="sample"/>, or NaN.</summary>
    public double Log2Intensity(int peptide, int sample) => _log2[peptide][sample];

    /// <summary>The state of peptide <paramref name="peptide"/> in sample <paramref name="sample"/>.</summary>
    public ObservationState State(int peptide, int sample) => _state[peptide][sample];

    /// <summary>Every cell, peptide by peptide, each peptide's samples in order: the long form of the table.</summary>
    public IEnumerable<Observation> Observations()
    {
        for (int p = 0; p < Peptides.Count; p++)
        for (int s = 0; s < Samples.Count; s++)
            yield return new Observation(Peptides[p], Samples[s], _log2[p][s], _state[p][s]);
    }

    /// <summary>The same table over a subset of its samples (one stratum, GR-5), with every peptide kept.</summary>
    /// <exception cref="ArgumentException">A sample the table does not have.</exception>
    public ObservationTable ForSamples(IEnumerable<string> sampleIds)
    {
        ArgumentNullException.ThrowIfNull(sampleIds);
        var keep = sampleIds.Distinct(StringComparer.Ordinal).Order(StringComparer.Ordinal).ToList();
        var index = new List<int>(keep.Count);
        foreach (string id in keep)
        {
            int s = IndexOf(Samples, id);
            if (s < 0)
                throw new ArgumentException($"The table has no sample '{id}'.", nameof(sampleIds));
            index.Add(s);
        }

        var keepSet = keep.ToHashSet(StringComparer.Ordinal);
        return new ObservationTable(Basis, Runs.Where(r => keepSet.Contains(r.SampleId)).ToList(), keep, Peptides,
            ProteinGroups, _log2.Select(row => index.Select(s => row[s]).ToArray()).ToArray(),
            _state.Select(row => index.Select(s => row[s]).ToArray()).ToArray(), Warnings);
    }

    /// <summary>The same table with each cell's log2 value replaced; states, peptides and samples unchanged.</summary>
    internal ObservationTable WithLog2(Func<int, int, double> log2) =>
        new(Basis, Runs, Samples, Peptides, ProteinGroups,
            Enumerable.Range(0, Peptides.Count)
                .Select(p => Enumerable.Range(0, Samples.Count).Select(s => log2(p, s)).ToArray()).ToArray(),
            _state, Warnings);

    /// <summary>
    /// Builds the table from a producer's per-run values.
    /// </summary>
    /// <param name="runs">The experimental design's runs.</param>
    /// <param name="sourceRunLabels">Every run label the producer wrote (its column labels), so a mismatch is refused.</param>
    /// <param name="peptides">Each peptide's values per run.</param>
    /// <param name="proteinGroups">The producer's protein groups.</param>
    /// <param name="basis">Which match-between-runs values enter.</param>
    /// <exception cref="ArgumentException">Input that cannot be read consistently; see the remarks on the class.</exception>
    public static ObservationTable Build(IEnumerable<ObservationRun> runs, IEnumerable<string> sourceRunLabels,
        IEnumerable<PeptideRunValues> peptides, IEnumerable<ProteinGroupInfo> proteinGroups, QuantBasis basis)
    {
        ArgumentNullException.ThrowIfNull(runs);
        ArgumentNullException.ThrowIfNull(sourceRunLabels);
        ArgumentNullException.ThrowIfNull(peptides);
        ArgumentNullException.ThrowIfNull(proteinGroups);

        var orderedRuns = ValidateRuns(runs);
        var labelToRun = MatchLabels(orderedRuns, sourceRunLabels);
        var groups = ValidateGroups(proteinGroups);

        var samples = orderedRuns.Select(r => r.SampleId).Distinct(StringComparer.Ordinal).Order(StringComparer.Ordinal).ToList();

        // Each sample's runs, replicate by replicate, fraction by fraction: the order every sum runs in.
        var sampleRuns = samples.Select(sample => orderedRuns
                .Where(r => string.Equals(r.SampleId, sample, StringComparison.Ordinal))
                .GroupBy(r => r.TechnicalReplicate)
                .OrderBy(g => g.Key)
                .Select(g => g.OrderBy(r => r.Fraction).Select(r => r.FileName).ToArray())
                .ToArray())
            .ToArray();

        var ordered = peptides.OrderBy(p => p.FullSequence, StringComparer.Ordinal).ToList();
        var observed = new List<ObservedPeptide>(ordered.Count);
        var log2 = new double[ordered.Count][];
        var state = new ObservationState[ordered.Count][];
        var undescribed = new SortedSet<string>(StringComparer.Ordinal);

        for (int p = 0; p < ordered.Count; p++)
        {
            PeptideRunValues peptide = ordered[p];
            if (p > 0 && string.Equals(peptide.FullSequence, ordered[p - 1].FullSequence, StringComparison.Ordinal))
                throw new ArgumentException($"Peptide '{peptide.FullSequence}' is listed twice.", nameof(peptides));

            var groupNames = ValidateGroupNames(peptide);
            foreach (string name in groupNames)
                if (!groups.ContainsKey(name) && !string.Equals(name, ObservedPeptide.UndefinedProteinGroup, StringComparison.Ordinal))
                    undescribed.Add(name);
            observed.Add(new ObservedPeptide(peptide.FullSequence, peptide.BaseSequence, groupNames));

            var byRun = RunValuesByFile(peptide, labelToRun);
            log2[p] = new double[samples.Count];
            state[p] = new ObservationState[samples.Count];
            for (int s = 0; s < samples.Count; s++)
                (log2[p][s], state[p][s]) = Combine(sampleRuns[s], byRun, basis);
        }

        var warnings = new List<string>();
        if (undescribed.Count > 0)
            warnings.Add($"{undescribed.Count} protein group(s) named by peptides are not in the protein-group table, so " +
                $"their q-value and decoy/contaminant status are unknown: {string.Join(", ", undescribed.Take(5))}" +
                (undescribed.Count > 5 ? ", ..." : "") + ".");

        return new ObservationTable(basis, orderedRuns, samples, observed, groups, log2, state, warnings);
    }

    /// <summary>
    /// One sample's cell: fractions summed per replicate, then the replicates' log2 values averaged. Under
    /// <see cref="QuantBasis.MsmsOnly"/> a transferred value is left out of every sum but still names the state.
    /// </summary>
    private static (double Log2, ObservationState State) Combine(string[][] replicates,
        IReadOnlyDictionary<string, RunValue> byRun, QuantBasis basis)
    {
        double log2Sum = 0;
        int replicatesWithValue = 0;
        bool anyTransfer = false;
        var best = ObservationState.NotDetected;

        foreach (string[] fractions in replicates)
        {
            double sum = 0;
            bool hasValue = false;
            foreach (string file in fractions)
            {
                RunValue value = byRun.TryGetValue(file, out var v) ? v : new RunValue(0, ObservationState.NotDetected);
                if (value.State < best)
                    best = value.State;

                bool enters = value.State == ObservationState.Quantified
                    || (value.State == ObservationState.MbrTransferred && basis == QuantBasis.MbrKept);
                if (!enters)
                    continue;

                sum += value.Intensity;
                hasValue = true;
                anyTransfer |= value.State == ObservationState.MbrTransferred;
            }

            if (hasValue)
            {
                log2Sum += Math.Log2(sum);
                replicatesWithValue++;
            }
        }

        if (replicatesWithValue == 0)
            return (double.NaN, best);

        return (log2Sum / replicatesWithValue, anyTransfer ? ObservationState.MbrTransferred : ObservationState.Quantified);
    }

    /// <summary>The peptide's values keyed by design file name, each checked.</summary>
    private static Dictionary<string, RunValue> RunValuesByFile(PeptideRunValues peptide,
        IReadOnlyDictionary<string, string> labelToRun)
    {
        var byFile = new Dictionary<string, RunValue>(StringComparer.Ordinal);
        foreach (var (label, value) in peptide.ByRun)
        {
            if (!labelToRun.TryGetValue(label, out string? file))
                throw new ArgumentException(
                    $"Peptide '{peptide.FullSequence}' has a value for run '{label}', which the producer did not list.",
                    "peptides");

            bool carriesValue = value.State is ObservationState.Quantified or ObservationState.MbrTransferred;
            if (carriesValue && !(double.IsFinite(value.Intensity) && value.Intensity > 0))
                throw new ArgumentException(
                    $"Peptide '{peptide.FullSequence}' in run '{label}' has intensity {value.Intensity} with state " +
                    $"{value.State}; a value must be finite and greater than 0.", "peptides");

            byFile[file] = value;
        }
        return byFile;
    }

    private static List<string> ValidateGroupNames(PeptideRunValues peptide)
    {
        foreach (string name in peptide.ProteinGroups)
            if (name.Contains(';'))
                throw new ArgumentException(
                    $"Peptide '{peptide.FullSequence}' names protein group '{name}'; a group name cannot contain ';', which " +
                    "separates groups in the peptide table.", "peptides");

        return peptide.ProteinGroups.Distinct(StringComparer.Ordinal).Order(StringComparer.Ordinal).ToList();
    }

    private static List<ObservationRun> ValidateRuns(IEnumerable<ObservationRun> runs)
    {
        var ordered = runs.OrderBy(r => r.FileName, StringComparer.Ordinal).ToList();
        for (int i = 1; i < ordered.Count; i++)
            if (string.Equals(ordered[i].FileName, ordered[i - 1].FileName, StringComparison.Ordinal))
                throw new ArgumentException($"Two runs have the file name '{ordered[i].FileName}'.", nameof(runs));

        var place = ordered.GroupBy(r => (r.SampleId, r.TechnicalReplicate, r.Fraction)).FirstOrDefault(g => g.Count() > 1);
        if (place != null)
            throw new ArgumentException(
                $"Runs {string.Join(" and ", place.Select(r => $"'{r.FileName}'"))} are both sample '{place.Key.SampleId}', " +
                $"technical replicate {place.Key.TechnicalReplicate}, fraction {place.Key.Fraction}.", nameof(runs));

        return ordered;
    }

    /// <summary>
    /// Matches each producer label to a design file: by name, or by name once a trailing <c>-calib</c> is removed. Every
    /// label must match one file and every file exactly one label.
    /// </summary>
    private static Dictionary<string, string> MatchLabels(IReadOnlyList<ObservationRun> runs, IEnumerable<string> labels)
    {
        var files = runs.Select(r => r.FileName).ToHashSet(StringComparer.Ordinal);
        var labelToRun = new Dictionary<string, string>(StringComparer.Ordinal);
        var runToLabel = new Dictionary<string, string>(StringComparer.Ordinal);

        foreach (string label in labels)
        {
            string? file = files.Contains(label) ? label
                : label.EndsWith(CalibratedSuffix, StringComparison.Ordinal) && files.Contains(label[..^CalibratedSuffix.Length])
                    ? label[..^CalibratedSuffix.Length]
                    : null;
            if (file == null)
                throw new ArgumentException($"The producer wrote run '{label}', which is not in the experimental design.",
                    nameof(labels));
            if (runToLabel.TryGetValue(file, out string? other))
                throw new ArgumentException($"The producer wrote both '{other}' and '{label}', which both match design run '{file}'.",
                    nameof(labels));
            runToLabel[file] = label;
            labelToRun[label] = file;
        }

        var unwritten = runs.Where(r => !runToLabel.ContainsKey(r.FileName)).Select(r => $"'{r.FileName}'").ToList();
        if (unwritten.Count > 0)
            throw new ArgumentException(
                $"The experimental design has run(s) the producer did not write: {string.Join(", ", unwritten)}.", nameof(labels));

        return labelToRun;
    }

    private static Dictionary<string, ProteinGroupInfo> ValidateGroups(IEnumerable<ProteinGroupInfo> proteinGroups)
    {
        var groups = new Dictionary<string, ProteinGroupInfo>(StringComparer.Ordinal);
        foreach (ProteinGroupInfo group in proteinGroups)
        {
            if (groups.TryGetValue(group.Name, out var existing) && existing != group)
                throw new ArgumentException(
                    $"Protein group '{group.Name}' is described twice, differently: {existing} and {group}.", nameof(proteinGroups));
            groups[group.Name] = group;
        }
        return groups;
    }

    private static int IndexOf(IReadOnlyList<string> list, string value)
    {
        for (int i = 0; i < list.Count; i++)
            if (string.Equals(list[i], value, StringComparison.Ordinal))
                return i;
        return -1;
    }
}
