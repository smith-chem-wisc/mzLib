using System;
using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;

namespace Quantification.Differential;

/// <summary>
/// One peptide-spectrum match's reporter-ion intensities as a search engine wrote them: the input to
/// <see cref="IsobaricObservationTable.Build"/>. It names no engine's format; an adapter maps an engine's PSM table to it,
/// after the engine's own filters (q-value, decoys).
/// </summary>
/// <param name="FileName">
/// The spectra file's name without its extension, as the engine writes it. A name that ends in <c>-calib</c> also matches
/// the design file without it, because MetaMorpheus quantifies the calibrated files.
/// </param>
/// <param name="ScanNumber">The MS2 scan. With the file it names the match in messages, and it breaks ties.</param>
/// <param name="FullSequence">The engine's full sequence.</param>
/// <param name="BaseSequence">The unmodified sequence.</param>
/// <param name="Charge">The precursor charge; at least 1.</param>
/// <param name="ProteinGroups">Every protein group the peptide's <c>Protein Groups</c> cell names (GR-23), in any order.</param>
/// <param name="ReporterIntensities">
/// One intensity per entry of the file in the experimental design, in the design's order: the order of
/// <see cref="IExperimentalDesign.FileNameSampleInfoDictionary"/>, which is the order of a match's intensities. 0 means
/// the reporter was not measured; a value must be finite and not negative.
/// </param>
public sealed record IsobaricPsm(string FileName, int ScanNumber, string FullSequence, string BaseSequence, int Charge,
    IReadOnlyList<string> ProteinGroups, IReadOnlyList<double> ReporterIntensities);

/// <summary>One row of an <see cref="IsobaricObservationTable"/>: a peptide ion in one plex run.</summary>
/// <param name="Peptide">The peptide and the protein groups it maps to.</param>
/// <param name="Charge">The ion's charge.</param>
/// <param name="Plex">The plex.</param>
/// <param name="TechnicalReplicate">The plex's technical replicate (its injection). Fractions are not a level of the table.</param>
/// <param name="FileName">The file of the match whose reporter intensities the row holds, as the match named it.</param>
/// <param name="ScanNumber">That match's scan.</param>
public sealed record IsobaricIonRow(ObservedPeptide Peptide, int Charge, int Plex, int TechnicalReplicate, string FileName,
    int ScanNumber);

/// <summary>One cell of an <see cref="IsobaricObservationTable"/>: a peptide ion in one plex run, in one channel's sample.</summary>
/// <param name="Row">The ion and its plex run.</param>
/// <param name="SampleId">The sample (<see cref="AnalysisSample.IdFor"/>).</param>
/// <param name="Log2Intensity">log2 of the reporter intensity, or NaN when the reporter was not measured. Never 0 for "missing".</param>
/// <param name="State">
/// <see cref="ObservationState.Quantified"/> with a value; <see cref="ObservationState.IdentifiedNotQuantified"/> without.
/// </param>
public readonly record struct IsobaricObservation(IsobaricIonRow Row, string SampleId, double Log2Intensity, ObservationState State);

/// <summary>
/// The input to differential statistics for an isobaric (TMT, iTRAQ) dataset (STAT1; REQ-QP-15 PLAN 6f): one row per
/// peptide ion (full sequence and charge) per plex run, with the log2 reporter intensity of every sample channel of that
/// plex, or NaN.
/// </summary>
/// <remarks>
/// <para>
/// <b>One value per ion per plex run</b> (MSstatsTMT 2.20's <c>summaryforMultipleRows = sum</c>): of an ion's matches in
/// one plex run, over all its fraction files, the one whose reporter intensities sum highest gives every channel's value,
/// so the ion's values come from its most intense fraction. A tie goes to the earlier file name (ordinal), then the
/// lower scan. A plex run is a plex and technical replicate; its fractions are its files.
/// </para>
/// <para>
/// <b>Few channels.</b> An ion whose chosen match has a reporter measured in at most 2 of its file's channels is dropped
/// (MSstatsTMT's <c>rmPSM_withfewMea_withinRun = TRUE</c>), counting every channel of the file: reference and unused
/// channels too. A match with no reporter at all goes the same way. <see cref="IonsWithFewChannels"/> counts them.
/// </para>
/// <para>
/// <b>Samples</b> are the design's (<see cref="AnalysisSample.FromExperimentalDesign"/>): one per sample channel of a plex,
/// keyed by <see cref="AnalysisSample.IdFor"/>, ordered by id. A pooled reference channel and a channel the plex does not
/// use hold no sample and no value here. A sample belongs to one plex, so a row has values only for the samples of its
/// plex. A plex with no sample channel (only pooled reference or unused channels, as a QC plex has) has no rows: its ions
/// are left out before the few-channel rule, and a warning names the plex and how many ions it held.
/// </para>
/// <para>
/// <b>Missing.</b> A reporter written as 0 was not measured: its cell is NaN with state
/// <see cref="ObservationState.IdentifiedNotQuantified"/>. Nothing is imputed.
/// </para>
/// <para>
/// <b>Refusals.</b> Input that cannot be read consistently throws <see cref="ArgumentException"/> naming the file or
/// match: a design entry that is not an isobaric channel, two files with one name, one file's channels in different
/// plexes, technical replicates or fractions, two files at one plex, technical replicate and fraction, a match whose file
/// is not in the design (or that names a design file two ways), whose intensity count differs from its file's channel
/// count, whose intensity is negative or not finite, whose charge is below 1, or whose peptide another match gives
/// different protein groups or a different base sequence, and a group name containing <c>;</c>.
/// </para>
/// <para>
/// <b>Determinism.</b> Rows are ordered by full sequence (ordinal), charge, plex and technical replicate, so the table does
/// not depend on input order.
/// </para>
/// </remarks>
public sealed class IsobaricObservationTable
{
    /// <summary>The suffix MetaMorpheus gives a calibrated file's name.</summary>
    private const string CalibratedSuffix = "-calib";

    // [row][position of the sample within its plex's samples]
    private readonly double[][] _log2;
    private readonly int[] _samplePlex;
    private readonly int[] _positionInPlex;
    private readonly IReadOnlyDictionary<int, int[]> _plexSamples;

    private IsobaricObservationTable(IReadOnlyList<string> samples, int[] samplePlex, IReadOnlyList<IsobaricIonRow> rows,
        double[][] log2, IReadOnlyDictionary<string, ProteinGroupInfo> proteinGroups, int ionsWithFewChannels,
        IReadOnlyList<string> warnings)
    {
        Samples = samples;
        _samplePlex = samplePlex;
        Rows = rows;
        _log2 = log2;
        ProteinGroups = proteinGroups;
        IonsWithFewChannels = ionsWithFewChannels;
        Warnings = warnings;

        var plexSamples = new Dictionary<int, int[]>();
        _positionInPlex = new int[samples.Count];
        foreach (var plex in Enumerable.Range(0, samples.Count).GroupBy(s => samplePlex[s]))
        {
            int[] members = plex.ToArray();
            plexSamples[plex.Key] = members;
            for (int k = 0; k < members.Length; k++)
                _positionInPlex[members[k]] = k;
        }
        _plexSamples = plexSamples;
    }

    /// <summary>The biological samples, ordered by id (ordinal); every sample of the design, with or without values.</summary>
    public IReadOnlyList<string> Samples { get; }

    /// <summary>The peptide ions in their plex runs, ordered by full sequence (ordinal), charge, plex and technical replicate.</summary>
    public IReadOnlyList<IsobaricIonRow> Rows { get; }

    /// <summary>The protein groups the producer described, by name.</summary>
    public IReadOnlyDictionary<string, ProteinGroupInfo> ProteinGroups { get; }

    /// <summary>
    /// How many ions in plex runs were dropped for a reporter measured in at most 2 channels. Ions of a plex with no sample
    /// channel are not counted here; a warning counts them.
    /// </summary>
    public int IonsWithFewChannels { get; }

    /// <summary>Things a reader should know that did not stop the build (e.g. groups the protein table does not describe).</summary>
    public IReadOnlyList<string> Warnings { get; }

    /// <summary>The plex of sample <paramref name="sample"/>.</summary>
    public int PlexOf(int sample) => _samplePlex[sample];

    /// <summary>The samples of plex <paramref name="plex"/>, as indices into <see cref="Samples"/>, in order; empty for a plex with none.</summary>
    public IReadOnlyList<int> SamplesOfPlex(int plex) => _plexSamples.TryGetValue(plex, out int[]? members) ? members : Array.Empty<int>();

    /// <summary>
    /// The log2 reporter intensity of row <paramref name="row"/> in sample <paramref name="sample"/>, or NaN: not measured,
    /// or the sample is in another plex.
    /// </summary>
    public double Log2Intensity(int row, int sample) =>
        _samplePlex[sample] == Rows[row].Plex ? _log2[row][_positionInPlex[sample]] : double.NaN;

    /// <summary>
    /// The state of row <paramref name="row"/> in sample <paramref name="sample"/>: <see cref="ObservationState.Quantified"/>
    /// with a value, <see cref="ObservationState.IdentifiedNotQuantified"/> without, and
    /// <see cref="ObservationState.NotDetected"/> for a sample of another plex.
    /// </summary>
    public ObservationState State(int row, int sample) =>
        _samplePlex[sample] != Rows[row].Plex ? ObservationState.NotDetected
        : double.IsNaN(_log2[row][_positionInPlex[sample]]) ? ObservationState.IdentifiedNotQuantified
        : ObservationState.Quantified;

    /// <summary>Every cell: row by row, each row's plex samples in order. The long form of the table.</summary>
    public IEnumerable<IsobaricObservation> Observations()
    {
        for (int r = 0; r < Rows.Count; r++)
            foreach (int s in _plexSamples[Rows[r].Plex])
                yield return new IsobaricObservation(Rows[r], Samples[s], Log2Intensity(r, s), State(r, s));
    }

    /// <summary>
    /// The sha256 of <see cref="ToCanonicalText"/> (lower-case hex): the <c>observation_table</c> input the STAT1 runner
    /// records (M5b ruling 10). The same table gives the same hash on any machine and culture; a different value, ion,
    /// plex run, group, group fact, few-channel count or warning gives another.
    /// </summary>
    internal string Sha256() => CanonicalText.Sha256(ToCanonicalText());

    /// <summary>
    /// The table as canonical text: a header (isobaric), the samples in order with their plex (<c>sample:plex</c>), each
    /// described protein group by name with its q-value, flags and genes, each row in order (peptide, charge, plex,
    /// technical replicate, file, scan) with its plex samples' log2 values in order, the few-channel count, then the
    /// warnings (<see cref="CanonicalText"/>). A cell's state follows from its value, so it is not written.
    /// </summary>
    internal string ToCanonicalText()
    {
        var text = new System.Text.StringBuilder();
        text.Append("observation_table\tisobaric\n");
        text.Append("samples");
        for (int s = 0; s < Samples.Count; s++)
            text.Append('\t').Append(CanonicalText.Escape(Samples[s])).Append(':').Append(CanonicalText.Number(_samplePlex[s]));
        text.Append('\n');
        CanonicalText.AppendGroups(text, ProteinGroups);
        for (int r = 0; r < Rows.Count; r++)
        {
            var row = Rows[r];
            text.Append("ion\t");
            CanonicalText.AppendPeptide(text, row.Peptide);
            text.Append('\t').Append(CanonicalText.Number(row.Charge)).Append('\t').Append(CanonicalText.Number(row.Plex))
                .Append('\t').Append(CanonicalText.Number(row.TechnicalReplicate)).Append('\t').Append(CanonicalText.Escape(row.FileName))
                .Append('\t').Append(CanonicalText.Number(row.ScanNumber));
            foreach (double v in _log2[r]) text.Append('\t').Append(CanonicalText.Number(v));
            text.Append('\n');
        }
        text.Append("ions_with_few_channels\t").Append(CanonicalText.Number(IonsWithFewChannels)).Append('\n');
        CanonicalText.AppendWarnings(text, Warnings);
        return text.ToString();
    }

    /// <summary>
    /// The same table with each cell of a row's plex replaced by <paramref name="log2"/>(row, sample); a missing cell must
    /// stay missing, because the state is read from the value.
    /// </summary>
    internal IsobaricObservationTable WithLog2(Func<int, int, double> log2)
    {
        var values = new double[Rows.Count][];
        for (int r = 0; r < Rows.Count; r++)
        {
            int[] members = _plexSamples[Rows[r].Plex];
            values[r] = new double[members.Length];
            for (int k = 0; k < members.Length; k++)
                values[r][k] = double.IsNaN(_log2[r][k]) ? double.NaN : log2(r, members[k]);
        }
        return new IsobaricObservationTable(Samples, _samplePlex, Rows, values, ProteinGroups, IonsWithFewChannels, Warnings);
    }

    /// <summary>One design file: where it sits, and the sample each of its channels holds (-1 for none).</summary>
    private sealed record DesignFile(string Name, int Plex, int TechnicalReplicate, int Fraction, int[] SampleIndex);

    /// <summary>An ion's best match so far in one plex run.</summary>
    private sealed record Candidate(IsobaricPsm Psm, DesignFile File, double Sum);

    /// <summary>Builds the table from a pipeline's isobaric experimental design and its matches. See the remarks on the class.</summary>
    /// <param name="design">The experimental design: every entry an isobaric channel.</param>
    /// <param name="psms">The matches to quantify, already filtered by the engine.</param>
    /// <param name="proteinGroups">The producer's protein groups.</param>
    /// <exception cref="ArgumentException">Input that cannot be read consistently; see the remarks on the class.</exception>
    public static IsobaricObservationTable Build(IExperimentalDesign design, IEnumerable<IsobaricPsm> psms,
        IEnumerable<ProteinGroupInfo> proteinGroups)
    {
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(design.FileNameSampleInfoDictionary, nameof(design));
        ArgumentNullException.ThrowIfNull(psms);
        ArgumentNullException.ThrowIfNull(proteinGroups);

        var groups = ObservationTable.ValidateGroups(proteinGroups);
        var samples = AnalysisSample.FromExperimentalDesign(design).Select(s => s.SampleId).ToList();
        var (files, samplePlex) = ReadDesign(design, samples);

        var best = new Dictionary<(string FullSequence, int Charge, int Plex, int TechnicalReplicate), Candidate>();
        var peptides = new Dictionary<string, ObservedPeptide>(StringComparer.Ordinal);
        var labelOfFile = new Dictionary<string, string>(StringComparer.Ordinal);
        foreach (IsobaricPsm psm in psms)
        {
            if (psm is null)
                throw new ArgumentException("A match is null.", nameof(psms));
            DesignFile file = CheckPsm(psm, files, labelOfFile);
            ObservedPeptide peptide = PeptideOf(psm, peptides);

            double sum = psm.ReporterIntensities.Sum();
            var key = (peptide.FullSequence, psm.Charge, file.Plex, file.TechnicalReplicate);
            if (!best.TryGetValue(key, out var held) || Precedes(psm, sum, held))
                best[key] = new Candidate(psm, file, sum);
        }

        var rows = new List<IsobaricIonRow>();
        var log2 = new List<double[]>();
        int fewChannels = 0;
        var plexPositions = PlexPositions(samplePlex);
        // A plex with no sample channel (a QC plex of pooled reference, say) has nowhere to put its ions' values.
        var plexesWithoutSamples = new SortedSet<int>(files.Values.Select(f => f.Plex).Where(p => !plexPositions.ContainsKey(p)));
        int ionsWithoutSamples = 0;
        foreach (var (key, chosen) in best
                     .OrderBy(b => b.Key.FullSequence, StringComparer.Ordinal)
                     .ThenBy(b => b.Key.Charge).ThenBy(b => b.Key.Plex).ThenBy(b => b.Key.TechnicalReplicate)
                     .Select(b => (b.Key, b.Value)))
        {
            if (!plexPositions.TryGetValue(key.Plex, out var plex))
            {
                ionsWithoutSamples++;
                continue;
            }

            var intensities = chosen.Psm.ReporterIntensities;
            if (intensities.Count(v => v > 0) <= 2)
            {
                fewChannels++;
                continue;
            }

            var (positions, count) = plex;
            var values = Enumerable.Repeat(double.NaN, count).ToArray();
            for (int c = 0; c < intensities.Count; c++)
            {
                int sample = chosen.File.SampleIndex[c];
                if (sample >= 0 && intensities[c] > 0)
                    values[positions[sample]] = Math.Log2(intensities[c]);
            }
            rows.Add(new IsobaricIonRow(peptides[key.FullSequence], key.Charge, key.Plex, key.TechnicalReplicate,
                chosen.Psm.FileName, chosen.Psm.ScanNumber));
            log2.Add(values);
        }

        var warnings = new List<string>();
        if (plexesWithoutSamples.Count > 0)
            warnings.Add($"{string.Join(", ", plexesWithoutSamples.Select(p => $"plex {p}"))}: no sample channel, only pooled " +
                $"reference or unused channels, so {ionsWithoutSamples} ion(s) there are not in the table.");
        var undescribed = new SortedSet<string>(rows.SelectMany(r => r.Peptide.ProteinGroups)
            .Where(g => !groups.ContainsKey(g) && !string.Equals(g, ObservedPeptide.UndefinedProteinGroup, StringComparison.Ordinal)),
            StringComparer.Ordinal);
        if (undescribed.Count > 0)
            warnings.Add($"{undescribed.Count} protein group(s) named by peptides are not in the protein-group table, so " +
                $"their q-value and decoy/contaminant status are unknown: {string.Join(", ", undescribed.Take(5))}" +
                (undescribed.Count > 5 ? ", ..." : "") + ".");

        return new IsobaricObservationTable(samples, samplePlex, rows, log2.ToArray(), groups, fewChannels, warnings);
    }

    /// <summary>A larger sum wins; a tie goes to the earlier file name (ordinal), then the lower scan.</summary>
    private static bool Precedes(IsobaricPsm psm, double sum, Candidate held)
    {
        if (sum != held.Sum)
            return sum > held.Sum;
        int byFile = string.CompareOrdinal(psm.FileName, held.Psm.FileName);
        return byFile != 0 ? byFile < 0 : psm.ScanNumber < held.Psm.ScanNumber;
    }

    /// <summary>For each plex, the position of each of its samples within the plex (indexed by sample) and the plex's sample count.</summary>
    private static Dictionary<int, (int[] Positions, int Count)> PlexPositions(int[] samplePlex)
    {
        var result = new Dictionary<int, (int[], int)>();
        foreach (var plex in Enumerable.Range(0, samplePlex.Length).GroupBy(s => samplePlex[s]))
        {
            var positions = Enumerable.Repeat(-1, samplePlex.Length).ToArray();
            int k = 0;
            foreach (int s in plex)
                positions[s] = k++;
            result[plex.Key] = (positions, k);
        }
        return result;
    }

    /// <summary>Each design file by its name without extension, and each sample's plex.</summary>
    private static (Dictionary<string, DesignFile> Files, int[] SamplePlex) ReadDesign(IExperimentalDesign design,
        List<string> samples)
    {
        var sampleIndex = samples.Select((id, i) => (id, i)).ToDictionary(t => t.id, t => t.i, StringComparer.Ordinal);
        var samplePlex = new int[samples.Count];
        var files = new Dictionary<string, DesignFile>(StringComparer.Ordinal);
        var places = new Dictionary<(int, int, int), string>();

        foreach (var (key, entries) in design.FileNameSampleInfoDictionary.OrderBy(e => e.Key, StringComparer.Ordinal))
        {
            if (entries is null || entries.Length == 0)
                throw new ArgumentException($"Design file '{key}' has no channels.", nameof(design));
            var channels = new IsobaricQuantSampleInfo[entries.Length];
            for (int c = 0; c < entries.Length; c++)
            {
                if (entries[c] is not IsobaricQuantSampleInfo channel)
                    throw new ArgumentException(
                        $"Design file '{key}' has an entry that is not an isobaric channel; the isobaric observation table " +
                        "reads isobaric designs only.", nameof(design));
                channels[c] = channel;
            }

            IsobaricQuantSampleInfo first = channels[0];
            string name = ((ISampleInfo)first).FilenameWithoutExtension;
            if (channels.Any(c => c.PlexId != first.PlexId || c.TechnicalReplicate != first.TechnicalReplicate
                                  || c.Fraction != first.Fraction))
                throw new ArgumentException(
                    $"The channels of design file '{name}' are not all in one plex, technical replicate and fraction.",
                    nameof(design));
            if (files.ContainsKey(name))
                throw new ArgumentException($"Two design files have the name '{name}'.", nameof(design));
            var place = (first.PlexId, first.TechnicalReplicate, first.Fraction);
            if (places.TryGetValue(place, out string? other))
                throw new ArgumentException(
                    $"Design files '{other}' and '{name}' are both plex {first.PlexId}, technical replicate " +
                    $"{first.TechnicalReplicate}, fraction {first.Fraction}.", nameof(design));
            places[place] = name;

            var index = new int[channels.Length];
            for (int c = 0; c < channels.Length; c++)
            {
                // A reference channel or an unused one holds no sample (AnalysisSample.FromExperimentalDesign).
                index[c] = !channels[c].IsReferenceChannel
                           && sampleIndex.TryGetValue(AnalysisSample.IdFor(channels[c]), out int s) ? s : -1;
                if (index[c] >= 0)
                    samplePlex[index[c]] = channels[c].PlexId;
            }
            files[name] = new DesignFile(name, first.PlexId, first.TechnicalReplicate, first.Fraction, index);
        }
        return (files, samplePlex);
    }

    /// <summary>The match's design file, once the match is checked.</summary>
    private static DesignFile CheckPsm(IsobaricPsm psm, Dictionary<string, DesignFile> files,
        Dictionary<string, string> labelOfFile)
    {
        string where = $"Match '{psm.FullSequence}' (file '{psm.FileName}', scan {psm.ScanNumber})";
        if (string.IsNullOrEmpty(psm.FullSequence))
            throw new ArgumentException($"{where} has no full sequence.", "psms");
        if (psm.Charge < 1)
            throw new ArgumentException($"{where} has charge {psm.Charge}; a charge is at least 1.", "psms");
        if (psm.ProteinGroups is null || psm.ReporterIntensities is null)
            throw new ArgumentException($"{where} has no protein groups or no reporter intensities.", "psms");

        string? name = psm.FileName is null ? null
            : files.ContainsKey(psm.FileName) ? psm.FileName
            : psm.FileName.EndsWith(CalibratedSuffix, StringComparison.Ordinal)
              && files.ContainsKey(psm.FileName[..^CalibratedSuffix.Length])
                ? psm.FileName[..^CalibratedSuffix.Length]
                : null;
        if (name is null)
            throw new ArgumentException($"{where} is from a file the experimental design does not have.", "psms");
        if (labelOfFile.TryGetValue(name, out string? label) && !string.Equals(label, psm.FileName, StringComparison.Ordinal))
            throw new ArgumentException(
                $"The matches name design file '{name}' both '{label}' and '{psm.FileName}'.", "psms");
        labelOfFile[name] = psm.FileName!;

        DesignFile file = files[name];
        if (psm.ReporterIntensities.Count != file.SampleIndex.Length)
            throw new ArgumentException(
                $"{where} has {psm.ReporterIntensities.Count} reporter intensities, but the design gives its file " +
                $"{file.SampleIndex.Length} channels.", "psms");
        foreach (double value in psm.ReporterIntensities)
            if (!double.IsFinite(value) || value < 0)
                throw new ArgumentException(
                    $"{where} has reporter intensity {value}; an intensity is finite and not negative, 0 when not measured.",
                    "psms");
        return file;
    }

    /// <summary>The match's peptide, the same object for every match of one full sequence, which must agree.</summary>
    private static ObservedPeptide PeptideOf(IsobaricPsm psm, Dictionary<string, ObservedPeptide> peptides)
    {
        foreach (string group in psm.ProteinGroups)
            if (group is null || group.Contains(';'))
                throw new ArgumentException(
                    $"Peptide '{psm.FullSequence}' names protein group '{group}'; a group name is not null and cannot " +
                    "contain ';', which separates groups in the peptide table.", "psms");
        var names = psm.ProteinGroups.Distinct(StringComparer.Ordinal).Order(StringComparer.Ordinal).ToList();

        if (!peptides.TryGetValue(psm.FullSequence, out var peptide))
        {
            peptide = new ObservedPeptide(psm.FullSequence, psm.BaseSequence, names);
            peptides[psm.FullSequence] = peptide;
            return peptide;
        }
        if (!string.Equals(peptide.BaseSequence, psm.BaseSequence, StringComparison.Ordinal)
            || !peptide.ProteinGroups.SequenceEqual(names, StringComparer.Ordinal))
            throw new ArgumentException(
                $"Peptide '{psm.FullSequence}' has base sequence '{peptide.BaseSequence}' and groups " +
                $"'{string.Join(";", peptide.ProteinGroups)}' in one match but '{psm.BaseSequence}' and " +
                $"'{string.Join(";", names)}' in another.", "psms");
        return peptide;
    }
}
