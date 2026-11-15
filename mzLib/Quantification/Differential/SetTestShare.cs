using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text.RegularExpressions;
using Omics;

namespace Quantification.Differential;

/// <summary>The three shapes of a <c>marker_region</c> token (<c>DEF-DIFF-SHARE-REGION</c>).</summary>
public enum ShareRegionKind
{
    /// <summary><c>a-b</c>: the residues a to b.</summary>
    Range,
    /// <summary><c>p</c>: one residue.</summary>
    Residue,
    /// <summary><c>junction a/b</c>: the junction between residues a and b.</summary>
    Junction,
}

/// <summary>
/// One region of an isoform, numbered on that accession's own sequence (QuantProject <c>DEF-DIFF-SHARE-REGION</c>).
/// </summary>
public sealed record ShareRegion
{
    /// <summary>Residues of the accession taken on each side of a region to make its signature.</summary>
    public const int SignatureFlank = 5;

    /// <summary>The accession the region is numbered on: the row's <c>isoform_accession</c>.</summary>
    public required string Accession { get; init; }
    /// <summary>The shape.</summary>
    public required ShareRegionKind Kind { get; init; }
    /// <summary>a for a range or a junction, p for a residue (1-based).</summary>
    public required int First { get; init; }
    /// <summary>b for a range or a junction, p for a residue (1-based).</summary>
    public required int Last { get; init; }

    /// <summary>The token as the contract writes it, without the accession: <c>574-627</c>, <c>994</c> or <c>junction 288/289</c>.</summary>
    public string Text => Kind switch
    {
        ShareRegionKind.Range => string.Create(CultureInfo.InvariantCulture, $"{First}-{Last}"),
        ShareRegionKind.Residue => First.ToString(CultureInfo.InvariantCulture),
        _ => string.Create(CultureInfo.InvariantCulture, $"junction {First}/{Last}"),
    };

    /// <summary>
    /// The region lies inside a sequence of <paramref name="length"/> residues: 1 ≤ a ≤ b ≤ length, and a junction joins
    /// two neighbouring residues (b = a + 1), as the accession's own numbering has them.
    /// </summary>
    public bool LiesIn(int length) => First >= 1 && Last <= length && Kind switch
    {
        ShareRegionKind.Range => First <= Last,
        ShareRegionKind.Residue => First == Last,
        _ => Last == First + 1,
    };

    /// <summary>
    /// A peptide at residues <paramref name="start"/> to <paramref name="end"/> of the accession overlaps the region: it
    /// covers at least one residue of a range or a residue, or residue a or b of a junction.
    /// </summary>
    public bool Overlaps(int start, int end) => Kind == ShareRegionKind.Junction
        ? (start <= First && First <= end) || (start <= Last && Last <= end)
        : start <= Last && end >= First;

    /// <summary>
    /// The region's signature in <paramref name="sequence"/>: its own residues with up to <see cref="SignatureFlank"/> of
    /// the accession on each side; for a junction, the 5 residues up to a and the 5 from b.
    /// </summary>
    /// <exception cref="ArgumentException">The region does not lie in <paramref name="sequence"/>.</exception>
    public string Signature(string sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (!LiesIn(sequence.Length))
            throw new ArgumentException($"{Accession}:{Text} does not lie in a sequence of {sequence.Length} residues.", nameof(sequence));
        // A range or residue: residues a - 5 to b + 5. A junction: a - 4 to a, then b to b + 4.
        int from = Kind == ShareRegionKind.Junction ? Math.Max(0, First - SignatureFlank) : Math.Max(0, First - 1 - SignatureFlank);
        int to = Kind == ShareRegionKind.Junction ? Math.Min(sequence.Length, Last - 1 + SignatureFlank) : Math.Min(sequence.Length, Last + SignatureFlank);
        return sequence[from..to];
    }
}

/// <summary>Reads a <c>marker_region</c> cell (<c>DEF-DIFF-SHARE-REGION</c>).</summary>
public static class ShareRegionReader
{
    private static readonly Regex Token = new(
        @"(?<![A-Za-z0-9_.\-])(?<acc>[A-Za-z0-9][A-Za-z0-9_.\-]*):\s*(?:junction\s+(?<ja>\d+)\s*/\s*(?<jb>\d+)|(?<a>\d+)\s*-\s*(?<b>\d+)|(?<p>\d+))",
        RegexOptions.CultureInvariant | RegexOptions.IgnoreCase);

    /// <summary>
    /// Every <c>&lt;isoform_accession&gt;:</c> token of <paramref name="markerRegion"/>, in the order written: <c>a-b</c>,
    /// <c>p</c> or <c>junction a/b</c>. A token of another accession, and any other text, is a note. A number too large to
    /// read gives a region that lies in no sequence.
    /// </summary>
    public static IReadOnlyList<ShareRegion> Read(string? markerRegion, string isoformAccession)
    {
        ArgumentNullException.ThrowIfNull(isoformAccession);
        if (string.IsNullOrWhiteSpace(markerRegion)) return Array.Empty<ShareRegion>();
        static int N(Group g) => int.TryParse(g.Value, NumberStyles.None, CultureInfo.InvariantCulture, out int v) ? v : -1;
        var regions = new List<ShareRegion>();
        foreach (Match m in Token.Matches(markerRegion))
        {
            if (!string.Equals(m.Groups["acc"].Value, isoformAccession, StringComparison.Ordinal)) continue;
            regions.Add(m.Groups["ja"].Success
                ? new ShareRegion { Accession = isoformAccession, Kind = ShareRegionKind.Junction, First = N(m.Groups["ja"]), Last = N(m.Groups["jb"]) }
                : m.Groups["a"].Success
                    ? new ShareRegion { Accession = isoformAccession, Kind = ShareRegionKind.Range, First = N(m.Groups["a"]), Last = N(m.Groups["b"]) }
                    : new ShareRegion { Accession = isoformAccession, Kind = ShareRegionKind.Residue, First = N(m.Groups["p"]), Last = N(m.Groups["p"]) });
        }
        return regions;
    }
}

/// <summary>
/// The sequences of a search's protein databases, for <c>DEF-DIFF-SHARE-PEPTIDES</c>: which entries contain a peptide.
/// Containment reads I and L as one residue (ruled by the user 2026-11-14), since a search cannot tell them apart.
/// </summary>
public sealed class ShareSequenceIndex
{
    private readonly IReadOnlyDictionary<string, string> _sequences;
    private readonly IReadOnlyDictionary<string, string> _foldedByAccession;
    private readonly IReadOnlyList<(string Accession, string Folded)> _folded;

    private ShareSequenceIndex(IReadOnlyDictionary<string, string> sequences)
    {
        _sequences = sequences;
        _foldedByAccession = sequences.ToDictionary(e => e.Key, e => Fold(e.Value), StringComparer.Ordinal);
        _folded = _foldedByAccession.OrderBy(e => e.Key, StringComparer.Ordinal).Select(e => (e.Key, e.Value)).ToList();
    }

    /// <summary>The index of a search's database entries, by <see cref="IBioPolymer.BaseSequence"/>.</summary>
    /// <exception cref="ArgumentException">One accession with two different sequences.</exception>
    public static ShareSequenceIndex FromBioPolymers(IEnumerable<IBioPolymer> entries)
    {
        ArgumentNullException.ThrowIfNull(entries);
        return FromSequences(entries.Select(e => (e.Accession, e.BaseSequence)));
    }

    /// <summary>
    /// The index from (accession, sequence) pairs. An accession given twice with the same sequence (a contaminant
    /// database repeating a target entry) is one entry.
    /// </summary>
    /// <exception cref="ArgumentException">An empty accession or sequence, or one accession with two different sequences.</exception>
    public static ShareSequenceIndex FromSequences(IEnumerable<(string Accession, string Sequence)> entries)
    {
        ArgumentNullException.ThrowIfNull(entries);
        var map = new Dictionary<string, string>(StringComparer.Ordinal);
        foreach (var (accession, sequence) in entries)
        {
            if (string.IsNullOrEmpty(accession) || string.IsNullOrEmpty(sequence))
                throw new ArgumentException("An entry has no accession or no sequence.", nameof(entries));
            if (map.TryGetValue(accession, out string? known) && !string.Equals(known, sequence, StringComparison.Ordinal))
                throw new ArgumentException($"{accession} is given with two different sequences, so its numbering is not one.", nameof(entries));
            map[accession] = sequence;
        }
        return new ShareSequenceIndex(map);
    }

    /// <summary>I and L read as one residue.</summary>
    public static string Fold(string sequence) => sequence.Replace('I', 'L');

    /// <summary>The databases hold this accession.</summary>
    public bool Contains(string accession) => _sequences.ContainsKey(accession);

    /// <summary>The accession's sequence, or null.</summary>
    public string? SequenceOf(string accession) => _sequences.TryGetValue(accession, out string? s) ? s : null;

    /// <summary>Every accession whose sequence contains <paramref name="peptide"/>, I and L read as one, ordinal-sorted.</summary>
    public IReadOnlyList<string> AccessionsContaining(string peptide)
    {
        ArgumentException.ThrowIfNullOrEmpty(peptide);
        string folded = Fold(peptide);
        return _folded.Where(e => e.Folded.Contains(folded, StringComparison.Ordinal)).Select(e => e.Accession).ToList();
    }

    /// <summary>The accession's sequence contains <paramref name="text"/>, I and L read as one.</summary>
    public bool Holds(string accession, string text) =>
        _foldedByAccession.TryGetValue(accession, out string? s) && s.Contains(Fold(text), StringComparison.Ordinal);
}

/// <summary>How one base sequence of a member's gene is used (<c>DEF-DIFF-SHARE-PEPTIDES</c>, rulings S2 and S3).</summary>
public enum SharePeptideClass
{
    /// <summary>Overlaps a region, and every entry containing it contains that region.</summary>
    Marker,
    /// <summary>Overlaps no region, is in every one of the member's entries, and in no other gene's entry.</summary>
    Reference,
    /// <summary>The gene's, but neither: counted, not used.</summary>
    Neither,
}

/// <summary>The peptides of a member's gene, by class.</summary>
public sealed record SharePeptides
{
    /// <summary>The marker peptides: every full sequence whose base sequence is a marker, in the table's order.</summary>
    public IReadOnlyList<ObservedPeptide> Markers { get; init; } = Array.Empty<ObservedPeptide>();
    /// <summary>The reference peptides, likewise.</summary>
    public IReadOnlyList<ObservedPeptide> References { get; init; } = Array.Empty<ObservedPeptide>();
    /// <summary>The gene's other peptides, likewise.</summary>
    public IReadOnlyList<ObservedPeptide> Neither { get; init; } = Array.Empty<ObservedPeptide>();

    /// <summary>The markers' base sequences, each once, ordinal-sorted (<c>marker_peptides</c>).</summary>
    public IReadOnlyList<string> MarkerBaseSequences => Markers.Select(p => p.BaseSequence).Distinct(StringComparer.Ordinal).Order(StringComparer.Ordinal).ToList();
}

/// <summary>Classes a member's peptides (<c>DEF-DIFF-SHARE-PEPTIDES</c>).</summary>
public static class SharePeptideClassifier
{
    /// <summary>
    /// The class of every peptide of the table that belongs to the member's gene (is in one of its entries); a peptide in
    /// none of them is not the gene's and is left out. Peptides are placed by base sequence in the isoform's sequence, at
    /// every place they occur, with no enzyme assumed; a localization-ambiguous row has its one base sequence.
    /// </summary>
    /// <param name="peptides">The dataset's peptide rows.</param>
    /// <param name="isoformSequence">The sequence of <c>isoform_accession</c>.</param>
    /// <param name="regions">Its regions, each lying in it.</param>
    /// <param name="memberAccessions">The member's entries (<c>DEF-DIFF-SET-MEMBER</c>'s gene rule).</param>
    /// <param name="memberGene">The member's gene when found by gene name; an entry naming any other gene (or none) is then
    /// another gene's. Null when found by <c>member_accessions</c>: any entry outside them is then another gene's.</param>
    /// <param name="sequences">The search's databases.</param>
    /// <param name="genes">The same databases' gene names.</param>
    public static SharePeptides Classify(IEnumerable<ObservedPeptide> peptides, string isoformSequence, IReadOnlyList<ShareRegion> regions,
        IReadOnlyCollection<string> memberAccessions, string? memberGene, ShareSequenceIndex sequences, SetGeneIndex genes)
    {
        ArgumentNullException.ThrowIfNull(peptides);
        ArgumentNullException.ThrowIfNull(isoformSequence);
        ArgumentNullException.ThrowIfNull(regions);
        ArgumentNullException.ThrowIfNull(memberAccessions);
        ArgumentNullException.ThrowIfNull(sequences);
        ArgumentNullException.ThrowIfNull(genes);
        if (regions.Any(r => !r.LiesIn(isoformSequence.Length)))
            throw new ArgumentException("Every region must lie in the isoform's sequence.", nameof(regions));

        var mine = new HashSet<string>(memberAccessions, StringComparer.Ordinal);
        var signatures = regions.Select(r => r.Signature(isoformSequence)).ToList();
        string foldedIsoform = ShareSequenceIndex.Fold(isoformSequence);
        var classOf = new Dictionary<string, SharePeptideClass?>(StringComparer.Ordinal);

        SharePeptideClass? ClassOf(string baseSequence)
        {
            if (!mine.Any(a => sequences.Holds(a, baseSequence))) return null;
            var holders = sequences.AccessionsContaining(baseSequence);

            string folded = ShareSequenceIndex.Fold(baseSequence);
            var overlapped = new HashSet<int>();
            for (int at = foldedIsoform.IndexOf(folded, StringComparison.Ordinal); at >= 0;
                 at = foldedIsoform.IndexOf(folded, at + 1, StringComparison.Ordinal))
                for (int r = 0; r < regions.Count; r++)
                    if (regions[r].Overlaps(at + 1, at + folded.Length)) overlapped.Add(r);

            if (overlapped.Count > 0)
                return overlapped.Any(r => holders.All(h => sequences.Holds(h, signatures[r]))) ? SharePeptideClass.Marker : SharePeptideClass.Neither;

            bool inEveryEntry = mine.All(a => sequences.Holds(a, baseSequence));
            bool inNoOtherGene = holders.All(h => mine.Contains(h) && (memberGene is null || genes.GenesOf(h) is { Count: 1 } g
                && string.Equals(g[0], memberGene, StringComparison.OrdinalIgnoreCase)));
            return inEveryEntry && inNoOtherGene ? SharePeptideClass.Reference : SharePeptideClass.Neither;
        }

        var markers = new List<ObservedPeptide>();
        var references = new List<ObservedPeptide>();
        var neither = new List<ObservedPeptide>();
        foreach (var peptide in peptides)
        {
            if (string.IsNullOrEmpty(peptide.BaseSequence)) continue;
            if (!classOf.TryGetValue(peptide.BaseSequence, out var c)) classOf[peptide.BaseSequence] = c = ClassOf(peptide.BaseSequence);
            (c switch { SharePeptideClass.Marker => markers, SharePeptideClass.Reference => references, SharePeptideClass.Neither => neither, _ => null })
                ?.Add(peptide);
        }
        return new SharePeptides { Markers = markers, References = references, Neither = neither };
    }
}

/// <summary>A share member's state and peptides on one dataset contrast (<c>DEF-DIFF-SHARE-STATE</c>).</summary>
public sealed record ShareMemberResult
{
    /// <summary>The member's resolution (<see cref="SetMemberResolver"/>), with its final <see cref="SetMemberAssignment.State"/>.</summary>
    public required SetMemberAssignment Assignment { get; init; }
    /// <summary>The member's row in the contrast: its <c>measure</c>, <c>isoform_accession</c> and <c>marker_region</c>.</summary>
    public required SetTestRow Row { get; init; }
    /// <summary>The regions read from <c>marker_region</c>.</summary>
    public IReadOnlyList<ShareRegion> Regions { get; init; } = Array.Empty<ShareRegion>();
    /// <summary>The gene's peptides by class; null when the member stopped before they were classed.</summary>
    public SharePeptides? Peptides { get; init; }
    /// <summary>Numerator-side samples in which a marker peptide has a value.</summary>
    public int MarkerSamplesNumerator { get; init; }
    /// <summary>Reference-side samples in which a marker peptide has a value.</summary>
    public int MarkerSamplesDenominator { get; init; }
    /// <summary>Numerator-side samples in which a reference peptide has a value.</summary>
    public int ReferenceSamplesNumerator { get; init; }
    /// <summary>Reference-side samples in which a reference peptide has a value.</summary>
    public int ReferenceSamplesDenominator { get; init; }
    /// <summary>Reference base sequences with a value in at least one of the contrast's samples.</summary>
    public int ReferencesWithValue { get; init; }
    /// <summary><c>listed_marker_used</c>: every <c>marker_peptide</c> of the row is among the markers found; null when none is listed or none was classed.</summary>
    public bool? ListedMarkerUsed { get; init; }
    /// <summary>The member enters the share model.</summary>
    public bool InEstimate => Assignment.State == SetMemberStateFromRow.Tested;
}

/// <summary>A share member's state, in the contract's order (<c>DEF-DIFF-SHARE-STATE</c>, ruling S4).</summary>
public static class ShareMemberStates
{
    /// <summary><c>isoform_peptide</c> only: <c>isoform_accession</c> is not an entry of the search's databases.</summary>
    public const string IsoformNotSearched = "isoform_not_searched";
    /// <summary>The row's <c>marker_region</c> has no token of its accession, a region lies outside the accession's sequence, or the accession is not searched.</summary>
    public const string UnreadableRegion = "unreadable_region";
    /// <summary>No marker peptide, or the markers have values in fewer than 2 samples on a side.</summary>
    public const string TooFewMarkers = "too_few:markers";
    /// <summary>Fewer than 2 reference peptides, or the reference has values in fewer than 2 samples on a side.</summary>
    public const string TooFewReference = "too_few:reference";
    /// <summary>The fewest reference peptides (base sequences with a value) a member needs.</summary>
    public const int MinimumReferencePeptides = 2;
    /// <summary>The fewest samples on each side in which the markers, and the reference, need a value.</summary>
    public const int MinimumSamplesPerSide = 2;
    /// <summary>The <c>measure</c> of a row tested on one form's peptides.</summary>
    public const string IsoformPeptide = "isoform_peptide";

    /// <summary>
    /// Each resolved member's final share state, the first that applies: <see cref="SetMemberResolver"/>'s states; then
    /// <see cref="IsoformNotSearched"/>, <see cref="UnreadableRegion"/>, <c>seen_one_side:*</c>, <see cref="TooFewMarkers"/>,
    /// <see cref="TooFewReference"/>, <c>incomplete</c> (no-change predictions only) and <c>tested</c>.
    /// </summary>
    /// <param name="contrast">The dataset contrast; each member's row gives its accession and region.</param>
    /// <param name="prediction">The part's prediction.</param>
    /// <param name="resolved">Its members, resolved (<see cref="SetMemberResolver.Resolve"/>).</param>
    /// <param name="setup">The contrast's setup; its samples on each side are counted.</param>
    /// <param name="normalized">The stratum's normalized observation table.</param>
    /// <param name="sequences">The search's databases' sequences.</param>
    /// <param name="genes">The same databases' gene names.</param>
    /// <returns>One result per member, in <paramref name="resolved"/>'s order.</returns>
    public static IReadOnlyList<ShareMemberResult> Apply(SetDatasetContrast contrast, SetPrediction prediction,
        IReadOnlyList<SetMemberAssignment> resolved, SetContrastSetup setup, ObservationTable normalized,
        ShareSequenceIndex sequences, SetGeneIndex genes)
    {
        ArgumentNullException.ThrowIfNull(contrast);
        ArgumentNullException.ThrowIfNull(resolved);
        ArgumentNullException.ThrowIfNull(setup);
        ArgumentNullException.ThrowIfNull(normalized);
        ArgumentNullException.ThrowIfNull(sequences);
        ArgumentNullException.ThrowIfNull(genes);
        var column = normalized.Samples.Select((s, i) => (s, i)).ToDictionary(t => t.s, t => t.i, StringComparer.Ordinal);
        int[] Columns(IReadOnlyList<string> samples) => samples.Select(s => column.TryGetValue(s, out int c) ? c
            : throw new ArgumentException($"The table has no sample '{s}'.", nameof(normalized))).ToArray();
        int[] numerator = Columns(setup.NumeratorSamples), denominator = Columns(setup.DenominatorSamples);
        var rowOf = new Dictionary<ObservedPeptide, int>(ReferenceEqualityComparer.Instance);
        for (int i = 0; i < normalized.Peptides.Count; i++) rowOf[normalized.Peptides[i]] = i;
        int SamplesWithValue(IReadOnlyList<ObservedPeptide> peptides, int[] side) =>
            side.Count(c => peptides.Any(p => double.IsFinite(normalized.Log2Intensity(rowOf[p], c))));

        var results = new List<ShareMemberResult>();
        foreach (var assignment in resolved)
        {
            var row = contrast.Rows.FirstOrDefault(r => string.Equals(r.Member, assignment.Member, StringComparison.Ordinal))
                ?? throw new ArgumentException($"{assignment.Member} is not a member of {contrast.Dataset}'s contrast.", nameof(resolved));
            var result = new ShareMemberResult { Assignment = assignment, Row = row };
            ShareMemberResult Stop(string state) => result with { Assignment = assignment with { State = state } };
            if (assignment.State is not null) { results.Add(result); continue; }

            string accession = row.IsoformAccession ?? "";
            if (string.Equals(row.Measure, IsoformPeptide, StringComparison.Ordinal) && !sequences.Contains(accession))
            {
                results.Add(Stop(IsoformNotSearched));
                continue;
            }
            var regions = ShareRegionReader.Read(row.MarkerRegion, accession);
            result = result with { Regions = regions };
            string? sequence = sequences.SequenceOf(accession);
            if (sequence is null || regions.Count == 0 || regions.Any(r => !r.LiesIn(sequence.Length)))
            {
                results.Add(Stop(UnreadableRegion));
                continue;
            }

            string? gene = row.MemberAccessions.Count > 0 ? null : row.Member;
            var accessions = row.MemberAccessions.Count > 0
                ? row.MemberAccessions.Distinct(StringComparer.Ordinal).ToList()
                : genes.AccessionsOf(row.Member).ToList();
            var classes = SharePeptideClassifier.Classify(normalized.Peptides, sequence, regions, accessions, gene, sequences, genes);
            var listed = row.MarkerPeptides;
            int mn = SamplesWithValue(classes.Markers, numerator), md = SamplesWithValue(classes.Markers, denominator);
            int rn = SamplesWithValue(classes.References, numerator), rd = SamplesWithValue(classes.References, denominator);
            int referencesWithValue = classes.References
                .Where(p => numerator.Concat(denominator).Any(c => double.IsFinite(normalized.Log2Intensity(rowOf[p], c))))
                .Select(p => p.BaseSequence).Distinct(StringComparer.Ordinal).Count();
            result = result with
            {
                Peptides = classes,
                MarkerSamplesNumerator = mn,
                MarkerSamplesDenominator = md,
                ReferenceSamplesNumerator = rn,
                ReferenceSamplesDenominator = rd,
                ReferencesWithValue = referencesWithValue,
                ListedMarkerUsed = listed.Count == 0 ? null : listed.All(classes.MarkerBaseSequences.Contains),
            };

            string state =
                mn > 0 && md == 0 ? SetMemberStateFromRow.SeenNumeratorOnly
                : md > 0 && mn == 0 ? SetMemberStateFromRow.SeenDenominatorOnly
                : mn < MinimumSamplesPerSide || md < MinimumSamplesPerSide ? TooFewMarkers
                : referencesWithValue < MinimumReferencePeptides || rn < MinimumSamplesPerSide || rd < MinimumSamplesPerSide ? TooFewReference
                : SetVerdict.IsNoChange(prediction) && (SetMemberStateFromRow.IncompleteSide(mn, numerator.Length)
                    || SetMemberStateFromRow.IncompleteSide(md, denominator.Length) || SetMemberStateFromRow.IncompleteSide(rn, numerator.Length)
                    || SetMemberStateFromRow.IncompleteSide(rd, denominator.Length)) ? SetMemberStateFromRow.Incomplete
                : SetMemberStateFromRow.Tested;
            results.Add(result with { Assignment = assignment with { State = state } });
        }
        return results;
    }
}
