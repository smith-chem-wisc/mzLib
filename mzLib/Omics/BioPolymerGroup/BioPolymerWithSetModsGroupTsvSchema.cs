namespace Omics.BioPolymerGroup;

/// <summary>
/// The tab-separated report schema for <see cref="BioPolymerWithSetModsGroup"/>.
///
/// Deliberately not the parent-level schema with different labels: a peptide group has no
/// sequence coverage, gene, or organism, and its occupancy positions are peptide-local. What the
/// two reports do share — the per-sample quantification block — comes from
/// <see cref="SampleGroupColumnBuilder"/>, so the schema is built from the whole dataset and every
/// row is the width of the header.
/// </summary>
public static class BioPolymerWithSetModsGroupTsvSchema
{
    /// <summary>
    /// Builds the schema for a set of groups that will be written to one file.
    /// </summary>
    /// <param name="groups">All groups destined for the file. Their sample group results are
    /// populated if they have not been already, since the quantification columns derive from them.</param>
    public static IReadOnlyList<TsvColumn<BioPolymerWithSetModsGroup>> For(
        IReadOnlyCollection<BioPolymerWithSetModsGroup> groups)
    {
        var columns = new List<TsvColumn<BioPolymerWithSetModsGroup>>
        {
            new("Base Sequence", g => g.BaseSequence),
            new("Peptidoforms", g => Truncate(string.Join("|", Peptidoforms(g)))),
            new("Number of Peptidoforms", g => Peptidoforms(g).Count().ToString()),
            new("Parent Accessions", g => Truncate(string.Join("|",
                Occurrences(g).Select(o => o.Accession)))),
            new("Number of Parents", g => g.ParentBioPolymers.Count.ToString()),
            new("Start Residue in Parent", g => Truncate(string.Join("|",
                Occurrences(g).Select(o => o.Start.ToString())))),
            new("End Residue in Parent", g => Truncate(string.Join("|",
                Occurrences(g).Select(o => o.End.ToString()))))
        };

        columns.AddRange(QuantificationColumns(groups));

        columns.AddRange(new TsvColumn<BioPolymerWithSetModsGroup>[]
        {
            new("Number of PSMs", g => g.AllPsmsBelowOnePercentFDR.Count.ToString()),
            new("Decoy/Contaminant/Target", TargetDecoyLabel),
            new("Cumulative Target", g => g.CumulativeTarget.ToString()),
            new("Cumulative Decoy", g => g.CumulativeDecoy.ToString()),
            new("QValue", g => g.QValue.ToString()),
            new("Best PSM Score", g => g.BestPsmScore.ToString()),
            new("Best PSM Notch QValue", g => g.BestPsmQValue.ToString())
        });

        return columns;
    }

    /// <summary>
    /// One block of columns per sample group, with occupancy rendered at peptide-local positions.
    /// </summary>
    private static IEnumerable<TsvColumn<BioPolymerWithSetModsGroup>> QuantificationColumns(
        IReadOnlyCollection<BioPolymerWithSetModsGroup> groups)
        => SampleGroupColumnBuilder.Build(
            groups,
            SampleGroupsOf,
            (g, result) => Truncate(result.FormatOccupancy([g.BaseSequence], proteinLevel: false, intensityBased: false)),
            (g, result) => Truncate(result.FormatOccupancy([g.BaseSequence], proteinLevel: false, intensityBased: true)));

    private static IReadOnlyList<SampleGroupResult> SampleGroupsOf(BioPolymerWithSetModsGroup group)
    {
        if (group.SampleGroupResults is null)
            group.PopulateSampleGroupResults();

        return group.SampleGroupResults!;
    }

    /// <summary>
    /// Distinct modified forms, ordered so output is stable across runs.
    /// </summary>
    private static IEnumerable<string> Peptidoforms(BioPolymerWithSetModsGroup group)
        => group.Peptidoforms.Select(p => p.FullSequence).Distinct().OrderBy(s => s, StringComparer.Ordinal);

    /// <summary>
    /// Every place this sequence was found: one entry per (parent, residue range), ordered by
    /// accession then position.
    /// </summary>
    /// <remarks>
    /// One entry per occurrence rather than per parent, because a sequence can occur more than once
    /// in the same parent. Reporting per parent silently dropped the later occurrences. The accession
    /// repeats for each occurrence so that Parent Accessions, Start Residue and End Residue remain
    /// readable as three parallel lists, one entry apiece.
    /// </remarks>
    private static IEnumerable<(string Accession, int Start, int End)> Occurrences(BioPolymerWithSetModsGroup group)
        => group.Peptidoforms
            .Where(p => p.Parent is not null)
            .Select(p => (Accession: p.Parent.Accession, Start: p.OneBasedStartResidue, End: p.OneBasedEndResidue))
            .Distinct()
            .OrderBy(o => o.Accession, StringComparer.Ordinal)
            .ThenBy(o => o.Start);

    /// <summary>
    /// Single-letter classification: entrapment decoy, entrapment target, decoy, contaminant, or target.
    /// </summary>
    private static string TargetDecoyLabel(BioPolymerWithSetModsGroup group)
    {
        if (group.IsEntrapment && group.IsDecoy) return "ED";
        if (group.IsEntrapment) return "ET";
        if (group.IsDecoy) return "D";
        if (group.IsContaminant) return "C";
        return "T";
    }

    private static string Truncate(string? input) => TsvWriter.Truncate(input);
}
