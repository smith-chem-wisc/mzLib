namespace Omics.BioPolymerGroup;

/// <summary>
/// The tab-separated report schema for <see cref="BioPolymerGroup"/>: which columns appear, in
/// what order, and how each value is formatted.
///
/// The schema is built from the whole dataset rather than from a single group, because the
/// quantification columns depend on the experimental design shared by the file. Every group is
/// then rendered against that one schema, so a group that was not quantified in some condition
/// contributes empty fields instead of silently emitting a narrower row.
/// </summary>
public static class BioPolymerGroupTsvSchema
{
    /// <summary>
    /// Builds the schema for a set of groups that will be written to one file.
    /// </summary>
    /// <param name="groups">All groups destined for the file. Their sample group results are
    /// populated if they have not been already, since the quantification columns derive from them.</param>
    public static IReadOnlyList<TsvColumn<BioPolymerGroup>> For(IReadOnlyCollection<BioPolymerGroup> groups)
    {
        var columns = new List<TsvColumn<BioPolymerGroup>>
        {
            new("BioPolymer Accession", g => g.BioPolymerGroupName),
            new("Gene", g => Truncate(string.Join("|",
                g.ListOfBioPolymersOrderedByAccession.Select(p => p.GeneNames.Select(x => x.Item2).FirstOrDefault())))),
            new("Organism", g => Truncate(string.Join("|",
                g.ListOfBioPolymersOrderedByAccession.Select(p => p.Organism).Distinct()))),
            new("BioPolymer Full Name", g => Truncate(string.Join("|",
                g.ListOfBioPolymersOrderedByAccession.Select(p => p.FullName).Distinct()))),
            new("BioPolymer Unmodified Mass", g => Truncate(string.Join("|", UnmodifiedMasses(g)))),
            new("Number of BioPolymers in Group", g => g.BioPolymers.Count.ToString()),
            new("Unique Sequences", g => Truncate(SequenceList(g, g.UniqueBioPolymersWithSetMods))),
            new("Shared Sequences", g => Truncate(SequenceList(g,
                g.AllBioPolymersWithSetMods.Except(g.UniqueBioPolymersWithSetMods)))),
            new("Number of Sequences", g => CountDistinct(g, g.AllBioPolymersWithSetMods).ToString()),
            new("Number of Unique Sequences", g => CountDistinct(g, g.UniqueBioPolymersWithSetMods).ToString()),
            new("Sequence Coverage Fraction", g => Truncate(string.Join("|",
                Coverage(g).SequenceCoverageFraction.Select(p => string.Format("{0:0.#####}", p))))),
            new("Sequence Coverage", g => Truncate(string.Join("|", Coverage(g).SequenceCoverageDisplayList))),
            new("Sequence Coverage with Mods", g => Truncate(string.Join("|", Coverage(g).SequenceCoverageDisplayListWithMods))),
            new("Fragment Sequence Coverage", g => Truncate(string.Join("|", Coverage(g).FragmentSequenceCoverageDisplayList)))
        };

        columns.AddRange(QuantificationColumns(groups));

        columns.AddRange(new TsvColumn<BioPolymerGroup>[]
        {
            new("Number of PSMs", g => g.AllPsmsBelowOnePercentFDR.Count.ToString()),
            new("BioPolymer Decoy/Contaminant/Target", TargetDecoyLabel),
            new("BioPolymer Cumulative Target", g => g.CumulativeTarget.ToString()),
            new("BioPolymer Cumulative Decoy", g => g.CumulativeDecoy.ToString()),
            new("BioPolymer QValue", g => g.QValue.ToString()),
            new("Best Sequence Score", g => g.BestBioPolymerWithSetModsScore.ToString()),
            new("Best Sequence Notch QValue", g => g.BestBioPolymerWithSetModsQValue.ToString())
        });

        return columns;
    }

    /// <summary>
    /// One block of columns per sample group in the dataset, with occupancy rendered in this
    /// group's coordinate space.
    /// </summary>
    private static IEnumerable<TsvColumn<BioPolymerGroup>> QuantificationColumns(
        IReadOnlyCollection<BioPolymerGroup> groups)
        => SampleGroupColumnBuilder.Build(
            groups,
            SampleGroupsOf,
            (g, result) => Truncate(result.FormatOccupancy(OccupancyKeys(g), IsParentLevel(g), intensityBased: false)),
            (g, result) => Truncate(result.FormatOccupancy(OccupancyKeys(g), IsParentLevel(g), intensityBased: true)));

    private static IReadOnlyList<SampleGroupResult> SampleGroupsOf(BioPolymerGroup group)
    {
        if (group.SampleGroupResults is null)
            group.PopulateSampleGroupResults();

        return group.SampleGroupResults!;
    }

    private static bool IsParentLevel(BioPolymerGroup group)
        => group.GroupType == BioPolymerGroupType.Parent;

    /// <summary>
    /// Entity keys the occupancy string is ordered by: accessions at parent level, base sequences
    /// at digestion-product level.
    /// </summary>
    private static List<string> OccupancyKeys(BioPolymerGroup group)
        => (IsParentLevel(group)
                ? group.ListOfBioPolymersOrderedByAccession.Select(p => p.Accession)
                : group.AllBioPolymersWithSetMods.Select(p => p.BaseSequence).Distinct().OrderBy(s => s))
            .ToList();

    /// <summary>
    /// Coverage already computed for this group, or an empty result. Deliberately does not trigger
    /// calculation — writing a report should not silently do the expensive work.
    /// </summary>
    /// <remarks>
    /// The empty result is allocated per call rather than shared: its list properties are get-only
    /// but the lists themselves are mutable, so a single shared instance would be one stray Add away
    /// from leaking coverage between unrelated groups.
    /// </remarks>
    private static BioPolymerGroup.SequenceCoverageResult Coverage(BioPolymerGroup group)
        => group.IsSequenceCoverageCalculated ? group.CoverageResult : new BioPolymerGroup.SequenceCoverageResult();

    /// <summary>
    /// Unmodified mass per distinct biopolymer sequence, in accession order.
    /// Yields NaN for a sequence with no matching identified form.
    /// </summary>
    private static IEnumerable<double> UnmodifiedMasses(BioPolymerGroup group)
        => group.ListOfBioPolymersOrderedByAccession
            .Select(p => p.BaseSequence)
            .Distinct()
            .Select(sequence => group.AllBioPolymersWithSetMods
                .FirstOrDefault(bpws => bpws.BaseSequence == sequence)?.MonoisotopicMass ?? double.NaN);

    private static string SequenceList(BioPolymerGroup group, IEnumerable<IBioPolymerWithSetMods> sequences)
        => string.Join("|", Keys(group, sequences).Distinct());

    private static int CountDistinct(BioPolymerGroup group, IEnumerable<IBioPolymerWithSetMods> sequences)
        => Keys(group, sequences).Distinct().Count();

    /// <summary>
    /// Projects sequences onto full or base sequence per <see cref="BioPolymerGroup.DisplayModsOnPeptides"/>.
    /// </summary>
    private static IEnumerable<string> Keys(BioPolymerGroup group, IEnumerable<IBioPolymerWithSetMods> sequences)
        => group.DisplayModsOnPeptides
            ? sequences.Select(p => p.FullSequence)
            : sequences.Select(p => p.BaseSequence);

    /// <summary>
    /// Single-letter classification: entrapment decoy, entrapment target, decoy, contaminant, or target.
    /// </summary>
    private static string TargetDecoyLabel(BioPolymerGroup group)
    {
        if (group.IsEntrapment && group.IsDecoy) return "ED";
        if (group.IsEntrapment) return "ET";
        if (group.IsDecoy) return "D";
        if (group.IsContaminant) return "C";
        return "T";
    }

    private static string Truncate(string? input) => TsvWriter.Truncate(input);
}
