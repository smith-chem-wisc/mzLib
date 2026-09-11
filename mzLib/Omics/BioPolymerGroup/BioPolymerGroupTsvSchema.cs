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
    /// Maximum length for string fields in output. Strings exceeding this length are truncated.
    ///
    /// Default is 32,000 characters, slightly below Excel's cell limit of 32,767, so output files
    /// open in Excel without truncation or corruption. Set to 0 or negative to disable truncation
    /// (useful for programmatic processing where Excel compatibility does not matter).
    /// </summary>
    /// <remarks>
    /// Excel specification: a cell can contain up to 32,767 characters.
    /// See: https://support.microsoft.com/en-us/office/excel-specifications-and-limits
    /// </remarks>
    public static int MaxStringLength { get; set; } = 32000;

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
    /// The quantification block, emitted as two declared lists rather than one.
    ///
    /// A spectral count and its count-based occupancy describe an acquired FILE; an intensity and
    /// its intensity-based occupancy describe a SAMPLE — which for an isobaric run is a channel,
    /// and one file carries many. Each sample group result declares which count section it belongs
    /// to (<see cref="SampleGroupResult.CountIdentity"/>), so the counting columns are emitted once
    /// per section and the intensity columns once per sample within it, instead of every column
    /// family being replicated across one shared key space.
    ///
    /// For a label-free design and for no design at all the two spaces coincide — a sample group
    /// is its own count section — so each section holds exactly one sample and the emitted columns
    /// are identical to before, name for name and in the same order. Only an isobaric table changes,
    /// and it changes by losing the eleven identical <c>SpectralCount_</c> values and eleven
    /// byte-identical <c>CountOccupancy_</c> strings an 11-plex row used to carry.
    /// </summary>
    private static IEnumerable<TsvColumn<BioPolymerGroup>> QuantificationColumns(
        IReadOnlyCollection<BioPolymerGroup> groups)
    {
        foreach (var group in groups)
        {
            if (group.SampleGroupResults is null)
                group.PopulateSampleGroupResults();
        }

        // Union of the dataset's count sections and of its sample groups, each in first-seen order.
        //
        // Both are matched across records by identity — never by Label and never by position.
        // Labels are not unique (SampleGroupBuilder names a sample group after its first file
        // whenever conditions are undefined or a file is missing from disk, so same-named files in
        // different directories collide), and position is not stable across records because a record
        // only has sample groups for the files it appears in.
        var countIdentities = new List<string>();
        var countLabelByIdentity = new Dictionary<string, (string Label, string? LabelSourcePath)>();
        var samplesInCountSection = new Dictionary<string, List<string>>();

        var identities = new List<string>();
        var seen = new HashSet<string>();
        var labelByIdentity = new Dictionary<string, (string Label, string? LabelSourcePath)>();

        foreach (var group in groups)
        {
            foreach (var result in group.SampleGroupResults!)
            {
                // A null CountIdentity means the result is its own count section, which is every
                // design but the isobaric one. Resolved here rather than behind the property so the
                // absent case reads the same way LabelSourcePath's does.
                string countIdentity = result.CountIdentity ?? result.Identity;

                if (!countLabelByIdentity.ContainsKey(countIdentity))
                {
                    countIdentities.Add(countIdentity);
                    countLabelByIdentity[countIdentity] = (result.CountLabel ?? result.Label, result.LabelSourcePath);
                    samplesInCountSection[countIdentity] = [];
                }

                if (seen.Add(result.Identity))
                {
                    identities.Add(result.Identity);
                    labelByIdentity[result.Identity] = (result.Label, result.LabelSourcePath);
                    samplesInCountSection[countIdentity].Add(result.Identity);
                }
            }
        }

        // Whether the intensity columns exist is a property of the run, not of a sample group:
        // an engine either populated these groups or it did not. Deliberately not
        // SampleGroupResult.HasIntensityData, which answers "did this sample group get a value" --
        // using that to choose columns made a group whose samples were all unobserved describe a
        // narrower table than its neighbour. Same predicate the group applied before rendering moved
        // out; both members are public, so the schema computes it rather than needing access.
        bool quantified = groups.Any(g => g.SamplesForQuantification is { Count: > 0 }
                                       && g.IntensitiesBySample is not null);

        // Disambiguated independently, because a count section and a sample group are different
        // things with different names: two files called sample.raw in different directories collide
        // as count sections even though their channel labels already separate them as samples.
        var countDisplayLabels = SampleGroupLabels.Disambiguate(countIdentities, countLabelByIdentity);
        var displayLabels = SampleGroupLabels.Disambiguate(identities, labelByIdentity);

        var index = new SampleGroupIndex();
        var columns = new List<TsvColumn<BioPolymerGroup>>();

        foreach (var countIdentity in countIdentities)
        {
            string thisCountIdentity = countIdentity;
            string countDisplay = countDisplayLabels[thisCountIdentity];
            var samples = samplesInCountSection[thisCountIdentity];

            columns.Add(new TsvColumn<BioPolymerGroup>($"SpectralCount_{countDisplay}",
                g => index.ResultForCountSection(g, thisCountIdentity)?.SpectralCount.ToString() ?? string.Empty));

            if (quantified)
                foreach (var identity in samples)
                {
                    string thisIdentity = identity;
                    columns.Add(new TsvColumn<BioPolymerGroup>($"Intensity_{displayLabels[thisIdentity]}",
                        g => index.ResultFor(g, thisIdentity) is { HasIntensityData: true } r ? r.Intensity.ToString() : string.Empty));
                }

            columns.Add(new TsvColumn<BioPolymerGroup>($"CountOccupancy_{countDisplay}",
                g => Truncate(index.ResultForCountSection(g, thisCountIdentity)?.FormatOccupancy(OccupancyKeys(g), IsParentLevel(g), intensityBased: false))));

            if (quantified)
                foreach (var identity in samples)
                {
                    string thisIdentity = identity;
                    columns.Add(new TsvColumn<BioPolymerGroup>($"IntensityOccupancy_{displayLabels[thisIdentity]}",
                        g => Truncate(index.ResultFor(g, thisIdentity)?.FormatOccupancy(OccupancyKeys(g), IsParentLevel(g), intensityBased: true))));
                }
        }

        return columns;
    }

    /// <summary>
    /// Finds a group's sample group by identity, in place of scanning its list once per cell.
    ///
    /// Decoupling a column from its position in the group's list is what fixes the ragged-row
    /// defect, but it also means a column can no longer read its value off the list in order. Left
    /// as a scan, each of the ~4 columns per sample group re-walks the whole list, so writing costs
    /// O(samples² × records) where rendering a row used to cost O(samples).
    ///
    /// A row is rendered one group at a time across every column, so remembering the most recent
    /// group's lookup restores O(samples × records).
    ///
    /// Keyed on the sample-group list instance rather than on the group: every invalidation path
    /// nulls SampleGroupResults, and repopulating assigns a fresh list, so a group whose
    /// quantification changed mid-write rebuilds its index instead of serving a stale one.
    /// Not thread-safe, in common with the lazy population it wraps.
    /// </summary>
    private sealed class SampleGroupIndex
    {
        private List<SampleGroupResult>? _source;
        private Dictionary<string, SampleGroupResult> _byIdentity = [];
        private Dictionary<string, SampleGroupResult> _byCountIdentity = [];

        public SampleGroupResult? ResultFor(BioPolymerGroup group, string identity)
            => Index(group)._byIdentity.GetValueOrDefault(identity);

        /// <summary>
        /// Any one result from a count section. Every result in a section reports the same
        /// <see cref="SampleGroupResult.SpectralCount"/> and the same count-based occupancy by
        /// construction — <see cref="SampleGroupBuilder"/> derives both from one file's PSM list and
        /// hands that same list to each of the file's channels — so which one answers does not
        /// matter, only that the section is reported once rather than once per channel.
        /// </summary>
        public SampleGroupResult? ResultForCountSection(BioPolymerGroup group, string countIdentity)
            => Index(group)._byCountIdentity.GetValueOrDefault(countIdentity);

        private SampleGroupIndex Index(BioPolymerGroup group)
        {
            if (group.SampleGroupResults is null)
                group.PopulateSampleGroupResults();

            var results = group.SampleGroupResults!;

            if (!ReferenceEquals(results, _source))
            {
                _source = results;
                _byIdentity = new Dictionary<string, SampleGroupResult>(results.Count);
                _byCountIdentity = new Dictionary<string, SampleGroupResult>(results.Count);

                // First wins, matching the FirstOrDefault this replaces — identities are unique
                // within a group by construction, but TryAdd keeps a hand-built list from throwing.
                // Count identities are deliberately NOT unique: a plex's channels share one.
                foreach (var result in results)
                {
                    _byIdentity.TryAdd(result.Identity, result);
                    _byCountIdentity.TryAdd(result.CountIdentity ?? result.Identity, result);
                }
            }

            return this;
        }
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

    /// <summary>
    /// Truncates to <see cref="MaxStringLength"/> so output stays within Excel's cell limit.
    /// </summary>
    private static string Truncate(string? input)
    {
        if (string.IsNullOrEmpty(input))
            return string.Empty;

        if (MaxStringLength <= 0 || input.Length <= MaxStringLength)
            return input;

        return input.Substring(0, MaxStringLength);
    }
}
