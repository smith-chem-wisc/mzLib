using System;
using System.Collections.Generic;
using System.Linq;
using Omics;

namespace Quantification.Differential;

/// <summary>
/// The UniProt primary gene names of every entry in a search's protein databases (<c>DEF-DIFF-SET-MEMBER</c>): every
/// <c>&lt;gene&gt;</c> element's primary name, so a multi-gene entry belongs to each gene it names. Synonyms and the
/// protein table's <c>Gene</c> cell are never used. Genes match case-insensitively.
/// </summary>
public sealed class SetGeneIndex
{
    private readonly IReadOnlyDictionary<string, IReadOnlyList<string>> _genesByAccession;
    private readonly IReadOnlyDictionary<string, IReadOnlyList<string>> _accessionsByGene;

    private SetGeneIndex(IReadOnlyDictionary<string, IReadOnlyList<string>> genesByAccession)
    {
        _genesByAccession = genesByAccession;
        _accessionsByGene = genesByAccession
            .SelectMany(e => e.Value.Select(g => (Gene: g, Accession: e.Key)))
            .GroupBy(x => x.Gene, StringComparer.OrdinalIgnoreCase)
            .ToDictionary(g => g.Key, g => (IReadOnlyList<string>)g.Select(x => x.Accession)
                .Distinct(StringComparer.Ordinal).Order(StringComparer.Ordinal).ToList(), StringComparer.OrdinalIgnoreCase);
    }

    /// <summary>Every accession in the databases, with or without a gene name.</summary>
    public IEnumerable<string> Accessions => _genesByAccession.Keys;

    /// <summary>
    /// The index of a search's database entries. An entry's primary names are its <see cref="IBioPolymer.GeneNames"/>
    /// whose type is <c>primary</c>. An accession seen in two databases keeps the names of both.
    /// </summary>
    public static SetGeneIndex FromBioPolymers(IEnumerable<IBioPolymer> entries)
    {
        ArgumentNullException.ThrowIfNull(entries);
        return FromPrimaryNames(entries.Select(e => (e.Accession,
            (IEnumerable<string>)(e.GeneNames ?? new List<Tuple<string, string>>())
                .Where(g => string.Equals(g.Item1, "primary", StringComparison.OrdinalIgnoreCase))
                .Select(g => g.Item2))));
    }

    /// <summary>The index from (accession, primary names) pairs; an entry with no name is given with none.</summary>
    public static SetGeneIndex FromPrimaryNames(IEnumerable<(string Accession, IEnumerable<string> PrimaryNames)> entries)
    {
        ArgumentNullException.ThrowIfNull(entries);
        var map = new Dictionary<string, HashSet<string>>(StringComparer.Ordinal);
        foreach (var (accession, names) in entries)
        {
            if (string.IsNullOrEmpty(accession))
                throw new ArgumentException("An entry has no accession.", nameof(entries));
            if (!map.TryGetValue(accession, out var set))
                map[accession] = set = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
            foreach (string name in names ?? Enumerable.Empty<string>())
                if (!string.IsNullOrWhiteSpace(name))
                    set.Add(name.Trim());
        }
        return new SetGeneIndex(map.ToDictionary(e => e.Key,
            e => (IReadOnlyList<string>)e.Value.Order(StringComparer.OrdinalIgnoreCase).ToList(), StringComparer.Ordinal));
    }

    /// <summary>The databases hold this accession.</summary>
    public bool Contains(string accession) => _genesByAccession.ContainsKey(accession);

    /// <summary>The accession's primary gene names; empty when it has none or is not in the databases.</summary>
    public IReadOnlyList<string> GenesOf(string accession) =>
        _genesByAccession.TryGetValue(accession, out var genes) ? genes : Array.Empty<string>();

    /// <summary>The accessions whose primary names include <paramref name="gene"/> (any case), ordinal-sorted.</summary>
    public IReadOnlyList<string> AccessionsOf(string gene) =>
        _accessionsByGene.TryGetValue(gene, out var accessions) ? accessions : Array.Empty<string>();
}

/// <summary>The <c>member_state</c> values decided from the databases and the protein groups (<c>DEF-DIFF-SET-MEMBER</c>).</summary>
public static class SetMemberState
{
    /// <summary>The member matches no gene name (or no accession) in the search's databases.</summary>
    public const string Unmapped = "unmapped";
    /// <summary>No protein group of the search holds the member's gene.</summary>
    public const string NoGroup = "no_group";
    /// <summary>Its only groups are outside STAT1's row set; followed by the reason.</summary>
    public const string NotInRowSetPrefix = "not_in_rowset:";
    /// <summary>The reason: a contaminant group.</summary>
    public const string ContaminantReason = "contaminant";
    /// <summary>The reason: protein q above 0.01.</summary>
    public const string AboveOnePercentQReason = "above_1pct_q";
    /// <summary>
    /// The reason: the group passes the protein filter but none of its peptides has a value in the stratum, so it is not
    /// in the stratum's table (<c>DEF-DIFF-ROWSET</c>). <b>Pending the user's ruling</b>: the contract lists only the two
    /// reasons above.
    /// </summary>
    public const string NoValueInStratumReason = "no_value_in_stratum";
    /// <summary>Its only groups also hold other genes, and the all-members exception does not apply.</summary>
    public const string SharedGroup = "shared_group";
}

/// <summary>Which protein group stands for one member in one dataset contrast, or why none does.</summary>
public sealed record SetMemberAssignment
{
    /// <summary>The member, as frozen.</summary>
    public required string Member { get; init; }
    /// <summary>
    /// A <see cref="SetMemberState"/> value when no group stands for the member; null when one does, and the member's
    /// state is then decided from STAT1's row for that group (<c>too_few</c>, <c>seen_one_side:*</c>, <c>incomplete</c>
    /// or <c>tested</c>).
    /// </summary>
    public string? State { get; init; }
    /// <summary>The chosen group's name, or null.</summary>
    public string? ProteinGroup { get; init; }
    /// <summary>Every gene of the chosen group, ordinal-sorted (case-insensitive order); empty without a group.</summary>
    public IReadOnlyList<string> GroupGenes { get; init; } = Array.Empty<string>();
    /// <summary>The other groups that could have stood for the member, not chosen; ordinal-sorted.</summary>
    public IReadOnlyList<string> OtherGroups { get; init; } = Array.Empty<string>();
    /// <summary>The accessions that define the member: its <c>member_accessions</c>, or those of its gene; ordinal-sorted.</summary>
    public IReadOnlyList<string> MemberAccessionsUsed { get; init; } = Array.Empty<string>();
    /// <summary>The group came in through the all-members exception, labelled with each of its genes.</summary>
    public bool ByAllMembersException { get; init; }
}

/// <summary>
/// Finds the protein group that stands for each member of one dataset contrast (<c>DEF-DIFF-SET-MEMBER</c>, up to
/// <c>shared_group</c>; the states read from STAT1's rows come after).
/// </summary>
/// <remarks>
/// <para><b>Which groups can stand for a member:</b> only groups in the stratum's row set (<c>DEF-DIFF-ROWSET</c>: target,
/// not a contaminant, protein q at most 0.01, and a value for one of its peptides in a sample of the stratum; decoys never),
/// and only a group every accession of which belongs to the member's gene and to no other: by gene name, the group's genes
/// are exactly the member's gene, and an accession with no gene name counts as another gene; with <c>member_accessions</c>,
/// every accession is one of them. So a multi-gene entry (one histone H4 entry encoded by many genes) stands for no single
/// member, and enters only through the exception.</para>
/// <para><b>The all-members exception:</b> a group whose genes are all members of the part enters once, labelled with each
/// of them. <b>Pending the user's ruling</b>, it enters only for members with no group of their own, and only when none of
/// its genes has a group of its own; otherwise a member left with it is <see cref="SetMemberState.SharedGroup"/>. So no
/// gene is counted twice.</para>
/// <para><b>One group per member, chosen without looking at any effect:</b> the most samples of the contrast in which one
/// of the group's unique peptides (GR-23) has a value; then the most such unique peptides with a value in one of the
/// contrast's samples; then the group name in ordinal order.</para>
/// </remarks>
public static class SetMemberResolver
{
    /// <summary>The protein q-value at or below which a target group is in the row set.</summary>
    public const double RowSetQValue = 0.01;

    private sealed record GroupFacts(string Name, IReadOnlyList<string> Accessions, IReadOnlyList<string> Genes,
        bool HasGenelessAccession, string? OutsideReason, int SamplesWithValue, int UniquePeptidesWithValue);

    /// <summary>Resolves every member of <paramref name="contrast"/>.</summary>
    /// <param name="contrast">One dataset contrast of a part.</param>
    /// <param name="partMembers">All members of the part, for the all-members exception.</param>
    /// <param name="stratum">The stratum's observation table (<see cref="ObservationTable.ForSamples"/>).</param>
    /// <param name="contrastSamples">The samples on either side of the contrast, after any <c>contexts</c> restriction.</param>
    /// <param name="genes">The search's databases.</param>
    /// <returns>One assignment per member of the contrast, in the contrast's order.</returns>
    public static IReadOnlyList<SetMemberAssignment> Resolve(SetDatasetContrast contrast, IReadOnlyList<string> partMembers,
        ObservationTable stratum, IReadOnlyCollection<string> contrastSamples, SetGeneIndex genes)
    {
        ArgumentNullException.ThrowIfNull(contrast);
        ArgumentNullException.ThrowIfNull(partMembers);
        ArgumentNullException.ThrowIfNull(stratum);
        ArgumentNullException.ThrowIfNull(contrastSamples);
        ArgumentNullException.ThrowIfNull(genes);

        var sampleIndex = new List<int>();
        foreach (string id in contrastSamples.Distinct(StringComparer.Ordinal))
        {
            int s = stratum.Samples.ToList().IndexOf(id);
            if (s < 0)
                throw new ArgumentException($"The stratum's table has no sample '{id}'.", nameof(contrastSamples));
            sampleIndex.Add(s);
        }

        var peptidesOf = new Dictionary<string, List<int>>(StringComparer.Ordinal);
        for (int p = 0; p < stratum.Peptides.Count; p++)
            foreach (string g in stratum.Peptides[p].ProteinGroups)
            {
                if (!peptidesOf.TryGetValue(g, out var list))
                    peptidesOf[g] = list = new List<int>();
                list.Add(p);
            }

        // Decoy groups never stand for a member and never hold its gene.
        var facts = stratum.ProteinGroups.Values
            .Where(g => !g.IsDecoy)
            .Select(g => Describe(g, peptidesOf.TryGetValue(g.Name, out var list) ? list : new List<int>(), stratum,
                sampleIndex, genes))
            .ToDictionary(f => f.Name, StringComparer.Ordinal);
        var members = contrast.Rows.GroupBy(r => r.Member, StringComparer.Ordinal).Select(g => g.First()).ToList();
        var partGenes = new HashSet<string>(partMembers, StringComparer.OrdinalIgnoreCase);

        // First pass: each member's own groups. By gene name, a group is the member's own only when its genes are exactly
        // the member's gene and no accession is gene-less; with member_accessions, when every accession is one of them.
        var own = new Dictionary<string, List<GroupFacts>>(StringComparer.Ordinal);
        var accessionsOf = new Dictionary<string, IReadOnlyList<string>>(StringComparer.Ordinal);
        foreach (var row in members)
        {
            var accessions = MemberAccessions(row, genes);
            accessionsOf[row.Member] = accessions;
            var mine = new HashSet<string>(accessions, StringComparer.Ordinal);
            bool byAccession = row.MemberAccessions.Count > 0;
            own[row.Member] = facts.Values
                .Where(f => f.OutsideReason is null && (byAccession
                    ? f.Accessions.All(mine.Contains)
                    : !f.HasGenelessAccession && f.Genes.Count == 1
                      && string.Equals(f.Genes[0], row.Member, StringComparison.OrdinalIgnoreCase)))
                .ToList();
        }
        var gotOwn = new HashSet<string>(own.Where(o => o.Value.Count > 0).Select(o => o.Key), StringComparer.OrdinalIgnoreCase);

        var result = new List<SetMemberAssignment>();
        foreach (var row in members)
        {
            var accessions = accessionsOf[row.Member];
            var assignment = new SetMemberAssignment { Member = row.Member, MemberAccessionsUsed = accessions };
            if (accessions.Count == 0 || !accessions.Any(genes.Contains))
            {
                result.Add(assignment with { State = SetMemberState.Unmapped });
                continue;
            }

            var mine = new HashSet<string>(accessions, StringComparer.Ordinal);
            var holding = facts.Values.Where(f => f.Accessions.Any(mine.Contains)).ToList();
            if (holding.Count == 0)
            {
                result.Add(assignment with { State = SetMemberState.NoGroup });
                continue;
            }
            var inRowSet = holding.Where(f => f.OutsideReason is null).ToList();
            if (inRowSet.Count == 0)
            {
                result.Add(assignment with { State = SetMemberState.NotInRowSetPrefix + ClosestReason(holding) });
                continue;
            }

            var candidates = own[row.Member];
            bool byException = false;
            if (candidates.Count == 0)
            {
                candidates = inRowSet
                    .Where(f => !f.HasGenelessAccession && f.Genes.Count > 1 && f.Genes.All(partGenes.Contains)
                                && !f.Genes.Any(gotOwn.Contains))
                    .ToList();
                byException = candidates.Count > 0;
            }
            if (candidates.Count == 0)
            {
                result.Add(assignment with { State = SetMemberState.SharedGroup });
                continue;
            }

            var ordered = candidates
                .OrderByDescending(f => f.SamplesWithValue)
                .ThenByDescending(f => f.UniquePeptidesWithValue)
                .ThenBy(f => f.Name, StringComparer.Ordinal)
                .ToList();
            result.Add(assignment with
            {
                ProteinGroup = ordered[0].Name,
                GroupGenes = ordered[0].Genes,
                OtherGroups = ordered.Skip(1).Select(f => f.Name).Order(StringComparer.Ordinal).ToList(),
                ByAllMembersException = byException,
            });
        }
        return result;
    }

    private static IReadOnlyList<string> MemberAccessions(SetTestRow row, SetGeneIndex genes) =>
        row.MemberAccessions.Count > 0
            ? row.MemberAccessions.Distinct(StringComparer.Ordinal).Order(StringComparer.Ordinal).ToList()
            : genes.AccessionsOf(row.Member);

    private static GroupFacts Describe(ProteinGroupInfo group, IReadOnlyList<int> peptides, ObservationTable stratum,
        IReadOnlyList<int> contrastSamples, SetGeneIndex genes)
    {
        var accessions = group.Accessions;
        var groupGenes = accessions.SelectMany(genes.GenesOf).Distinct(StringComparer.OrdinalIgnoreCase)
            .Order(StringComparer.OrdinalIgnoreCase).ToList();
        bool geneless = accessions.Any(a => genes.GenesOf(a).Count == 0);

        bool valueInStratum = peptides.Any(p => Enumerable.Range(0, stratum.Samples.Count)
            .Any(s => double.IsFinite(stratum.Log2Intensity(p, s))));

        string? outside = group.IsContaminant ? SetMemberState.ContaminantReason
            : !(group.QValue <= RowSetQValue) ? SetMemberState.AboveOnePercentQReason
            : !valueInStratum ? SetMemberState.NoValueInStratumReason
            : null;

        var unique = peptides.Where(p => string.Equals(stratum.Peptides[p].UniqueProteinGroup, group.Name, StringComparison.Ordinal))
            .ToList();
        int samplesWithValue = contrastSamples.Count(s => unique.Any(p => double.IsFinite(stratum.Log2Intensity(p, s))));
        int uniqueWithValue = unique.Count(p => contrastSamples.Any(s => double.IsFinite(stratum.Log2Intensity(p, s))));
        return new GroupFacts(group.Name, accessions, groupGenes, geneless, outside, samplesWithValue, uniqueWithValue);
    }

    /// <summary>
    /// The reason of the member's group that came closest to the row set: no value in the stratum, then above 1% q, then
    /// contaminant.
    /// </summary>
    private static string ClosestReason(IEnumerable<GroupFacts> groups)
    {
        var reasons = groups.Select(g => g.OutsideReason!).ToHashSet(StringComparer.Ordinal);
        return reasons.Contains(SetMemberState.NoValueInStratumReason) ? SetMemberState.NoValueInStratumReason
            : reasons.Contains(SetMemberState.AboveOnePercentQReason) ? SetMemberState.AboveOnePercentQReason
            : SetMemberState.ContaminantReason;
    }
}
