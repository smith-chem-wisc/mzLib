#nullable enable
using MzLibUtil;
using Omics;
using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;

namespace UsefulProteomicsDatabases.EntrapmentGeneration;

/// <summary>
/// A set of modifications a search engine cannot tell apart on a precursor, and the residues they
/// can sit on.
/// </summary>
/// <remarks>
/// <para><b>Why a group and not a modification.</b> The localization candidate set a search actually
/// explores is defined by <b>mass shift</b>, not by modification name. A precursor carries a mass;
/// every modification of that mass competes for it. Measured on a real GPTMD database by the first
/// consumer of this code: grouping by modification description, 32.7% of (peptide, modification)
/// pairs have two or more candidate sites; grouping by <b>mass</b>, 71.3% do. <c>+15.995</c> alone
/// spans seventeen modification IDs -- Hydroxylation on K/N/P and Oxidation on C/D/E/F/H/I/L/Q/R/S/T/V/W/Y.
/// An invariant stated per modification is therefore not the invariant the search engine sees, and
/// the two sentences are easy to confuse: <i>"the companion carries the same modifications"</i> and
/// <i>"the companion has the same mass groups"</i> are different claims.</para>
/// <para><b>Residues, not names.</b> Membership is by <c>(residue, mass)</c>. Names cannot be relied
/// on: the modification dictionary has no residue-specific N-terminal acetyl record, so a wildcard
/// <c>Acetylation on X</c> is kept in real configurations, and being "explicit" would mean
/// hand-authoring modifications -- exactly the sort of free choice this project avoids.</para>
/// </remarks>
public sealed class MassGroup
{
    internal MassGroup(double representativeMass, double lowestMass, double highestMass,
        IReadOnlySet<char> residues, IReadOnlyList<string> modifications)
    {
        RepresentativeMass = representativeMass;
        LowestMass = lowestMass;
        HighestMass = highestMass;
        Residues = residues;
        Modifications = modifications;
    }

    /// <summary>The group's lowest member mass -- a stable label, not an average.</summary>
    /// <remarks>
    /// An average would move when a member is added, so two reports of the same database under
    /// different modification lists would name the same group differently and could not be diffed.
    /// </remarks>
    public double RepresentativeMass { get; }

    public double LowestMass { get; }

    public double HighestMass { get; }

    /// <summary>Residues that some member of this group can occupy.</summary>
    public IReadOnlySet<char> Residues { get; }

    /// <summary>The member modifications, by <c>IdWithMotif</c>, for a reader rather than for logic.</summary>
    public IReadOnlyList<string> Modifications { get; }

    /// <summary>How wide the group actually is, in daltons.</summary>
    public double Span => HighestMass - LowestMass;

    public override string ToString() =>
        RepresentativeMass.ToString("F5", CultureInfo.InvariantCulture)
        + " {" + string.Concat(Residues.OrderBy(c => c)) + "}";
}

/// <summary>
/// The mass groups of a set of modifications, and the counts a companion peptide must match.
/// </summary>
/// <remarks>
/// <para><b>The tolerance is a parameter and deliberately has no clever default.</b> Where the
/// boundary sits is the whole question -- it decides whether near-isobars such as phospho and sulfo
/// (9.5 mDa apart) are one candidate set or two -- and picking it here would be choosing a consumer's
/// analysis for them. Zero groups only exactly-equal masses, which is already the right unit and is
/// the conservative answer.</para>
/// <para><b>Grouping is single-linkage</b>, so a chain of members each within tolerance of the next
/// forms one group even when its ends are further apart than the tolerance. That is the honest
/// reading of "cannot be told apart" transitively, and <see cref="MassGroup.Span"/> is published so a
/// reader can see when it has happened rather than discover it later.</para>
/// <para><b>What this does not model.</b> A modification whose composition is an integer multiple of
/// a smaller one acquires a combinatorial competitor once several modifications are allowed on one
/// peptide -- trimethylation is exactly isobaric with three methylations. That competitor is a
/// property of the search's settings rather than of the database, and composition preservation
/// already guarantees the companion offers the same residues to it, so it is documented here rather
/// than counted.</para>
/// </remarks>
public sealed class MassGroupIndex
{
    private readonly List<MassGroup> _groups;

    /// <param name="modifications">The modifications in play -- normally every one annotated in the
    /// target database.</param>
    /// <param name="toleranceDaltons">Masses within this of one another join the same group.
    /// Zero groups only exactly-equal masses.</param>
    public MassGroupIndex(IEnumerable<Modification> modifications, double toleranceDaltons = 0.0)
    {
        if (modifications is null)
        {
            throw new MzLibException("Mass groups need the modifications they are grouping.");
        }
        if (toleranceDaltons < 0 || double.IsNaN(toleranceDaltons))
        {
            throw new MzLibException(
                $"Mass-group tolerance must be zero or more, but was {toleranceDaltons}.");
        }

        Tolerance = toleranceDaltons;

        // A modification with no mass cannot define a mass group and is dropped rather than treated
        // as zero -- a zero-mass group would swallow every other unmassed modification and read as a
        // real candidate set.
        var usable = modifications
            .Where(m => m is not null && m.MonoisotopicMass.HasValue
                        && double.IsFinite(m.MonoisotopicMass.Value))
            .OrderBy(m => m.MonoisotopicMass!.Value)
            .ToList();

        _groups = new List<MassGroup>();
        var pending = new List<Modification>();
        foreach (Modification m in usable)
        {
            if (pending.Count > 0
                && m.MonoisotopicMass!.Value - pending[^1].MonoisotopicMass!.Value > toleranceDaltons)
            {
                _groups.Add(Close(pending));
                pending.Clear();
            }
            pending.Add(m);
        }
        if (pending.Count > 0)
        {
            _groups.Add(Close(pending));
        }
    }

    public double Tolerance { get; }

    public IReadOnlyList<MassGroup> Groups => _groups;

    private static MassGroup Close(List<Modification> members)
    {
        var residues = new HashSet<char>();
        var names = new List<string>();
        foreach (Modification m in members)
        {
            foreach (char residue in ResiduesOf(m))
            {
                residues.Add(residue);
            }
            string name = m.IdWithMotif ?? m.OriginalId ?? "unnamed";
            if (!names.Contains(name))
            {
                names.Add(name);
            }
        }

        return new MassGroup(members[0].MonoisotopicMass!.Value, members[0].MonoisotopicMass!.Value,
            members[^1].MonoisotopicMass!.Value, residues, names);
    }

    /// <summary>
    /// The residues a modification's motif can sit on -- the capital letters of the motif.
    /// </summary>
    /// <remarks>
    /// A motif is written with the modified residue capitalised and its context in lower case
    /// ("Nxs" for an N-glycosylation sequon), so the capitals are the positions that carry the mass.
    /// A wildcard "X" is expanded to every amino acid rather than kept as a literal, because a
    /// wildcard modification really can sit anywhere and treating "X" as a residue would make a
    /// peptide's capacity for it zero.
    /// </remarks>
    private static IEnumerable<char> ResiduesOf(Modification modification)
    {
        string motif = modification.Target?.ToString() ?? string.Empty;
        foreach (char c in motif)
        {
            if (!char.IsUpper(c))
            {
                continue;
            }
            if (c == 'X')
            {
                foreach (char any in "ACDEFGHIKLMNPQRSTVWY")
                {
                    yield return any;
                }
            }
            else
            {
                yield return c;
            }
        }
    }

    /// <summary>The group holding <paramref name="mass"/>, or null when none does.</summary>
    public MassGroup? GroupFor(double mass)
    {
        foreach (MassGroup g in _groups)
        {
            if (mass >= g.LowestMass - Tolerance && mass <= g.HighestMass + Tolerance)
            {
                return g;
            }
        }
        return null;
    }

    /// <summary>
    /// How many residues of <paramref name="peptide"/> could carry <paramref name="group"/>'s mass.
    /// </summary>
    public static int Capacity(string peptide, MassGroup group)
    {
        if (string.IsNullOrEmpty(peptide) || group is null)
        {
            return 0;
        }

        int n = 0;
        foreach (char c in peptide)
        {
            if (group.Residues.Contains(c))
            {
                n++;
            }
        }
        return n;
    }

    /// <summary>
    /// Annotated sites per mass group on one digestion product, counting only annotations a search
    /// could actually apply there.
    /// </summary>
    /// <param name="polymerSequence">The whole protein the peptide came from.</param>
    /// <param name="peptideStart">Zero-based offset of the peptide within it.</param>
    /// <param name="peptideLength">Length of the peptide.</param>
    /// <param name="annotations">The protein's modifications, keyed one-based.</param>
    /// <remarks>
    /// <b>The fit test is the point of this method.</b> An annotation is counted only where
    /// <c>ModificationLocalization.ModFits</c> accepts it at its position <i>within the peptide</i>,
    /// which is where a peptide-level restriction is finally judged. Counting annotations without it
    /// reproduces the failure this exists to detect: a modification restricted to a peptide terminus
    /// is transported, written, and counted by every protein-level tally, and then quietly does not
    /// exist as a hypothesis because it no longer sits at the terminus.
    /// </remarks>
    public Dictionary<MassGroup, int> AnnotatedSites(string polymerSequence, int peptideStart,
        int peptideLength, IDictionary<int, List<Modification>> annotations)
    {
        var counts = new Dictionary<MassGroup, int>();
        if (annotations is null || string.IsNullOrEmpty(polymerSequence))
        {
            return counts;
        }

        for (int offset = 0; offset < peptideLength; offset++)
        {
            int oneBasedInPolymer = peptideStart + offset + 1;
            if (!annotations.TryGetValue(oneBasedInPolymer, out List<Modification>? mods))
            {
                continue;
            }

            foreach (Modification m in mods)
            {
                if (m?.MonoisotopicMass is null || !double.IsFinite(m.MonoisotopicMass.Value))
                {
                    continue;
                }
                if (!ModificationLocalization.ModFits(m, polymerSequence, offset + 1, peptideLength,
                        oneBasedInPolymer))
                {
                    continue;   // transported, written, and unusable -- exactly what we are counting
                }

                MassGroup? g = GroupFor(m.MonoisotopicMass.Value);
                if (g is null)
                {
                    continue;
                }
                counts[g] = counts.GetValueOrDefault(g) + 1;
            }
        }

        return counts;
    }
}

/// <summary>Target-versus-companion totals for one mass group, over a whole database.</summary>
public sealed class MassGroupTally
{
    internal MassGroupTally(MassGroup group)
    {
        Group = group;
    }

    public MassGroup Group { get; }

    /// <summary>Residues able to carry this mass, summed over target peptides.</summary>
    public long TargetCapacity { get; internal set; }

    /// <summary>The same, over their companions.</summary>
    public long CompanionCapacity { get; internal set; }

    /// <summary>Annotated sites a search could apply, summed over target peptides.</summary>
    public long TargetAnnotatedSites { get; internal set; }

    /// <summary>The same, over their companions.</summary>
    public long CompanionAnnotatedSites { get; internal set; }

    /// <summary>Peptides whose companion offers a different number of usable residues.</summary>
    public int PeptidesWithCapacityMismatch { get; internal set; }

    /// <summary>Peptides whose companion offers a different number of usable annotated sites.</summary>
    public int PeptidesWithAnnotatedSiteMismatch { get; internal set; }

    public bool Holds => TargetCapacity == CompanionCapacity
                         && TargetAnnotatedSites == CompanionAnnotatedSites
                         && PeptidesWithCapacityMismatch == 0
                         && PeptidesWithAnnotatedSiteMismatch == 0;
}

/// <summary>
/// Accumulates the mass-group invariant across a database, peptide by peptide.
/// </summary>
/// <remarks>
/// <para><b>The invariant, in the units a search engine sees:</b> for each (peptide, mass group), the
/// number of residues in the companion that can carry that mass equals the number in its target, and
/// the number of <i>annotated</i> sites for that group is equal too.</para>
/// <para>Under composition-preserving permutation with each modification transported to the residue
/// it was on, this already holds. It is asserted rather than argued because the two clauses fail for
/// different reasons and only the second one has ever failed here: a modification restricted to a
/// peptide terminus survives every protein-level count and then does not fit its new position. That
/// is invisible to a modification tally and visible here.</para>
/// <para><b>Report the distribution, not only the per-peptide equality.</b> A per-peptide check is
/// pointwise and is satisfied vacuously by a peptide that was excised and no longer exists; only the
/// totals catch a whole class going missing. Both are kept, which is why there are four counters per
/// group rather than one flag.</para>
/// </remarks>
public sealed class MassGroupComparison
{
    private readonly MassGroupIndex _index;
    private readonly Dictionary<MassGroup, MassGroupTally> _tallies = new();

    public MassGroupComparison(MassGroupIndex index)
    {
        _index = index ?? throw new MzLibException("A mass-group comparison needs an index.");
        foreach (MassGroup g in index.Groups)
        {
            _tallies[g] = new MassGroupTally(g);
        }
    }

    public int PeptidesCompared { get; private set; }

    public IReadOnlyList<MassGroupTally> Tallies =>
        _tallies.Values.OrderBy(t => t.Group.RepresentativeMass).ToList();

    /// <summary>Whether every group's totals and per-peptide counts agree.</summary>
    public bool Holds => _tallies.Values.All(t => t.Holds);

    /// <summary>
    /// Compares one target protein against its companion, piece by piece.
    /// </summary>
    /// <param name="target">The target polymer, with its annotations.</param>
    /// <param name="companion">Its entrapment partner, with the transported annotations.</param>
    /// <param name="assembly">What happened to each piece.</param>
    /// <param name="minPeptideLength">Pieces shorter than this are skipped: no search reports them on
    /// their own, so counting them would put unidentifiable peptides into an invariant about what a
    /// search can see.</param>
    public void Add(IBioPolymer target, IBioPolymer companion, EntrapmentAssembly assembly,
        int minPeptideLength)
    {
        if (target is null || companion is null || assembly is null)
        {
            throw new MzLibException("A mass-group comparison needs a target, a companion and an assembly.");
        }

        foreach (EntrapmentPiece piece in assembly.Pieces)
        {
            // An excised piece has no companion peptide at all. It is not a mismatch -- there is
            // nothing to mismatch with -- and the achieved fold ratio is where that loss is reported.
            if (piece.Outcome == PieceOutcome.Excised || piece.EntrapmentPiece_ is null)
            {
                continue;
            }
            if (piece.TargetPiece.Length < minPeptideLength)
            {
                continue;
            }

            PeptidesCompared++;

            Dictionary<MassGroup, int> targetSites = _index.AnnotatedSites(target.BaseSequence,
                piece.TargetStart, piece.TargetPiece.Length, target.OneBasedPossibleLocalizedModifications);
            Dictionary<MassGroup, int> companionSites = _index.AnnotatedSites(companion.BaseSequence,
                piece.EntrapmentStart, piece.EntrapmentPiece_.Length,
                companion.OneBasedPossibleLocalizedModifications);

            foreach (MassGroupTally tally in _tallies.Values)
            {
                int targetCapacity = MassGroupIndex.Capacity(piece.TargetPiece, tally.Group);
                int companionCapacity = MassGroupIndex.Capacity(piece.EntrapmentPiece_, tally.Group);
                int targetAnnotated = targetSites.GetValueOrDefault(tally.Group);
                int companionAnnotated = companionSites.GetValueOrDefault(tally.Group);

                tally.TargetCapacity += targetCapacity;
                tally.CompanionCapacity += companionCapacity;
                tally.TargetAnnotatedSites += targetAnnotated;
                tally.CompanionAnnotatedSites += companionAnnotated;
                if (targetCapacity != companionCapacity)
                {
                    tally.PeptidesWithCapacityMismatch++;
                }
                if (targetAnnotated != companionAnnotated)
                {
                    tally.PeptidesWithAnnotatedSiteMismatch++;
                }
            }
        }
    }

    /// <summary>
    /// The comparison as a tab-separated table, one row per mass group plus a total.
    /// </summary>
    /// <remarks>
    /// Every group is emitted, including those with no annotation anywhere, because "this mass had no
    /// annotated site on either side" and "this mass was not looked at" are different statements and
    /// only the table can distinguish them.
    /// </remarks>
    public string ToTabSeparated()
    {
        var text = new StringBuilder();
        text.AppendLine("# massGroupToleranceDaltons\t"
                        + _index.Tolerance.ToString("G17", CultureInfo.InvariantCulture));
        text.AppendLine("# peptidesCompared\t" + PeptidesCompared.ToString(CultureInfo.InvariantCulture));
        text.AppendLine("# invariantHolds\t" + (Holds ? "true" : "false"));
        text.AppendLine(string.Join("\t", "massGroup", "spanDaltons", "residues", "modifications",
            "targetCapacity", "companionCapacity", "targetAnnotatedSites", "companionAnnotatedSites",
            "peptidesWithCapacityMismatch", "peptidesWithAnnotatedSiteMismatch"));

        foreach (MassGroupTally t in Tallies)
        {
            text.AppendLine(string.Join("\t",
                t.Group.RepresentativeMass.ToString("F5", CultureInfo.InvariantCulture),
                t.Group.Span.ToString("F5", CultureInfo.InvariantCulture),
                string.Concat(t.Group.Residues.OrderBy(c => c)),
                string.Join("; ", t.Group.Modifications),
                t.TargetCapacity.ToString(CultureInfo.InvariantCulture),
                t.CompanionCapacity.ToString(CultureInfo.InvariantCulture),
                t.TargetAnnotatedSites.ToString(CultureInfo.InvariantCulture),
                t.CompanionAnnotatedSites.ToString(CultureInfo.InvariantCulture),
                t.PeptidesWithCapacityMismatch.ToString(CultureInfo.InvariantCulture),
                t.PeptidesWithAnnotatedSiteMismatch.ToString(CultureInfo.InvariantCulture)));
        }

        return text.ToString();
    }
}
