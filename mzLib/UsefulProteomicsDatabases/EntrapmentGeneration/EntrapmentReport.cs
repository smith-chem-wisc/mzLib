#nullable enable
using MzLibUtil;
using Omics.Digestion;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;

namespace UsefulProteomicsDatabases.EntrapmentGeneration;

/// <summary>How the database was made, recorded alongside what came out of it.</summary>
/// <remarks>
/// Without these an entrapment database cannot be regenerated or compared against another, and the
/// numbers below cannot be interpreted -- the peptide population, and therefore every distribution
/// in the report, depends on the enzyme and the digestion settings as much as on the method.
/// </remarks>
public sealed class EntrapmentProvenance
{
    internal EntrapmentProvenance(string method, string enzyme, int seed, int foldCount,
        int maxMissedCleavages, int minPeptideLength, int maxPeptideLength)
    {
        Method = method;
        Enzyme = enzyme;
        Seed = seed;
        FoldCount = foldCount;
        MaxMissedCleavages = maxMissedCleavages;
        MinPeptideLength = minPeptideLength;
        MaxPeptideLength = maxPeptideLength;
    }

    public string Method { get; }
    public string Enzyme { get; }
    public int Seed { get; }
    public int FoldCount { get; }
    public int MaxMissedCleavages { get; }
    public int MinPeptideLength { get; }
    public int MaxPeptideLength { get; }
}

/// <summary>Counts for one stratum -- one value of whatever the report is stratified by.</summary>
public sealed class EntrapmentStratum
{
    internal EntrapmentStratum(int siteCount) => SiteCount = siteCount;

    /// <summary>The value being stratified by, e.g. the number of candidate sites in a peptide.</summary>
    public int SiteCount { get; }

    /// <summary>Distinct target peptides in this stratum, counted once however many folds ran.</summary>
    /// <remarks>
    /// <b>Base pieces</b> -- the unit the generator permutes -- not what a search reports. For that,
    /// see <see cref="SearchSpacePeptides"/>. The two differ by roughly five-fold at two missed
    /// cleavages, so a rate computed against the wrong one is meaningless.
    /// </remarks>
    public int TargetPeptides { get; internal set; }

    /// <summary>
    /// Distinct peptides in this stratum a search could report: runs of up to
    /// <c>MaxMissedCleavages + 1</c> base pieces, within the length bounds.
    /// </summary>
    /// <remarks>
    /// This is the population an FDP estimator's <c>r</c> is over, and the denominator
    /// <see cref="Ambiguous"/> belongs to -- <see cref="Ambiguous"/> is counted over missed-cleavage
    /// peptides too, so dividing it by <see cref="TargetPeptides"/> divides two different
    /// populations. On the reviewed human proteome that mistake turned 0.24% into "1.19%".
    /// </remarks>
    public int SearchSpacePeptides { get; internal set; }

    /// <summary>Entrapment peptides produced for them, across every fold.</summary>
    public int EntrapmentPeptides { get; internal set; }

    /// <summary>Exactly one arrangement exists. Arithmetic: no seed and no fold count can help.</summary>
    public int UnpairableNoPermutationExists { get; internal set; }

    /// <summary>Arrangements exist but every one available to the fold was already spoken for.</summary>
    public int UnpairableAllPermutationsTaken { get; internal set; }

    /// <summary>The space is real but too small to give every fold its own partner.</summary>
    public int UnpairableSpaceTooSmallForFoldCount { get; internal set; }

    /// <summary>
    /// Target peptides sharing a composition-and-pinning key with another peptide of the same
    /// protein, so a discovery cannot be traced back to one of them.
    /// </summary>
    public int Ambiguous { get; internal set; }

    /// <summary>
    /// Entrapment missed-cleavage peptides whose pieces were not adjacent in the target, so they
    /// have no isomeric counterpart. The only peptides in the database that break that invariant.
    /// </summary>
    public int MissedCleavagePeptidesSpanningAnExcision { get; internal set; }

    /// <summary>
    /// Missed-cleavage peptides of the entrapment database that equal a real target peptide and
    /// could not be permuted away from it, because the run's final piece has only one arrangement.
    /// </summary>
    /// <remarks>
    /// A search matching one of these counts a <i>true</i> peptide as an entrapment discovery,
    /// inflating an FDP estimate directly. They are counted rather than repaired -- repair would
    /// mean backtracking into a placed piece or excising a good one -- so a consumer that needs an
    /// exact count must exclude them. See <see cref="EntrapmentAssembly.UnrepairableRunCollisions"/>.
    /// </remarks>
    public int UnrepairableRunCollisions { get; internal set; }

    public int Unpairable => UnpairableNoPermutationExists + UnpairableAllPermutationsTaken
                             + UnpairableSpaceTooSmallForFoldCount;

    /// <summary>
    /// Partners actually delivered per target peptide. Compare against the requested fold count:
    /// the shortfall is where the database could not honour what was asked of it.
    /// </summary>
    public double AchievedFoldRatio =>
        TargetPeptides == 0 ? 0d : EntrapmentPeptides / (double)TargetPeptides;
}

/// <summary>
/// What an entrapment database actually contains, stratified, next to how it was made.
/// </summary>
/// <remarks>
/// <para>The report emits the <b>distribution</b>, not a summary of it. There is no universal
/// distribution of peptides to validate against -- it moves with the proteome, the enzyme and the
/// digestion settings -- so only this database's own numbers mean anything, and a single headline
/// figure would hide exactly the stratum-by-stratum differences worth looking at.</para>
/// <para>Failures are kept apart by cause as well as by stratum, because the remedies differ:
/// nothing helps a peptide with one arrangement, a smaller fold count helps one whose space is too
/// small, and a different target database may help one whose space is merely crowded.</para>
/// </remarks>
public sealed class EntrapmentReport
{
    internal EntrapmentReport(EntrapmentProvenance provenance, IReadOnlyList<EntrapmentStratum> strata,
        EntrapmentStratum total,
        IReadOnlyDictionary<string, IReadOnlyCollection<string>> ambiguousByAccession,
        IReadOnlyDictionary<string, IReadOnlyCollection<string>> unrepairableByAccession)
    {
        Provenance = provenance;
        Strata = strata;
        Total = total;
        AmbiguousPeptidesByAccession = ambiguousByAccession;
        UnrepairableRunCollisionsByAccession = unrepairableByAccession;
    }

    /// <summary>
    /// Target peptides that cannot be traced back to one target, by accession -- two peptides of the
    /// same protein sharing a composition-and-pinning key.
    /// </summary>
    /// <remarks>
    /// A consumer computing the paired FDP estimator has to exclude the <i>same</i> peptides this
    /// generator excluded, or its <c>r = 1</c> assumption fails silently -- no error, just a wrong
    /// number. That is only possible if it can see which ones they are, so the list is the
    /// deliverable and the count is the summary.
    /// </remarks>
    public IReadOnlyDictionary<string, IReadOnlyCollection<string>> AmbiguousPeptidesByAccession { get; }

    /// <summary>
    /// Entrapment missed-cleavage peptides that are also real target peptides, by target accession.
    /// </summary>
    /// <remarks>
    /// A search matching one of these counts a <i>true</i> peptide as an entrapment discovery. They
    /// could not be permuted away because the run's final piece has only one arrangement, and this
    /// project counts collisions rather than repairing them by backtracking or excising a good
    /// piece. Excluding them is the consumer's call, and needs the sequences.
    /// </remarks>
    public IReadOnlyDictionary<string, IReadOnlyCollection<string>> UnrepairableRunCollisionsByAccession { get; }

    /// <summary>
    /// Peptides of the foreign-species arm that are also target peptides, by the foreign protein's
    /// own accession.
    /// </summary>
    /// <remarks>
    /// Homology, not a defect. A conserved protein shares peptides with its ortholog, and a shared
    /// peptide is a real target peptide sitting in the entrapment database. Nothing can be permuted
    /// away -- the foreign sequence is what it is -- so the arm reports them and the caller decides
    /// whether to exclude the peptides or drop the proteins.
    /// </remarks>
    public IReadOnlyDictionary<string, IReadOnlyCollection<string>> ForeignPeptidesSharedWithTarget { get; internal set; }
        = new Dictionary<string, IReadOnlyCollection<string>>();

    /// <summary>Entries contributed by the foreign-species arm.</summary>
    public int ForeignEntries { get; internal set; }

    /// <summary>
    /// The two exclusion lists as a tab-separated table: <c>accession, peptide, reason</c>.
    /// </summary>
    /// <remarks>
    /// A sidecar rather than more columns on the stratified table, because these are per-peptide
    /// facts and that table is per-stratum. Empty but for its header when nothing is excluded, which
    /// is a meaningful answer rather than a missing file.
    /// </remarks>
    public string ExclusionsToTabSeparated()
    {
        var text = new StringBuilder();
        text.AppendLine(string.Join("\t", "accession", "peptide", "reason"));

        foreach ((string accession, IReadOnlyCollection<string> peptides) in
                 AmbiguousPeptidesByAccession.OrderBy(kv => kv.Key, StringComparer.Ordinal))
        {
            foreach (string peptide in peptides.OrderBy(p => p, StringComparer.Ordinal))
            {
                text.AppendLine(string.Join("\t", accession, peptide, "ambiguous"));
            }
        }

        foreach ((string accession, IReadOnlyCollection<string> peptides) in
                 UnrepairableRunCollisionsByAccession.OrderBy(kv => kv.Key, StringComparer.Ordinal))
        {
            foreach (string peptide in peptides.OrderBy(p => p, StringComparer.Ordinal))
            {
                text.AppendLine(string.Join("\t", accession, peptide, "unrepairableRunCollision"));
            }
        }

        foreach ((string accession, IReadOnlyCollection<string> peptides) in
                 ForeignPeptidesSharedWithTarget.OrderBy(kv => kv.Key, StringComparer.Ordinal))
        {
            foreach (string peptide in peptides.OrderBy(p => p, StringComparer.Ordinal))
            {
                text.AppendLine(string.Join("\t", accession, peptide, "sharedWithTarget"));
            }
        }

        return text.ToString();
    }

    public EntrapmentProvenance Provenance { get; }

    /// <summary>One entry per observed stratum, ascending.</summary>
    public IReadOnlyList<EntrapmentStratum> Strata { get; }

    /// <summary>The same counts across every stratum.</summary>
    public EntrapmentStratum Total { get; }

    /// <summary>Counts occurrences of any of <paramref name="residues"/> in a peptide.</summary>
    /// <remarks>O-glycosylation localization stratifies on "ST": the count of candidate sites in a
    /// peptide is what sets the difficulty of placing a modification on it.</remarks>
    public static Func<string, int> CountResidues(string residues)
    {
        if (string.IsNullOrEmpty(residues))
        {
            throw new MzLibException("CountResidues needs at least one residue to count.");
        }

        var wanted = new HashSet<char>(residues);
        return peptide => peptide.Count(wanted.Contains);
    }

    /// <summary>Counts N-X-S/T sequons, where X is any residue but proline.</summary>
    /// <remarks>
    /// A rearrangement destroys and creates these, so unlike a plain residue count this one is not
    /// preserved. Reporting it is the point: drift is acceptable, drift nobody can see is not.
    /// </remarks>
    public static int CountNGlycoSequons(string peptide)
    {
        int count = 0;
        for (int i = 0; i + 2 < peptide.Length; i++)
        {
            if (peptide[i] == 'N' && peptide[i + 1] != 'P' && (peptide[i + 2] == 'S' || peptide[i + 2] == 'T'))
            {
                count++;
            }
        }

        return count;
    }

    /// <summary>The report as a tab-separated table, provenance in leading comment lines.</summary>
    public string ToTabSeparated()
    {
        var text = new StringBuilder();
        text.AppendLine($"# method\t{Provenance.Method}");
        text.AppendLine($"# enzyme\t{Provenance.Enzyme}");
        text.AppendLine($"# seed\t{Provenance.Seed}");
        text.AppendLine($"# foldCount\t{Provenance.FoldCount}");
        text.AppendLine($"# maxMissedCleavages\t{Provenance.MaxMissedCleavages}");
        text.AppendLine($"# minPeptideLength\t{Provenance.MinPeptideLength}");
        text.AppendLine($"# maxPeptideLength\t{Provenance.MaxPeptideLength}");
        // targetPeptides/entrapmentPeptides count BASE PIECES; searchSpacePeptides counts what a
        // search reports, missed cleavages included. `ambiguous` belongs to the latter -- dividing
        // it by the former divides two different populations, which is how 0.24% got written
        // as "1.19%".
        text.AppendLine(string.Join("\t", "siteCount", "targetPeptides", "entrapmentPeptides",
            "achievedFoldRatio", "unpairableNoPermutationExists", "unpairableAllPermutationsTaken",
            "unpairableSpaceTooSmallForFoldCount", "searchSpacePeptides", "ambiguous",
            "mcSpanningAnExcision", "unrepairableRunCollisions"));

        foreach (EntrapmentStratum stratum in Strata.Append(Total))
        {
            text.AppendLine(string.Join("\t",
                ReferenceEquals(stratum, Total) ? "all" : stratum.SiteCount.ToString(CultureInfo.InvariantCulture),
                stratum.TargetPeptides.ToString(CultureInfo.InvariantCulture),
                stratum.EntrapmentPeptides.ToString(CultureInfo.InvariantCulture),
                stratum.AchievedFoldRatio.ToString("0.0000", CultureInfo.InvariantCulture),
                stratum.UnpairableNoPermutationExists.ToString(CultureInfo.InvariantCulture),
                stratum.UnpairableAllPermutationsTaken.ToString(CultureInfo.InvariantCulture),
                stratum.UnpairableSpaceTooSmallForFoldCount.ToString(CultureInfo.InvariantCulture),
                stratum.SearchSpacePeptides.ToString(CultureInfo.InvariantCulture),
                stratum.Ambiguous.ToString(CultureInfo.InvariantCulture),
                stratum.MissedCleavagePeptidesSpanningAnExcision.ToString(CultureInfo.InvariantCulture),
                stratum.UnrepairableRunCollisions.ToString(CultureInfo.InvariantCulture)));
        }

        return text.ToString();
    }
}

/// <summary>
/// Accumulates a report while a database is generated, one target protein and fold at a time.
/// </summary>
/// <remarks>
/// Counts are over the peptides a search could actually report -- base pieces at least
/// <see cref="IDigestionParams.MinLength"/> long. Shorter pieces are never identified on their own,
/// so counting them would flatter every ratio here.
/// </remarks>
public sealed class EntrapmentReportBuilder
{
    private readonly IDigestionParams _digestionParams;
    private readonly int _foldCount;
    private readonly int _seed;
    private readonly Func<string, int> _siteCounter;
    private readonly Dictionary<int, EntrapmentStratum> _strata = new();
    private readonly HashSet<string> _countedTargetPieces = new();
    private readonly Dictionary<string, HashSet<string>> _ambiguousByAccession = new();
    private readonly Dictionary<string, List<string>> _unrepairableByAccession = new();
    private readonly Dictionary<string, IReadOnlyCollection<string>> _foreignShared = new();
    private int _foreignEntries;

    /// <param name="siteCounter">What to stratify by, e.g.
    /// <see cref="EntrapmentReport.CountResidues"/>("ST"). Null puts everything in one stratum.</param>
    public EntrapmentReportBuilder(IDigestionParams digestionParams, int foldCount, int seed,
        Func<string, int>? siteCounter = null)
    {
        if (digestionParams is null)
        {
            throw new MzLibException("A report needs the digestion parameters it was generated under.");
        }

        _digestionParams = digestionParams;
        _foldCount = foldCount;
        _seed = seed;
        _siteCounter = siteCounter ?? (_ => 0);
    }

    /// <summary>Records one target protein's assembly for one fold.</summary>
    public void Add(Protein target, int fold, EntrapmentAssembly assembly)
    {
        if (target is null || assembly is null)
        {
            throw new MzLibException("A report entry needs both a target protein and its assembly.");
        }

        if (!_ambiguousByAccession.TryGetValue(target.Accession, out HashSet<string>? ambiguous))
        {
            var pairing = new EntrapmentPairing(target, _digestionParams);
            ambiguous = pairing.AmbiguousPeptides.ToHashSet();
            _ambiguousByAccession[target.Accession] = ambiguous;

            foreach (string peptide in ambiguous)
            {
                StratumFor(peptide).Ambiguous++;
            }

            // Not attributable to a stratum: the count is over the protein's whole search space,
            // and stratifying it would mean re-enumerating every run to classify each one.
            StratumFor(string.Empty).SearchSpacePeptides += pairing.SearchablePeptideCount;
        }

        foreach (EntrapmentPiece piece in assembly.Pieces)
        {
            if (piece.TargetPiece.Length < _digestionParams.MinLength)
            {
                continue;   // never identified on its own, so not part of any ratio here
            }

            EntrapmentStratum stratum = StratumFor(piece.TargetPiece);

            // Each target peptide is counted once however many folds run, so the ratio below is
            // partners per target -- the achieved r -- rather than a fraction of what was asked.
            if (_countedTargetPieces.Add(target.Accession + " " + piece.Index))
            {
                stratum.TargetPeptides++;
            }

            switch (piece.Outcome)
            {
                case PieceOutcome.Permuted:
                case PieceOutcome.KeptVerbatimTooShort:
                    stratum.EntrapmentPeptides++;
                    break;
                case PieceOutcome.Excised when piece.Failure == EntrapmentFailure.NoPermutationExists:
                    stratum.UnpairableNoPermutationExists++;
                    break;
                case PieceOutcome.Excised when piece.Failure == EntrapmentFailure.SpaceTooSmallForFoldCount:
                    stratum.UnpairableSpaceTooSmallForFoldCount++;
                    break;
                case PieceOutcome.Excised:
                    stratum.UnpairableAllPermutationsTaken++;
                    break;
            }
        }

        // Not attributable to any one stratum: a broken run spans several peptides at once.
        StratumFor(string.Empty).MissedCleavagePeptidesSpanningAnExcision +=
            assembly.MissedCleavagePeptidesSpanningAnExcision;
        StratumFor(string.Empty).UnrepairableRunCollisions += assembly.UnrepairableRunCollisions;
        if (assembly.UnrepairableRunCollisionPeptides.Count > 0)
        {
            if (!_unrepairableByAccession.TryGetValue(target.Accession, out List<string>? collisions))
            {
                collisions = new List<string>();
                _unrepairableByAccession[target.Accession] = collisions;
            }

            // Distinct peptides, not placements. The same run can collide at two points in one
            // low-complexity protein, and an exclusion list is read as "these peptides" -- a
            // consumer subtracting a duplicate twice would over-correct.
            foreach (string collision in assembly.UnrepairableRunCollisionPeptides)
            {
                if (!collisions.Contains(collision))
                {
                    collisions.Add(collision);
                }
            }
        }
    }

    /// <summary>
    /// Records the foreign-species arm: how many entries it contributed, and which of its peptides
    /// are also target peptides.
    /// </summary>
    /// <remarks>
    /// Separate from <see cref="Add"/> because a foreign entry has no assembly and no target -- it
    /// was relabelled, not rearranged -- so none of the per-stratum permutation figures apply to it.
    /// Folding it into those would put entries in the ratio that were never at risk of failing.
    /// </remarks>
    public void AddForeign(int entryCount,
        IReadOnlyDictionary<string, IReadOnlyCollection<string>> sharedWithTarget)
    {
        if (entryCount < 0)
        {
            throw new MzLibException($"Foreign entry count must not be negative, but was {entryCount}.");
        }

        _foreignEntries += entryCount;
        foreach ((string accession, IReadOnlyCollection<string> peptides) in
                 sharedWithTarget ?? new Dictionary<string, IReadOnlyCollection<string>>())
        {
            _foreignShared[accession] = peptides;
        }
    }

    /// <summary>The finished report.</summary>
    public EntrapmentReport Build()
    {
        var total = new EntrapmentStratum(-1);
        foreach (EntrapmentStratum stratum in _strata.Values)
        {
            total.TargetPeptides += stratum.TargetPeptides;
            total.EntrapmentPeptides += stratum.EntrapmentPeptides;
            total.UnpairableNoPermutationExists += stratum.UnpairableNoPermutationExists;
            total.UnpairableAllPermutationsTaken += stratum.UnpairableAllPermutationsTaken;
            total.UnpairableSpaceTooSmallForFoldCount += stratum.UnpairableSpaceTooSmallForFoldCount;
            total.Ambiguous += stratum.Ambiguous;
            total.MissedCleavagePeptidesSpanningAnExcision += stratum.MissedCleavagePeptidesSpanningAnExcision;
            total.SearchSpacePeptides += stratum.SearchSpacePeptides;
            total.UnrepairableRunCollisions += stratum.UnrepairableRunCollisions;
        }

        var provenance = new EntrapmentProvenance(
            method: "composition-preserving permutation (deterministic unranking)",
            enzyme: _digestionParams.DigestionAgent.Name,
            seed: _seed,
            foldCount: _foldCount,
            maxMissedCleavages: _digestionParams.MaxMissedCleavages,
            minPeptideLength: _digestionParams.MinLength,
            maxPeptideLength: _digestionParams.MaxLength);

        return new EntrapmentReport(provenance,
            _strata.Values.OrderBy(s => s.SiteCount).ToList(),
            total,
            _ambiguousByAccession.Where(kv => kv.Value.Count > 0)
                .ToDictionary(kv => kv.Key, kv => (IReadOnlyCollection<string>)kv.Value),
            _unrepairableByAccession.Where(kv => kv.Value.Count > 0)
                .ToDictionary(kv => kv.Key, kv => (IReadOnlyCollection<string>)kv.Value))
        {
            ForeignEntries = _foreignEntries,
            ForeignPeptidesSharedWithTarget = _foreignShared,
        };
    }

    private EntrapmentStratum StratumFor(string peptide)
    {
        int siteCount = peptide.Length == 0 ? 0 : _siteCounter(peptide);
        if (!_strata.TryGetValue(siteCount, out EntrapmentStratum? stratum))
        {
            stratum = new EntrapmentStratum(siteCount);
            _strata[siteCount] = stratum;
        }

        return stratum;
    }
}
