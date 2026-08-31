using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases.EntrapmentGeneration;

namespace Test.DatabaseTests;

[TestFixture]
[ExcludeFromCodeCoverage]
public class EntrapmentReportTests
{
    private static readonly HashSet<string> NothingForbidden = new();

    private static IDigestionParams Tryptic =>
        new DigestionParams("trypsin", minPeptideLength: 7, maxMissedCleavages: 2);

    private static readonly string[] Proteins =
    {
        "SYKALADQMNLLLSKGGVDTTPFAWENDRQISTLGGYK",
        "AKMTTSTTPAPTTTKQEEEKKKSSSSSSRGGVDTTPFAWENDR",
        "MPSSTSRLLGVDKAAATPTGGGSQSCCCKYPERDNRSTTTK"
    };

    private static EntrapmentReport BuildReport(int foldCount = 1, Func<string, int> siteCounter = null,
        IEnumerable<string> sequences = null)
    {
        IDigestionParams digestion = Tryptic;
        var builder = new EntrapmentReportBuilder(digestion, foldCount, seed: 1,
            siteCounter ?? EntrapmentReport.CountResidues("ST"));

        foreach (string sequence in sequences ?? Proteins)
        {
            var target = new Protein(sequence, "P" + sequence.Length + sequence[0]);
            for (int fold = 0; fold < foldCount; fold++)
            {
                EntrapmentProteinGenerator.Create(target, digestion, NothingForbidden,
                    out EntrapmentAssembly assembly, fold: fold, foldCount: foldCount);
                builder.Add(target, fold, assembly);
            }
        }

        return builder.Build();
    }

    // ---- the distribution itself -------------------------------------------

    [Test]
    public void Report_EmitsEveryStratumRatherThanASummary()
    {
        EntrapmentReport report = BuildReport();

        Assert.That(report.Strata.Count, Is.GreaterThan(1),
            "a single number would hide exactly the differences worth looking at");
        Assert.That(report.Strata.Select(s => s.SiteCount), Is.Ordered,
            "strata are emitted in ascending order so two reports can be compared line for line");
        Assert.That(report.Strata.Sum(s => s.TargetPeptides), Is.EqualTo(report.Total.TargetPeptides),
            "the strata must add up to the total");
    }

    [Test]
    public void Report_CountsOnlyPeptidesASearchCouldReport()
    {
        IDigestionParams digestion = Tryptic;
        EntrapmentReport report = BuildReport();

        // "AK" and other sub-length pieces exist in the fixtures but can never be identified alone,
        // so counting them would flatter every ratio in the report.
        int searchable = 0;
        foreach (string sequence in Proteins)
        {
            List<int> sites = digestion.DigestionAgent.GetDigestionSiteIndices(sequence);
            sites.Sort();
            searchable += Enumerable.Range(0, sites.Count - 1)
                .Count(i => sites[i + 1] - sites[i] >= digestion.MinLength);
        }

        Assert.That(report.Total.TargetPeptides, Is.EqualTo(searchable));
    }

    // ---- failures kept apart -----------------------------------------------

    [Test]
    public void Report_SeparatesFailuresByCauseAndByStratum()
    {
        EntrapmentReport report = BuildReport();

        // "SSSSSSR" carries six candidate sites and has exactly one arrangement, so it lands in the
        // no-permutation bucket of a high stratum -- not in a single undifferentiated "failed".
        EntrapmentStratum sixSites = report.Strata.Single(s => s.SiteCount == 6);
        Assert.That(sixSites.UnpairableNoPermutationExists, Is.EqualTo(1));
        Assert.That(sixSites.UnpairableAllPermutationsTaken, Is.Zero);
        Assert.That(sixSites.UnpairableSpaceTooSmallForFoldCount, Is.Zero);

        // Against a literal, not against the same expression the property computes -- restating
        // the implementation proves only that it equals itself.
        Assert.That(report.Total.UnpairableNoPermutationExists, Is.EqualTo(1));
        Assert.That(report.Total.UnpairableAllPermutationsTaken, Is.Zero);
        Assert.That(report.Total.UnpairableSpaceTooSmallForFoldCount, Is.Zero);
        Assert.That(report.Total.Unpairable, Is.EqualTo(1));
    }

    [Test]
    public void Report_CountsMissedCleavagePeptidesThatBreakTheIsomerismInvariant()
    {
        EntrapmentReport report = BuildReport();

        // The unpairable tract sits mid-protein on purpose: excising a TERMINAL piece breaks no
        // runs at all, so a fixture with it at the end would assert nothing.
        Assert.That(report.Total.MissedCleavagePeptidesSpanningAnExcision, Is.GreaterThan(0),
            "an excision inside a protein necessarily breaks the runs that spanned it");
    }

    // ---- achieved r --------------------------------------------------------

    [Test]
    public void Report_ReportsPartnersDeliveredPerTargetNotAFractionOfWhatWasAsked()
    {
        EntrapmentReport oneFold = BuildReport(foldCount: 1);
        EntrapmentReport fourFolds = BuildReport(foldCount: 4);

        // A target peptide is counted once however many folds run, so the ratio is the achieved r.
        Assert.That(fourFolds.Total.TargetPeptides, Is.EqualTo(oneFold.Total.TargetPeptides));
        Assert.That(fourFolds.Total.AchievedFoldRatio,
            Is.GreaterThan(oneFold.Total.AchievedFoldRatio * 3),
            "four folds must deliver roughly four partners per target, not one");
        Assert.That(fourFolds.Total.AchievedFoldRatio, Is.LessThanOrEqualTo(4d));
    }

    [Test]
    public void Report_ShowsWhereTheDatabaseFellShortOfWhatWasAsked()
    {
        EntrapmentReport report = BuildReport(foldCount: 4);

        // The peptide with one arrangement can never be paired, so its stratum cannot reach 4.
        EntrapmentStratum sixSites = report.Strata.Single(s => s.SiteCount == 6);
        Assert.That(sixSites.AchievedFoldRatio, Is.LessThan(4d));
        Assert.That(sixSites.UnpairableNoPermutationExists, Is.GreaterThan(0));
    }

    // ---- provenance --------------------------------------------------------

    [Test]
    public void Report_RecordsHowTheDatabaseWasMade()
    {
        EntrapmentReport report = BuildReport(foldCount: 3);

        // Without these the numbers above cannot be interpreted or reproduced: the peptide
        // population depends on the enzyme and digestion settings as much as on the method.
        Assert.That(report.Provenance.Enzyme, Is.EqualTo("trypsin"));
        Assert.That(report.Provenance.Seed, Is.EqualTo(1));
        Assert.That(report.Provenance.FoldCount, Is.EqualTo(3));
        Assert.That(report.Provenance.MaxMissedCleavages, Is.EqualTo(2));
        Assert.That(report.Provenance.MinPeptideLength, Is.EqualTo(7));
        Assert.That(report.Provenance.Method, Does.Contain("composition-preserving"));
    }

    // ---- stratifying by something else -------------------------------------

    [Test]
    public void Report_CanStratifyByAnythingTheCallerCounts()
    {
        EntrapmentReport bySites = BuildReport(siteCounter: EntrapmentReport.CountResidues("ST"));
        EntrapmentReport bySequons = BuildReport(siteCounter: EntrapmentReport.CountNGlycoSequons);

        Assert.That(bySites.Total.TargetPeptides, Is.EqualTo(bySequons.Total.TargetPeptides));
        Assert.That(bySites.Strata.Select(s => s.SiteCount),
            Is.Not.EqualTo(bySequons.Strata.Select(s => s.SiteCount)),
            "counting different things must produce different strata");
    }

    [Test]
    public void CountNGlycoSequons_MatchesTheMotifAndExcludesProline()
    {
        Assert.That(EntrapmentReport.CountNGlycoSequons("NAS"), Is.EqualTo(1));
        Assert.That(EntrapmentReport.CountNGlycoSequons("NAT"), Is.EqualTo(1));
        Assert.That(EntrapmentReport.CountNGlycoSequons("NPS"), Is.Zero, "X must not be proline");
        Assert.That(EntrapmentReport.CountNGlycoSequons("NAA"), Is.Zero);
        Assert.That(EntrapmentReport.CountNGlycoSequons("NASNAT"), Is.EqualTo(2));
        Assert.That(EntrapmentReport.CountNGlycoSequons("NA"), Is.Zero, "too short to hold a sequon");
    }

    [Test]
    public void CountResidues_RefusesAnEmptyResidueSet()
    {
        Assert.Throws<MzLibUtil.MzLibException>(() => EntrapmentReport.CountResidues(""));
        Assert.Throws<MzLibUtil.MzLibException>(() => EntrapmentReport.CountResidues(null));
    }

    // ---- the emitted table --------------------------------------------------

    [Test]
    public void Report_WritesProvenanceAndOneRowPerStratumPlusATotal()
    {
        EntrapmentReport report = BuildReport(foldCount: 2);
        string[] lines = report.ToTabSeparated()
            .Split('\n').Select(l => l.TrimEnd('\r')).Where(l => l.Length > 0).ToArray();

        string[] comments = lines.Where(l => l.StartsWith("#")).ToArray();
        Assert.That(comments.Any(l => l.Contains("trypsin")), Is.True);
        Assert.That(comments.Any(l => l.StartsWith("# seed")), Is.True);

        string header = lines.First(l => !l.StartsWith("#"));
        Assert.That(header.Split('\t'), Does.Contain("achievedFoldRatio"));

        string[] rows = lines.Where(l => !l.StartsWith("#")).Skip(1).ToArray();
        Assert.That(rows.Length, Is.EqualTo(report.Strata.Count + 1), "one row per stratum, plus the total");
        Assert.That(rows.Last(), Does.StartWith("all"));
        Assert.That(rows.All(r => r.Split('\t').Length == header.Split('\t').Length), Is.True);
    }

    [Test]
    public void Report_NumbersReconcileWithTheDatabaseActuallyWritten()
    {
        // The report is only worth reading if it describes the database that came out, so count the
        // entrapment peptides independently and compare.
        IDigestionParams digestion = Tryptic;
        var builder = new EntrapmentReportBuilder(digestion, foldCount: 1, seed: 1,
            EntrapmentReport.CountResidues("ST"));
        int actualEntrapmentPeptides = 0;

        foreach (string sequence in Proteins)
        {
            var target = new Protein(sequence, "P" + sequence.Length + sequence[0]);
            Protein entrapment = EntrapmentProteinGenerator.Create(target, digestion, NothingForbidden,
                out EntrapmentAssembly assembly);
            builder.Add(target, 0, assembly);

            List<int> sites = digestion.DigestionAgent.GetDigestionSiteIndices(entrapment.BaseSequence);
            sites.Sort();
            actualEntrapmentPeptides += Enumerable.Range(0, sites.Count - 1)
                .Count(i => sites[i + 1] - sites[i] >= digestion.MinLength);
        }

        Assert.That(builder.Build().Total.EntrapmentPeptides, Is.EqualTo(actualEntrapmentPeptides));
    }

    [Test]
    public void Report_SeesTheOtherTwoFailureCauses()
    {
        // The fixtures above only ever produce "no permutation exists". These two paths need their
        // own shapes: a space too small for the folds asked of it, and a space fully occupied.
        IDigestionParams digestion = Tryptic;
        // "EEEEEEQK" is one base piece (trypsin cuts after every K, so QEEEKKK would split into
        // three) and its seven residues before the pinned K admit exactly seven arrangements.
        const string sequence = "EEEEEEQKGGVDTTPFAWENDR";
        var target = new Protein(sequence, "P1");

        var tooSmall = new EntrapmentReportBuilder(digestion, foldCount: 9, seed: 1,
            EntrapmentReport.CountResidues("ST"));
        EntrapmentProteinGenerator.Create(target, digestion, NothingForbidden,
            out EntrapmentAssembly smallAssembly, fold: 0, foldCount: 9);
        tooSmall.Add(target, 0, smallAssembly);
        Assert.That(tooSmall.Build().Total.UnpairableSpaceTooSmallForFoldCount, Is.EqualTo(1));

        var everyArrangement = new HashSet<string>();
        System.Numerics.BigInteger size = UsefulProteomicsDatabases.DecoySequenceValidator
            .PermutationSpaceSize("EEEEEEQK", digestion.DigestionAgent.DigestionMotifs);
        for (System.Numerics.BigInteger i = 0; i < size; i++)
        {
            everyArrangement.Add(UsefulProteomicsDatabases.DecoySequenceValidator
                .UnrankPermutation("EEEEEEQK", digestion.DigestionAgent.DigestionMotifs, i, out _));
        }

        var taken = new EntrapmentReportBuilder(digestion, foldCount: 1, seed: 1,
            EntrapmentReport.CountResidues("ST"));
        EntrapmentProteinGenerator.Create(target, digestion, everyArrangement,
            out EntrapmentAssembly takenAssembly);
        taken.Add(target, 0, takenAssembly);
        Assert.That(taken.Build().Total.UnpairableAllPermutationsTaken, Is.EqualTo(1));
    }

    [Test]
    public void Report_TotalsAddUpFieldByField()
    {
        EntrapmentReport report = BuildReport(foldCount: 2);

        // Every field, not just the peptide count: the total is accumulated field by field, so a
        // slip in one of them would otherwise go unseen.
        Assert.That(report.Total.EntrapmentPeptides, Is.EqualTo(report.Strata.Sum(s => s.EntrapmentPeptides)));
        Assert.That(report.Total.UnpairableNoPermutationExists,
            Is.EqualTo(report.Strata.Sum(s => s.UnpairableNoPermutationExists)));
        Assert.That(report.Total.UnpairableAllPermutationsTaken,
            Is.EqualTo(report.Strata.Sum(s => s.UnpairableAllPermutationsTaken)));
        Assert.That(report.Total.UnpairableSpaceTooSmallForFoldCount,
            Is.EqualTo(report.Strata.Sum(s => s.UnpairableSpaceTooSmallForFoldCount)));
        Assert.That(report.Total.Ambiguous, Is.EqualTo(report.Strata.Sum(s => s.Ambiguous)));
        Assert.That(report.Total.MissedCleavagePeptidesSpanningAnExcision,
            Is.EqualTo(report.Strata.Sum(s => s.MissedCleavagePeptidesSpanningAnExcision)));
        Assert.That(report.Total.EntrapmentPeptides, Is.GreaterThan(0), "the fixture must be non-trivial");
    }

    [Test]
    public void Report_CountsAnAmbiguousPeptideOnceHoweverManyFoldsRun()
    {
        // Ambiguity is a property of the target protein, not of a fold, so it is worked out once per
        // accession. Counting it again per fold would inflate it in step with the fold count.
        var pairing = new EntrapmentPairing(
            new Protein("LIHTGVKAAADEFGHIKLIHTVGKMMWWYYPPQQR", "P57071"), Tryptic);
        Assert.That(pairing.AmbiguousPeptides.Count, Is.GreaterThan(0), "the fixture must be ambiguous");

        int[] counts = new[] { 1, 4 }.Select(folds =>
            BuildReport(foldCount: folds,
                sequences: new[] { "LIHTGVKAAADEFGHIKLIHTVGKMMWWYYPPQQR" }).Total.Ambiguous).ToArray();

        Assert.That(counts[0], Is.GreaterThan(0));
        Assert.That(counts[1], Is.EqualTo(counts[0]), "four folds must not quadruple the ambiguity count");
    }

    [Test]
    public void Report_HandlesHavingBeenGivenNothing()
    {
        EntrapmentReport empty = new EntrapmentReportBuilder(Tryptic, 1, 1).Build();

        Assert.That(empty.Strata, Is.Empty);
        Assert.That(empty.Total.TargetPeptides, Is.Zero);
        Assert.That(empty.Total.AchievedFoldRatio, Is.Zero, "no targets means no ratio, not a divide by zero");
        Assert.That(empty.ToTabSeparated(), Does.Contain("# enzyme"));
    }

    [Test]
    public void ReportBuilder_RefusesIncompleteInput()
    {
        Assert.Throws<MzLibUtil.MzLibException>(() => new EntrapmentReportBuilder(null, 1, 1));

        var builder = new EntrapmentReportBuilder(Tryptic, 1, 1);
        var target = new Protein(Proteins[0], "P1");
        EntrapmentProteinGenerator.Create(target, Tryptic, NothingForbidden, out EntrapmentAssembly assembly);

        Assert.Throws<MzLibUtil.MzLibException>(() => builder.Add(null, 0, assembly));
        Assert.Throws<MzLibUtil.MzLibException>(() => builder.Add(target, 0, null));
    }

    [Test]
    public void Report_TableCarriesTheNumbersAndNotJustTheShape()
    {
        EntrapmentReport report = BuildReport();
        string[] lines = report.ToTabSeparated()
            .Split('\n').Select(l => l.TrimEnd('\r')).Where(l => l.Length > 0).ToArray();

        // Provenance values, not merely the presence of a label.
        Assert.That(lines, Does.Contain("# enzyme\ttrypsin"));
        Assert.That(lines, Does.Contain("# seed\t1"));
        Assert.That(lines, Does.Contain("# maxMissedCleavages\t2"));
        Assert.That(lines, Does.Contain("# minPeptideLength\t7"));

        // And the total row's figures must match the object the table was rendered from.
        string[] total = lines.Last().Split('\t');
        Assert.That(total[0], Is.EqualTo("all"));
        Assert.That(total[1], Is.EqualTo(report.Total.TargetPeptides.ToString()));
        Assert.That(total[2], Is.EqualTo(report.Total.EntrapmentPeptides.ToString()));
        Assert.That(total[4], Is.EqualTo(report.Total.UnpairableNoPermutationExists.ToString()));
        Assert.That(total[7], Is.EqualTo(report.Total.SearchSpacePeptides.ToString()));
        Assert.That(total[8], Is.EqualTo(report.Total.Ambiguous.ToString()));
        Assert.That(total[10], Is.EqualTo(report.Total.UnrepairableRunCollisions.ToString()));
    }

    /// <summary>
    /// A protein whose short pieces collide by composition while its searchable peptides do not.
    /// "AK" and "KA" cannot both occur, but two 2-residue pieces sharing a key are easy to build,
    /// and before the length bound they were counted as ambiguous despite being unsearchable.
    /// </summary>
    private const string ShortPieceCollisions = "AKAKMSTQAEVDLNSGWKALADQMNLLLSK";

    [Test]
    public void PairingIgnoresPeptidesTooShortToBeReported()
    {
        // The index must hold only what a search could report. A key begins with the peptide's
        // length, so short peptides never collide with long ones and dropping them cannot change
        // how a searchable peptide resolves -- it only stops them being counted as ambiguous.
        var lenient = new DigestionParams("trypsin", minPeptideLength: 1, maxMissedCleavages: 2);
        var realistic = new DigestionParams("trypsin", minPeptideLength: 7, maxMissedCleavages: 2);
        var protein = new Protein(ShortPieceCollisions, "P00001");

        var withShort = new EntrapmentPairing(protein, lenient);
        var withoutShort = new EntrapmentPairing(protein, realistic);

        Assert.That(withShort.SearchablePeptideCount, Is.GreaterThan(withoutShort.SearchablePeptideCount),
            "the lenient bound must admit more peptides, or the fixture proves nothing");
        Assert.That(withoutShort.AmbiguousPeptides.Any(pep => pep.Length < 7), Is.False,
            "no peptide below the minimum length may be reported as ambiguous");
    }

    [Test]
    public void SearchablePeptideCountSpansMissedCleavagesNotJustBasePieces()
    {
        // The denominator an ambiguity rate and an FDP estimator's r both need. It must exceed the
        // base-piece count, because a run of adjacent pieces is a peptide a search reports too.
        var digestion = new DigestionParams("trypsin", minPeptideLength: 7, maxMissedCleavages: 2);
        var protein = new Protein("MSTQAEVDLNSGWKALADQMNLLLSKGGVDTTPFAWENDRQISTLGGYK", "P00001");

        var pairing = new EntrapmentPairing(protein, digestion);

        // four base pieces, all >= 7, plus runs of two and of three
        Assert.That(pairing.SearchablePeptideCount, Is.EqualTo(4 + 3 + 2));
    }

    [Test]
    public void ReportCarriesTheSearchSpaceDenominatorBesideTheBasePieceCounts()
    {
        var digestion = Tryptic;
        var builder = new EntrapmentReportBuilder(digestion, 1, 1, EntrapmentReport.CountResidues("ST"));
        var protein = new Protein("MSTQAEVDLNSGWKALADQMNLLLSKGGVDTTPFAWENDRQISTLGGYK", "P00001");

        Protein _ = EntrapmentProteinGenerator.Create(protein, digestion, NothingForbidden,
            out EntrapmentAssembly assembly);
        builder.Add(protein, 0, assembly);
        EntrapmentReport report = builder.Build();

        Assert.That(report.Total.SearchSpacePeptides, Is.EqualTo(9),
            "runs of one, two and three base pieces");
        Assert.That(report.Total.SearchSpacePeptides, Is.GreaterThan(report.Total.TargetPeptides),
            "the search space must exceed the base pieces it is built from");

        string tsv = report.ToTabSeparated();
        Assert.That(tsv, Does.Contain("searchSpacePeptides"));
        Assert.That(tsv, Does.Contain("unrepairableRunCollisions"));
    }

    [Test]
    public void ReportCountsUnrepairableRunCollisions()
    {
        // A run ending in a piece with no alternative arrangement cannot be permuted away from a
        // real peptide, so it is counted rather than repaired -- and the count has to reach the
        // report, or a consumer cannot exclude the peptides it describes.
        var digestion = Tryptic;
        const string target = "MSTQAEVDLNSGWKAAR";

        EntrapmentAssembly free = EntrapmentAssembler.Assemble(target, digestion, NothingForbidden);
        string run = free.Pieces[0].EntrapmentPiece_ + free.Pieces[1].EntrapmentPiece_;

        var protein = new Protein(target, "P00001");
        var builder = new EntrapmentReportBuilder(digestion, 1, 1);
        Protein _ = EntrapmentProteinGenerator.Create(protein, digestion,
            new HashSet<string> { run }, out EntrapmentAssembly guarded);
        builder.Add(protein, 0, guarded);

        Assert.That(guarded.UnrepairableRunCollisions, Is.EqualTo(1));
        Assert.That(builder.Build().Total.UnrepairableRunCollisions, Is.EqualTo(1));
    }
}
