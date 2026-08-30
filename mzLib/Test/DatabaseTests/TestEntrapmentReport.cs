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

        Assert.That(report.Total.Unpairable,
            Is.EqualTo(report.Total.UnpairableNoPermutationExists
                       + report.Total.UnpairableAllPermutationsTaken
                       + report.Total.UnpairableSpaceTooSmallForFoldCount));
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
}
