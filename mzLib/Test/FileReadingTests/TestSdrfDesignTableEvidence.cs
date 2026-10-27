using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MzLibUtil;
using NUnit.Framework;
using Readers;
using UsefulProteomicsDatabases;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Design tables a deposit holds beside its runs (sdrf MASSIVE-SDRF.md, MV3): Spectronaut's condition setup, an
    /// MSstats annotation and MaxQuant's experimental design, read as Certain claims from the deposit, and the drafter
    /// taking conditions, replicates, fractions and sample names from them. The main case is MSV000096226 (PXD057243):
    /// its condition setup names 10 of its 11 brain runs, and the deposit also holds 48 runs of another study's hearts.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrfDesignTableEvidence
    {
        private static SupplementTable Table(string file, string[] header, params string[][] rows) =>
            new(file, "", "", header, rows.Select(r => (IReadOnlyList<string>)r).ToList(), Enumerable.Range(2, rows.Length).ToList(), false);

        private const string SetupFile = "metadata/220424_AM_MouseBrain_DR_WholeProteome_ConditionSetup.tsv";

        /// <summary>MSV000096226's condition setup, verbatim but for the colour column.</summary>
        private static SupplementTable Setup() => Table(SetupFile,
            new[] { "#", "Reference", "Run Label", "Condition", "Fraction", "Replicate", "Quantity Correction Factor", "Label", "Color", "File Name" },
            new[] { "1", "True", "240418_AM_O1.raw", "Old", "NA", "1", "1", "Old", "#B6EC13", "240418_AM_O1" },
            new[] { "2", "True", "240418_AM_O2.raw", "Old", "NA", "2", "1", "Old", "#B6EC13", "240418_AM_O2" },
            new[] { "3", "True", "240418_AM_O3.raw", "Old", "NA", "3", "1", "Old", "#B6EC13", "240418_AM_O3" },
            new[] { "4", "True", "240418_AM_O4.raw", "Old", "NA", "4", "1", "Old", "#B6EC13", "240418_AM_O4" },
            new[] { "5", "False", "240418_AM_ODR_2.raw", "RF", "NA", "1", "1", "RF", "#1349EC", "240418_AM_ODR_2" },
            new[] { "6", "False", "240418_AM_ODR_3_20240420123730.raw", "RF", "NA", "2", "1", "RF", "#1349EC", "240418_AM_ODR_3_20240420123730" },
            new[] { "7", "False", "240418_AM_ODR_4.raw", "RF", "NA", "3", "1", "RF", "#1349EC", "240418_AM_ODR_4" },
            new[] { "8", "False", "240418_AM_ODR_RF_2.raw", "RF", "NA", "4", "1", "RF", "#1349EC", "240418_AM_ODR_RF_2" },
            new[] { "9", "False", "240418_AM_ODR_RF_3.raw", "RF", "NA", "5", "1", "RF", "#1349EC", "240418_AM_ODR_RF_3" },
            new[] { "10", "False", "240418_AM_ODR_RF_4.raw", "RF", "NA", "6", "1", "RF", "#1349EC", "240418_AM_ODR_RF_4" });

        private static readonly string[] Brain =
        {
            "240418_AM_O1.raw", "240418_AM_O2.raw", "240418_AM_O3.raw", "240418_AM_O4.raw", "240418_AM_ODR_2.raw", "240418_AM_ODR_3.raw",
            "240418_AM_ODR_3_20240420123730.raw", "240418_AM_ODR_4.raw", "240418_AM_ODR_RF_2.raw", "240418_AM_ODR_RF_3.raw", "240418_AM_ODR_RF_4.raw",
        };

        /// <summary>Another study's heart runs in the same deposit (pride 025, trap 1): two arms, high-pH fractions.</summary>
        private static readonly string[] Heart = Enumerable.Range(1, 3)
            .SelectMany(i => new[] { $"190418_NR_Hanna_DR_HpH_fr{i}.raw", $"190419_NR_Hanna_NF_HpH_fr{i}.raw" }).ToArray();

        private static string[] Runs => Brain.Concat(Heart).ToArray();

        private static PrideProject Project() => new()
        {
            Title = "Old mouse brain proteome under calorie restriction and refeeding",
            Organisms = { new CvParam("NEWT", "NEWT:10090", "Mus musculus (mouse)", "") },
            Instruments = { new CvParam("MS", "MS:1003028", "Orbitrap Exploris 480", "") },
        };

        private static SdrfEvidence One(IEnumerable<SdrfEvidence> e, string file, string column, string label = "") =>
            e.Single(x => x.DataFile == file && x.Column == column && x.Label == label);

        private static SdrfDraftRow Row(SdrfDraft d, string file) => d.Rows.Single(r => r.DataFile == file);

        private static string Factor(SdrfDraft d, string file, string column) =>
            Row(d, file).Factors[d.FactorColumns.ToList().IndexOf(column)].Value;

        // ---------------- reading the tables ----------------

        [Test]
        public void ASpectronautConditionSetupGivesEachRunItsConditionAndReplicate()
        {
            var e = SampleEvidenceExtractor.ExtractDeposited(new[] { Setup() }, Runs);

            var condition = One(e, "240418_AM_ODR_RF_4.raw", "factor value[condition]");
            Assert.That(condition.Value, Is.EqualTo("RF"));
            Assert.That((condition.Method, condition.Source, condition.Confidence, condition.Locator),
                Is.EqualTo(("design-table", "deposit", SdrfEvidenceConfidence.Certain, SetupFile + "!R11C4")));
            Assert.That(One(e, "240418_AM_O3.raw", "characteristics[biological replicate]").Value, Is.EqualTo("3"));
            Assert.That(e.Select(x => x.DataFile).Distinct(), Has.Exactly(10).Items, "the setup names 10 runs");
            Assert.That(e.Where(x => x.DataFile == "240418_AM_ODR_3.raw" || x.DataFile.Contains("Hanna")), Is.Empty,
                "a run the setup does not name gets nothing");
        }

        [Test]
        public void ASpectronautLabelColourAndMissingFractionAreNotClaims()
        {
            var e = SampleEvidenceExtractor.ExtractDeposited(new[] { Setup() }, Runs);
            Assume.That(e.Where(x => x.Method == "design-table"), Is.Not.Empty);

            // Spectronaut's Label repeats (or truncates) the condition; it is not a chemical label.
            Assert.That(e.Select(x => x.Column), Has.None.EqualTo("comment[label]"));
            Assert.That(e.Select(x => x.Label), Has.All.Empty);
            Assert.That(e.Select(x => x.Column), Has.None.EqualTo("comment[fraction identifier]"), "NA is no fraction");
            Assert.That(e.Select(x => x.Value), Has.None.StartsWith("#"), "a colour is no characteristic");
        }

        [Test]
        public void SpectronautsNotDefinedIsNoCondition()
        {
            // MSV000090777 keeps 27 setups; a library-only analysis writes "Not Defined" where a comparison writes the arm.
            string[] header = { "#", "Is Reference", "Run Label", "Condition", "Fraction", "Replicate" };
            var library = Table("quant/Non_Library_ConditionSetup.tsv", header, new[] { "1", "True", "run_inj_27.wiff", "Not Defined", "NA", "1" });
            var comparison = Table("quant/InjvsNon_ConditionSetup.tsv", header, new[] { "2", "False", "run_inj_27.wiff", "Injured", "NA", "1" });
            var runs = new[] { "run_inj_27.wiff", "run_non_23.wiff" };

            var e = SampleEvidenceExtractor.ExtractDeposited(new[] { library, comparison }, runs);
            Assert.That(e.Where(x => x.Column == "factor value[condition]").Select(x => x.Value), Is.EqualTo(new[] { "Injured" }));

            var d = SdrfDrafter.Draft(Project(), runs, e);
            Assert.That(Factor(d, "run_inj_27.wiff", "factor value[condition]"), Is.EqualTo("Injured"));
        }

        [Test]
        public void ARunListedInTwoFormatsGetsTheClaimsOnBoth()
        {
            // MSV000090777 lists every SCIEX run as .wiff (with its .wiff.scan) and as .wiff2; the setup names the .wiff.
            var t = Table("quant/ConditionSetup.tsv", new[] { "#", "Is Reference", "Run Label", "Condition", "Fraction", "Replicate" },
                new[] { "1", "True", "22_0505_002_SSL2_ND_Non_23_01.wiff", "Non-injured", "NA", "1" });
            var runs = new[] { "22_0505_002_SSL2_ND_Non_23_01.wiff", "22_0505_002_SSL2_ND_Non_23_01.wiff.scan", "22_0505_002_SSL2_ND_Non_23_01.wiff2" };

            var e = SampleEvidenceExtractor.ExtractDeposited(new[] { t }, runs);

            Assert.That(e.Where(x => x.Column == "factor value[condition]").Select(x => x.DataFile),
                Is.EquivalentTo(new[] { "22_0505_002_SSL2_ND_Non_23_01.wiff", "22_0505_002_SSL2_ND_Non_23_01.wiff2" }));
        }

        [Test]
        public void TheSameSetupReadAsASupplementSaysSupplement()
        {
            var e = SampleEvidenceExtractor.Extract(new[] { Setup() }, Runs);

            Assert.That(e.Where(x => x.Method == "design-table"), Is.Not.Empty);
            Assert.That(e.Where(x => x.Method == "design-table").Select(x => x.Source), Has.All.EqualTo("supplement"));
        }

        [Test]
        public void AnMsstatsAnnotationGivesConditionSampleFractionAndTechnicalReplicate()
        {
            // MassIVE.quant's label-free format: a sample identifier in BioReplicate.
            var t = Table("metadata/run_annotation.csv", new[] { "Run", "Condition", "BioReplicate", "Fraction", "TechReplicate" },
                new[] { "ctrl_a_f1.raw", "Control", "patient_A", "1", "1" },
                new[] { "ctrl_a_f2.raw", "Control", "patient_A", "2", "1" },
                new[] { "case_b_f1.raw", "Case", "patient_B", "1", "1" },
                new[] { "case_b_f1_rep.raw", "Case", "patient_B", "1", "2" });
            var runs = new[] { "ctrl_a_f1.raw", "ctrl_a_f2.raw", "case_b_f1.raw", "case_b_f1_rep.raw" };

            var e = SampleEvidenceExtractor.ExtractDeposited(new[] { t }, runs);

            Assert.That(One(e, "ctrl_a_f2.raw", "factor value[condition]").Value, Is.EqualTo("Control"));
            Assert.That(One(e, "ctrl_a_f2.raw", "source name").Value, Is.EqualTo("patient_A"));
            Assert.That(One(e, "ctrl_a_f2.raw", "comment[fraction identifier]").Value, Is.EqualTo("2"));
            Assert.That(One(e, "case_b_f1_rep.raw", "comment[technical replicate]").Value, Is.EqualTo("2"));
            Assert.That(e.Where(x => x.DataFile == "ctrl_a_f1.raw").Select(x => x.Column), Has.None.EqualTo("characteristics[biological replicate]"),
                "a sample identifier is not a replicate number");
        }

        [Test]
        public void AnMsstatsNumericBioReplicateIsAReplicateNumber()
        {
            var t = Table("annotation.csv", new[] { "Run", "Condition", "BioReplicate" },
                new[] { "c1.raw", "Control", "1" }, new[] { "c2.raw", "Control", "2" }, new[] { "t1.raw", "Treated", "1" });

            var e = SampleEvidenceExtractor.ExtractDeposited(new[] { t }, new[] { "c1.raw", "c2.raw", "t1.raw" });

            var bio = One(e, "c2.raw", "characteristics[biological replicate]");
            Assert.That((bio.Value, bio.Method, bio.Confidence), Is.EqualTo(("2", "design-table", SdrfEvidenceConfidence.Certain)));
            Assert.That(e.Select(x => x.Column), Has.None.EqualTo("source name"));
        }

        [Test]
        public void AnMsstatsTmtAnnotationGivesEachChannelItsOwnClaims()
        {
            // MSV000084266's SpikeIn5mix_MS2_PD_annotation.csv, in shape.
            var t = Table("metadata/SpikeIn5mix_MS2_PD_annotation.csv", new[] { "Run", "Fraction", "TechRepMixture", "Channel", "Condition", "Mixture", "BioReplicate" },
                new[] { "161122_TMT10_Mixture1_01.raw", "1", "1", "126", "Norm", "Mixture1", "Norm" },
                new[] { "161122_TMT10_Mixture1_01.raw", "1", "1", "127N", "0.667", "Mixture1", "0.667" });

            var e = SampleEvidenceExtractor.ExtractDeposited(new[] { t }, new[] { "161122_TMT10_Mixture1_01.raw" });

            Assert.That(One(e, "161122_TMT10_Mixture1_01.raw", "factor value[condition]", "TMT127N").Value, Is.EqualTo("0.667"));
            Assert.That(One(e, "161122_TMT10_Mixture1_01.raw", "factor value[condition]", "TMT126").Value, Is.EqualTo("Norm"));
        }

        [Test]
        public void AMaxQuantDesignGroupsFractionsIntoTheirExperiment()
        {
            // MSV000086944's experimental.design.txt, in shape: Name has no extension.
            var t = Table("metadata/experimental.design.txt", new[] { "Name", "Fraction", "Experiment" },
                new[] { "20151114_alr_GM12878_ETD_AspN_fr10", "10", "GM12878_AspN_ETD" },
                new[] { "20151114_alr_GM12878_ETD_AspN_fr11", "11", "GM12878_AspN_ETD" },
                new[] { "20151114_alr_GM12878_HCD_Trypsin_fr1", "1", "GM12878_Trypsin_HCD" });
            var runs = new[] { "20151114_alr_GM12878_ETD_AspN_fr10.raw", "20151114_alr_GM12878_ETD_AspN_fr11.raw", "20151114_alr_GM12878_HCD_Trypsin_fr1.raw" };

            var e = SampleEvidenceExtractor.ExtractDeposited(new[] { t }, runs);

            Assert.That(One(e, runs[1], "source name").Value, Is.EqualTo("GM12878_AspN_ETD"));
            Assert.That(One(e, runs[1], "comment[fraction identifier]").Value, Is.EqualTo("11"));
            Assert.That(One(e, runs[2], "source name").Value, Is.EqualTo("GM12878_Trypsin_HCD"));
        }

        [Test]
        public void AMaxQuantDesignWithoutFractionsIsKnownByItsPtmColumn()
        {
            var t = Table("metadata/experimentalDesignTemplate.txt", new[] { "Name", "Experiment", "PTM" },
                new[] { "alpha", "S1", "False" }, new[] { "beta", "S2", "False" });

            var e = SampleEvidenceExtractor.ExtractDeposited(new[] { t }, new[] { "alpha.raw", "beta.raw" });

            Assert.That(One(e, "beta.raw", "source name").Value, Is.EqualTo("S2"));
        }

        [Test]
        public void ASpectronautRunLabelThatNamesNoRunFallsBackToTheFileNameAndABlankConditionIsNoClaim()
        {
            var t = Table("quant/ConditionSetup.tsv", new[] { "#", "Run Label", "Condition", "Replicate", "File Name" },
                new[] { "1", "Old 1", "Old", "1", "240418_AM_O1" },
                new[] { "2", "Old 2", "", "2", "240418_AM_O2" });

            var e = SampleEvidenceExtractor.ExtractDeposited(new[] { t }, Runs);

            Assert.That(One(e, "240418_AM_O1.raw", "factor value[condition]").Value, Is.EqualTo("Old"));
            Assert.That(One(e, "240418_AM_O2.raw", "characteristics[biological replicate]").Value, Is.EqualTo("2"));
            Assert.That(e.Where(x => x.DataFile == "240418_AM_O2.raw").Select(x => x.Column), Has.None.EqualTo("factor value[condition]"));
        }

        // ---------------- drafting with them ----------------

        [Test]
        public void ADesignTableGivesTheRunsItNamesTheirConditionAndLeavesTheOthersToTheirNames()
        {
            var plain = SdrfDrafter.Draft(Project(), Runs);
            var d = SdrfDrafter.Draft(Project(), Runs, SampleEvidenceExtractor.ExtractDeposited(new[] { Setup() }, Runs));

            Assert.That(d.FactorColumns[0], Is.EqualTo("factor value[condition]"));
            Assert.That(Factor(d, "240418_AM_O2.raw", "factor value[condition]"), Is.EqualTo("Old"));
            Assert.That(Factor(d, "240418_AM_ODR_RF_2.raw", "factor value[condition]"), Is.EqualTo("RF"));
            var cell = Row(d, "240418_AM_O2.raw").Factors[0];
            Assert.That((cell.Source, cell.Reference, cell.Method), Is.EqualTo((SdrfDraftSource.Deposit, SetupFile + "!R3C4", "rules")));

            // The heart runs and the brain run the setup does not name keep what their names said, column for column.
            foreach (var run in Heart.Append("240418_AM_ODR_3.raw"))
                foreach (var column in plain.FactorColumns)
                    Assert.That(Factor(d, run, column), Is.EqualTo(Factor(plain, run, column)), $"{run} {column}");
            // A run the setup describes takes none of the names' other factors.
            foreach (var column in d.FactorColumns.Skip(1))
                Assert.That(Factor(d, "240418_AM_O2.raw", column), Is.EqualTo(SdrfReserved.NotAvailable), column);
        }

        [Test]
        public void ARunTheTableNamesTakesNoneOfTheNamesOtherFactors()
        {
            // The names give two conditions (genotype, treatment); the table states only one, for two of the runs.
            var runs = new[] { "WT_ctrl_1.raw", "WT_drug_1.raw", "KO_ctrl_1.raw", "KO_drug_1.raw", "WT_ctrl_2.raw", "KO_drug_2.raw" };
            var plain = SdrfDrafter.Draft(Project(), runs);
            Assume.That(plain.FactorColumns, Has.Count.EqualTo(2));
            var t = Table("metadata/ConditionSetup.tsv", new[] { "Run Label", "Condition", "Replicate" },
                new[] { "WT_ctrl_1.raw", "baseline", "1" }, new[] { "KO_drug_1.raw", "treated", "1" });

            var d = SdrfDrafter.Draft(Project(), runs, SampleEvidenceExtractor.ExtractDeposited(new[] { t }, runs));

            Assert.That(d.FactorColumns, Has.Count.EqualTo(2), "the names' second column is still the others'");
            Assert.That(Row(d, "WT_ctrl_1.raw").Factors.Select(f => f.Value), Is.EqualTo(new[] { "baseline", SdrfReserved.NotAvailable }));
            Assert.That(Row(d, "WT_drug_1.raw").Factors.Select(f => f.Value), Is.EqualTo(Row(plain, "WT_drug_1.raw").Factors.Select(f => f.Value)));
        }

        [Test]
        public void ADesignTableReplicateReplacesTheDraftsRankAndALikelyClaimDoesNot()
        {
            var evidence = SampleEvidenceExtractor.ExtractDeposited(new[] { Setup() }, Runs).ToList();
            var d = SdrfDrafter.Draft(Project(), Runs, evidence);

            var bio = Row(d, "240418_AM_ODR_RF_4.raw").BiologicalReplicate;
            Assert.That((bio.Value, bio.Source), Is.EqualTo(("6", SdrfDraftSource.Deposit)));

            // The same claim at Likely confidence fills only what the draft does not know.
            var likely = evidence.Select(e => e with { Confidence = SdrfEvidenceConfidence.Likely }).ToList();
            var plain = SdrfDrafter.Draft(Project(), Runs);
            var weak = SdrfDrafter.Draft(Project(), Runs, likely);
            var before = Row(plain, "240418_AM_ODR_RF_4.raw").BiologicalReplicate;
            if (before.Source != SdrfDraftSource.Default)
                Assert.That(Row(weak, "240418_AM_ODR_RF_4.raw").BiologicalReplicate.Value, Is.EqualTo(before.Value), "a reading stands against a Likely claim");
            Assert.That(Factor(weak, "240418_AM_O2.raw", "factor value[condition]"), Is.EqualTo(Factor(plain, "240418_AM_O2.raw", "factor value[condition]")),
                "a Likely condition is not applied");
            Assert.That(weak.EvidenceNotes.Select(n => n.Column), Has.Member("factor value[condition]"), "and says so");
        }

        [Test]
        public void AFractionAndTechnicalReplicateFromADesignTableFillTheDefaultsAndSayDeposited()
        {
            var t = Table("metadata/experimental.design.txt", new[] { "Name", "Fraction", "Experiment" },
                new[] { "alpha", "3", "S1" }, new[] { "beta", "4", "S1" });
            var runs = new[] { "alpha.raw", "beta.raw" };
            Assume.That(SdrfDrafter.Draft(Project(), runs).Rows[0].Fraction.Source, Is.EqualTo(SdrfDraftSource.Default));

            var d = SdrfDrafter.Draft(Project(), runs, SampleEvidenceExtractor.ExtractDeposited(new[] { t }, runs));

            Assert.That((Row(d, "beta.raw").Fraction.Value, Row(d, "beta.raw").Fraction.Source), Is.EqualTo(("4", SdrfDraftSource.Deposit)));
            Assert.That(Row(d, "alpha.raw").SourceName.Value, Is.EqualTo("S1"));
            Assert.That(Row(d, "beta.raw").SourceName.Value, Is.EqualTo("S1"), "both fractions are one sample");

            var row = SdrfDrafter.ToDocument(d).Results.Single(r => r["comment[data file]"] == "beta.raw");
            Assert.That(row["comment[fraction identifier]"], Is.EqualTo("4"));
            Assert.That(row["comment[fraction identifier source]"], Is.EqualTo("deposited"));
            Assert.That(row["comment[fraction identifier source reference]"], Is.EqualTo("metadata/experimental.design.txt!R3C2"));
            Assert.That(row["comment[fraction identifier source method]"], Is.EqualTo("rules"));
        }

        [Test]
        public void ACharacteristicFromADesignTableIsWrittenDepositedWithItsReference()
        {
            var d = SdrfDrafter.Draft(Project(), Runs, SampleEvidenceExtractor.ExtractDeposited(new[] { Setup() }, Runs));

            var doc = SdrfDrafter.ToDocument(d);
            var row = doc.Results.Single(r => r["comment[data file]"] == "240418_AM_O3.raw");
            Assert.That(row["characteristics[biological replicate]"], Is.EqualTo("3"));
            Assert.That(row["comment[biological replicate source]"], Is.EqualTo("deposited"));
            Assert.That(row["comment[biological replicate source reference]"], Is.EqualTo(SetupFile + "!R4C6"));
            Assert.That(SdrfValidator.Validate(doc).Errors, Is.Empty);
        }

        [Test]
        public void TheImproverFillsAGapFromADesignTableAndSaysDeposited()
        {
            string[] columns = { "source name", "characteristics[organism]", "characteristics[biological replicate]", "assay name", "comment[label]", "comment[data file]" };
            var header = new SdrfHeader(columns);
            var deposited = new SdrfDocument(header, new[]
            {
                new SdrfRow(header, new[] { "s1", "mus musculus", "not available", "run 1", "label free sample", "240418_AM_O3.raw" }),
            });
            var draft = SdrfDrafter.Draft(Project(), Runs, SampleEvidenceExtractor.ExtractDeposited(new[] { Setup() }, Runs));

            var row = SdrfImprover.Improve(deposited, draft).Document.Results.Single(r => r["comment[data file]"] == "240418_AM_O3.raw");

            Assert.That(row["characteristics[biological replicate]"], Is.EqualTo("3"));
            Assert.That(row["comment[biological replicate source]"], Is.EqualTo("deposited"));
            Assert.That(row["comment[biological replicate source reference]"], Is.EqualTo(SetupFile + "!R4C6"));
        }

        // ---------------- claims the drafter declines ----------------

        private static readonly string[] TwoRuns = { "alpha.raw", "beta.raw" };

        private static SdrfEvidence Claim(string file, string column, string value) =>
            new(file, "", column, value, "deposit", "metadata/design.tsv!R2C2", "design-table", SdrfEvidenceConfidence.Certain);

        [Test]
        public void AClaimAboutAColumnTheDrafterDoesNotTakeIsOnlyNoted()
        {
            var d = SdrfDrafter.Draft(Project(), TwoRuns, new[] { Claim("alpha.raw", "comment[instrument]", "Q Exactive") });

            Assert.That(d.EvidenceNotes.Single(n => n.Column == "comment[instrument]").Why, Does.StartWith("only characteristics"));
            Assert.That(Row(d, "alpha.raw").Characteristics.Keys, Has.None.EqualTo("comment[instrument]"));
        }

        [Test]
        public void ATechnicalReplicateFromEvidenceFillsTheDefaultAndAFractionThatIsNoNumberIsNoted()
        {
            var d = SdrfDrafter.Draft(Project(), TwoRuns, new[]
            {
                Claim("alpha.raw", "comment[technical replicate]", "2"),
                Claim("beta.raw", "comment[fraction identifier]", "first"),
            });

            var technical = Row(d, "alpha.raw").TechnicalReplicate;
            Assert.That((technical.Value, technical.Source), Is.EqualTo(("2", SdrfDraftSource.Deposit)));
            Assert.That(Row(d, "beta.raw").Fraction.Source, Is.Not.EqualTo(SdrfDraftSource.Deposit));
            Assert.That(d.EvidenceNotes.Single(n => n.DataFile == "beta.raw").Why, Does.Contain("whole number from 1"));
        }

        [Test]
        public void AConditionTheWholeDepositSharesIsNotAFactor()
        {
            var d = SdrfDrafter.Draft(Project(), TwoRuns, new[] { Claim("", "factor value[condition]", "fed") });

            Assert.That(d.FactorColumns, Has.None.EqualTo("factor value[condition]"));
            Assert.That(d.EvidenceNotes.Single(n => n.Column == "factor value[condition]").Why, Does.Contain("whole deposit shares"));
        }

        [Test]
        public void ARunTakesEveryConditionTheTableStatesAndNoneWhereItsClaimsDisagree()
        {
            var d = SdrfDrafter.Draft(Project(), TwoRuns, new[]
            {
                Claim("alpha.raw", "factor value[condition]", "Old"),
                Claim("alpha.raw", "factor value[diet]", "fed"),
                Claim("beta.raw", "factor value[condition]", "Old"),
                Claim("beta.raw", "factor value[condition]", "Young"),
            });

            Assert.That(Factor(d, "alpha.raw", "factor value[condition]"), Is.EqualTo("Old"));
            Assert.That(Factor(d, "alpha.raw", "factor value[diet]"), Is.EqualTo("fed"));
            Assert.That(Factor(d, "beta.raw", "factor value[condition]"), Is.Not.EqualTo("Old").And.Not.EqualTo("Young"));
            Assert.That(d.EvidenceNotes.Single(n => n.DataFile == "beta.raw").Why, Does.Contain("disagree"));
        }
    }
}
