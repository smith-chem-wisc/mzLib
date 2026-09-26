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
    /// Tests for drafting an SDRF for a PRIDE deposit from its project record and raw-file names. Each
    /// rule was measured by blind grading before it was written here (sdrf project, results/benchmark,
    /// 2026-09-23); the deposits below are shaped after ones the graders quoted.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrfDrafter
    {
        private static CvParam Term(string label, string accession, string name) => new(label, accession, name, "");

        private static PrideProject Covid() => new()
        {
            Accession = "PXD020394",
            Title = "Proteomes of oro- and naso-pharyngeal swabs",
            ProjectDescription = "We recovered proteins from SARS-CoV-2 positive and negative swabs.",
            Organisms = { Term("NEWT", "NEWT:1431340", "Sars bat coronavirus"), Term("NEWT", "NEWT:9606", "Homo sapiens (human)") },
            OrganismParts = { Term("CL", "UBERON:0000355", "Pharyngeal mucosa") },
            Diseases = { Term("DOID", "DOID:0080600", "Covid-19") },
            Instruments = { Term("MS", "MS:1001911", "Q Exactive") },
        };

        private static List<string> CovidFiles() => new[] { "NEG", "POS" }
            .SelectMany(arm => Enumerable.Range(1, 5).SelectMany(i => new[] { $"{arm}{i}.raw", $"{arm}{i}rep.raw" })).ToList();

        private static SdrfDraftRow Row(SdrfDraft d, string file) => d.Rows.Single(r => r.DataFile == file);

        [Test]
        public void EveryRawFileGetsExactlyOneRow()
        {
            var files = CovidFiles().Append("Blank.raw").ToList();

            var d = SdrfDrafter.Draft(Covid(), files);

            Assert.That(d.Rows.Select(r => r.DataFile), Is.EquivalentTo(files));
        }

        [Test]
        public void ReinjectionsShareTheirSampleAndCountAsTechnicalReplicates()
        {
            var d = SdrfDrafter.Draft(Covid(), CovidFiles());

            var first = Row(d, "NEG2.raw");
            var again = Row(d, "NEG2rep.raw");
            Assert.That(again.SourceName.Value, Is.EqualTo(first.SourceName.Value));
            Assert.That((first.TechnicalReplicate.Value, again.TechnicalReplicate.Value), Is.EqualTo(("1", "2")));
            Assert.That(d.Rows.Select(r => r.SourceName.Value).Distinct().Count(), Is.EqualTo(10));
        }

        [Test]
        public void BiologicalReplicatesAreNumberedWithinTheirCondition()
        {
            var d = SdrfDrafter.Draft(Covid(), CovidFiles());

            Assert.That(Row(d, "POS4.raw").BiologicalReplicate.Value, Is.EqualTo("4"));
            Assert.That(Row(d, "NEG4.raw").BiologicalReplicate.Value, Is.EqualTo("4"), "each arm counts from 1");
        }

        [Test]
        public void ACaseControlSplitGivesTheControlArmNormalAndTheCasesTheProjectDisease()
        {
            var d = SdrfDrafter.Draft(Covid(), CovidFiles());

            Assert.That(Row(d, "NEG1.raw").Disease.Value, Is.EqualTo("normal"));
            Assert.That(Row(d, "POS1.raw").Disease.Value, Is.EqualTo("Covid-19"));
            Assert.That(Row(d, "POS1.raw").Disease.Source, Is.EqualTo(SdrfDraftSource.Inferred));
        }

        [Test]
        public void AProjectFactIsTakenOnlyWhenPrideListsOneValue()
        {
            var d = SdrfDrafter.Draft(Covid(), CovidFiles());
            var r = Row(d, "POS1.raw");

            Assert.That(r.Organism.Source, Is.EqualTo(SdrfDraftSource.NotAvailable), "PRIDE lists two organisms");
            Assert.That(r.OrganismPart.Value, Is.EqualTo("Pharyngeal mucosa"));
            Assert.That(r.OrganismPart.Source, Is.EqualTo(SdrfDraftSource.PrideProjectRecord));
            Assert.That(r.OrganismPart.Term!.Accession, Is.EqualTo("UBERON:0000355"));
            Assert.That(r.Instrument.Term!.Accession, Is.EqualTo("MS:1001911"));
        }

        [Test]
        public void ANewtOrganismIsWrittenAsNcbiTaxon()
        {
            var p = Covid();
            p.Organisms.RemoveAt(0);

            var r = Row(SdrfDrafter.Draft(p, CovidFiles()), "POS1.raw");

            Assert.That(r.Organism.Term!.Accession, Is.EqualTo("NCBITaxon:9606"));
            Assert.That(r.Organism.Value, Is.EqualTo("homo sapiens"));
        }

        /// <summary>
        /// G31 / pride 006 (SDRF-P6): 17 of 692 PRIDE part and disease terms carry the wrong cvLabel, so the
        /// ontology is read from the accession's prefix. pride's named examples.
        /// </summary>
        [TestCase("UBERON:0002048", "Lung", "UBERON")]    // PXD006179, labelled CL
        [TestCase("BTO:0004102", "Kidney cell", "BTO")]    // PXD008795, labelled CL
        [TestCase("UBERON:0000355", "Pharyngeal mucosa", "UBERON")]
        public void ATermsOntologyIsItsAccessionPrefixNotPridesLabel(string accession, string name, string ontology)
        {
            var p = Covid();
            p.OrganismParts.Clear();
            p.OrganismParts.Add(Term("CL", accession, name));

            var part = Row(SdrfDrafter.Draft(p, CovidFiles()), "POS1.raw").OrganismPart;

            Assert.That(part.Term!.CvLabel, Is.EqualTo(ontology));
            Assert.That(part.Term.Accession, Is.EqualTo(accession));
            Assert.That(part.Value, Is.EqualTo(name));
        }

        /// <summary>
        /// A term from an ontology of the wrong KIND is not written in the column: an anatomy or sequence term
        /// is not a disease (PXD031485 lists UBERON:0000474 as its disease; pride names SO:0000337 in
        /// PXD010991), and a cellular component or a chemical is not an organism part (PXD013322 GO:0005634).
        /// </summary>
        [TestCase("disease", "UBERON", "UBERON:0000474", "Female reproductive system")]
        [TestCase("disease", "DOID", "SO:0000337", "Rnai_reagent")]
        [TestCase("organism part", "CL", "GO:0005634", "Nucleus")]
        [TestCase("organism part", "CL", "CHEBI:15377", "Water")]
        public void ATermFromTheWrongKindOfOntologyIsNotAvailableAndSaysWhy(string column, string label, string accession, string name)
        {
            var p = Covid();
            var list = column == "disease" ? p.Diseases : p.OrganismParts;
            list.Clear();
            list.Add(Term(label, accession, name));

            var row = Row(SdrfDrafter.Draft(p, CovidFiles()), "POS1.raw");
            var cell = column == "disease" ? row.Disease : row.OrganismPart;

            Assert.That(cell.Source, Is.EqualTo(SdrfDraftSource.NotAvailable));
            Assert.That(cell.Evidence, Does.Contain(accession).And.Contain(accession.Split(':')[0]));
        }

        [TestCase("SYMP:0000292", "Heart failure")]    // PXD007171: a symptom is written as the disease
        [TestCase("HP:0001627", "Abnormal heart morphology")]
        public void APhenotypeOrSymptomTermIsStillADisease(string accession, string name)
        {
            var p = Covid();
            p.Diseases.Clear();
            p.Diseases.Add(Term("DOID", accession, name));

            var disease = Row(SdrfDrafter.Draft(p, CovidFiles()), "POS1.raw").Disease;

            Assert.That(disease.Value, Is.EqualTo(name));
            Assert.That(disease.Term!.CvLabel, Is.EqualTo(accession.Split(':')[0]));
        }

        [Test]
        public void ADiseaseFreeProjectIsNormal()
        {
            var p = Covid();
            p.Diseases.Clear();
            p.Diseases.Add(Term("DOID", "DOID:4", "Disease free"));

            var r = Row(SdrfDrafter.Draft(p, new[] { "A_1.raw", "A_2.raw" }), "A_1.raw");

            Assert.That(r.Disease.Value, Is.EqualTo("normal"));
        }

        /// <summary>
        /// A file with no known condition is not ranked. Its only company in "no condition" is every other
        /// unexplained file, so ranking would count across the deposit, which a draft never does.
        /// PXD009495 numbered unrelated samples 1..5 this way.
        /// </summary>
        [Test]
        public void AFileWithNoKnownConditionIsNotRankedWithTheOtherUnknowns()
        {
            var files = new[] { "WT_1.raw", "WT_2.raw", "KO_1.raw", "KO_2.raw", "QC_pool_x.raw", "QC_pool_y.raw", "QC_pool_z.raw" };

            var d = SdrfDrafter.Draft(Covid(), files);

            Assert.That(d.FactorColumns, Is.Not.Empty, "WT and KO are a condition");
            Assert.That(Row(d, "KO_2.raw").BiologicalReplicate.Value, Is.EqualTo("2"));
            foreach (var qc in files.Where(f => f.StartsWith("QC", StringComparison.Ordinal)))
                Assert.That(Row(d, qc).BiologicalReplicate.Value, Is.EqualTo("1"),
                    qc + ": " + Row(d, qc).BiologicalReplicate.Evidence);
        }

        [Test]
        public void TheFactorIsTheConditionTheNamesAndRecordShare()
        {
            var d = SdrfDrafter.Draft(Covid(), CovidFiles());

            Assert.That(d.FactorColumns, Is.EqualTo(new[] { "factor value[condition]" }));
            Assert.That(Row(d, "NEG3.raw").Factors[0].Value, Is.EqualTo("NEG").IgnoreCase);
            Assert.That(Row(d, "NEG3.raw").Factors[0].Evidence, Is.Not.Empty);
        }

        [Test]
        public void TheRecordCanMakeAMarkerTechnical()
        {
            var p = Covid();
            p.OrganismParts.Clear();
            p.SampleProcessingProtocol = "A single HeLa lysate was analyzed in triplicate.";
            var files = new[] { "HeLa_1.raw", "HeLa_2.raw", "HeLa_3.raw" };

            var d = SdrfDrafter.Draft(p, files);

            Assert.That(d.Rows.Select(r => r.SourceName.Value).Distinct().Count(), Is.EqualTo(1));
            Assert.That(d.Rows.Select(r => r.TechnicalReplicate.Value), Is.EqualTo(new[] { "1", "2", "3" }));
            Assert.That(d.Rows.Select(r => r.BiologicalReplicate.Value), Is.All.EqualTo("1"));
        }

        [Test]
        public void WithNoStructureEveryFileIsItsOwnSampleAndReplicateOne()
        {
            var p = Covid();
            p.ProjectDescription = "";
            p.Title = "";
            var files = new[] { "alpha.raw", "gamma.raw", "omega.raw" };

            var d = SdrfDrafter.Draft(p, files);

            Assert.That(d.Rows.Select(r => r.SourceName.Value).Distinct().Count(), Is.EqualTo(3));
            Assert.That(d.Rows.Select(r => r.BiologicalReplicate.Value), Is.All.EqualTo("1"), "never a count across the deposit");
            Assert.That(d.FactorColumns, Is.Empty);
        }

        /// <summary>
        /// Found by improving the whole curated corpus: two families of names with different numbers of
        /// factors gave rows fewer factor cells than the draft had factor columns.
        /// </summary>
        [Test]
        public void EveryRowHasOneCellPerFactorColumnWhenFamiliesDiffer()
        {
            var p = Covid();
            p.ProjectDescription = "";
            var files = new List<string> { "Blank_std.raw", "Blank_std2.raw" };
            files.AddRange(new[] { "WT", "KO" }.SelectMany(g => new[] { "A", "B" }.SelectMany(t => Enumerable.Range(1, 2).Select(i => $"{g}_{t}_{i}.raw"))));
            files.AddRange(new[] { "X", "Y" }.SelectMany(g => Enumerable.Range(1, 2).Select(i => $"20200101_run_{g}{i}.raw")));

            var d = SdrfDrafter.Draft(p, files);

            Assert.That(d.Rows.All(r => r.Factors.Count == d.FactorColumns.Count));
            var wt = Row(d, "WT_A_1.raw");
            Assert.That(wt.Factors.Count(f => f.Source != SdrfDraftSource.NotAvailable), Is.EqualTo(2), "its own family's two factors");
        }

        [Test]
        public void ASidecarFileIsNotARow()
        {
            var d = SdrfDrafter.Draft(Covid(), new[] { "A_1.wiff", "A_1.wiff.scan", "A_2.wiff", "A_2.wiff.scan" });

            Assert.That(d.Rows.Select(r => r.DataFile), Is.EquivalentTo(new[] { "A_1.wiff", "A_2.wiff" }));
        }

        [Test]
        public void AValueWrittenOnlyBecauseSdrfNeedsOneIsMarkedDefault()
        {
            var d = SdrfDrafter.Draft(Covid(), CovidFiles().Append("Blank.raw"));

            Assert.That(Row(d, "POS1.raw").Fraction.Source, Is.EqualTo(SdrfDraftSource.Default), "no fraction marker: 1 is a default, not a reading");
            Assert.That(Row(d, "POS1.raw").TechnicalReplicate.Source, Is.EqualTo(SdrfDraftSource.Inferred), "POS1rep exists: POS1 is the first injection");
            Assert.That(Row(d, "Blank.raw").TechnicalReplicate.Source, Is.EqualTo(SdrfDraftSource.Default), "no re-injection: 1 is a default");
        }

        [Test]
        public void ABiologicalReplicateNothingMarkedIsLabelledDefaultNotTheRowSource()
        {
            // D39: a replicate 1 written only because an SDRF needs one is `default`. Under D31's grain the row's
            // comment[characteristics source] covers every characteristic without its own override, so without
            // one a default replicate would read as "pride project record".
            var draft = SdrfDrafter.Draft(Covid(), CovidFiles().Append("Blank.raw"));
            var doc = SdrfDrafter.ToDocument(draft, "PXD020394");

            string Effective(SdrfRow r) => doc.Header.Contains("comment[biological replicate source]")
                && r["comment[biological replicate source]"] is { } o && o != "not applicable" ? o : r["comment[characteristics source]"];
            foreach (var row in draft.Rows)
            {
                var written = doc.Results.Single(r => r["comment[data file]"] == row.DataFile);
                string expected = row.BiologicalReplicate.Source == SdrfDraftSource.Default ? "default" : "inferred";
                Assert.That(Effective(written), Is.EqualTo(expected), row.DataFile);
            }
            Assert.That(draft.Rows.Any(r => r.BiologicalReplicate.Source == SdrfDraftSource.Default), "the fixture must carry a default");
            Assert.That(doc.Header.Count(h => h == "comment[characteristics source]"), Is.EqualTo(1));
        }

        [Test]
        public void EveryRowSaysWhereItsFractionAndTechnicalReplicateCameFrom()
        {
            // dataRepo SDRF-DR10: they pool runs by fraction and technical replicate, and a defaulted 1 must say so.
            var draft = SdrfDrafter.Draft(Covid(), CovidFiles().Append("Blank.raw"));
            var doc = SdrfDrafter.ToDocument(draft, "PXD020394");

            foreach (var row in draft.Rows)
            {
                var written = doc.Results.Single(r => r["comment[data file]"] == row.DataFile);
                Assert.That(written["comment[fraction identifier source]"], Is.EqualTo(row.Fraction.Source == SdrfDraftSource.Default ? "default" : "inferred"), row.DataFile);
                Assert.That(written["comment[technical replicate source]"], Is.EqualTo(row.TechnicalReplicate.Source == SdrfDraftSource.Default ? "default" : "inferred"), row.DataFile);
            }
            var blank = doc.Results.Single(r => r["comment[data file]"] == "Blank.raw");
            Assert.That((blank["comment[fraction identifier source]"], blank["comment[technical replicate source]"]), Is.EqualTo(("default", "default")));
            var pos1 = doc.Results.Single(r => r["comment[data file]"] == "POS1.raw");
            Assert.That(pos1["comment[technical replicate source]"], Is.EqualTo("inferred"), "POS1rep exists: the injection was read");
        }

        [Test]
        public void EveryInferredCellSaysWhy()
        {
            var d = SdrfDrafter.Draft(Covid(), CovidFiles());

            var cells = d.Rows.SelectMany(r => new[] { r.SourceName, r.BiologicalReplicate, r.TechnicalReplicate, r.Fraction, r.Disease }.Concat(r.Factors));
            Assert.That(cells.Where(c => c.Source == SdrfDraftSource.Inferred).All(c => c.Evidence.Length > 0));
        }

        // ---- writing the draft as an SDRF ----

        [Test]
        public void ADraftBecomesAValidSdrfWithOneRowPerFile()
        {
            var doc = SdrfDrafter.ToDocument(SdrfDrafter.Draft(Covid(), CovidFiles()), "PXD020394");

            Assert.That(doc.Results.Count, Is.EqualTo(20));
            Assert.That(doc.Results.Select(r => r["comment[data file]"]), Is.EquivalentTo(CovidFiles()));
            Assert.That(doc.Header, Does.Contain("factor value[condition]"));
            var v = SdrfValidator.Validate(doc);
            Assert.That(v.Errors, Is.Empty, string.Join("\n", v.Errors.Select(e => e.ToString())));
        }

        /// <summary>
        /// G37: a file the drafter could not place has an UNKNOWN condition -- `not available` -- not a condition
        /// that does not apply to it. PXD032040's shape: IgG / ctrl IPs beside a batch of runs with no condition token.
        /// </summary>
        [Test]
        public void AFactorTheDrafterCouldNotPlaceIsNotAvailableNotNotApplicable()
        {
            var files = new[] { "IP_IgG1.raw", "IP_IgG2.raw", "IP_IgG3.raw", "IP_ctrl1.raw", "IP_ctrl2.raw", "IP_ctrl3.raw",
                                "Q1_ColID_292_10.raw", "Q1_ColID_292_11.raw", "Q1_ColID_292_12.raw" };
            var draft = SdrfDrafter.Draft(Covid(), files);
            Assume.That(draft.FactorColumns, Is.Not.Empty, "the fixture must yield a condition");
            var doc = SdrfDrafter.ToDocument(draft, "PXD032040");

            var unplaced = doc.Results.Single(r => r["comment[data file]"] == "Q1_ColID_292_10.raw");
            var placed = doc.Results.Single(r => r["comment[data file]"] == "IP_IgG1.raw");

            Assert.That(unplaced[draft.FactorColumns[0]], Is.EqualTo("not available"));
            Assert.That(placed[draft.FactorColumns[0]], Is.EqualTo("IgG"));
            Assert.That(SdrfValidator.Validate(doc).Errors, Is.Empty);
        }

        [Test]
        public void ProvenanceIsARowDefaultWithAnOverrideOnlyWhereACellDiffers()
        {
            var doc = SdrfDrafter.ToDocument(SdrfDrafter.Draft(Covid(), CovidFiles()), "PXD020394");
            var neg = doc.Results.Single(r => r["comment[data file]"] == "NEG1.raw");

            Assert.That(neg["comment[characteristics source]"], Is.EqualTo("pride project record"));
            Assert.That(neg["comment[disease source]"], Is.EqualTo("inferred"), "the control arm's 'normal' was inferred");
            Assert.That(neg["characteristics[disease]"], Does.Contain("AC=PATO:0000461"));
            var h = doc.Header.ToList();
            Assert.That(h.IndexOf("comment[characteristics source]"), Is.GreaterThan(h.IndexOf("assay name")), "comment columns after assay name (N10)");
        }

        [Test]
        public void ACellNothingStatesIsNotAvailableAndCarriesNoSource()
        {
            var doc = SdrfDrafter.ToDocument(SdrfDrafter.Draft(Covid(), CovidFiles()), "PXD020394");
            var row = doc.Results.Single(r => r["comment[data file]"] == "POS1.raw");

            Assert.That(row["characteristics[organism]"], Is.EqualTo("not available"), "PRIDE lists two organisms");
            Assert.That(row["assay name"], Is.EqualTo("run POS1"));
        }

        [Test]
        public void MalformedArgumentsThrow()
        {
            Assert.Throws<ArgumentNullException>(() => SdrfDrafter.Draft(null!, new[] { "a.raw" }));
            Assert.Throws<ArgumentNullException>(() => SdrfDrafter.Draft(Covid(), null!));
        }
    }
}
