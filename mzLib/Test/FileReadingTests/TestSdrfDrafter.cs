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

        [Test]
        public void ADiseaseFreeProjectIsNormal()
        {
            var p = Covid();
            p.Diseases.Clear();
            p.Diseases.Add(Term("DOID", "DOID:4", "Disease free"));

            var r = Row(SdrfDrafter.Draft(p, new[] { "A_1.raw", "A_2.raw" }), "A_1.raw");

            Assert.That(r.Disease.Value, Is.EqualTo("normal"));
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
