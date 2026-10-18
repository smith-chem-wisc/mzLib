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
    /// Tests for drafting from a repository-neutral deposit record (sdrf MASSIVE-SDRF.md, MV1): a MassIVE dataset record
    /// drafts the way a PRIDE project record does, from what it has, and its provenance names MassIVE, never PRIDE. The
    /// MassIVE records below are shaped after MSV000096226 (PXD057243): a title, a summary, keywords, one species and one
    /// instrument, and no protocols, organism parts or diseases.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrfDepositRecord
    {
        private static CvParam Term(string label, string accession, string name) => new(label, accession, name, "");

        private static readonly CvParam Mouse = Term("NCBITaxon", "NCBITaxon:10090", "Mus musculus");
        private static readonly CvParam Human = Term("NCBITaxon", "NCBITaxon:9606", "Homo sapiens");
        private static readonly CvParam Exploris = Term("MS", "MS:1003028", "Orbitrap Exploris 480");

        private static SdrfDepositRecord Massive(params CvParam[] organisms) => new(SdrfRepository.Massive)
        {
            Title = "Old mouse brain proteome under calorie restriction and refeeding",
            Description = "Brains of old mice, ad libitum (Old) against refed after restriction (RF).",
            Keywords = new[] { "brain", "aging", "dietary restriction" },
            Organisms = organisms,
            Instruments = new[] { Exploris },
        };

        private static readonly string[] Runs = { "Old_1.raw", "Old_2.raw", "RF_1.raw", "RF_2.raw" };

        private static SdrfDraftRow Row(SdrfDraft d, string file) => d.Rows.Single(r => r.DataFile == file);

        /// <summary>Every cell a draft row carries, with the characteristics beyond the fixed ones.</summary>
        private static IEnumerable<SdrfDraftCell> Cells(SdrfDraftRow r) =>
            new[] { r.SourceName, r.Organism, r.OrganismPart, r.Disease, r.Instrument, r.BiologicalReplicate, r.TechnicalReplicate, r.Fraction }
                .Concat(r.Factors).Concat(r.Characteristics.Values);

        [Test]
        public void AMassiveRecordDraftsOneRowPerRunFromWhatItHas()
        {
            var d = SdrfDrafter.Draft(Massive(Mouse), Runs);

            Assert.That(d.Rows.Select(r => r.DataFile), Is.EquivalentTo(Runs));
            var old1 = Row(d, "Old_1.raw");
            Assert.That(old1.Organism.Value, Is.EqualTo("mus musculus"));
            Assert.That(old1.Organism.Term!.Accession, Is.EqualTo("NCBITaxon:10090"));
            Assert.That(old1.Instrument.Value, Is.EqualTo("Orbitrap Exploris 480"));
            // MassIVE's record has neither, so neither is claimed.
            Assert.That(old1.OrganismPart.Source, Is.EqualTo(SdrfDraftSource.NotAvailable));
            Assert.That(old1.Disease.Source, Is.EqualTo(SdrfDraftSource.NotAvailable));
            // The names still carry the design.
            Assert.That(d.FactorColumns, Has.Count.EqualTo(1));
            Assert.That(Row(d, "RF_2.raw").Factors.Single().Value, Is.EqualTo("RF"));
        }

        [Test]
        public void AValueFromAMassiveRecordSaysSo()
        {
            var d = SdrfDrafter.Draft(Massive(Mouse), Runs);

            var old1 = Row(d, "Old_1.raw");
            Assert.That(old1.Organism.Source, Is.EqualTo(SdrfDraftSource.MassiveDatasetRecord));
            Assert.That(old1.Instrument.Source, Is.EqualTo(SdrfDraftSource.MassiveDatasetRecord));
            Assert.That(old1.Organism.Evidence, Is.EqualTo("the one organism MassIVE's dataset record lists"));

            var doc = SdrfDrafter.ToDocument(d, "PXD057243");
            foreach (var r in doc.Results)
            {
                Assert.That(r["comment[characteristics source]"], Is.EqualTo("massive dataset record"), r["comment[data file]"]);
                Assert.That(r.Cells, Has.None.Contain("pride"), "no written cell names PRIDE");
            }
        }

        /// <summary>Each case reaches one of the drafter's reasons that name the record.</summary>
        private static IEnumerable<TestCaseData> RecordReasons()
        {
            yield return new TestCaseData(Massive(Mouse), Runs).SetName("one organism");
            yield return new TestCaseData(Massive(), Runs).SetName("no organism");
            yield return new TestCaseData(Massive(Mouse, Human), Runs).SetName("several organisms");
            yield return new TestCaseData(Massive(Mouse, Human), new[] { "Mouse_1.raw", "Mouse_2.raw", "Human_1.raw", "Human_2.raw" })
                .SetName("an organism the name gives");
            yield return new TestCaseData(Massive(Mouse) with { Diseases = new[] { Term("UBERON", "UBERON:0000955", "brain") } }, Runs)
                .SetName("a disease of the wrong kind");
            yield return new TestCaseData(Massive(Mouse) with { Diseases = new[] { Term("DOID", "DOID:4", "Disease free") } }, Runs)
                .SetName("disease free");
        }

        [TestCaseSource(nameof(RecordReasons))]
        public void NoReasonInAMassiveDraftNamesPride(object record, string[] runs)
        {
            // object: a public test method cannot take the internal record type.
            var d = SdrfDrafter.Draft((SdrfDepositRecord)record, runs);

            var cells = d.Rows.SelectMany(Cells).ToList();
            Assert.That(cells.Select(c => c.Source), Has.None.EqualTo(SdrfDraftSource.PrideProjectRecord));
            Assert.That(cells.Select(c => c.Evidence), Has.None.Contain("PRIDE"));
        }

        [Test]
        public void ARowDefaultTiedBetweenANameAndTheMassiveRecordGoesToTheRecord()
        {
            // Organism read off the name (inferred), organism part from the record: one each. A tie goes to the record,
            // as it does for PRIDE's, and the organism gets its own word.
            var record = Massive(Mouse, Human) with { OrganismParts = new[] { Term("UBERON", "UBERON:0000955", "brain") } };
            var d = SdrfDrafter.Draft(record, new[] { "Mouse_1.raw", "Mouse_2.raw", "Human_1.raw", "Human_2.raw" });

            var mouse = SdrfDrafter.ToDocument(d).Results.Single(r => r["comment[data file]"] == "Mouse_1.raw");
            Assert.That(mouse["comment[characteristics source]"], Is.EqualTo("massive dataset record"));
            Assert.That(mouse["comment[organism source]"], Is.EqualTo("inferred"));
        }

        [Test]
        public void AMassiveRecordWithNoSpeciesLeavesTheOrganismNotAvailableAndSaysWhy()
        {
            // Old MassIVE records give "taxonomy: common name = null" (MSV000078494): no species at all.
            var d = SdrfDrafter.Draft(Massive(), Runs);

            var organism = Row(d, "Old_1.raw").Organism;
            Assert.That(organism.Value, Is.EqualTo(SdrfReserved.NotAvailable));
            Assert.That(organism.Evidence, Is.EqualTo("MassIVE's dataset record lists no organism"));
        }

        [Test]
        public void AMassiveRecordFactFillsAGapButNeverDisputesAStatedValue()
        {
            string[] columns = { "source name", "characteristics[organism]", "assay name", "comment[label]", "comment[data file]", "comment[instrument]" };
            var header = new SdrfHeader(columns);
            var deposited = new SdrfDocument(header, new[]
            {
                new SdrfRow(header, new[] { "o1", "not available", "run 1", "label free sample", "Old_1.raw", "NT=Q Exactive;AC=MS:1001911" }),
                new SdrfRow(header, new[] { "o2", "not available", "run 2", "label free sample", "Old_2.raw", "NT=Q Exactive;AC=MS:1001911" }),
            });

            var i = SdrfImprover.Improve(deposited, SdrfDrafter.Draft(Massive(Mouse), Runs));

            // A deposit-wide summary is weaker evidence than a per-file statement: the instrument stands, undisputed.
            Assert.That(i.Disagreements.Where(x => x.Column == "comment[instrument]"), Is.Empty);
            var old1 = i.Document.Results.Single(r => r["comment[data file]"] == "Old_1.raw");
            Assert.That(old1["comment[instrument]"], Does.Contain("MS:1001911"));
            // The gap is filled, and its source is MassIVE's record.
            Assert.That(old1["characteristics[organism]"], Does.Contain("10090"));
            Assert.That(old1["comment[organism source]"], Is.EqualTo("massive dataset record"));
        }

        [Test]
        public void APrideProjectDraftsExactlyAsItsRecordDoes()
        {
            var project = new PrideProject
            {
                Accession = "PXD020394",
                Title = "Proteomes of oro- and naso-pharyngeal swabs",
                ProjectDescription = "We recovered proteins from SARS-CoV-2 positive and negative swabs.",
                SampleProcessingProtocol = "Three biological replicates per arm.",
                DataProcessingProtocol = "Searched with MetaMorpheus.",
                Keywords = { "COVID-19" },
                ProjectTags = { "Biomedical" },
                Organisms = { Term("NEWT", "NEWT:9606", "Homo sapiens (human)") },
                OrganismParts = { Term("CL", "UBERON:0000355", "Pharyngeal mucosa") },
                Diseases = { Term("DOID", "DOID:0080600", "Covid-19") },
                Instruments = { Term("MS", "MS:1001911", "Q Exactive") },
            };
            var files = new[] { "NEG1.raw", "NEG2.raw", "POS1.raw", "POS2.raw" };

            var record = SdrfDepositRecord.FromPride(project);
            Assert.That(record.Repository, Is.EqualTo(SdrfRepository.Pride));
            Assert.That((record.Title, record.Description, record.SampleProcessingProtocol, record.DataProcessingProtocol),
                Is.EqualTo((project.Title, project.ProjectDescription, project.SampleProcessingProtocol, project.DataProcessingProtocol)));
            Assert.That(record.Keywords, Is.EqualTo(project.Keywords));
            Assert.That(record.ProjectTags, Is.EqualTo(project.ProjectTags));
            Assert.That(record.Organisms, Is.EqualTo(project.Organisms));
            Assert.That(record.OrganismParts, Is.EqualTo(project.OrganismParts));
            Assert.That(record.Diseases, Is.EqualTo(project.Diseases));
            Assert.That(record.Instruments, Is.EqualTo(project.Instruments));

            static string Written(SdrfDraft d)
            {
                var path = System.IO.Path.Combine(System.IO.Path.GetTempPath(), Guid.NewGuid() + ".sdrf.tsv");
                try { SdrfDrafter.ToDocument(d).WriteResults(path); return System.IO.File.ReadAllText(path); }
                finally { System.IO.File.Delete(path); }
            }
            Assert.That(Written(SdrfDrafter.Draft(project, files)), Is.EqualTo(Written(SdrfDrafter.Draft(record, files))));
            Assert.That(Row(SdrfDrafter.Draft(project, files), "POS1.raw").Organism.Source, Is.EqualTo(SdrfDraftSource.PrideProjectRecord));
        }

        [Test]
        public void MalformedArgumentsThrow()
        {
            Assert.Throws<ArgumentNullException>(() => SdrfDrafter.Draft((SdrfDepositRecord)null!, Runs));
            Assert.Throws<ArgumentNullException>(() => SdrfDrafter.Draft(Massive(Mouse), null!));
            Assert.Throws<ArgumentNullException>(() => SdrfDepositRecord.FromPride(null!));
        }
    }
}
