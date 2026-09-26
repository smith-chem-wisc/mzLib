using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MzLibUtil;
using NUnit.Framework;
using Readers;
using UsefulProteomicsDatabases;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Publication evidence in the drafter (sdrf design SAMPLE-EVIDENCE.md, E1; D41 all in mzLib, D42 rules first).
    /// Evidence fills what the drafter left unknown or defaulted and adds characteristics it never writes; it never
    /// overrides a reading, and every filled cell says it came from the publication and where.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrfEvidence
    {
        private static CvParam Term(string label, string accession, string name) => new(label, accession, name, "");

        /// <summary>One organism and one instrument from PRIDE; no organism part, no disease.</summary>
        private static PrideProject Project() => new()
        {
            Accession = "PXD060431",
            Title = "Human heart proteomes",
            Organisms = { Term("NEWT", "NEWT:9606", "Homo sapiens (human)") },
            Instruments = { Term("MS", "MS:1001911", "Q Exactive") },
        };

        private static readonly string[] Files = { "HumanHFpEF_1.raw", "HumanHFpEF_2.raw", "HumanControl_1.raw" };

        private static SdrfEvidence Claim(string file, string column, string value, string label = "",
            SdrfEvidenceConfidence confidence = SdrfEvidenceConfidence.Likely, string locator = "mmc1.xlsx!DatasetS1!R2C4") =>
            new(file, label, column, value, "supplement", locator, "file-key", confidence);

        private static SdrfDraftRow Row(SdrfDraft d, string file) => d.Rows.Single(r => r.DataFile == file);

        private static SdrfRow Written(SdrfDocument doc, string file) => doc.Results.Single(r => r["comment[data file]"] == file);

        /// <summary>The source that applies to one characteristic: its own override, else the row default (D31).</summary>
        private static string SourceOf(SdrfDocument doc, SdrfRow row, string name) =>
            doc.Header.Contains($"comment[{name} source]") && row[$"comment[{name} source]"] is { } own && own != "not applicable"
                ? own : row["comment[characteristics source]"];

        [Test]
        public void EvidenceAddsTheCharacteristicsTheDrafterNeverWrites()
        {
            var evidence = new[]
            {
                Claim("HumanHFpEF_1.raw", "characteristics[age]", "77Y"),
                Claim("HumanHFpEF_1.raw", "characteristics[sex]", "female"),
                Claim("HumanHFpEF_2.raw", "characteristics[age]", "73Y"),
                Claim("HumanHFpEF_2.raw", "characteristics[sex]", "male"),
            };
            var draft = SdrfDrafter.Draft(Project(), Files, evidence);

            var cell = Row(draft, "HumanHFpEF_1.raw").Characteristics!["characteristics[age]"];
            Assert.That((cell.Value, cell.Source), Is.EqualTo(("77Y", SdrfDraftSource.Publication)));
            Assert.That(cell.Reference, Is.EqualTo("mmc1.xlsx!DatasetS1!R2C4"));
            Assert.That(cell.Evidence, Does.Contain("supplement"));

            var doc = SdrfDrafter.ToDocument(draft, "PXD060431");
            var one = Written(doc, "HumanHFpEF_1.raw");
            Assert.That(one["characteristics[age]"], Is.EqualTo("77Y"));
            Assert.That(one["characteristics[sex]"], Is.EqualTo("female"));
            Assert.That(SourceOf(doc, one, "age"), Is.EqualTo("publication"));
            Assert.That(SourceOf(doc, one, "organism"), Is.EqualTo("pride project record"), "the organism still came from PRIDE");
            Assert.That(one["comment[age source reference]"], Is.EqualTo("mmc1.xlsx!DatasetS1!R2C4"));
            Assert.That(one["comment[age source method]"], Is.EqualTo("rules"), "aging 038: how it was read is its own column");
            Assert.That(Written(doc, "HumanControl_1.raw")["characteristics[age]"], Is.EqualTo("not available"), "no claim, no value");
            Assert.That(SdrfValidator.Validate(doc).Errors, Is.Empty);
        }

        [TestCase("model", "model")]
        [TestCase("curator", "curator")]
        [TestCase("channel-map", "rules")]
        [TestCase("isa-tab", "rules")]
        public void TheSourceMethodSaysWhetherRulesACuratorOrAModelReadIt(string method, string word)
        {
            var claim = new SdrfEvidence("HumanHFpEF_1.raw", "", "characteristics[age]", "77Y", "paper", "PMC1/sec:methods", method, SdrfEvidenceConfidence.Likely);
            var doc = SdrfDrafter.ToDocument(SdrfDrafter.Draft(Project(), Files, new[] { claim }), "PXD060431");

            Assert.That(Written(doc, "HumanHFpEF_1.raw")["comment[age source method]"], Is.EqualTo(word));
        }

        [Test]
        public void EvidenceFillsAnUnknownCellButNeverOverridesAReading()
        {
            var p = Project();
            p.Diseases.Add(Term("DOID", "DOID:0050700", "Cardiomyopathy"));
            var evidence = new[]
            {
                Claim("HumanHFpEF_1.raw", "characteristics[organism part]", "heart left ventricle"),
                Claim("HumanHFpEF_1.raw", "characteristics[disease]", "heart failure with preserved ejection fraction"),
            };
            var draft = SdrfDrafter.Draft(p, Files, evidence);
            var row = Row(draft, "HumanHFpEF_1.raw");

            Assert.That((row.OrganismPart.Value, row.OrganismPart.Source), Is.EqualTo(("heart left ventricle", SdrfDraftSource.Publication)));
            Assert.That(row.Disease.Source, Is.EqualTo(SdrfDraftSource.PrideProjectRecord), "a reading is never overridden");
            var note = draft.EvidenceNotes.Single(n => n.Column == "characteristics[disease]");
            Assert.That((note.DataFile, note.EvidenceValue), Is.EqualTo(("HumanHFpEF_1.raw", "heart failure with preserved ejection fraction")));

            var doc = SdrfDrafter.ToDocument(draft, "PXD060431");
            Assert.That(Written(doc, "HumanHFpEF_1.raw")["characteristics[organism part]"], Is.EqualTo("heart left ventricle"));
            Assert.That(SdrfValidator.Validate(doc).Errors, Is.Empty);
        }

        [Test]
        public void ADefaultedReplicateTakesEvidenceAndAReadOneKeepsItsReading()
        {
            // No structure in these names, so the drafter defaults replicate 1 (D39).
            var bare = new[] { "alpha.raw", "beta.raw" };
            Assume.That(SdrfDrafter.Draft(Project(), bare).Rows[0].BiologicalReplicate.Source, Is.EqualTo(SdrfDraftSource.Default));
            var bio = Row(SdrfDrafter.Draft(Project(), bare, new[] { Claim("alpha.raw", "characteristics[biological replicate]", "3") }),
                "alpha.raw").BiologicalReplicate;
            Assert.That((bio.Value, bio.Source), Is.EqualTo(("3", SdrfDraftSource.Publication)));

            // HumanControl_1's replicate was read off its name: evidence does not override it, and says so.
            var draft = SdrfDrafter.Draft(Project(), Files, new[] { Claim("HumanControl_1.raw", "characteristics[biological replicate]", "3") });
            Assert.That(Row(draft, "HumanControl_1.raw").BiologicalReplicate.Source, Is.EqualTo(SdrfDraftSource.Inferred));
            Assert.That(draft.EvidenceNotes.Single().Column, Is.EqualTo("characteristics[biological replicate]"));
        }

        [Test]
        public void AFileClaimBeatsADepositClaimAndADepositClaimReachesEveryFile()
        {
            var evidence = new[]
            {
                Claim("", "characteristics[sex]", "female"),
                Claim("HumanHFpEF_2.raw", "characteristics[sex]", "male"),
            };
            var draft = SdrfDrafter.Draft(Project(), Files, evidence);

            Assert.That(Files.Select(f => Row(draft, f).Characteristics!["characteristics[sex]"].Value),
                Is.EqualTo(new[] { "female", "male", "female" }));
        }

        [Test]
        public void GuessesConflictsAndChannelClaimsAreNotApplied()
        {
            var evidence = new[]
            {
                Claim("HumanHFpEF_1.raw", "characteristics[age]", "70Y", confidence: SdrfEvidenceConfidence.Guess),
                Claim("HumanHFpEF_2.raw", "characteristics[sex]", "male"),
                Claim("HumanHFpEF_2.raw", "characteristics[sex]", "female", locator: "mmc2.xlsx!S1!R9C2"),
                Claim("HumanControl_1.raw", "source name", "patient 7", label: "TMT127N"),
            };
            var draft = SdrfDrafter.Draft(Project(), Files, evidence);

            Assert.That(Row(draft, "HumanHFpEF_1.raw").Characteristics, Is.Empty, "a guess is for review only");
            Assert.That(Row(draft, "HumanHFpEF_2.raw").Characteristics, Is.Empty, "two claims that disagree fill nothing");
            Assert.That(draft.EvidenceNotes.Select(n => n.Column),
                Is.SupersetOf(new[] { "characteristics[sex]", "source name" }), "the conflict and the channel claim are reported");
        }

        [Test]
        public void WithoutEvidenceTheDraftIsUnchanged()
        {
            var plain = SdrfDrafter.Draft(Project(), Files);
            var empty = SdrfDrafter.Draft(Project(), Files, Array.Empty<SdrfEvidence>());

            Assert.That(empty.Rows.Select(r => (r.DataFile, r.SourceName, r.OrganismPart, r.BiologicalReplicate)),
                Is.EqualTo(plain.Rows.Select(r => (r.DataFile, r.SourceName, r.OrganismPart, r.BiologicalReplicate))));
            Assert.That(plain.EvidenceNotes, Is.Empty);
        }

        [Test]
        public void TheImproverFillsAGapFromEvidenceAndSaysPublication()
        {
            var header = new SdrfHeader(new[] { "source name", "characteristics[organism]", "characteristics[organism part]",
                "characteristics[biological replicate]", "assay name", "comment[label]", "comment[data file]", "comment[technical replicate]" });
            var deposited = new SdrfDocument(header, new[]
            {
                new SdrfRow(header, new[] { "s1", "homo sapiens", "not available", "1", "run 1", "label free sample", "HumanHFpEF_1.raw", "1" })
            });
            var draft = SdrfDrafter.Draft(Project(), Files,
                new[] { Claim("HumanHFpEF_1.raw", "characteristics[organism part]", "heart left ventricle") });

            var improved = SdrfImprover.Improve(deposited, draft).Document.Results.Single(r => r["comment[data file]"] == "HumanHFpEF_1.raw");

            Assert.That(improved["characteristics[organism part]"], Is.EqualTo("heart left ventricle"));
            Assert.That(improved["comment[organism part source]"], Is.EqualTo("publication"));
            Assert.That(improved["comment[organism part source reference]"], Is.EqualTo("mmc1.xlsx!DatasetS1!R2C4"));
            Assert.That(improved["comment[organism part source method]"], Is.EqualTo("rules"));
        }

        [Test]
        public void TheEvidenceFileRoundTrips()
        {
            var evidence = new[]
            {
                Claim("HumanHFpEF_1.raw", "characteristics[age]", "77Y"),
                new SdrfEvidence("", "TMT127N", "source name", "SCC070", "supplement", "mmc6.xlsx!TMT Metadata!R2", "channel-map",
                    SdrfEvidenceConfidence.Certain),
            };
            string path = Path.Combine(TestContext.CurrentContext.WorkDirectory, "evidence-roundtrip.tsv");

            SdrfEvidenceFile.Write(path, evidence);
            var back = SdrfEvidenceFile.Read(path);

            Assert.That(back, Is.EqualTo(evidence));
            Assert.That(File.ReadAllLines(path)[0],
                Is.EqualTo("data file\tlabel\tcolumn\tvalue\tsource\tlocator\tmethod\tconfidence"));
            File.Delete(path);
        }

        [Test]
        public void AnEvidenceFileWithoutItsColumnsIsRefused()
        {
            string path = Path.Combine(TestContext.CurrentContext.WorkDirectory, "evidence-bad.tsv");
            File.WriteAllText(path, "data file\tcolumn\tvalue\nA.raw\tcharacteristics[age]\t40Y\n");

            var e = Assert.Throws<MzLibException>(() => SdrfEvidenceFile.Read(path));
            Assert.That(e!.Message, Does.Contain("confidence"));
            File.Delete(path);
        }
    }
}
