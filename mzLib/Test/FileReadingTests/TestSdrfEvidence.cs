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
        [TestCase("sdrf", "rules")]
        [TestCase("file-key", "rules")]
        [TestCase("text", "rules")]
        public void TheSourceMethodSaysWhetherRulesACuratorOrAModelReadIt(string method, string word)
        {
            var claim = new SdrfEvidence("HumanHFpEF_1.raw", "", "characteristics[age]", "77Y", "paper", "PMC1/sec:methods", method, SdrfEvidenceConfidence.Likely);
            var doc = SdrfDrafter.ToDocument(SdrfDrafter.Draft(Project(), Files, new[] { claim }), "PXD060431");

            Assert.That(Written(doc, "HumanHFpEF_1.raw")["comment[age source method]"], Is.EqualTo(word));
        }

        /// <summary>
        /// G42 (dataRepo 024/025): where a publication is also the row default, option B (D31) wrote no
        /// column source, so the builder filled it <c>not applicable</c> beside a reference and a method
        /// (PXD016662's 65 replicate rows), which reads as a contradiction. A column that carries a source
        /// reference or method now always carries its source word too. A cell with no evidence still has none.
        /// </summary>
        [Test]
        public void AColumnWithASourceReferenceAlwaysSaysItsSource()
        {
            var bare = new[] { "alpha.raw", "beta.raw" };
            var evidence = new[]
            {
                Claim("alpha.raw", "characteristics[organism part]", "heart"),
                Claim("alpha.raw", "characteristics[age]", "77Y"),
                Claim("alpha.raw", "characteristics[sex]", "female"),
                Claim("alpha.raw", "characteristics[biological replicate]", "3"),
            };
            var doc = SdrfDrafter.ToDocument(SdrfDrafter.Draft(Project(), bare, evidence), "PXD060431");
            var alpha = Written(doc, "alpha.raw");
            Assume.That(alpha["comment[characteristics source]"], Is.EqualTo("publication"), "publication is the row default here");

            foreach (var name in new[] { "organism part", "age", "sex", "biological replicate" })
            {
                Assert.That(alpha[$"comment[{name} source reference]"], Is.EqualTo("mmc1.xlsx!DatasetS1!R2C4"), name);
                Assert.That(alpha[$"comment[{name} source]"], Is.EqualTo("publication"),
                    $"{name}: a reference is never written beside 'not applicable'");
            }
            Assert.That(alpha["comment[organism source]"], Is.EqualTo("pride project record"), "an override is still an override");

            var beta = Written(doc, "beta.raw");
            Assert.That(beta["comment[age source reference]"], Is.EqualTo("not applicable"));
            Assert.That(beta["comment[age source]"], Is.EqualTo("not applicable"), "no evidence, no word of its own");
            Assert.That(SdrfValidator.Validate(doc).Errors, Is.Empty);
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

        /// <summary>A row the improver drafts for a file the deposit does not list says "publication" too.</summary>
        [Test]
        public void ADraftedRowsPublicationCellSaysPublication()
        {
            var header = new SdrfHeader(new[] { "source name", "characteristics[organism]", "characteristics[organism part]",
                "characteristics[biological replicate]", "assay name", "comment[label]", "comment[data file]", "comment[technical replicate]" });
            var deposited = new SdrfDocument(header, new[]
            {
                new SdrfRow(header, new[] { "s1", "homo sapiens", "heart", "1", "run 1", "label free sample", "HumanHFpEF_1.raw", "1" })
            });
            var draft = SdrfDrafter.Draft(Project(), Files,
                new[] { Claim("HumanControl_1.raw", "characteristics[organism part]", "heart left ventricle") });

            var doc = SdrfImprover.Improve(deposited, draft).Document;
            var added = Written(doc, "HumanControl_1.raw");

            Assert.That(added["characteristics[organism part]"], Is.EqualTo("heart left ventricle"));
            Assert.That(SourceOf(doc, added, "organism part"), Is.EqualTo("publication"));
            Assert.That(added["comment[organism part source reference]"], Is.EqualTo("mmc1.xlsx!DatasetS1!R2C4"));
            Assert.That(added["comment[organism part source method]"], Is.EqualTo("rules"));
        }

        [Test]
        public void TheEvidenceFileRoundTrips()
        {
            var evidence = new[]
            {
                Claim("HumanHFpEF_1.raw", "characteristics[age]", "77Y"),
                new SdrfEvidence("", "TMT127N", "source name", "SCC070", "supplement", "mmc6.xlsx!TMT Metadata!R2", "channel-map",
                    SdrfEvidenceConfidence.Certain),
                new SdrfEvidence("", "", "characteristics[organism part]", "flower", "paper", "figure: Fig1.jpg", "model",
                    SdrfEvidenceConfidence.Likely, "*TMT6*"),
            };
            string path = Path.Combine(TestContext.CurrentContext.WorkDirectory, "evidence-roundtrip.tsv");

            SdrfEvidenceFile.Write(path, evidence);
            var back = SdrfEvidenceFile.Read(path);

            Assert.That(back, Is.EqualTo(evidence));
            Assert.That(File.ReadAllLines(path)[0],
                Is.EqualTo("data file\tlabel\tcolumn\tvalue\tsource\tlocator\tmethod\tconfidence\tdata file pattern"));
            File.Delete(path);
        }

        [Test]
        public void AnEvidenceFileWrittenBeforeThePatternColumnStillReads()
        {
            string path = Path.Combine(TestContext.CurrentContext.WorkDirectory, "evidence-v1.tsv");
            File.WriteAllText(path, "data file\tlabel\tcolumn\tvalue\tsource\tlocator\tmethod\tconfidence\n" +
                "A.raw\t\tcharacteristics[age]\t40Y\tsupplement\tx\tfile-key\tlikely\n");

            Assert.That(SdrfEvidenceFile.Read(path).Single().DataFilePattern, Is.Empty);
            File.Delete(path);
        }

        [Test]
        public void AClaimAboutASetOfFilesFillsEachFileInItAndNoOther()
        {
            // G38 pilot: "the flower mix is the TMT6 files" is neither one file nor the whole deposit.
            string[] files = { "P020835_TMT6_F01.raw", "P020835_TMT6_F02.raw", "P017712_TMT10_F01.raw" };
            var evidence = new[]
            {
                new SdrfEvidence("", "", "characteristics[organism part]", "flower", "paper", "figure: Fig1.jpg", "model",
                    SdrfEvidenceConfidence.Likely, "*_tmt6_*"),
                new SdrfEvidence("", "", "characteristics[organism part]", "root", "paper", "x", "model",
                    SdrfEvidenceConfidence.Likely, "*_iTRAQ_*"),
            };

            var draft = SdrfDrafter.Draft(Project(), files, evidence);

            Assert.That(Row(draft, "P020835_TMT6_F01.raw").OrganismPart.Value, Is.EqualTo("flower"), "the glob ignores case");
            Assert.That(Row(draft, "P020835_TMT6_F02.raw").OrganismPart.Value, Is.EqualTo("flower"));
            Assert.That(Row(draft, "P017712_TMT10_F01.raw").OrganismPart.Source, Is.EqualTo(SdrfDraftSource.NotAvailable));
            Assert.That(draft.EvidenceNotes.Single().Why, Does.Contain("matches none"));
        }

        [TestCase("*TMT6*", "P020835_TMT6_F01.raw", true)]
        [TestCase("P0208??_*", "P020835_TMT6_F01.raw", true)]
        [TestCase("*TMT6*", "P017712_TMT10_F01.raw", false)]
        [TestCase("*.raw", "a.mzML", false)]
        [TestCase("a+b*", "a+b_1.raw", true)]
        public void AFilePatternIsAGlobOverTheWholeName(string pattern, string file, bool matches) =>
            Assert.That(SdrfEvidence.GlobMatches(pattern, file), Is.EqualTo(matches));

        [Test]
        public void AnEvidenceFileSkipsBlankLinesReadsShortRowsAndRefusesAnUnknownConfidence()
        {
            string path = Path.Combine(TestContext.CurrentContext.WorkDirectory, "evidence-edges.tsv");
            const string header = "confidence\tdata file\tlabel\tcolumn\tvalue\tsource\tlocator\tmethod\n";
            File.WriteAllText(path, header + "\nlikely\tA.raw\t\tcharacteristics[age]\n");

            var claim = SdrfEvidenceFile.Read(path).Single();
            Assert.That((claim.Confidence, claim.Column, claim.Value, claim.Method), Is.EqualTo((SdrfEvidenceConfidence.Likely, "characteristics[age]", "", "")),
                "columns are found by name, and cells past a short row's end are empty");

            File.WriteAllText(path, header + "stated\tA.raw\t\tcharacteristics[age]\t40Y\tpaper\tx\trules\n");
            var e = Assert.Throws<MzLibException>(() => SdrfEvidenceFile.Read(path));
            Assert.That(e!.Message, Does.Contain("line 2").And.Contain("stated"));

            File.WriteAllText(path, "");
            Assert.Throws<MzLibException>(() => SdrfEvidenceFile.Read(path));
            File.Delete(path);
        }

        [TestCase("0")]
        [TestCase("2")]
        [TestCase("7")]
        public void AnEvidenceFileRefusesANumericConfidence(string confidence)
        {
            string path = Path.Combine(TestContext.CurrentContext.WorkDirectory, $"evidence-numeric-{confidence}.tsv");
            File.WriteAllText(path, "data file\tlabel\tcolumn\tvalue\tsource\tlocator\tmethod\tconfidence\n"
                + $"A.raw\t\tcharacteristics[age]\t40Y\tpaper\tx\tfile-key\t{confidence}\n");

            var e = Assert.Throws<MzLibException>(() => SdrfEvidenceFile.Read(path));
            Assert.That(e!.Message, Does.Contain("line 2"));
            File.Delete(path);
        }

        [TestCase("llm")]
        [TestCase("curated")]
        [TestCase("")]
        public void AClaimReadByAnUnknownMethodIsReportedNotWrittenAsRules(string method)
        {
            var claim = Claim("HumanHFpEF_1.raw", "characteristics[age]", "77Y") with { Method = method };
            var draft = SdrfDrafter.Draft(Project(), Files, new[] { claim });

            Assert.That(Row(draft, "HumanHFpEF_1.raw").Characteristics!.ContainsKey("characteristics[age]"), Is.False);
            Assert.That(draft.EvidenceNotes.Single(n => n.Column == "characteristics[age]").Why, Does.Contain("method"));
        }

        [Test]
        public void EveryClaimTheDrafterDoesNotApplyIsReported()
        {
            var evidence = new[]
            {
                Claim("HumanHFpEF_1_typo.raw", "characteristics[age]", "50Y"),                                  // names no file
                Claim("HumanHFpEF_2.raw", "characteristics[age]", "60Y", confidence: SdrfEvidenceConfidence.Guess), // a guess
                Claim("HumanHFpEF_2.raw", "characteristics[sex]", " "),                                          // no value
                Claim("", "characteristics[strain]", "B6"),                                                      // deposit-wide ...
                Claim("HumanControl_1.raw", "characteristics[strain]", "C57"),                                   // ... overridden here
            };
            var draft = SdrfDrafter.Draft(Project(), Files, evidence);
            var notes = draft.EvidenceNotes;

            Assert.That(notes.Any(n => n.DataFile == "HumanHFpEF_1_typo.raw" && n.EvidenceValue == "50Y"), "a claim naming no raw file");
            Assert.That(notes.Any(n => n.DataFile == "HumanHFpEF_2.raw" && n.EvidenceValue == "60Y"), "a guess");
            Assert.That(notes.Any(n => n.DataFile == "HumanHFpEF_2.raw" && n.Column == "characteristics[sex]"), "a claim with no value");
            Assert.That(notes.Any(n => n.DataFile == "HumanControl_1.raw" && n.EvidenceValue == "B6"), "a deposit-wide claim a file claim overrode");
            Assert.That(Row(draft, "HumanControl_1.raw").Characteristics!["characteristics[strain]"].Value, Is.EqualTo("C57"));
            Assert.That(notes.Count(n => n.EvidenceValue == "B6"), Is.EqualTo(1), "only the overridden file gets the note");
        }

        [Test]
        public void AFilePatternWithManyWildcardsThatMatchesNothingReturnsQuickly()
        {
            // A backtracking glob took 36 s on this pattern: exponential in the wildcards (review of #1377).
            string pattern = string.Concat(Enumerable.Repeat("*0", 10)) + "*ZZZ.raw";
            string file = "20200101_QE_HF_" + new string('0', 30) + "_x.raw";
            var watch = System.Diagnostics.Stopwatch.StartNew();
            Assert.That(SdrfEvidence.GlobMatches(pattern, file), Is.False);
            Assert.That(watch.Elapsed.TotalSeconds, Is.LessThan(2));
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
