using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MzLibUtil;
using NUnit.Framework;
using Omics.Modifications;
using Readers;
using UsefulProteomicsDatabases;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Drafting an isobaric (TMT, iTRAQ) deposit as what it is (sdrf PXR-F15 P1). Before this, every deposit was drafted as
    /// label free with one row per file, so every TMT draft was wrong in kind. A drafted isobaric SDRF states the kit only
    /// from evidence or the record's own words, one row per file and channel, and says <c>not available</c> for everything
    /// a design needs that nothing states: the plex, the fraction, the technical replicate, and each channel's sample.
    /// QuantProject's readers then refuse the design instead of guessing it (their QP-S48, 2026-11-16).
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrfDrafterIsobaric
    {
        private static CvParam Term(string label, string accession, string name) => new(label, accession, name, "");

        private static PrideProject Project(string protocol) => new()
        {
            Accession = "PXD034564",
            Title = "Reduced polymerase pausing in the aging liver",
            ProjectDescription = "Liver proteomes of young and old mice.",
            SampleProcessingProtocol = protocol,
            Organisms = { Term("NEWT", "NEWT:10090", "Mus musculus (mouse)") },
            OrganismParts = { Term("BTO", "BTO:0000759", "liver") },
            Instruments = { Term("MS", "MS:1002416", "Orbitrap Fusion") },
        };

        private const string Tmt10 = "Peptides were labelled with TMT10plex reagents and fractionated by high-pH reversed phase.";

        private static readonly string[] Files = { "Plex1_F01.raw", "Plex1_F02.raw" };

        private static SdrfEvidence Kit(string file, string value) =>
            new(file, "", "comment[label]", value, "deposit", "02b_qc/qc_report.json", "spectra", SdrfEvidenceConfidence.Certain);

        private static IEnumerable<SdrfDraftRow> Of(SdrfDraft d, string file) => d.Rows.Where(r => r.DataFile == file);

        private static IEnumerable<SdrfRow> Written(SdrfDocument doc, string file) => doc.Results.Where(r => r["comment[data file]"] == file);

        [Test]
        public void ARecordThatNamesItsKitDraftsOneRowPerFileAndChannel()
        {
            var d = SdrfDrafter.Draft(Project(Tmt10), Files);

            var expected = new[] { "TMT126", "TMT127N", "TMT127C", "TMT128N", "TMT128C", "TMT129N", "TMT129C", "TMT130N", "TMT130C", "TMT131" };
            foreach (var file in Files)
                Assert.That(Of(d, file).Select(r => r.Label!.Value), Is.EqualTo(expected), file);
            Assert.That(d.Rows.All(r => r.Label!.Source == SdrfDraftSource.PrideProjectRecord));
            Assert.That(d.Rows[0].Label!.Evidence, Does.Contain("TMT10plex"), "the label cites the record's words");
        }

        /// <summary>The kit's channels come from mzLib's own table (IsobaricMassTag, SdrfIsobaricLabels), so a draft writes
        /// what MetaMorpheus writes: TMT10's tenth channel is TMT131, TMT11's are TMT131N and TMT131C, iTRAQ 8 has no 120.</summary>
        [TestCase("TMT11plex", 11, "TMT131N", "TMT131C")]
        [TestCase("TMTpro 18plex", 18, "TMT134C", "TMT135N")]
        [TestCase("TMTpro 16plex", 16, "TMT133C", "TMT134N")]
        [TestCase("iTRAQ 8-plex", 8, "ITRAQ119", "ITRAQ121")]
        [TestCase("iTRAQ 4-plex", 4, "ITRAQ116", "ITRAQ117")]
        public void EachKitDraftsItsOwnChannels(string kit, int channels, string secondLast, string last)
        {
            var d = SdrfDrafter.Draft(Project($"Labelled with {kit} reagents."), new[] { "S1.raw" });

            var labels = d.Rows.Select(r => r.Label!.Value).ToList();
            Assert.That(labels, Has.Count.EqualTo(channels));
            Assert.That(labels.Skip(channels - 2), Is.EqualTo(new[] { secondLast, last }));
            Assert.That(labels, Has.None.EqualTo("ITRAQ120"));
        }

        [Test]
        public void TheDocumentWritesBareLabelsAndOneAssayNamePerRun()
        {
            var doc = SdrfDrafter.ToDocument(SdrfDrafter.Draft(Project(Tmt10), Files), "PXD034564");

            var run = Written(doc, "Plex1_F01.raw").ToList();
            Assert.That(run, Has.Count.EqualTo(10));
            Assert.That(run[0]["comment[label]"], Is.EqualTo("TMT126"), "D25: a label is written bare");
            Assert.That(run.Select(r => r["assay name"]).Distinct(), Is.EqualTo(new[] { "run Plex1_F01" }),
                "the channels of one run share its assay name");
            var v = SdrfValidator.Validate(doc);
            Assert.That(v.Errors, Is.Empty, string.Join("\n", v.Errors.Select(e => e.ToString())));
        }

        /// <summary>PXR-F15: identity from evidence, never from numbering. <c>F01</c> in a name is a number, not a statement
        /// that the file is a fraction of the plex another file belongs to.</summary>
        [Test]
        public void APlexFractionAndTechnicalReplicateNothingStatesAreNotAvailable()
        {
            var d = SdrfDrafter.Draft(Project(Tmt10), Files);

            Assert.That(d.Rows.All(r => r.Plex!.Source == SdrfDraftSource.NotAvailable));
            Assert.That(d.Rows.All(r => r.Fraction.Source == SdrfDraftSource.NotAvailable), "read from F01/F02 in the names");
            Assert.That(d.Rows.All(r => r.TechnicalReplicate.Source == SdrfDraftSource.NotAvailable));

            var doc = SdrfDrafter.ToDocument(d);
            foreach (var row in doc.Results)
            {
                Assert.That(row["comment[plex]"], Is.EqualTo("not available"), "QP-S48: exactly `not available`, which the reader refuses");
                Assert.That(row["comment[fraction identifier]"], Is.EqualTo("not available"));
                Assert.That(row["comment[technical replicate]"], Is.EqualTo("not available"));
                if (doc.Header.Contains("comment[fraction identifier source]"))
                    Assert.That(row["comment[fraction identifier source]"], Is.Not.EqualTo("inferred"), "a value nobody states has no source");
            }
        }

        /// <summary>QP-S48 (b): QuantProject's reader reads any source name as a real sample, and refuses a channel only by
        /// QP-S23's marker. So a channel nothing names gets a name scoped to its run and keeps <c>default</c>.</summary>
        [Test]
        public void AChannelNothingNamesIsNamedForItsRunAndKeepsTheDefaultMarker()
        {
            var d = SdrfDrafter.Draft(Project(Tmt10), Files);

            var row = Of(d, "Plex1_F01.raw").First();
            Assert.That(row.SourceName.Value, Is.EqualTo("Plex1_F01 TMT126"));
            Assert.That(row.SourceName.Source, Is.EqualTo(SdrfDraftSource.Default));
            Assert.That((row.BiologicalReplicate.Value, row.BiologicalReplicate.Source), Is.EqualTo(("1", SdrfDraftSource.Default)));
            Assert.That(d.Rows.Select(r => r.SourceName.Value).Distinct().Count(), Is.EqualTo(20));

            var doc = SdrfDrafter.ToDocument(d);
            Assert.That(doc.Results.All(r => r["comment[biological replicate source]"] == "default"));
        }

        /// <summary>QP-S48 (b): a made-up name another row already uses would make two channels one sample to every reader.
        /// It is numbered, as MetaMorpheus's writer numbers one (IsobaricSourceNames).</summary>
        [Test]
        public void AMadeUpNameAlreadyTakenIsNumbered()
        {
            var d = SdrfDrafter.Draft(Project(Tmt10), new[] { "S1.raw", "S1.mzML" });

            var first = d.Rows.Where(r => r.Label!.Value == "TMT126").Select(r => r.SourceName.Value).ToList();
            Assert.That(first, Is.EquivalentTo(new[] { "S1 TMT126", "S1 TMT126 (2)" }));
            Assert.That(d.Rows.Select(r => r.SourceName.Value).Distinct().Count(), Is.EqualTo(d.Rows.Count));
        }

        /// <summary>A multiplexed run's name describes the run (its plex, its fraction), never one channel's sample: a
        /// condition read from <c>WT_1</c> would put every channel of the run in WT. The record's deposit-wide facts stay.</summary>
        [Test]
        public void CellsReadFromAFileNameAreNotWrittenOnItsChannels()
        {
            var names = new[] { "WT_1.raw", "WT_2.raw", "KO_1.raw", "KO_2.raw" };
            var labelFree = SdrfDrafter.Draft(Project("Label-free quantification."), names);
            Assume.That(labelFree.Rows.Any(r => r.Factors.Any(f => f.Source == SdrfDraftSource.Inferred)),
                "the fixture must carry a condition read from the names");

            var d = SdrfDrafter.Draft(Project(Tmt10), names);

            Assert.That(d.Rows.SelectMany(r => r.Factors).All(f => f.Source == SdrfDraftSource.NotAvailable));
            Assert.That(d.Rows.All(r => r.Organism.Source == SdrfDraftSource.PrideProjectRecord));
            Assert.That(d.Rows.All(r => r.OrganismPart.Value == "liver"));
        }

        [Test]
        public void ARecordWhoseOnlyKitWordIsTmt6plexStatesNoKit()
        {
            var d = SdrfDrafter.Draft(Project("Peptides were labelled with TMT6plex reagents."), Files);

            Assert.That(d.Rows, Has.Count.EqualTo(2), "one row per file: no kit, so no channels");
            Assert.That(d.Rows.All(r => r.Label!.Source == SdrfDraftSource.NotAvailable));
            Assert.That(d.Rows[0].Label!.Evidence, Does.Contain("TMT6plex"));
            var doc = SdrfDrafter.ToDocument(d);
            Assert.That(doc.Results.Select(r => r["comment[label]"]), Is.All.EqualTo("not available"),
                "not `label free`: the record says the runs are labelled");
        }

        [Test]
        public void ADepositThatNamesNoIsobaricLabelDraftsExactlyAsBefore()
        {
            var d = SdrfDrafter.Draft(Project("Label-free quantification with MaxLFQ."), Files);

            Assert.That(d.Rows, Has.Count.EqualTo(2));
            Assert.That(d.Rows.All(r => r.Label == null && r.Plex == null));
            Assert.That(d.Rows.Select(r => r.Fraction.Value), Is.EqualTo(new[] { "1", "2" }), "the names' fractions, as before");
            var doc = SdrfDrafter.ToDocument(d);
            Assert.That(doc.Results.Select(r => r["comment[label]"]), Is.All.EqualTo("NT=label free sample;AC=MS:1002038"));
            Assert.That(doc.Header, Has.None.EqualTo("comment[plex]"));
        }

        /// <summary>A kit read from the spectra (PXReprise's QC) outranks the record's words, as a measured instrument
        /// outranks text (PXR-F10): PXD025618's record says TMT10 and its spectra carry TMT11's 131C.</summary>
        [Test]
        public void AKitReadFromTheSpectraOutranksTheRecord()
        {
            var d = SdrfDrafter.Draft(Project(Tmt10), Files, new[] { Kit("", "TMT11plex") });

            Assert.That(Of(d, "Plex1_F01.raw").Select(r => r.Label!.Value).Last(), Is.EqualTo("TMT131C"));
            Assert.That(d.Rows.All(r => r.Label!.Source == SdrfDraftSource.Deposit));
            var doc = SdrfDrafter.ToDocument(d);
            var row = doc.Results.First();
            Assert.That(row["comment[label source]"], Is.EqualTo("deposited"));
            Assert.That(row["comment[label source reference]"], Is.EqualTo("02b_qc/qc_report.json"));
            Assert.That(row["comment[label source method]"], Is.EqualTo("rules"));
            Assert.That(d.EvidenceNotes.Select(n => n.Why), Has.None.Contain("only characteristics"),
                "a kit claim is the label's evidence, not a column evidence cannot fill");
        }

        [Test]
        public void KitClaimsThatDisagreeStateNoKit()
        {
            var d = SdrfDrafter.Draft(Project(Tmt10), Files, new[] { Kit("", "TMT11plex"), Kit("", "TMTpro 16plex") });

            Assert.That(d.Rows, Has.Count.EqualTo(2));
            Assert.That(d.Rows.All(r => r.Label!.Source == SdrfDraftSource.NotAvailable));
            Assert.That(d.EvidenceNotes.Where(n => n.Column == "comment[label]").Select(n => n.Why), Has.Some.Contain("disagree"));
        }

        /// <summary>One deposit can hold labelled and label-free runs (PXD014565 has label-free, SILAC and TMT runs). A kit
        /// claimed for one file makes only that file channels.</summary>
        [Test]
        public void AKitClaimedForOneFileMakesOnlyThatFileChannels()
        {
            var d = SdrfDrafter.Draft(Project("Label-free and TMT quantification."), new[] { "lf.raw", "tmt.raw" },
                new[] { Kit("tmt.raw", "TMT10plex") });

            Assert.That(Of(d, "tmt.raw").Count(), Is.EqualTo(10));
            Assert.That(Of(d, "lf.raw").Single().Label!.Source, Is.EqualTo(SdrfDraftSource.NotAvailable),
                "the record says TMT for the deposit, so a file no claim describes is labelled, with no kit known");
        }

        [Test]
        public void AKitClaimedForOneFileOfALabelFreeRecordLeavesTheOthersLabelFree()
        {
            var d = SdrfDrafter.Draft(Project("Label-free quantification."), new[] { "lf.raw", "tmt.raw" },
                new[] { Kit("tmt.raw", "TMT10plex") });

            Assert.That(Of(d, "tmt.raw").Count(), Is.EqualTo(10));
            Assert.That(Of(d, "lf.raw").Single().Label, Is.Null);
            var doc = SdrfDrafter.ToDocument(d);
            Assert.That(Written(doc, "lf.raw").Single()["comment[label]"], Is.EqualTo("label free sample"),
                "bare in a document that holds channels");
        }

        /// <summary>A claim about a whole multiplexed run describes the run: its characteristics (organism part) hold for every
        /// channel, but a condition or a sample name would put every channel in one condition or make them one sample.</summary>
        [Test]
        public void AClaimAboutAWholeRunNamesNoChannelsConditionOrSample()
        {
            var claims = new[]
            {
                new SdrfEvidence("Plex1_F01.raw", "", "characteristics[strain]", "C57BL/6J", "publication", "Methods", "text", SdrfEvidenceConfidence.Certain),
                new SdrfEvidence("Plex1_F01.raw", "", "factor value[age]", "old", "deposit", "design.xlsx!A2", "design-table", SdrfEvidenceConfidence.Certain),
                new SdrfEvidence("Plex1_F01.raw", "", "source name", "mouse 7", "publication", "Table S1", "file-key", SdrfEvidenceConfidence.Certain),
                new SdrfEvidence("Plex1_F01.raw", "", "characteristics[biological replicate]", "3", "publication", "Table S1", "file-key", SdrfEvidenceConfidence.Certain),
            };

            var d = SdrfDrafter.Draft(Project(Tmt10), Files, claims);

            var run = Of(d, "Plex1_F01.raw").ToList();
            Assert.That(run.All(r => r.Characteristics["characteristics[strain]"].Value == "C57BL/6J"));
            Assert.That(run.SelectMany(r => r.Factors).All(f => f.Source == SdrfDraftSource.NotAvailable));
            Assert.That(run.Select(r => r.SourceName.Value), Has.None.EqualTo("mouse 7"));
            Assert.That(run.All(r => r.BiologicalReplicate.Source == SdrfDraftSource.Default), "one replicate number for a whole run would number every channel alike");
            Assert.That(d.EvidenceNotes.Where(n => n.DataFile == "Plex1_F01.raw").Select(n => n.Why),
                Has.Some.Contain("channel"));
        }

        /// <summary>The point of the whole change: a drafted TMT SDRF goes into QuantProject's reader and is refused, naming
        /// what is missing, instead of being searched with an invented design (PXR-F15; their QP-S48).</summary>
        [Test]
        public void ADraftedIsobaricSdrfIsRefusedByTheDesignReaderNotGuessed()
        {
            Assume.That(IsobaricMassTag.TryGetIsobaricMassTag(IsobaricMassTagType.TMT10, out var tag), "TMT10 must be known");
            var doc = SdrfDrafter.ToDocument(SdrfDrafter.Draft(Project(Tmt10), Files));

            var design = SdrfIsobaricDesign.Read(doc, new SdrfIsobaricDesignOptions { Tag = tag });

            Assert.That(design.IsValid, Is.False);
        }

        [Test]
        public void AnImprovedCuratedSdrfStillJoinsAnIsobaricDraft()
        {
            var draft = SdrfDrafter.Draft(Project(Tmt10), Files);
            var curated = SdrfDrafter.ToDocument(SdrfDrafter.Draft(Project("Label-free quantification."), Files));

            Assert.DoesNotThrow(() => SdrfImprover.Improve(curated, draft));
        }
    }
}
