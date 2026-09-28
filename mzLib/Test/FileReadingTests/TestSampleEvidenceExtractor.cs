using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Rule-based evidence from supplement tables (sdrf design SAMPLE-EVIDENCE.md, E3). Each case is shaped like a real
    /// deposit the rules were measured on (results/benchmark, 2026-09-26).
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSampleEvidenceExtractor
    {
        private static SupplementTable Table(string file, string sheet, string[] header, params string[][] rows) =>
            new(file, sheet, "", header, rows.Select(r => (IReadOnlyList<string>)r).ToList(),
                Enumerable.Range(2, rows.Length).ToList(), false);

        private static SdrfEvidence One(IEnumerable<SdrfEvidence> e, string file, string column, string label = "") =>
            e.Single(x => x.DataFile == file && x.Column == column && x.Label == label);

        [Test]
        public void ATableKeyedByFileNameGivesPerFileCharacteristicsDecodedFromTheHeader()
        {
            // PXD060431: the ID is the file stem; sex is coded in the header; age has no unit.
            var t = Table("JAH3-s001.xlsx", "Data", new[] { "ID", "Gender (1F,0M)", "Age", "Well_Location", "Tissue enriched" },
                new[] { "HumanHFpEF_1", "1", "77", "P1:D3", "Fat" },
                new[] { "HumanHFpEF_2", "0", "73", "P1:D4", "Fat" },
                new[] { "HumanControl_1", "0", "72", "P1:D5", "Fat" });
            var files = new[] { "HumanHFpEF_1.raw", "HumanHFpEF_2.raw", "HumanControl_1.raw", "HumanHFrEF_1.raw" };

            var e = SampleEvidenceExtractor.Extract(new[] { t }, files);

            Assert.That(One(e, "HumanHFpEF_1.raw", "characteristics[sex]").Value, Is.EqualTo("female"));
            Assert.That(One(e, "HumanHFpEF_2.raw", "characteristics[sex]").Value, Is.EqualTo("male"));
            Assert.That(One(e, "HumanHFpEF_1.raw", "characteristics[age]").Value, Is.EqualTo("77Y"));
            var c = One(e, "HumanHFpEF_1.raw", "characteristics[age]");
            Assert.That((c.Method, c.Source, c.Confidence, c.Locator), Is.EqualTo(("file-key", "supplement", SdrfEvidenceConfidence.Likely, "JAH3-s001.xlsx!Data!R2C3")));
            Assert.That(e.Select(x => x.Column), Has.No.Member("characteristics[organism part]"),
                "a well location and a protein-atlas 'tissue enriched' are not the sample's organism part");
            Assert.That(e.Where(x => x.DataFile == "HumanHFrEF_1.raw"), Is.Empty, "a file the table does not name gets nothing");
        }

        [Test]
        public void AChannelMapInOneCellGivesEveryFractionOfThePlexItsPatient()
        {
            // PXD010429: plex and channel share a cell; the plex token is in the file names.
            var t = Table("mmc6.xlsx", "TMT Metadata", new[] { "sample_name", "TMT.Sample.Name", "Proteomic Subtype" },
                new[] { "SCC070", "TMT05_TMT-129", "Redox B" },
                new[] { "SCC031", "TMT05_TMT-127", "Inflamed A" },
                new[] { "SCC080", "TMT08_TMT-131", "Redox B" },
                new[] { "SCC012", "TMT14_TMT-129", "Inflamed A" });
            var files = new[] { "Haura-SCC-PG_TMT05-2-Fx05-run1.raw", "Haura-SCC-PG_TMT05-2-Fx11-run1.raw", "Haura-SCC-PG_TMT08-2-Fx01-run1.raw" };

            var e = SampleEvidenceExtractor.Extract(new[] { t }, files);

            Assert.That(One(e, "Haura-SCC-PG_TMT05-2-Fx05-run1.raw", "source name", "TMT129").Value, Is.EqualTo("SCC070"));
            Assert.That(One(e, "Haura-SCC-PG_TMT05-2-Fx11-run1.raw", "source name", "TMT127").Value, Is.EqualTo("SCC031"));
            Assert.That(One(e, "Haura-SCC-PG_TMT08-2-Fx01-run1.raw", "source name", "TMT131").Value, Is.EqualTo("SCC080"));
            Assert.That(One(e, "Haura-SCC-PG_TMT05-2-Fx05-run1.raw", "characteristics[phenotype]", "TMT129").Value, Is.EqualTo("Redox B"),
                "a subtype is a phenotype, not a disease");
            Assert.That(e.Where(x => x.Value == "SCC012"), Is.Empty, "a plex with no file here is skipped");
            Assert.That(e.All(x => x.Method == "channel-map"));
        }

        [Test]
        public void AChannelMapWithItsOwnPlexColumn()
        {
            var t = Table("s3.xlsx", "Labels", new[] { "Plex", "Channel", "Sample" },
                new[] { "Plex1", "126", "Ctrl-1" }, new[] { "Plex1", "127N", "Ctrl-2" },
                new[] { "Plex1", "127C", "Drug-1" }, new[] { "Plex1", "128N", "Drug-2" });
            var e = SampleEvidenceExtractor.Extract(new[] { t }, new[] { "exp_Plex1_F01.raw", "exp_Plex1_F02.raw" });

            Assert.That(One(e, "exp_Plex1_F02.raw", "source name", "TMT127C").Value, Is.EqualTo("Drug-1"));
            Assert.That(e.Count(x => x.Column == "source name"), Is.EqualTo(8), "4 channels x 2 fractions");
        }

        [Test]
        public void FigureSourceDataIsNotAChannelMap()
        {
            // PXD020586: a figure's data sheet with a channel column and a numeric "sample".
            var t = Table("MOESM7.xlsx", "Supp. Fig. 6", new[] { "Plex", "Channel", "Sample", "time" },
                new[] { "P5", "127N", "1", "0.08" }, new[] { "P5", "128N", "1", "0.91" },
                new[] { "P5", "128C", "1", "1.78" }, new[] { "P5", "129N", "1", "2.11" });

            Assert.That(SampleEvidenceExtractor.Extract(new[] { t }, new[] { "run_P5_S1.raw" }), Is.Empty);
        }

        [Test]
        public void AnIsaTabStudyJoinsItsAssayToTheRawFiles()
        {
            var study = Table("isa1__s_study.txt", "", new[] { "Source Name", "Characteristics[organism]", "Characteristics[organism part]", "Sample Name" },
                new[] { "mouse 1", "Mus musculus", "heart", "M1_heart" },
                new[] { "mouse 2", "Mus musculus", "liver", "M2_liver" });
            var assay = Table("isa1__a_assay.txt", "", new[] { "Sample Name", "Assay Name", "Raw Data File" },
                new[] { "M1_heart", "a1", "m1_heart.raw" },
                new[] { "M2_liver", "a2", "m2_liver.raw" });

            var e = SampleEvidenceExtractor.Extract(new[] { study, assay }, new[] { "m1_heart.raw", "m2_liver.raw" });

            Assert.That(One(e, "m2_liver.raw", "characteristics[organism part]").Value, Is.EqualTo("liver"));
            Assert.That(One(e, "m1_heart.raw", "source name").Value, Is.EqualTo("mouse 1"));
            Assert.That(e.All(x => x.Method == "isa-tab" && x.Confidence == SdrfEvidenceConfidence.Certain));
        }

        [Test]
        public void AnSdrfShippedAsASupplementIsReadForItsOwnFilesOnly()
        {
            // PXD046357 ships an SDRF covering sibling deposits too: only rows naming this deposit's files count.
            var t = Table("MOESM8.xlsx", "Experiment overview",
                new[] { "source name", "characteristics[organism]", "characteristics[cell type]", "comment[data file]" },
                new[] { "hela_1", "homo sapiens", "HeLa", "HeLa_1cell_01.raw" },
                new[] { "hek_1", "homo sapiens", "HEK293", "HeK293T_01.raw" });

            var e = SampleEvidenceExtractor.Extract(new[] { t }, new[] { "HeLa_1cell_01.raw" });

            Assert.That(One(e, "HeLa_1cell_01.raw", "characteristics[cell type]").Value, Is.EqualTo("HeLa"));
            Assert.That(e.Where(x => x.Value == "HEK293"), Is.Empty);
            Assert.That(e.All(x => x.Method == "sdrf"));
        }

        [Test]
        public void AnIsaTabMassSpectrometryAssayJoinsThroughItsRawSpectralDataFile()
        {
            // The ISA-Tab spec names an MS assay's file column "Raw Spectral Data File", not "Raw Data File".
            var study = Table("isa1__s_study.txt", "", new[] { "Source Name", "Characteristics[organism part]", "Sample Name" },
                new[] { "mouse 1", "heart", "M1_heart" });
            var assay = Table("isa1__a_ms.txt", "", new[] { "Sample Name", "Assay Name", "Raw Spectral Data File" },
                new[] { "M1_heart", "a1", "m1_heart.raw" });

            var e = SampleEvidenceExtractor.Extract(new[] { study, assay }, new[] { "m1_heart.raw" });

            Assert.That(One(e, "m1_heart.raw", "characteristics[organism part]").Value, Is.EqualTo("heart"));
        }

        [Test]
        public void AMultiplexedSdrfSupplementGivesEachChannelItsOwnClaims()
        {
            // One raw file, three TMT channels, three patients: no channel's sex or age is the file's.
            var t = Table("sdrf.tsv", "", new[] { "source name", "characteristics[sex]", "characteristics[age]", "comment[label]", "comment[data file]" },
                new[] { "p1", "male", "70Y", "TMT126", "plex1_F1.raw" },
                new[] { "p2", "female", "40Y", "TMT127N", "plex1_F1.raw" },
                new[] { "p3", "female", "55Y", "TMT127C", "plex1_F1.raw" });

            var e = SampleEvidenceExtractor.Extract(new[] { t }, new[] { "plex1_F1.raw" });

            Assert.That(e.Where(x => x.Label == ""), Is.Empty, "every row is one channel, so no claim is about the whole file");
            Assert.That(One(e, "plex1_F1.raw", "characteristics[sex]", "TMT127N").Value, Is.EqualTo("female"));
            Assert.That(One(e, "plex1_F1.raw", "characteristics[age]", "TMT126").Value, Is.EqualTo("70Y"));
            Assert.That(e.Select(x => x.Label).Distinct(), Is.EquivalentTo(new[] { "TMT126", "TMT127N", "TMT127C" }));
            Assert.That(e.Select(x => x.Column), Has.No.Member("comment[label]"), "the label is the claim's channel, not a claim");
        }

        [Test]
        public void AMultiplexedIsaTabAssayGivesEachChannelItsOwnClaims()
        {
            var study = Table("isa1__s_study.txt", "", new[] { "Source Name", "Characteristics[sex]", "Sample Name" },
                new[] { "p1", "male", "S1" }, new[] { "p2", "female", "S2" });
            var assay = Table("isa1__a_ms.txt", "", new[] { "Sample Name", "Label", "Raw Spectral Data File" },
                new[] { "S1", "TMT 126", "plex1.raw" }, new[] { "S2", "TMT 127N", "plex1.raw" });

            var e = SampleEvidenceExtractor.Extract(new[] { study, assay }, new[] { "plex1.raw" });

            Assert.That(e.Where(x => x.Label == ""), Is.Empty);
            Assert.That(One(e, "plex1.raw", "characteristics[sex]", "TMT127N").Value, Is.EqualTo("female"));
        }

        [Test]
        public void TwoTablesThatDisagreeAboutAFileBothReachTheDrafter()
        {
            var files = new[] { "A1.raw", "A2.raw" };
            var a = Table("a.xlsx", "S1", new[] { "File", "Sex" }, new[] { "A1", "male" }, new[] { "A2", "male" });
            var b = Table("b.xlsx", "S1", new[] { "File", "Sex" }, new[] { "A1", "female" }, new[] { "A2", "male" });

            var e = SampleEvidenceExtractor.Extract(new[] { a, b }, files);

            Assert.That(e.Where(x => x.DataFile == "A1.raw" && x.Column == "characteristics[sex]").Select(x => x.Value),
                Is.EquivalentTo(new[] { "male", "female" }), "the disagreement is kept so the drafter can report it");
            Assert.That(e.Count(x => x.DataFile == "A2.raw" && x.Column == "characteristics[sex]"), Is.EqualTo(1), "an agreeing repeat collapses");
            var draft = SdrfDrafter.Draft(new UsefulProteomicsDatabases.PrideProject { Accession = "PXD000001" }, files, e);
            Assert.That(draft.Rows.Single(r => r.DataFile == "A1.raw").Characteristics!.ContainsKey("characteristics[sex]"), Is.False);
            Assert.That(draft.EvidenceNotes.Any(n => n.DataFile == "A1.raw" && n.Why.Contains("disagree")));
        }

        [Test]
        public void AStrongerRulesClaimStillHidesAWeakerRulesRepeat()
        {
            // The SDRF supplement (Certain) and a file-keyed table (Likely) both speak about A1's sex: the SDRF wins.
            var sdrf = Table("sdrf.tsv", "", new[] { "source name", "characteristics[sex]", "comment[data file]" },
                new[] { "p1", "male", "A1.raw" }, new[] { "p2", "male", "A2.raw" });
            var key = Table("a.xlsx", "S1", new[] { "File", "Sex" }, new[] { "A1", "female" }, new[] { "A2", "male" });

            var e = SampleEvidenceExtractor.Extract(new[] { sdrf, key }, new[] { "A1.raw", "A2.raw" });

            Assert.That(One(e, "A1.raw", "characteristics[sex]").Method, Is.EqualTo("sdrf"));
        }

        [TestCase("Organism", "Not correct", null)]
        [TestCase("Organism", "Escherichia coli", "Escherichia coli")]
        [TestCase("Biological replicate", "lung_1w_F", null)]
        [TestCase("Biological replicate", "2", "2")]
        [TestCase("Age (months)", "9", "9M")]
        [TestCase("Age", "0.15 years", "0.15Y")]
        [TestCase("Sex", "M", "male")]
        [TestCase("Sex", "subject", null)]
        [TestCase("Sex (0=male, 1=female)", "1", "female")]
        [TestCase("Sex (0=male, 1=female)", "0", "male")]
        [TestCase("Gender (1F,0M)", "0", "male")]
        [TestCase("Sex (1=F, 1=M)", "1", null)]
        [TestCase("Sex", "Women", "female")]
        [TestCase("Sex", "n/a", null)]
        [TestCase("Age (weeks)", "12", "12W")]
        [TestCase("Age (days)", "3", "3D")]
        [TestCase("Age", "2.50", "2.5Y")]
        [TestCase("Age", "0", null)]
        [TestCase("Age", "adult", null)]
        [TestCase("Species", "E. coli", "E. coli")]
        [TestCase("Organism", "human", null)]
        [TestCase("Biological replicate", "0", null)]
        public void AValueIsDecodedByItsHeaderOrDropped(string header, string value, string? expected)
        {
            var t = Table("t.xlsx", "S1", new[] { "File", header }, new[] { "A1", value }, new[] { "A2", value });
            var e = SampleEvidenceExtractor.Extract(new[] { t }, new[] { "A1.raw", "A2.raw" });

            if (expected == null) Assert.That(e, Is.Empty);
            else Assert.That(e.Where(x => x.DataFile == "A1.raw").Select(x => x.Value).Single(), Is.EqualTo(expected));
        }

        [Test]
        public void AResultTableIsNeverReadForSamples()
        {
            var t = Table("mmc2.xlsx", "Proteins", new[] { "Raw file", "Protein IDs", "Gene names", "Intensity", "Sex" },
                new[] { "A1", "P12345", "ALB", "1e9", "female" }, new[] { "A2", "Q99999", "APOA1", "2e9", "male" });

            Assert.That(SampleEvidenceExtractor.Extract(new[] { t }, new[] { "A1.raw", "A2.raw" }), Is.Empty);
        }
    }
}
