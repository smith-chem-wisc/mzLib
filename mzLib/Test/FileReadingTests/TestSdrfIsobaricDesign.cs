using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Modifications;
using Readers;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Covers <see cref="SdrfIsobaricDesign"/>, the SDRF → isobaric design projection (QuantProject M7,
    /// TMT half). The contract is the label-free one: a design MetaMorpheus accepts, or a refusal that
    /// lists every reason. The two real fixtures are corpus subsets chosen for the plex: PXD008841 names
    /// it only in a file-name token, and PXD061609's batch column is not a plex. Everything else is built
    /// in memory so the one defect under test is the only thing wrong.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrfIsobaricDesign
    {
        private const string Batch = "comment[sample preparation batch]";
        private const string PlexColumn = "comment[plex]";
        private const string Pxd008841Pattern = @"TMT_?pool(?<plex>\d+)";

        private static IsobaricMassTag Tag(IsobaricMassTagType type)
        {
            Assert.That(IsobaricMassTag.TryGetIsobaricMassTag(type, out var tag), Is.True, type.ToString());
            return tag!;
        }

        private static string Fixture(string name) =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, "FileReadingTests",
                "ExternalFileTypes", "SdrfDesign", name);

        private sealed record Row(string File, string Label, string Sample, string Condition = "A", string Biorep = "1",
            string Plex = "P1", string Fraction = "1", string Techrep = "1", string ReplicateSource = "deposited");

        /// <summary>An isobaric SDRF in memory, one row per (file, channel), with the plex in <c>comment[plex]</c>.</summary>
        private static SdrfDocument Document(params Row[] rows)
        {
            var header = new SdrfHeader(new[]
            {
                "source name", "characteristics[biological replicate]", "assay name", "comment[label]",
                "comment[fraction identifier]", "comment[technical replicate]", "comment[data file]",
                PlexColumn, SdrfIsobaricDesign.BiologicalReplicateSourceColumn, "factor value[condition]"
            });
            return new SdrfDocument(header, rows.Select(r => new SdrfRow(header, new[]
            {
                r.Sample, r.Biorep, "run " + r.File, r.Label, r.Fraction, r.Techrep, r.File, r.Plex,
                r.ReplicateSource, r.Condition
            })));
        }

        private static SdrfIsobaricDesignOptions ByColumn(IsobaricMassTagType type = IsobaricMassTagType.TMT10) =>
            new() { Tag = Tag(type), PlexColumn = PlexColumn };

        /// <summary>A valid two-channel plex: file a.raw, 126 = S1 (A), 127N = S2 (B).</summary>
        private static Row[] TwoChannels(string file = "a.raw", string fraction = "1") => new[]
        {
            new Row(file, "TMT126", "S1", "A", Fraction: fraction),
            new Row(file, "TMT127N", "S2", "B", Fraction: fraction),
        };

        // ---------------------------------------------------------------- PXD008841

        [Test]
        public void Pxd008841WithTheFileNamePatternGivesTwoPlexesOfTenChannels()
        {
            var design = SdrfIsobaricDesign.Read(Fixture("PXD008841.subset.sdrf.tsv"),
                new SdrfIsobaricDesignOptions { Tag = Tag(IsobaricMassTagType.TMT10), PlexFileNamePattern = Pxd008841Pattern });

            Assert.That(design.Refusals, Is.Empty, design.Report());
            Assert.That(design.Notes, Is.Empty);
            Assert.That(design.FileKeyColumn, Is.EqualTo("comment[data file]"));
            Assert.That(design.ConditionColumns, Is.EqualTo(new[] { "factor value[disease]" }));
            Assert.That(design.Files, Has.Count.EqualTo(4));

            // Both spellings of the token give the plex, and the ids are MetaMorpheus's: name order.
            Assert.That(design.Files.Select(f => f.Plex).Distinct(), Is.EquivalentTo(new[] { "1", "5" }));
            Assert.That(design.PlexIds, Is.EquivalentTo(new Dictionary<string, int> { ["1"] = 1, ["5"] = 2 }));

            var pool1 = design.Files.First(f => f.Plex == "1");
            Assert.That(pool1.Channels.Select(c => c.Label), Is.EqualTo(Tag(IsobaricMassTagType.TMT10).ChannelLabels),
                "channels are in reporter m/z order, so 127N comes before 127C");
            Assert.That(pool1.Channels[0].SampleName, Is.EqualTo("OSL.53E"));
            Assert.That(pool1.Channels[0].Condition, Is.EqualTo("basal-like breast carcinoma"));
            Assert.That(pool1.Channels.Last().SampleName, Is.EqualTo("pool"));

            Assert.That(design.Files.Select(f => f.Fraction), Is.EquivalentTo(new[] { 1, 2, 1, 2 }));
            Assert.That(design.Files.All(f => f.TechnicalReplicate == 1));
        }

        [Test]
        public void Pxd008841WithNoPlexSourceIsRefusedRatherThanGuessed()
        {
            var design = SdrfIsobaricDesign.Read(Fixture("PXD008841.subset.sdrf.tsv"),
                new SdrfIsobaricDesignOptions { Tag = Tag(IsobaricMassTagType.TMT10) });

            Assert.That(design.IsValid, Is.False);
            Assert.That(design.Refusals.Single(), Does.Contain("No plex source was declared"));
            Assert.That(design.Files, Is.Empty);
            Assert.Throws<InvalidOperationException>(() => design.ToExperimentalDesign());
        }

        [Test]
        public void Pxd008841ReadAsTmt11GivesAnEmpty131CInEveryFile()
        {
            var tmt11 = Tag(IsobaricMassTagType.TMT11);
            var design = SdrfIsobaricDesign.Read(Fixture("PXD008841.subset.sdrf.tsv"),
                new SdrfIsobaricDesignOptions { Tag = tmt11, PlexFileNamePattern = Pxd008841Pattern });
            Assert.That(design.Refusals, Is.Empty, design.Report());

            var experimentalDesign = design.ToExperimentalDesign();
            Assert.That(experimentalDesign.FileNameSampleInfoDictionary, Has.Count.EqualTo(4));
            foreach (var samples in experimentalDesign.FileNameSampleInfoDictionary.Values)
            {
                var channels = samples.Cast<IsobaricQuantSampleInfo>().ToList();
                Assert.That(channels.Select(c => c.ChannelLabel), Is.EqualTo(tmt11.ChannelLabels));
                Assert.That(channels.Select(c => c.ReporterIonMz), Is.EqualTo(tmt11.ReporterIonMzs));

                var empty = channels.Last();
                Assert.That(empty.ChannelLabel, Is.EqualTo("131C"));
                Assert.That(empty.SampleName, Is.Null);
                Assert.That(empty.Condition, Is.Empty);
                Assert.That(empty.BiologicalReplicate, Is.EqualTo(0));
                Assert.That(channels.Take(10).All(c => c.SampleName != null && c.BiologicalReplicate == 1));
                Assert.That(channels.All(c => !c.IsReferenceChannel), "sample type is not read until M4 (MAP-08)");
            }
        }

        // ---------------------------------------------------------------- PXD061609

        [Test]
        public void Pxd061609AsOnePlexIsRefusedBecauseEachSaxArmReusesTheFractionNumbers()
        {
            var design = SdrfIsobaricDesign.Read(Fixture("PXD061609.subset.sdrf.tsv"),
                new SdrfIsobaricDesignOptions { Tag = Tag(IsobaricMassTagType.TMT18), SinglePlex = "tube1" });

            Assert.That(design.IsValid, Is.False);
            // 18 samples x 2 fractions, each named by three files (one per arm): MetaMorpheus's
            // (sample, biorep, fraction, techrep) rule, reported in full.
            Assert.That(design.Refusals, Has.Count.EqualTo(36), design.Report());
            Assert.That(design.Refusals, Has.All.Contains("more than once"));
            Assert.That(design.Refusals.First(), Does.Contain("'ea12578.raw'").And.Contain("'ea12603.raw'").And.Contain("'ea12628.raw'"));
        }

        [Test]
        public void Pxd061609WithItsBatchColumnDeclaredAsThePlexGivesThreePlexes()
        {
            var design = SdrfIsobaricDesign.Read(Fixture("PXD061609.subset.sdrf.tsv"),
                new SdrfIsobaricDesignOptions { Tag = Tag(IsobaricMassTagType.TMT18), PlexColumn = Batch });

            Assert.That(design.Refusals, Is.Empty, design.Report());
            Assert.That(design.PlexIds.Keys, Is.EquivalentTo(new[] { "HighSalt", "LowSalt", "NoSAX" }));
            Assert.That(design.PlexSource, Does.Contain(Batch).And.Contain("declared"));

            // Biological replicates are passed through as written (N3), not ranked: brain_male is 2.
            var noSax = design.Files.First(f => f.Plex == "NoSAX");
            Assert.That(noSax.Channels, Has.Count.EqualTo(18));
            Assert.That(noSax.Channels.Single(c => c.Label == "127N").SampleName, Is.EqualTo("brain_male"));
            Assert.That(noSax.Channels.Single(c => c.Label == "127N").BiologicalReplicate, Is.EqualTo(2));
            Assert.That(noSax.Channels.Single(c => c.Label == "127N").Condition, Is.EqualTo("NT=brain;AC=UBERON:0000955"));
        }

        [Test]
        public void Pxd061609ReadAsTmt10IsRefusedForEveryChannelTmt10DoesNotHave()
        {
            var design = SdrfIsobaricDesign.Read(Fixture("PXD061609.subset.sdrf.tsv"),
                new SdrfIsobaricDesignOptions { Tag = Tag(IsobaricMassTagType.TMT10), PlexColumn = Batch });

            // 131C and 132N..135N: 8 channels x 6 files.
            Assert.That(design.Refusals, Has.Count.EqualTo(48), design.Report());
            Assert.That(design.Refusals, Has.All.Contains("not a channel of TMT10"));
        }

        // ---------------------------------------------------------------- options

        [Test]
        public void MissingTagAndPlexSourceAreBothReported()
        {
            var design = SdrfIsobaricDesign.Read(Document(TwoChannels()), new SdrfIsobaricDesignOptions());

            Assert.That(design.Refusals, Has.Count.EqualTo(2));
            Assert.That(design.Refusals[0], Does.Contain("No isobaric tag"));
            Assert.That(design.Refusals[1], Does.Contain("No plex source"));
            Assert.That(design.Report(), Does.Contain("REFUSED: 2 reason(s)").And.Contain("Plex from: (not declared)"));
        }

        [Test]
        public void TwoPlexSourcesAreRefused()
        {
            var design = SdrfIsobaricDesign.Read(Document(TwoChannels()),
                new SdrfIsobaricDesignOptions { Tag = Tag(IsobaricMassTagType.TMT10), PlexColumn = PlexColumn, SinglePlex = "P1" });

            Assert.That(design.Refusals.Single(), Does.Contain("2 plex sources were declared"));
        }

        [TestCase("comment[no such column]", null, "is not in the SDRF")]
        [TestCase(null, "pool(", "is not a valid regular expression")]
        [TestCase(null, @"pool\d+", "has no group")]
        [TestCase(null, @"(zzz)", "does not capture a plex from 'a.raw'")]
        public void AnUnusablePlexSourceIsRefused(string? column, string? pattern, string expected)
        {
            var design = SdrfIsobaricDesign.Read(Document(TwoChannels()), new SdrfIsobaricDesignOptions
            {
                Tag = Tag(IsobaricMassTagType.TMT10),
                PlexColumn = column,
                PlexFileNamePattern = pattern,
            });

            Assert.That(design.IsValid, Is.False);
            Assert.That(design.Refusals, Has.Some.Contains(expected), design.Report());
        }

        [Test]
        public void AFileNamePatternUsesItsFirstGroupWhenNoneIsNamedPlex()
        {
            var design = SdrfIsobaricDesign.Read(Document(TwoChannels("run_set7_f1.raw")),
                new SdrfIsobaricDesignOptions { Tag = Tag(IsobaricMassTagType.TMT10), PlexFileNamePattern = @"SET(\d+)" });

            Assert.That(design.Refusals, Is.Empty, design.Report());
            Assert.That(design.Files.Single().Plex, Is.EqualTo("7"), "the pattern ignores case");
        }

        [TestCase("")]
        [TestCase("not available")]
        public void APlexColumnCellThatNamesNoPlexIsRefused(string plex)
        {
            var design = SdrfIsobaricDesign.Read(Document(
                new Row("a.raw", "TMT126", "S1", Plex: plex),
                new Row("a.raw", "TMT127N", "S2", Plex: "P1")), ByColumn());

            Assert.That(design.Refusals, Has.Some.Contains($"the plex column '{PlexColumn}' is '{plex}'"), design.Report());
        }

        [Test]
        public void MissingRequiredColumnsAreAllReported()
        {
            var header = new SdrfHeader(new[] { "assay name", "comment[data file]", PlexColumn, "factor value[condition]" });
            var sdrf = new SdrfDocument(header, new[] { new SdrfRow(header, new[] { "run a", "a.raw", "P1", "A" }) });

            var design = SdrfIsobaricDesign.Read(sdrf, ByColumn());

            Assert.That(design.Refusals, Has.Count.EqualTo(3), design.Report());
            Assert.That(design.Refusals, Has.Some.Contains("'comment[label]'"));
            Assert.That(design.Refusals, Has.Some.Contains("'source name'"));
            Assert.That(design.Refusals, Has.Some.Contains("'characteristics[biological replicate]'"));
        }

        // ---------------------------------------------------------------- labels (MAP-06)

        [TestCase("TMT127N")]
        [TestCase("tmt127n")]
        [TestCase("127N")]
        [TestCase("NT=TMT127N;AC=PRIDE:0000519")]
        [TestCase("AC=PRIDE:0000519;NT=TMT127N")]
        public void EveryLabelFormReadsAsTheKitsChannel(string label)
        {
            var design = SdrfIsobaricDesign.Read(Document(
                new Row("a.raw", "TMT126", "S1"),
                new Row("a.raw", label, "S2")), ByColumn());

            Assert.That(design.Refusals, Is.Empty, design.Report());
            Assert.That(design.Files.Single().Channels.Select(c => c.Label), Is.EqualTo(new[] { "126", "127N" }));
        }

        [Test]
        public void TmtProLabelsReadForTmtProKitsOnly()
        {
            var rows = new[] { new Row("a.raw", "TMTpro126", "S1") };
            Assert.That(SdrfIsobaricDesign.Read(Document(rows), ByColumn(IsobaricMassTagType.TMT18)).Refusals, Is.Empty);
            Assert.That(SdrfIsobaricDesign.Read(Document(rows), ByColumn(IsobaricMassTagType.TMT10)).Refusals.Single(),
                Does.Contain("not a channel of TMT10"));
        }

        [TestCase("ITRAQ115", IsobaricMassTagType.diLeu4, "not a channel of diLeu4")]
        [TestCase("TMT131C", IsobaricMassTagType.TMT10, "not a channel of TMT10")]
        [TestCase("label free sample", IsobaricMassTagType.TMT10, "a label-free row in an isobaric design")]
        [TestCase("TMT10plex", IsobaricMassTagType.TMT10, "not a channel of TMT10")]
        [TestCase("", IsobaricMassTagType.TMT10, "'comment[label]' is empty")]
        public void ALabelThatIsNotAChannelOfTheKitIsRefused(string label, IsobaricMassTagType type, string expected)
        {
            // 115 is a channel of iTRAQ4 and of DiLeu4; only the family prefix tells them apart.
            var design = SdrfIsobaricDesign.Read(Document(new Row("a.raw", label, "S1")), ByColumn(type));

            Assert.That(design.Refusals.Single(), Does.Contain(expected));
        }

        // ---------------------------------------------------------------- QP-S23

        [Test]
        public void DraftedChannelsWithUnknownSamplesAreEachRefusedByName()
        {
            var design = SdrfIsobaricDesign.Read(Document(
                new Row("a.raw", "TMT126", "S1"),
                new Row("a.raw", "TMT127N", "not available"),
                new Row("a.raw", "TMT127C", "a 127C", ReplicateSource: "default"),
                new Row("a.raw", "TMT128N", "", Condition: "B")), ByColumn());

            Assert.That(design.Refusals, Has.Count.EqualTo(3), design.Report());
            Assert.That(design.Refusals[0], Does.Contain("channel 127N").And.Contain("'not available'"));
            Assert.That(design.Refusals[1], Does.Contain("channel 127C").And.Contain("'default'"));
            Assert.That(design.Refusals[2], Does.Contain("'source name' is empty"));
        }

        // ---------------------------------------------------------------- MetaMorpheus's checks

        [Test]
        public void AFileGivenTwoFractionsOrTwoPlexesIsRefused()
        {
            var design = SdrfIsobaricDesign.Read(Document(
                new Row("a.raw", "TMT126", "S1", Fraction: "1"),
                new Row("a.raw", "TMT127N", "S2", Fraction: "2"),
                new Row("b.raw", "TMT126", "S1", Plex: "P1"),
                new Row("b.raw", "TMT127N", "S2", Plex: "P2")), ByColumn());

            Assert.That(design.Refusals, Has.Some.Matches<string>(r => r.Contains("'a.raw' is given 2 different") && r.Contains("fraction 2")), design.Report());
            Assert.That(design.Refusals, Has.Some.Matches<string>(r => r.Contains("'b.raw' is given 2 different") && r.Contains("plex 'P2'")));
        }

        [Test]
        public void AChannelNamedTwiceInOneFileIsRefused()
        {
            var design = SdrfIsobaricDesign.Read(Document(
                new Row("a.raw", "TMT126", "S1"),
                new Row("a.raw", "NT=TMT126;AC=PRIDE:0000516", "S1")), ByColumn());

            Assert.That(design.Refusals.Single(), Does.Contain("names channel 126 on 2 rows (lines 2, 3)"));
        }

        [Test]
        public void APlexThatDescribesAChannelTwoWaysIsRefused()
        {
            // Fraction 2 of the same plex says 127N holds another sample.
            var design = SdrfIsobaricDesign.Read(Document(TwoChannels("f1.raw", "1")
                .Concat(new[]
                {
                    new Row("f2.raw", "TMT126", "S1", "A", Fraction: "2"),
                    new Row("f2.raw", "TMT127N", "S3", "B", Fraction: "2"),
                }).ToArray()), ByColumn());

            Assert.That(design.Refusals.Single(), Does.Contain("Plex 'P1' channel 127N is described 2 ways")
                .And.Contain("'S2'").And.Contain("'S3'"));
        }

        [Test]
        public void TheSameSampleTwiceInAPlexAtOneFractionIsRefused()
        {
            var design = SdrfIsobaricDesign.Read(Document(
                new Row("a.raw", "TMT126", "S1"),
                new Row("a.raw", "TMT127N", "S1")), ByColumn());

            Assert.That(design.Refusals.Single(), Does.Contain("sample 'S1' biorep 1 fraction 1 techrep 1 more than once"));
        }

        [Test]
        public void ABridgeSampleInTwoPlexesIsNotARepeat()
        {
            var design = SdrfIsobaricDesign.Read(Document(
                new Row("a.raw", "TMT126", "S1"), new Row("a.raw", "TMT131N", "bridge", "pool", Plex: "P1"),
                new Row("b.raw", "TMT126", "S2", Plex: "P2"), new Row("b.raw", "TMT131N", "bridge", "pool", Plex: "P2")),
                ByColumn());

            Assert.That(design.Refusals, Is.Empty, design.Report());
            Assert.That(design.PlexIds, Is.EquivalentTo(new Dictionary<string, int> { ["P1"] = 1, ["P2"] = 2 }));
        }

        [Test]
        public void AFileMissingAChannelRowTakesItsPlexsAnnotationAndSaysSo()
        {
            var design = SdrfIsobaricDesign.Read(Document(TwoChannels("f1.raw", "1")
                .Append(new Row("f2.raw", "TMT126", "S1", "A", Fraction: "2")).ToArray()), ByColumn());

            Assert.That(design.Refusals, Is.Empty, design.Report());
            Assert.That(design.Files[1].Channels.Select(c => c.SampleName), Is.EqualTo(new[] { "S1", "S2" }));
            Assert.That(design.Notes.Single(), Does.Contain("'f2.raw' has no row for channel(s) 127N"));
        }

        [Test]
        public void ConditionsFollowTheLabelFreeRules()
        {
            var design = SdrfIsobaricDesign.Read(Document(
                new Row("a.raw", "TMT126", "S1", "not applicable"),
                new Row("a.raw", "TMT127N", "S2", "")), ByColumn());

            Assert.That(design.Refusals, Has.Count.EqualTo(2), design.Report());
            Assert.That(design.Refusals[0], Does.Contain("'factor value[condition]' is 'not applicable'"));
            Assert.That(design.Refusals[1], Does.Contain("'factor value[condition]' is empty"));
        }

        [TestCase("0")]
        [TestCase("two")]
        public void ABiologicalReplicateThatIsNotPositiveIsRefused(string biorep)
        {
            var design = SdrfIsobaricDesign.Read(Document(new Row("a.raw", "TMT126", "S1", Biorep: biorep)), ByColumn());

            Assert.That(design.Refusals.Single(), Does.Contain($"'characteristics[biological replicate]' is '{biorep}'"));
        }

        // ---------------------------------------------------------------- searched files

        [Test]
        public void SearchedFilesDropOtherRowsAndGiveTheSearchedPath()
        {
            string searched = Path.Combine("data", "a.raw");
            var design = SdrfIsobaricDesign.Read(Document(TwoChannels("a.raw")
                .Concat(new[] { new Row("b.raw", "TMT126", "not available") }).ToArray()),
                new SdrfIsobaricDesignOptions { Tag = Tag(IsobaricMassTagType.TMT10), PlexColumn = PlexColumn, SearchedFiles = new[] { searched } });

            Assert.That(design.Refusals, Is.Empty, "a problem in a row nothing reads must not refuse the design");
            Assert.That(design.Notes.Single(), Does.Contain("Line 4 ('b.raw') dropped"));
            Assert.That(design.Files.Single().FilePath, Is.EqualTo(searched));
        }

        [Test]
        public void ASearchedFileWithNoRowIsRefused()
        {
            var design = SdrfIsobaricDesign.Read(Document(TwoChannels("a.raw")),
                new SdrfIsobaricDesignOptions { Tag = Tag(IsobaricMassTagType.TMT10), PlexColumn = PlexColumn, SearchedFiles = new[] { "a.raw", "a-calib.mzML" } });

            Assert.That(design.Refusals.Single(), Does.Contain("Searched file 'a-calib.mzML' has no SDRF row"));
        }

        [Test]
        public void AnEmptySearchedListIsRefused()
        {
            var design = SdrfIsobaricDesign.Read(Document(TwoChannels()),
                new SdrfIsobaricDesignOptions { Tag = Tag(IsobaricMassTagType.TMT10), PlexColumn = PlexColumn, SearchedFiles = new[] { " " } });

            Assert.That(design.Refusals.Single(), Does.Contain("names no file"));
        }

        // ---------------------------------------------------------------- output

        [Test]
        public void WriteTmtDesignWritesOneRowPerFileAndChannelInReporterOrder()
        {
            var design = SdrfIsobaricDesign.Read(Document(
                new Row("a.raw", "TMT127C", "S3", "B", "2", Fraction: "3", Techrep: "2"),
                new Row("a.raw", "TMT126", "S1", "A", "1", Fraction: "3", Techrep: "2"),
                new Row("a.raw", "TMT127N", "S2", "A", "2", Fraction: "3", Techrep: "2")), ByColumn());
            Assert.That(design.Refusals, Is.Empty, design.Report());

            string path = Path.Combine(TestContext.CurrentContext.WorkDirectory, "TestSdrfIsobaricDesign_TmtDesign.txt");
            try
            {
                design.WriteTmtDesign(path);
                Assert.That(File.ReadAllText(path), Is.EqualTo(
                    SdrfIsobaricDesign.TmtDesignHeader + "\n" +
                    "a.raw\tP1\tS1\t126\tA\t1\t3\t2\n" +
                    "a.raw\tP1\tS2\t127N\tA\t2\t3\t2\n" +
                    "a.raw\tP1\tS3\t127C\tB\t2\t3\t2\n"));
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public void ARefusedDesignWritesNothing()
        {
            var design = SdrfIsobaricDesign.Read(Document(new Row("a.raw", "TMT126", "not available")), ByColumn());
            string path = Path.Combine(TestContext.CurrentContext.WorkDirectory, "TestSdrfIsobaricDesign_Refused.txt");

            Assert.Throws<InvalidOperationException>(() => design.WriteTmtDesign(path));
            Assert.That(File.Exists(path), Is.False);
        }

        [Test]
        public void TheSampleNameReachesTheColumnLabel()
        {
            var design = SdrfIsobaricDesign.Read(Document(TwoChannels()), ByColumn());
            var samples = design.ToExperimentalDesign().FileNameSampleInfoDictionary["a.raw"].Cast<IsobaricQuantSampleInfo>().ToList();

            Assert.That(samples, Has.Count.EqualTo(10));
            Assert.That(samples[0].SampleName, Is.EqualTo("S1"));
            Assert.That(samples[0].PlexId, Is.EqualTo(1));
            Assert.That(samples[0].Fraction, Is.EqualTo(1), "1-based, passed through (N3)");
            Assert.That(samples[2].SampleName, Is.Null, "127C is not annotated");
        }

        [Test]
        public void TheReportNamesItsSources()
        {
            var report = SdrfIsobaricDesign.Read(Document(TwoChannels()), ByColumn()).Report();

            Assert.That(report, Does.Contain("1 file(s), 1 plex(es), kit TMT10"));
            Assert.That(report, Does.Contain($"Plex from: column '{PlexColumn}' (declared)"));
            Assert.That(report, Does.Contain("Condition from: factor value[condition]"));
            Assert.That(report, Does.Contain("Sample type: not read (MAP-08)"));
        }

        [Test]
        public void ALabelFreeSdrfReadAsIsobaricPointsAtTheLabelFreeReader()
        {
            var design = SdrfIsobaricDesign.Read(Fixture("PXD067622.sdrf.tsv"),
                new SdrfIsobaricDesignOptions { Tag = Tag(IsobaricMassTagType.TMT10), SinglePlex = "P", ConditionColumns = new[] { "factor value[genotype]" } });

            Assert.That(design.Refusals, Has.All.Contains(nameof(SdrfLabelFreeDesign)));
        }
    }
}
