using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Omics.Modifications;
using Readers;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Covers <see cref="SdrfSearchDesign"/>: the one call a pipeline makes before a search, and MetaMorpheus makes when
    /// it writes a design, to learn whether an SDRF gives a design for exactly the searched files. It adds no rule of
    /// its own to the readers, so most tests compare it with them; the rest cover what it does that they do not: the
    /// join to the searched files, and refusing what would otherwise throw.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrfSearchDesign
    {
        private const string Genotype = "factor value[genotype]";
        private const string Treatment = "factor value[treatment]";
        private const string Pxd008841Pattern = @"TMT_?pool(?<plex>\d+)";

        private string _directory = null!;

        [SetUp]
        public void SetUp()
        {
            _directory = Path.Combine(Path.GetTempPath(), "SdrfSearchDesign_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(_directory);
        }

        [TearDown]
        public void TearDown()
        {
            if (Directory.Exists(_directory))
                Directory.Delete(_directory, true);
        }

        private static string Fixture(string name) =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, "FileReadingTests",
                "ExternalFileTypes", "SdrfDesign", name);

        private static List<string> DataFiles(string fixture) =>
            new SdrfDocument(Fixture(fixture)).Results.Select(r => r["comment[data file]"]!).Distinct().ToList();

        /// <summary>What a calibrated search reads for an acquired file.</summary>
        private static string Calibrated(string rawFile) => Path.GetFileNameWithoutExtension(rawFile) + "-calib.mzML";

        private static IsobaricMassTag Tmt10()
        {
            Assert.That(IsobaricMassTag.TryGetIsobaricMassTag(IsobaricMassTagType.TMT10, out var tag), Is.True);
            return tag!;
        }

        // ---------------------------------------------------------------- the readers' own answer

        /// <summary>
        /// A calibrated search of every file of PXD067622: the rows name X.raw and the search reads X-calib.mzML. The
        /// design is the label-free reader's own for the restricted rows, file for file, named as the search reads them.
        /// </summary>
        [Test]
        public void ACalibratedSearchGetsTheLabelFreeReadersOwnDesign()
        {
            var searched = DataFiles("PXD067622.sdrf.tsv").Select(Calibrated).ToList();

            var design = SdrfSearchDesign.Read(Fixture("PXD067622.sdrf.tsv"), searched,
                new SdrfSearchDesignOptions { ConditionColumns = new[] { Genotype, Treatment } });

            Assert.That(design.Refusals, Is.Empty, design.Report());
            Assert.That(design.IsValid, Is.True);
            Assert.That(design.Isobaric, Is.Null);
            Assert.That(design.LabelFree, Is.Not.Null);
            Assert.That(design.Notes, Is.Empty, "every row is searched and the SDRF numbers its replicates");

            var direct = SdrfLabelFreeDesign.Read(
                SdrfSearchScope.Restrict(new SdrfDocument(Fixture("PXD067622.sdrf.tsv")), searched).Document,
                new SdrfLabelFreeDesignOptions { ConditionColumns = new[] { Genotype, Treatment }, SearchedFiles = searched });
            Assert.That(design.LabelFree!.Files.Select(f => (Path.GetFileName(f.FullFilePathWithExtension), f.Condition,
                    f.BiologicalReplicate, f.Fraction, f.TechnicalReplicate)),
                Is.EqualTo(direct.Files.Select(f => (Path.GetFileName(f.FullFilePathWithExtension), f.Condition,
                    f.BiologicalReplicate, f.Fraction, f.TechnicalReplicate))));
            Assert.That(design.LabelFree.Files.Select(f => Path.GetFileName(f.FullFilePathWithExtension)),
                Is.EquivalentTo(searched), "named as the search reads them");
            Assert.That(design.LabelFree.Files.Select(f => f.Condition).Distinct().Count(), Is.EqualTo(8));
            Assert.That(design.Report(), Does.Contain(design.LabelFree.Report()));
        }

        [Test]
        public void AnIsobaricSearchGetsTheIsobaricReadersOwnDesign()
        {
            var searched = DataFiles("PXD008841.subset.sdrf.tsv");
            var options = new SdrfSearchDesignOptions { Tag = Tmt10(), PlexFileNamePattern = Pxd008841Pattern };

            var design = SdrfSearchDesign.Read(Fixture("PXD008841.subset.sdrf.tsv"), searched, options);

            Assert.That(design.Refusals, Is.Empty, design.Report());
            Assert.That(design.LabelFree, Is.Null);
            Assert.That(design.Isobaric, Is.Not.Null);
            Assert.That(design.Isobaric!.Files, Has.Count.EqualTo(4));
            Assert.That(design.Isobaric.PlexIds.Values.Distinct().Count(), Is.EqualTo(2));

            string written = Path.Combine(_directory, "TmtDesign.txt");
            string direct = Path.Combine(_directory, "direct.txt");
            design.WriteDesign(written);
            design.Isobaric.WriteTmtDesign(direct);
            Assert.That(File.ReadAllBytes(written), Is.EqualTo(File.ReadAllBytes(direct)));
        }

        [Test]
        public void AnIsobaricSearchWithNoPlexSourceIsRefusedByTheReader()
        {
            var design = SdrfSearchDesign.Read(Fixture("PXD008841.subset.sdrf.tsv"), DataFiles("PXD008841.subset.sdrf.tsv"),
                new SdrfSearchDesignOptions { Tag = Tmt10() });

            Assert.That(design.IsValid, Is.False);
            Assert.That(design.Isobaric, Is.Not.Null, "the reader ran and said why");
            Assert.That(design.Refusals, Is.EqualTo(design.Isobaric!.Refusals));
            Assert.That(design.Report(), Does.Contain(design.Isobaric.Report()));
        }

        [Test]
        public void TheReadersRefusalsAreTheDesignsRefusals()
        {
            var design = SdrfSearchDesign.Read(Fixture("PXD067622.sdrf.tsv"), DataFiles("PXD067622.sdrf.tsv"));

            Assert.That(design.IsValid, Is.False);
            Assert.That(design.LabelFree, Is.Not.Null);
            Assert.That(design.Refusals.Single(), Does.Contain("Several factor value columns and none declared"));
            Assert.That(design.Refusals, Is.EqualTo(design.LabelFree!.Refusals));
        }

        // ---------------------------------------------------------------- condition columns

        [TestCase("genotype", "treatment")]
        [TestCase("Factor Value[genotype]", "factor value[treatment]")]
        [TestCase(" genotype ", "FACTOR VALUE[treatment]")]
        public void AConditionColumnMayBeNamedBareOrInAnyCase(string first, string second)
        {
            var design = SdrfSearchDesign.Read(Fixture("PXD067622.sdrf.tsv"), DataFiles("PXD067622.sdrf.tsv"),
                new SdrfSearchDesignOptions { ConditionColumns = new[] { first, second } });

            Assert.That(design.Refusals, Is.Empty, design.Report());
            Assert.That(design.LabelFree!.ConditionColumns, Is.EqualTo(new[] { Genotype, Treatment }));
        }

        [Test]
        public void AConditionColumnThatIsNotAFactorIsUsedAsGiven()
        {
            var design = SdrfSearchDesign.Read(Fixture("PXD067622.sdrf.tsv"), DataFiles("PXD067622.sdrf.tsv"),
                new SdrfSearchDesignOptions { ConditionColumns = new[] { "characteristics[organism]" } });

            Assert.That(design.LabelFree!.ConditionColumns, Is.EqualTo(new[] { "characteristics[organism]" }));
        }

        // ---------------------------------------------------------------- the join to the searched files

        /// <summary>
        /// A search of half the files: the other rows are set aside before the reader sees them, and the join report says
        /// so. The reader reads the join itself (QP-S38), so no note has to warn that its line numbers or its "every row is
        /// searched" describe only the kept rows.
        /// </summary>
        [Test]
        public void ASearchOfSomeFilesSetsTheOtherRowsAsideAndTheJoinReportSaysSo()
        {
            var all = DataFiles("PXD067622.sdrf.tsv");
            var searched = all.Where(f => f.Contains("DMSO")).ToList();
            Assume.That(searched.Count, Is.InRange(1, all.Count - 1));

            var design = SdrfSearchDesign.Read(Fixture("PXD067622.sdrf.tsv"), searched,
                new SdrfSearchDesignOptions { ConditionColumns = new[] { Genotype, Treatment } });

            Assert.That(design.Refusals, Is.Empty, design.Report());
            Assert.That(design.LabelFree!.Files, Has.Count.EqualTo(searched.Count));
            Assert.That(design.Scoping!.DroppedDataFiles, Has.Count.EqualTo(all.Count - searched.Count));
            Assert.That(design.Notes, Is.EqualTo(design.LabelFree.Notes), "only the reader's own notes");
            Assert.That(design.Report(), Does.Contain($"SDRF joined to the searched files: {searched.Count} of {all.Count} SDRF rows kept.")
                .And.Contain("DroppedDataFiles (the search does not read them): ")
                .And.Not.Contain("count the kept rows"));
        }

        /// <summary>
        /// The readers read the join itself (QP-S38): a refusal names the SDRF file's line (b2.raw's bad number is on line
        /// 5, the third kept row), for a label-free and an isobaric search alike.
        /// </summary>
        [Test]
        public void TheReadersLineNumbersAreTheSdrfFilesLines()
        {
            var header = new SdrfHeader(new[]
            {
                "source name", "characteristics[biological replicate]", "comment[label]", "comment[fraction identifier]",
                "comment[technical replicate]", "comment[data file]", "comment[plex]", "factor value[condition]"
            });
            SdrfRow Row(string file, string label, string biorep, string condition, string plex = "P1") =>
                new(header, new[] { $"{condition} {biorep}", biorep, label, "1", "1", file, plex, condition });

            var labelFree = new SdrfDocument(header, new[]
            {
                Row("a1.raw", "label free sample", "1", "A"), Row("a2.raw", "label free sample", "2", "A"),
                Row("b1.raw", "label free sample", "1", "B"), Row("b2.raw", "label free sample", "x", "B"),
            });
            var design = SdrfSearchDesign.Read(labelFree, new[] { "a1-calib.mzML", "b1-calib.mzML", "b2-calib.mzML" });
            Assert.That(design.Refusals.Single(), Does.StartWith("Line 5 (b2-calib.mzML):"), design.Report());

            var isobaric = new SdrfDocument(header, new[]
            {
                Row("a.raw", "TMT126", "1", "A"), Row("c.raw", "TMT126", "1", "C", "P3"), Row("b.raw", "TMT131C", "1", "B", "P2"),
            });
            var tmt = SdrfSearchDesign.Read(isobaric, new[] { "a.raw", "b.raw" },
                new SdrfSearchDesignOptions { Tag = Tmt10(), PlexColumn = "comment[plex]" });
            Assert.That(tmt.Refusals.Single(), Does.StartWith("Line 4 (b.raw):").And.Contain("not a channel of TMT10"), tmt.Report());
        }

        /// <summary>A searched file whose name the suffix rule cannot join, matched by the acquired name the caller knows.</summary>
        [Test]
        public void AKnownAcquiredNameJoinsAFileTheSuffixRuleCannot()
        {
            var all = DataFiles("PXD067622.sdrf.tsv");
            var searched = all.Select((f, i) => i == 0 ? "renamed.mzML" : f).ToList();

            var design = SdrfSearchDesign.Read(Fixture("PXD067622.sdrf.tsv"), searched, new SdrfSearchDesignOptions
            {
                ConditionColumns = new[] { Genotype, Treatment },
                AcquiredNameOf = new Dictionary<string, string> { ["renamed.mzML"] = all[0] },
            });

            Assert.That(design.Refusals, Is.Empty, design.Report());
            Assert.That(design.LabelFree!.Files.Select(f => Path.GetFileName(f.FullFilePathWithExtension)), Does.Contain("renamed.mzML"));
        }

        [Test]
        public void ASearchOfEveryFileReportsNoJoin()
        {
            var design = SdrfSearchDesign.Read(Fixture("PXD067622.sdrf.tsv"), DataFiles("PXD067622.sdrf.tsv"),
                new SdrfSearchDesignOptions { ConditionColumns = new[] { Genotype, Treatment } });

            Assert.That(design.Report(), Does.Not.Contain("SDRF joined to the searched files"));
            Assert.That(design.Report(), Does.StartWith(design.LabelFree!.Report()));
        }

        [Test]
        public void NoRowNamingASearchedFileIsRefusedWithoutRunningTheReader()
        {
            var design = SdrfSearchDesign.Read(Fixture("PXD067622.sdrf.tsv"), new[] { "elsewhere.raw" },
                new SdrfSearchDesignOptions { ConditionColumns = new[] { Genotype, Treatment } });

            Assert.That(design.IsValid, Is.False);
            Assert.That(design.LabelFree, Is.Null);
            Assert.That(design.Isobaric, Is.Null);
            Assert.That(design.Refusals.Single(), Does.StartWith("No SDRF row names a searched file"));
            Assert.That(design.Scoping!.SearchedWithoutRow, Is.EqualTo(new[] { "elsewhere.raw" }));
            Assert.That(design.Report(), Does.StartWith("Design REFUSED: 1 reason(s).")
                .And.Contain("SearchedWithoutRow: elsewhere.raw")
                .And.Contain("  refused: No SDRF row names a searched file"));
        }

        /// <summary>An SDRF with no comment[data file] cannot be joined, so the reader keys on comment[searched data file].</summary>
        [Test]
        public void AnSdrfWithNoDataFileColumnIsReadAsGiven()
        {
            var header = new SdrfHeader(new[]
            {
                "source name", "characteristics[biological replicate]", "comment[label]", "comment[fraction identifier]",
                "comment[technical replicate]", "comment[searched data file]", "factor value[condition]"
            });
            var sdrf = new SdrfDocument(header, new[]
            {
                new SdrfRow(header, new[] { "A 1", "1", "label free sample", "1", "1", "a.mzML", "A" }),
                new SdrfRow(header, new[] { "B 1", "1", "label free sample", "1", "1", "b.mzML", "B" }),
            });

            var design = SdrfSearchDesign.Read(sdrf, new[] { "a.mzML", "b.mzML" });

            Assert.That(design.Refusals, Is.Empty, design.Report());
            Assert.That(design.Scoping, Is.Null);
            Assert.That(design.LabelFree!.FileKeyColumn, Is.EqualTo("comment[searched data file]"));
            Assert.That(design.LabelFree.Files, Has.Count.EqualTo(2));
        }

        // ---------------------------------------------------------------- refused, never thrown

        [Test]
        public void AMissingFileIsRefusedNotThrown()
        {
            string missing = Path.Combine(_directory, "nothing.sdrf.tsv");

            var design = SdrfSearchDesign.Read(missing, new[] { "a.raw" });

            Assert.That(design.IsValid, Is.False);
            Assert.That(design.Refusals.Single(), Is.EqualTo("The SDRF file does not exist: " + missing));
            Assert.That(design.Report(), Does.StartWith("Design REFUSED: 1 reason(s)."));
        }

        [Test]
        public void AnEmptyFileIsRefusedNotThrown()
        {
            string empty = Path.Combine(_directory, "empty.sdrf.tsv");
            File.WriteAllText(empty, "");

            var design = SdrfSearchDesign.Read(empty, new[] { "a.raw" });

            Assert.That(design.IsValid, Is.False);
            Assert.That(design.Refusals.Single(), Does.StartWith("The SDRF file could not be read: " + empty + ".")
                .And.Contain("SDRF file is empty"));
        }

        [Test]
        public void NoSearchedFileIsRefused()
        {
            var design = SdrfSearchDesign.Read(Fixture("PXD067622.sdrf.tsv"), new[] { " ", "" });

            Assert.That(design.Refusals.Single(), Does.StartWith("No searched file was given"));
            Assert.That(design.LabelFree, Is.Null);
        }

        [TestCase("column")]
        [TestCase("pattern")]
        [TestCase("single")]
        public void APlexSourceWithoutAKitIsRefused(string source)
        {
            var options = source switch
            {
                "column" => new SdrfSearchDesignOptions { PlexColumn = "comment[plex]" },
                "pattern" => new SdrfSearchDesignOptions { PlexFileNamePattern = Pxd008841Pattern },
                _ => new SdrfSearchDesignOptions { SinglePlex = "P1" },
            };

            var design = SdrfSearchDesign.Read(Fixture("PXD008841.subset.sdrf.tsv"), DataFiles("PXD008841.subset.sdrf.tsv"), options);

            Assert.That(design.Refusals.Single(), Does.StartWith("A plex source was given but no isobaric kit"));
            Assert.That(design.LabelFree, Is.Null);
            Assert.That(design.Isobaric, Is.Null);
        }

        [Test]
        public void ANullArgumentThrows()
        {
            Assert.Throws<ArgumentNullException>(() => SdrfSearchDesign.Read((string)null!, new[] { "a.raw" }));
            Assert.Throws<ArgumentNullException>(() => SdrfSearchDesign.Read(Fixture("PXD067622.sdrf.tsv"), null!));
            Assert.Throws<ArgumentNullException>(() => SdrfSearchDesign.Read((SdrfDocument)null!, new[] { "a.raw" }));
            Assert.Throws<ArgumentNullException>(() => SdrfSearchDesign.Read(new SdrfDocument(Fixture("PXD067622.sdrf.tsv")), null!));
        }

        // ---------------------------------------------------------------- writing

        [Test]
        public void AValidLabelFreeDesignWritesTheReadersFile()
        {
            var design = SdrfSearchDesign.Read(Fixture("PXD067622.sdrf.tsv"), DataFiles("PXD067622.sdrf.tsv"),
                new SdrfSearchDesignOptions { ConditionColumns = new[] { Genotype, Treatment } });
            string written = Path.Combine(_directory, "ExperimentalDesign.tsv");
            string direct = Path.Combine(_directory, "direct.tsv");

            design.WriteDesign(written);
            design.LabelFree!.WriteExperimentalDesignTsv(direct);

            Assert.That(File.ReadAllBytes(written), Is.EqualTo(File.ReadAllBytes(direct)));
            Assert.That(design.ToExperimentalDesign().FileNameSampleInfoDictionary, Has.Count.EqualTo(24));
        }

        [Test]
        public void ARefusedDesignWritesNothingAndBuildsNothing()
        {
            string path = Path.Combine(_directory, "ExperimentalDesign.tsv");
            var refusedByReader = SdrfSearchDesign.Read(Fixture("PXD067622.sdrf.tsv"), DataFiles("PXD067622.sdrf.tsv"));
            var refusedHere = SdrfSearchDesign.Read(Fixture("PXD067622.sdrf.tsv"), new[] { "elsewhere.raw" });

            foreach (var design in new[] { refusedByReader, refusedHere })
            {
                var thrown = Assert.Throws<InvalidOperationException>(() => design.WriteDesign(path));
                Assert.That(thrown!.Message, Does.Contain(design.Refusals.First()));
                Assert.Throws<InvalidOperationException>(() => design.ToExperimentalDesign());
                Assert.That(File.Exists(path), Is.False);
            }
            Assert.Throws<ArgumentNullException>(() => refusedByReader.WriteDesign(null!));
        }
    }
}
