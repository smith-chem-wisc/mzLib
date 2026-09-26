using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using MassSpectrometry;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Covers <see cref="SdrfLabelFreeDesign"/>, the SDRF → label-free design projection (QuantProject
    /// M7). Its whole contract is "a design MetaMorpheus accepts, or a refusal that says why", because
    /// an invalid ExperimentalDesign.tsv makes MetaMorpheus skip quantification with only a warning.
    /// The two real fixtures are aging's hand-written SDRFs; everything else is built in memory so the
    /// one defect under test is the only thing wrong.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrfLabelFreeDesign
    {
        private const string Genotype = "factor value[genotype]";
        private const string Treatment = "factor value[treatment]";

        private static string Fixture(string name) =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, "FileReadingTests",
                "ExternalFileTypes", "SdrfDesign", name);

        private static SdrfLabelFreeDesignOptions Declared(params string[] columns) =>
            new() { ConditionColumns = columns };

        /// <summary>
        /// A label-free SDRF in memory: one row per (file, condition, biorep, fraction, techrep).
        /// </summary>
        private static SdrfDocument Document(params (string File, string Condition, string Biorep, string Fraction, string Techrep)[] rows)
        {
            var header = new SdrfHeader(new[]
            {
                "source name", "characteristics[biological replicate]", "assay name", "comment[label]",
                "comment[fraction identifier]", "comment[technical replicate]", "comment[data file]",
                "factor value[condition]"
            });
            return new SdrfDocument(header, rows.Select(r => new SdrfRow(header, new[]
            {
                $"{r.Condition} {r.Biorep}", r.Biorep, "run " + r.File, "label free sample",
                r.Fraction, r.Techrep, r.File, r.Condition
            })));
        }

        /// <summary>Loads a fixture and rewrites one column's cells, for a variant of a real file.</summary>
        private static SdrfDocument Rewritten(string fixture, string column, Func<SdrfRow, string> value)
        {
            var original = new SdrfDocument(Fixture(fixture));
            var header = original.Header;
            int index = header.IndexOf(column);
            Assert.That(index, Is.GreaterThanOrEqualTo(0), $"{fixture} has no '{column}'");
            return new SdrfDocument(header, original.Results.Select(row =>
            {
                var cells = row.Cells.ToArray();
                cells[index] = value(row);
                return new SdrfRow(header, cells);
            }).ToList());
        }

        // ---------------------------------------------------------------- PXD067622

        [Test]
        public void Pxd067622WithBothFactorsDeclaredGivesEightConditionsOfThree()
        {
            var design = SdrfLabelFreeDesign.Read(Fixture("PXD067622.sdrf.tsv"), Declared(Genotype, Treatment));

            Assert.That(design.Refusals, Is.Empty, design.Report());
            Assert.That(design.Notes, Is.Empty, "the SDRF already numbers replicates 1..3 per condition");
            Assert.That(design.FileKeyColumn, Is.EqualTo("comment[data file]"));
            Assert.That(design.Files, Has.Count.EqualTo(24));

            var conditions = design.Files.GroupBy(f => f.Condition).ToList();
            Assert.That(conditions, Has.Count.EqualTo(8));
            Assert.That(conditions.Select(c => c.Key),
                Does.Contain("SPRTN-TurboID WT_DMSO (vehicle)"), "the declared columns joined with '_', in declared order");
            foreach (var condition in conditions)
            {
                // The model is 0-based: replicates 1..3 are 0..2.
                Assert.That(condition.Select(f => f.BiologicalReplicate), Is.EquivalentTo(new[] { 0, 1, 2 }), condition.Key);
                Assert.That(condition.Select(f => f.Fraction), Is.All.EqualTo(0));
                Assert.That(condition.Select(f => f.TechnicalReplicate), Is.All.EqualTo(0));
            }
        }

        [Test]
        public void Pxd067622WithNoDeclarationIsRefusedBecauseItHasTwoFactorColumns()
        {
            var design = SdrfLabelFreeDesign.Read(Fixture("PXD067622.sdrf.tsv"));

            Assert.That(design.IsValid, Is.False);
            Assert.That(design.Files, Is.Empty);
            Assert.That(design.Refusals.Single(), Does.Contain("Several factor value columns and none declared")
                .And.Contain(Genotype).And.Contain(Treatment));
        }

        /// <summary>
        /// MAP-33. A drafted SDRF that copies the study-wide index from the file names (WT_DMSO1-3,
        /// CA_DMSO4-6, ... CA_FA22-24) is ranked back to 1..3 within each condition, and the mapping
        /// is reported, giving exactly the design the hand-written numbering gives.
        /// </summary>
        [Test]
        public void StudyWideReplicateNumbersAreRankedWithinEachConditionAndReported()
        {
            var studyWide = Rewritten("PXD067622.sdrf.tsv", "characteristics[biological replicate]",
                row => Regex.Match(row["comment[data file]"]!, @"(\d+)\.raw$").Groups[1].Value);

            var design = SdrfLabelFreeDesign.Read(studyWide, Declared(Genotype, Treatment));
            var reference = SdrfLabelFreeDesign.Read(Fixture("PXD067622.sdrf.tsv"), Declared(Genotype, Treatment));

            Assert.That(design.Refusals, Is.Empty, design.Report());
            Assert.That(design.Files, Is.EqualTo(reference.Files), "same files, conditions and replicates as the per-condition numbering");
            Assert.That(design.Notes, Has.Count.EqualTo(7), "every condition but WT_DMSO, which is already 1..3");
            Assert.That(design.Notes, Does.Contain(
                "Condition 'SPRTN-TurboID CA_formaldehyde 1 mM, 1 h': biological replicates renumbered 22 -> 1, 23 -> 2, 24 -> 3."));
            Assert.That(design.Report(), Does.Contain("22 -> 1"), "the mapping is printed, not only stored");
        }

        // ---------------------------------------------------------------- PXD049018

        [Test]
        public void Pxd049018IsRefusedWhenTheUnknownTreatmentIsDeclared()
        {
            var design = SdrfLabelFreeDesign.Read(Fixture("PXD049018.sdrf.tsv"), Declared(Genotype, Treatment));

            Assert.That(design.IsValid, Is.False);
            Assert.That(design.Refusals, Has.Count.EqualTo(20), "one per row, each naming its file");
            Assert.That(design.Refusals[0], Does.Contain("MSB67868ABand_01.raw")
                .And.Contain($"'{Treatment}' is 'not available'"));
        }

        /// <summary>
        /// Leaving the unknown treatment out does NOT make the two pulldowns replicates of each
        /// other: the SDRF gives both biological replicate 1, so each band is claimed by two files.
        /// MetaMorpheus would reject that as a duplicate, so it is refused here.
        /// </summary>
        [Test]
        public void Pxd049018WithOnlyGenotypeDeclaredIsRefusedBecauseBothSamplesAreReplicateOne()
        {
            var design = SdrfLabelFreeDesign.Read(Fixture("PXD049018.sdrf.tsv"), Declared(Genotype));

            Assert.That(design.IsValid, Is.False);
            Assert.That(design.Refusals, Has.Count.EqualTo(10), "one per band");
            Assert.That(design.Refusals[0], Does.Contain("biorep 1 fraction 1 techrep 1 is named by 2 files")
                .And.Contain("MSB67868ABand_01.raw").And.Contain("MSB67869ABand_01.raw"));
        }

        /// <summary>
        /// With the treatment known, the fixture is two conditions of one sample each, ten fractions
        /// copied verbatim (MAP-34).
        /// </summary>
        [Test]
        public void Pxd049018WithTheTreatmentKnownGivesTwoConditionsOfTenFractions()
        {
            var known = Rewritten("PXD049018.sdrf.tsv", Treatment,
                row => row["source name"]!.EndsWith("68A") ? "cGAMP" : "mock");

            var design = SdrfLabelFreeDesign.Read(known, Declared(Treatment));

            Assert.That(design.Refusals, Is.Empty, design.Report());
            Assert.That(design.Files.GroupBy(f => f.Condition).Select(g => g.Key), Is.EquivalentTo(new[] { "cGAMP", "mock" }));
            var band10 = design.Files.Single(f => f.FullFilePathWithExtension == "MSB67869ABand_10.raw");
            Assert.That(band10.Fraction, Is.EqualTo(9), "fraction 10, 0-based");
            Assert.That(band10.BiologicalReplicate, Is.EqualTo(0));
        }

        // ---------------------------------------------------------------- what MetaMorpheus would reject

        [Test]
        public void AFractionGapIsRefusedButAMissingLastFractionIsNot()
        {
            var gap = SdrfLabelFreeDesign.Read(Document(
                ("a1.raw", "A", "1", "1", "1"), ("a3.raw", "A", "1", "3", "1"),
                ("b1.raw", "A", "2", "1", "1"), ("b2.raw", "A", "2", "2", "1"), ("b3.raw", "A", "2", "3", "1")));
            Assert.That(gap.Refusals.Single(), Does.StartWith("Condition 'A' biorep 1 fraction 2 is missing"));

            var missingLast = SdrfLabelFreeDesign.Read(Document(
                ("a1.raw", "A", "1", "1", "1"), ("a2.raw", "A", "1", "2", "1"),
                ("b1.raw", "A", "2", "1", "1"), ("b2.raw", "A", "2", "2", "1"), ("b3.raw", "A", "2", "3", "1")));
            Assert.That(missingLast.Refusals, Is.Empty, missingLast.Report());
        }

        [Test]
        public void ATechnicalReplicateGapIsRefused()
        {
            var design = SdrfLabelFreeDesign.Read(Document(
                ("a1.raw", "A", "1", "1", "1"), ("a3.raw", "A", "1", "1", "3")));

            Assert.That(design.Refusals.Single(), Is.EqualTo("Condition 'A' biorep 1 fraction 1 techrep 2 is missing."));
        }

        [TestCase("0")]
        [TestCase("-1")]
        [TestCase("1.5")]
        [TestCase("F1")]
        [TestCase("not available")]
        [TestCase("")]
        public void AFractionThatIsNotAPositiveIntegerIsRefused(string fraction)
        {
            var design = SdrfLabelFreeDesign.Read(Document(("a1.raw", "A", "1", fraction, "1")));

            Assert.That(design.Refusals.Single(), Is.EqualTo(
                $"Line 2 (a1.raw): 'comment[fraction identifier]' is '{fraction}', not an integer of 1 or more."));
        }

        [Test]
        public void EveryReasonIsReportedNotOnlyTheFirst()
        {
            var design = SdrfLabelFreeDesign.Read(Document(
                ("a1.raw", "A", "x", "1", "1"), ("a2.raw", "A", "1", "0", "1"), ("a3.raw", "A", "1", "1", "")));

            Assert.That(design.Refusals, Has.Count.EqualTo(3), design.Report());
        }

        /// <summary>
        /// A document with no column naming a file, no biological replicate column and no rows is
        /// refused for all three at once, and builds nothing.
        /// </summary>
        [Test]
        public void ADocumentMissingWhatEveryDesignNeedsIsRefusedForEachReason()
        {
            var header = new SdrfHeader(new[] { "source name", "comment[label]", "factor value[condition]" });
            var design = SdrfLabelFreeDesign.Read(new SdrfDocument(header, new List<SdrfRow>()));

            Assert.That(design.Refusals, Is.EqualTo(new[]
            {
                "The SDRF has neither 'comment[searched data file]' nor 'comment[data file]', so no row names a file.",
                "The SDRF has no 'characteristics[biological replicate]' column. MetaMorpheus needs a biological replicate for every file.",
                "The SDRF has no rows."
            }), design.Report());
            Assert.That(design.FileKeyColumn, Is.Null);
            Assert.That(design.Files, Is.Empty);
        }

        [Test]
        public void AFileNamedByTwoRowsIsRefused()
        {
            var design = SdrfLabelFreeDesign.Read(Document(
                ("a1.raw", "A", "1", "1", "1"), ("A1.raw", "B", "1", "1", "1")));

            Assert.That(design.Refusals.Single(), Does.Contain("'a1.raw' is named by 2 rows (lines 2, 3)"));
        }

        // ---------------------------------------------------------------- the condition (MAP-07, QP-S14)

        [Test]
        public void TheOnlyFactorColumnIsUsedWhenNoneIsDeclared()
        {
            var design = SdrfLabelFreeDesign.Read(Document(("a1.raw", "A", "1", "1", "1")));

            Assert.That(design.ConditionColumns, Is.EqualTo(new[] { "factor value[condition]" }));
            Assert.That(design.Files.Single().Condition, Is.EqualTo("A"));
        }

        /// <summary>
        /// SDRF-P1 (sdrf 012): the only factor column is refused when it is unknown in even one row,
        /// and since nothing was declared, the remedy is to declare, not to leave a column out.
        /// </summary>
        [Test]
        public void AnUnknownFactorInTheOnlyFactorColumnIsRefusedAndTheRemedyIsToDeclare()
        {
            var design = SdrfLabelFreeDesign.Read(Document(
                ("a1.raw", "A", "1", "1", "1"), ("a2.raw", "not available", "1", "1", "1")));

            Assert.That(design.Refusals.Single(), Does.StartWith("Line 3 (a2.raw): 'factor value[condition]' is 'not available'.")
                .And.Contain("it is the only factor value column and none was declared")
                .And.Contain("declare the column(s) the condition should be built from"));
        }

        [Test]
        public void AnUnknownFactorInADeclaredColumnNamesLeavingItOutAsTheRemedy()
        {
            var design = SdrfLabelFreeDesign.Read(Document(("a1.raw", "not applicable", "1", "1", "1")),
                Declared("factor value[condition]"));

            Assert.That(design.Refusals.Single(), Does.Contain(
                "fill it in, or leave the column out of the declared condition columns to pool these rows."));
        }

        [Test]
        public void ADeclaredColumnThatIsNotInTheSdrfIsRefusedNamingTheOnesThatAre()
        {
            var design = SdrfLabelFreeDesign.Read(Document(("a1.raw", "A", "1", "1", "1")), Declared("factor value[Condition]"));

            Assert.That(design.Refusals.Single(), Does.Contain("'factor value[Condition]' is not in the SDRF")
                .And.Contain("'factor value[condition]'"), "column names match case-sensitively");
        }

        [Test]
        public void AConditionColumnDeclaredTwiceIsRefused()
        {
            var design = SdrfLabelFreeDesign.Read(Document(("a1.raw", "A", "1", "1", "1")),
                Declared("factor value[condition]", "factor value[condition]"));

            Assert.That(design.Refusals.Single(), Is.EqualTo("The condition column 'factor value[condition]' is declared more than once."));
        }

        [Test]
        public void TwoFactorCombinationsThatJoinToOneConditionAreRefused()
        {
            var header = new SdrfHeader(new[] { "comment[data file]", "characteristics[biological replicate]", "factor value[x]", "factor value[y]" });
            var sdrf = new SdrfDocument(header, new[]
            {
                new SdrfRow(header, new[] { "a.raw", "1", "a_b", "c" }),
                new SdrfRow(header, new[] { "b.raw", "1", "a", "b_c" }),
            });

            var design = SdrfLabelFreeDesign.Read(sdrf, Declared("factor value[x]", "factor value[y]"));

            Assert.That(design.Refusals.Single(), Does.StartWith("Condition 'a_b_c' comes from 2 different sets of factor values"));
        }

        [Test]
        public void OneConditionSpelledTwoWaysIsRefused()
        {
            var design = SdrfLabelFreeDesign.Read(Document(("a1.raw", "DMSO", "1", "1", "1"), ("a2.raw", "dmso", "2", "1", "1")));

            Assert.That(design.Refusals.Single(), Does.Contain("comes from 2 different sets of factor values"));
        }

        // ---------------------------------------------------------------- files (MAP-12) and the search

        [Test]
        public void TheSearchedDataFileIsTheKeyWhenPresent()
        {
            var header = new SdrfHeader(new[]
            {
                "characteristics[biological replicate]", "comment[data file]", "comment[searched data file]", "factor value[condition]"
            });
            var sdrf = new SdrfDocument(header, new[] { new SdrfRow(header, new[] { "1", "a1.raw", "a1-calib.mzML", "A" }) });

            var design = SdrfLabelFreeDesign.Read(sdrf);

            Assert.That(design.FileKeyColumn, Is.EqualTo("comment[searched data file]"));
            Assert.That(design.Files.Single().FullFilePathWithExtension, Is.EqualTo("a1-calib.mzML"));
        }

        /// <summary>
        /// A row the search does not read is dropped and reported, and the ranking happens AFTER the
        /// drop: a partial download of replicates 2 and 3 is a design of replicates 1 and 2.
        /// </summary>
        [Test]
        public void RowsForFilesNotSearchedAreDroppedBeforeRanking()
        {
            var sdrf = Document(("a1.raw", "A", "1", "1", "1"), ("a2.raw", "A", "2", "1", "1"), ("a3.raw", "A", "3", "1", "1"));
            var searched = new[] { @"C:\data\a2.raw", @"C:\data\a3.raw" };

            var design = SdrfLabelFreeDesign.Read(sdrf, new SdrfLabelFreeDesignOptions { SearchedFiles = searched });

            Assert.That(design.Refusals, Is.Empty, design.Report());
            Assert.That(design.Notes, Does.Contain("Line 2 ('a1.raw') dropped: the search does not read that file."));
            Assert.That(design.Files.Select(f => (f.FullFilePathWithExtension, f.BiologicalReplicate)),
                Is.EqualTo(new[] { (@"C:\data\a2.raw", 0), (@"C:\data\a3.raw", 1) }), "the searched paths, ranked 1..2");
        }

        /// <summary>
        /// MetaMorpheus skips rows for files it does not search before it validates anything, so a
        /// problem confined to such a row must not refuse the design.
        /// </summary>
        [Test]
        public void AProblemInARowTheSearchDoesNotReadDoesNotRefuse()
        {
            var searched = new SdrfLabelFreeDesignOptions { SearchedFiles = new[] { "a2.raw", "a3.raw" } };

            var badFraction = SdrfLabelFreeDesign.Read(Document(
                ("a1.raw", "A", "1", "not available", "1"), ("a2.raw", "A", "1", "1", "1"), ("a3.raw", "A", "2", "1", "1")), searched);
            Assert.That(badFraction.Refusals, Is.Empty, badFraction.Report());
            Assert.That(badFraction.Notes, Does.Contain("Line 2 ('a1.raw') dropped: the search does not read that file."));
            Assert.That(badFraction.Files, Has.Count.EqualTo(2));

            var repeated = SdrfLabelFreeDesign.Read(Document(
                ("a1.raw", "A", "1", "1", "1"), ("a1.raw", "A", "2", "1", "1"), ("a2.raw", "A", "1", "1", "1"), ("a3.raw", "A", "2", "1", "1")), searched);
            Assert.That(repeated.Refusals, Is.Empty, repeated.Report());

            var caseCollision = SdrfLabelFreeDesign.Read(Document(
                ("a1.raw", "a", "1", "1", "1"), ("a2.raw", "A", "1", "1", "1"), ("a3.raw", "A", "2", "1", "1")), searched);
            Assert.That(caseCollision.Refusals, Is.Empty, caseCollision.Report());

            var emptyFile = SdrfLabelFreeDesign.Read(Document(
                ("", "A", "1", "1", "1"), ("a2.raw", "A", "1", "1", "1"), ("a3.raw", "A", "2", "1", "1")), searched);
            Assert.That(emptyFile.Refusals, Is.Empty, emptyFile.Report());
            Assert.That(emptyFile.Notes, Does.Contain("Line 2 dropped: 'comment[data file]' is empty, so it names no searched file."));
        }

        /// <summary>
        /// The same problems in a searched row still refuse: dropping first must not hide them.
        /// </summary>
        [Test]
        public void AProblemInASearchedRowStillRefuses()
        {
            var searched = new SdrfLabelFreeDesignOptions { SearchedFiles = new[] { "a1.raw", "a2.raw" } };
            var design = SdrfLabelFreeDesign.Read(Document(
                ("a1.raw", "A", "1", "not available", "1"), ("a2.raw", "A", "2", "1", "1")), searched);

            Assert.That(design.Refusals.Single(), Does.Contain("Line 2 (a1.raw)"));
        }

        /// <summary>
        /// MetaMorpheus matches each design row to the FIRST searched path with that name, so a second
        /// same-named file is never defined and quantification is skipped.
        /// </summary>
        [Test]
        public void TwoSearchedFilesSharingANameAreRefused()
        {
            var sdrf = Document(("x.raw", "A", "1", "1", "1"), ("y.raw", "A", "2", "1", "1"));

            var design = SdrfLabelFreeDesign.Read(sdrf, new SdrfLabelFreeDesignOptions
            {
                SearchedFiles = new[] { @"C:\dir1\x.raw", @"C:\dir2\x.raw", @"C:\dir1\y.raw" }
            });
            Assert.That(design.Refusals.Single(), Does.Contain(@"'C:\dir1\x.raw', 'C:\dir2\x.raw' share the name 'x.raw'"));
            Assert.That(design.Files, Is.Empty);

            var one = SdrfLabelFreeDesign.Read(sdrf, new SdrfLabelFreeDesignOptions
            {
                SearchedFiles = new[] { @"C:\dir1\x.raw", @"C:\dir1\x.raw", @"C:\dir1\y.raw" }
            });
            Assert.That(one.Refusals, Is.Empty, "the same path listed twice is one file");
        }

        [Test]
        public void ASearchedFileListNamingNoFileIsRefused()
        {
            foreach (var searched in new[] { Array.Empty<string>(), new[] { " " } })
            {
                var design = SdrfLabelFreeDesign.Read(Document(("a1.raw", "A", "1", "1", "1")),
                    new SdrfLabelFreeDesignOptions { SearchedFiles = searched });

                Assert.That(design.IsValid, Is.False);
                Assert.That(design.Refusals.Single(), Is.EqualTo("SearchedFiles was given but names no file, so no row could be kept."));
                Assert.That(design.Files, Is.Empty);
            }
        }

        [Test]
        public void ASearchedFileWithNoRowIsRefused()
        {
            var design = SdrfLabelFreeDesign.Read(Document(("a1.raw", "A", "1", "1", "1")),
                new SdrfLabelFreeDesignOptions { SearchedFiles = new[] { "a1.raw", "a2.raw" } });

            Assert.That(design.Refusals.Single(), Does.StartWith("Searched file 'a2.raw' has no SDRF row."));
        }

        /// <summary>
        /// MetaMorpheus matches names exactly, so a stem match is not a match. The refusal points to
        /// SdrfSearchScope, which owns the MAP-12 join, rather than joining here.
        /// </summary>
        [Test]
        public void ASearchedFileMatchedOnlyByStemIsRefusedWithAPointerToTheJoin()
        {
            var design = SdrfLabelFreeDesign.Read(Document(("a1.raw", "A", "1", "1", "1")),
                new SdrfLabelFreeDesignOptions { SearchedFiles = new[] { "A1-calib.mzML" } });

            Assert.That(design.Refusals, Has.Count.EqualTo(1), design.Report());

            var exact = SdrfLabelFreeDesign.Read(Document(("a1.raw", "A", "1", "1", "1")),
                new SdrfLabelFreeDesignOptions { SearchedFiles = new[] { "a1.mzML" } });
            Assert.That(exact.Refusals.Single(), Does.Contain("line 2 names 'a1.raw'").And.Contain("SdrfSearchScope"));
        }

        [Test]
        public void AnIsobaricSdrfIsRefused()
        {
            var header = new SdrfHeader(new[] { "characteristics[biological replicate]", "comment[label]", "comment[data file]", "factor value[condition]" });
            var sdrf = new SdrfDocument(header, new[] { new SdrfRow(header, new[] { "1", "TMT126", "plex1.raw", "A" }) });

            var design = SdrfLabelFreeDesign.Read(sdrf);

            Assert.That(design.Refusals.Single(), Does.Contain("'TMT126', not label free"));
        }

        [Test]
        public void TheAccessionedLabelFreeFormIsAccepted()
        {
            var header = new SdrfHeader(new[] { "characteristics[biological replicate]", "comment[label]", "comment[data file]", "factor value[condition]" });
            var sdrf = new SdrfDocument(header, new[] { new SdrfRow(header, new[] { "1", "AC=MS:1002038;NT=label free sample", "a.raw", "A" }) });

            Assert.That(SdrfLabelFreeDesign.Read(sdrf).Refusals, Is.Empty);
        }

        // ---------------------------------------------------------------- the written file

        /// <summary>
        /// The model is 0-based and ExperimentalDesign.tsv is 1-based; pinned in both directions,
        /// byte for byte, in MetaMorpheus's column order (FileName, Condition, Biorep, Fraction, Techrep).
        /// </summary>
        [Test]
        public void TheWrittenFileIsOneBasedInMetaMorpheusColumnOrder()
        {
            var design = SdrfLabelFreeDesign.Read(Document(
                ("a1.raw", "A", "5", "1", "1"), ("a2.raw", "A", "5", "2", "1"), ("a3.raw", "A", "5", "2", "2"), ("b1.raw", "B", "1", "1", "1")));
            Assert.That(design.Files[2], Is.EqualTo(new SpectraFileInfo("a3.raw", "A", 0, 1, 1)), "0-based, and 5 ranked to 1");

            string path = Path.Combine(TestContext.CurrentContext.WorkDirectory, $"{Guid.NewGuid():N}_ExperimentalDesign.tsv");
            try
            {
                design.WriteExperimentalDesignTsv(path);

                Assert.That(File.ReadAllText(path), Is.EqualTo(
                    "FileName\tCondition\tBiorep\tFraction\tTechrep\n" +
                    "a1.raw\tA\t1\t1\t1\n" +
                    "a2.raw\tA\t1\t2\t1\n" +
                    "a3.raw\tA\t1\t2\t2\n" +
                    "b1.raw\tB\t1\t1\t1\n"));
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public void ARefusedDesignWritesNothingAndBuildsNothing()
        {
            var design = SdrfLabelFreeDesign.Read(Document(("a1.raw", "A", "1", "0", "1")));
            string path = Path.Combine(TestContext.CurrentContext.WorkDirectory, $"{Guid.NewGuid():N}_ExperimentalDesign.tsv");

            Assert.That(() => design.WriteExperimentalDesignTsv(path), Throws.InvalidOperationException);
            Assert.That(File.Exists(path), Is.False);
            Assert.That(() => design.ToExperimentalDesign(), Throws.InvalidOperationException);
        }

        [Test]
        public void AValidDesignBecomesAnExperimentalDesignKeyedByFileName()
        {
            var design = SdrfLabelFreeDesign.Read(Document(("a1.raw", "A", "1", "1", "1"), ("b1.raw", "B", "1", "1", "1")))
                .ToExperimentalDesign();

            Assert.That(design.FileNameSampleInfoDictionary.Keys, Is.EquivalentTo(new[] { "a1.raw", "b1.raw" }));
            Assert.That(design.FileNameSampleInfoDictionary["A1.RAW"].Single().Condition, Is.EqualTo("A"));
        }
    }
}
