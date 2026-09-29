using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Tests for lifting the sample half of a document out of it, keyed by source name.
    ///
    /// The fixtures are deliberately awkward in the ways real documents are awkward: one sample over
    /// several rows, mis-cased names and values, a repeated column, a column whose rows disagree, and
    /// a blank source name. What is under test is never "does it parse" -- SdrfDocument does that --
    /// but whether the right cells end up under the right sample, and whether the cases that cannot
    /// be resolved are REPORTED rather than resolved by guessing.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrfSampleBlock
    {
        private static readonly string[] Columns =
        {
            "source name", "characteristics[organism]", "characteristics[organism part]",
            "characteristics[biological replicate]", "assay name", "comment[data file]",
            "comment[label]", "factor value[disease]"
        };

        private static SdrfDocument Doc(string[] columns, params string[][] rows)
        {
            var header = new SdrfHeader(columns);
            return new SdrfDocument(header, rows.Select(r => new SdrfRow(header, r)));
        }

        /// <summary>
        /// One sample measured in three runs is ONE sample. That is the whole point of keying on
        /// source name alone (D27): the block is what a search matches a row it is writing against.
        /// </summary>
        [Test]
        public void RowsOfOneSampleBecomeOneBlock()
        {
            var blocks = SdrfSampleBlock.BySourceName(Doc(Columns,
                new[] { "S1", "Homo sapiens", "liver", "1", "run 1", "a.raw", "TMT126", "normal" },
                new[] { "S1", "Homo sapiens", "liver", "1", "run 2", "b.raw", "TMT127N", "normal" },
                new[] { "S2", "Homo sapiens", "liver", "2", "run 3", "c.raw", "TMT126", "carcinoma" }),
                out var problems);

            Assert.That(problems, Is.Empty);
            Assert.That(blocks.Keys, Is.EquivalentTo(new[] { "S1", "S2" }));
            Assert.That(blocks["S1"].RowCount, Is.EqualTo(2));
            Assert.That(blocks["S1"]["characteristics[organism part]"], Is.EqualTo("liver"));
            Assert.That(blocks["S2"]["factor value[disease]"], Is.EqualTo("carcinoma"));
        }

        /// <summary>
        /// Only the sample half travels. The assay columns are the search's own business -- it knows
        /// its enzyme and its instrument better than a deposited file does -- and copying them
        /// through would let someone else's assay overwrite the one actually run.
        /// </summary>
        [Test]
        public void AssayColumnsAreNotPartOfASampleBlock()
        {
            var block = SdrfSampleBlock.BySourceName(Doc(Columns,
                new[] { "S1", "Homo sapiens", "liver", "1", "run 1", "a.raw", "TMT126", "normal" }),
                out _)["S1"];

            Assert.That(block.Cells.Keys, Is.EquivalentTo(new[]
            {
                "source name", "characteristics[organism]", "characteristics[organism part]",
                "characteristics[biological replicate]", "factor value[disease]"
            }));
            Assert.That(block["comment[data file]"], Is.Null, "An assay column is not a sample fact.");
            Assert.That(block["assay name"], Is.Null);
            Assert.That(block.CharacteristicColumns.Count(), Is.EqualTo(3));
            Assert.That(block.FactorValueColumns.Single(), Is.EqualTo("factor value[disease]"));
        }

        /// <summary>
        /// The two matching rules are different rules, and this pins both at once: a source name is a
        /// VALUE, so "S1" and "s1" are one sample; a column name is a NAME, so
        /// characteristics[age] and characteristics[Age] are two columns.
        /// </summary>
        [Test]
        public void SourceNamesMatchCaseInsensitively_ButColumnNamesDoNot()
        {
            string[] columns = { "source name", "characteristics[age]", "characteristics[Age]" };

            var blocks = SdrfSampleBlock.BySourceName(Doc(columns,
                new[] { "S1", "58Y", "58Y" },
                new[] { "s1", "58Y", "58Y" }),
                out _);

            Assert.That(blocks, Has.Count.EqualTo(1), "One sample, spelled two ways.");
            Assert.That(blocks["S1"].RowCount, Is.EqualTo(2));
            Assert.That(blocks["s1"].SourceName, Is.EqualTo("S1"),
                "The block keeps the spelling of the first row that named the sample.");

            var block = blocks["S1"];
            Assert.That(block.Cells.Keys, Is.EquivalentTo(new[]
                { "source name", "characteristics[age]", "characteristics[Age]" }),
                "Two spellings of a column are two columns: unifying them is a curation decision, " +
                "and this is not the layer that gets to make it.");
        }

        /// <summary>
        /// Values agree case-insensitively -- PXD023158 writes tmt126 and means TMT126 -- but the
        /// value kept is whichever the first row wrote. A block is a copy, not a normalisation.
        /// </summary>
        [Test]
        public void ValuesAgreeIgnoringCaseAndSurroundingSpace_AndAreKeptVerbatim()
        {
            var block = SdrfSampleBlock.BySourceName(Doc(Columns,
                new[] { "S1", "Homo sapiens", "Liver", "1", "run 1", "a.raw", "TMT126", "normal" },
                new[] { "S1", "homo sapiens", " liver ", "1", "run 2", "b.raw", "TMT127N", "NORMAL" }),
                out _)["S1"];

            Assert.That(block.ConflictingColumns, Is.Empty,
                "Case and padding are not disagreement; treating them as such would withhold most " +
                "of a real document.");
            Assert.That(block["characteristics[organism part]"], Is.EqualTo("Liver"),
                "Kept exactly as the first row wrote it, capital and all.");
        }

        /// <summary>
        /// A column whose rows genuinely disagree is WITHHELD and named, never resolved by taking the
        /// first. "First row wins" is how a sample silently acquires the wrong organism part, and
        /// nothing downstream could ever see that it happened.
        /// </summary>
        [Test]
        public void AColumnItsRowsDisagreeAboutIsWithheldAndNamed()
        {
            var block = SdrfSampleBlock.BySourceName(Doc(Columns,
                new[] { "S1", "Homo sapiens", "liver", "1", "run 1", "a.raw", "TMT126", "normal" },
                new[] { "S1", "Homo sapiens", "brain", "1", "run 2", "b.raw", "TMT127N", "normal" }),
                out _)["S1"];

            Assert.That(block.ConflictingColumns, Is.EqualTo(new[] { "characteristics[organism part]" }));
            Assert.That(block["characteristics[organism part]"], Is.Null,
                "Neither liver nor brain: the document does not know, and neither does this.");
            Assert.That(block["characteristics[organism]"], Is.EqualTo("Homo sapiens"),
                "One conflicted column does not cost the sample every other cell.");
        }

        /// <summary>
        /// SDRF columns may repeat -- nine corpus files repeat characteristics[organism part] -- so a
        /// block holds a list per column. Joining them would invent a value no row contains.
        /// </summary>
        [Test]
        public void ARepeatedColumnKeepsEveryValue()
        {
            string[] columns =
                { "source name", "characteristics[organism part]", "characteristics[organism part]" };

            var block = SdrfSampleBlock.BySourceName(Doc(columns,
                new[] { "S1", "liver", "left lobe" }),
                out _)["S1"];

            Assert.That(block.All("characteristics[organism part]"),
                Is.EqualTo(new[] { "liver", "left lobe" }));
            Assert.That(block["characteristics[organism part]"], Is.EqualTo("liver"),
                "The indexer gives the first, exactly as SdrfRow's does.");
        }

        /// <summary>
        /// SourceName is the join key, so it is trimmed; the source name CELL is a copy like every
        /// other cell, so it is not.
        /// </summary>
        [Test]
        public void APaddedSourceNameIsTrimmedAsTheKeyButKeptVerbatimAsACell()
        {
            var blocks = SdrfSampleBlock.BySourceName(Doc(Columns,
                new[] { " S1 ", "Homo sapiens", "liver", "1", "run 1", "a.raw", "TMT126", "normal" }),
                out _);

            Assert.That(blocks["S1"].SourceName, Is.EqualTo("S1"));
            Assert.That(blocks["S1"]["source name"], Is.EqualTo(" S1 "));
        }

        /// <summary>
        /// The column lists come back in the order the header wrote them, whatever order a
        /// dictionary happens to enumerate in.
        /// </summary>
        [Test]
        public void ColumnsComeBackInHeaderOrder()
        {
            string[] columns =
            {
                "source name", "characteristics[sex]", "factor value[treatment]",
                "characteristics[age]", "assay name", "factor value[disease]"
            };

            var block = SdrfSampleBlock.BySourceName(Doc(columns,
                new[] { "S1", "female", "drug", "58Y", "run 1", "normal" }),
                out _)["S1"];

            Assert.That(block.CharacteristicColumns,
                Is.EqualTo(new[] { "characteristics[sex]", "characteristics[age]" }));
            Assert.That(block.FactorValueColumns,
                Is.EqualTo(new[] { "factor value[treatment]", "factor value[disease]" }));
        }

        /// <summary>
        /// A row with no source name cannot be keyed to a sample, so it is reported rather than
        /// dropped quietly or filed under the empty string.
        /// </summary>
        [Test]
        public void ARowWithNoSourceNameIsReported()
        {
            var blocks = SdrfSampleBlock.BySourceName(Doc(Columns,
                new[] { "S1", "Homo sapiens", "liver", "1", "run 1", "a.raw", "TMT126", "normal" },
                new[] { "   ", "Homo sapiens", "liver", "1", "run 2", "b.raw", "TMT127N", "normal" }),
                out var problems);

            Assert.That(blocks, Has.Count.EqualTo(1));
            Assert.That(problems, Has.Count.EqualTo(1));
            Assert.That(problems.Single(), Does.Contain("Row 2").And.Contain("source name"));
        }

        /// <summary>
        /// A document with no source name column describes no samples, and says so rather than
        /// returning an empty dictionary that looks like a document with no rows.
        /// </summary>
        [Test]
        public void ADocumentWithNoSourceNameColumnSaysSo()
        {
            var blocks = SdrfSampleBlock.BySourceName(
                Doc(new[] { "assay name", "comment[data file]" }, new[] { "run 1", "a.raw" }),
                out var problems);

            Assert.That(blocks, Is.Empty);
            Assert.That(problems.Single(), Does.Contain("source name"));
        }

        /// <summary>
        /// Reserved words travel through untouched. An SDRF built from a config design may carry
        /// "not available" for facts PRIDE does not hold (D27), and that word is the file's one
        /// honest statement about those cells -- rewriting or dropping it would erase it.
        /// </summary>
        [Test]
        public void ReservedWordsAreCarriedThroughVerbatim()
        {
            var block = SdrfSampleBlock.BySourceName(Doc(Columns,
                new[] { "S1", "Homo sapiens", "not available", "1", "run 1", "a.raw", "TMT126", "not applicable" }),
                out _)["S1"];

            Assert.That(block["characteristics[organism part]"], Is.EqualTo("not available"));
            Assert.That(block["factor value[disease]"], Is.EqualTo("not applicable"));
        }

        /// <summary>
        /// Two corpus files spell the prefix "Factor Value[", and missing a factor column loses the
        /// one cell that says what the study varied -- so the prefix match ignores case even though
        /// exact column names do not.
        /// </summary>
        [Test]
        public void AMisCasedColumnPrefixIsStillASampleColumn()
        {
            var block = SdrfSampleBlock.BySourceName(
                Doc(new[] { "source name", "Factor Value[disease]", "Characteristics[organism]" },
                    new[] { "S1", "normal", "Homo sapiens" }),
                out _)["S1"];

            Assert.That(block["Factor Value[disease]"], Is.EqualTo("normal"),
                "Found by the header's own spelling, which is what SdrfRow.All matches ordinally.");
            Assert.That(block.FactorValueColumns.Single(), Is.EqualTo("Factor Value[disease]"));
            Assert.That(block.CharacteristicColumns.Single(), Is.EqualTo("Characteristics[organism]"));
        }

        /// <summary>
        /// A real curated document, so the shape is exercised against something nobody wrote for this
        /// test. PXD000070 is the committed fixture the round-trip tests use.
        /// </summary>
        [Test]
        public void ACuratedDocumentSplitsIntoItsSamples()
        {
            var document = new SdrfDocument(Path.Combine(
                TestContext.CurrentContext.TestDirectory, "FileReadingTests", "ExternalFileTypes",
                "PXD000070.sdrf.tsv"));

            var blocks = SdrfSampleBlock.BySourceName(document, out var problems);

            Assert.That(problems, Is.Empty, string.Join(" | ", problems));
            Assert.That(blocks, Is.Not.Empty);
            Assert.That(blocks.Values.Sum(b => b.RowCount), Is.EqualTo(document.Results.Count),
                "Every row belongs to exactly one sample.");
            Assert.That(blocks.Values.All(b => b.Cells.ContainsKey("source name")), Is.True);
            Assert.That(blocks.Values.All(b => b.Cells.Keys.All(
                    c => !c.StartsWith("comment[") && c != "assay name" && c != "technology type")), Is.True,
                "No assay column may appear in a sample block.");
        }
    }
}
