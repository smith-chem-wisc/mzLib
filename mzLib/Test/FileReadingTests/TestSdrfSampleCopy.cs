using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using Readers;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Tests for carrying an input SDRF row's sample columns into a built SDRF cell for cell (SdrfSampleCopy,
    /// SdrfSample.CopiedFrom; QuantProject 037, SDRF-D1). What these pin is that nothing on the way reads, normalises
    /// or reorders a copied cell, and that a copy is the row's only source of sample columns.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrfSampleCopy
    {
        private const string OrganismPart = "NT=liver;AC=UBERON:0002107;TA=x";

        private static readonly string[] Columns =
        {
            "source name", "characteristics[organism]", "characteristics[organism part]", "characteristics[age]",
            "characteristics[disease]", "characteristics[biological replicate]", "assay name", "comment[label]",
            "comment[cleavage agent details]", "comment[data file]", "comment[sample preparation batch]",
            "factor value[disease]", "factor value[treatment]"
        };

        private static SdrfDocument Input(params string[][] rows)
        {
            var header = new SdrfHeader(Columns);
            return new SdrfDocument(header, rows.Select(r => new SdrfRow(header, r)));
        }

        private static string[] InputRow(string source = "Mouse 1", string file = "a.raw", string replicate = "1") => new[]
        {
            source, "NT=Mus musculus;AC=NCBITaxon:10090", OrganismPart, "12M", "not available", replicate, "run " + file,
            "label free sample", "NT=Glu-C;AC=MS:1001917", file, "batch_20220601", "normal", "vehicle"
        };

        private static SdrfSample Copied(SdrfDocument input, int row = 0)
        {
            var copy = SdrfSampleCopy.FromRow(input, input.Results[row]);
            return new SdrfSample
            {
                SourceName = copy.SourceName,
                CopiedFrom = copy,
                Label = new CvParam("", "", "label free sample", "")
            };
        }

        private static SdrfAssay Assay(string file = "a.raw") => new()
        {
            DataFileName = file,
            AssayName = "run " + file,
            Instrument = new CvParam("MS", "MS:1001911", "Q Exactive", ""),
            PrecursorMassTolerance = new PpmTolerance(5),
            ProductMassTolerance = new PpmTolerance(20),
            CleavageAgent = ProteaseDictionary.Dictionary["trypsin"],
            DissociationType = DissociationType.HCD,
            AcquisitionMethod = new CvParam("PRIDE", "PRIDE:0000627", "Data-dependent acquisition", ""),
            TechnicalReplicate = 1,
            Fraction = 1
        };

        private static SdrfDocument Build(params SdrfRowInput[] rows) => SdrfBuilder.Build(rows);

        /// <summary>
        /// The test the copy exists for. Every sample cell comes out byte-equal, in the input's order, including a key
        /// the term parser does not know: the same organism part through CvParam loses its TA=x.
        /// </summary>
        [Test]
        public void EveryCopiedSampleCellIsWrittenByteForByteInTheInputsOrder()
        {
            var input = Input(InputRow());
            var built = Build(new SdrfRowInput(Copied(input), Assay()));

            var copiedColumns = Columns.Where(SdrfSampleCopy.IsSampleColumn).ToList();
            var header = built.Header.ToList();
            Assert.That(header.Where(SdrfSampleCopy.IsSampleColumn), Is.EqualTo(copiedColumns), "same columns, same order");
            foreach (var column in copiedColumns)
                Assert.That(built.Results[0][column], Is.EqualTo(input.Results[0][column]), column);

            Assert.That(SdrfCell.TryParseTerm(OrganismPart, out var term), Is.True);
            Assert.That(SdrfCell.ToCell(term), Is.Not.EqualTo(OrganismPart), "a round trip through CvParam would have changed it");
        }

        /// <summary>source name first, then characteristics, the assay block, the copied comments, and factors last.</summary>
        [Test]
        public void CopiedColumnsAreGroupedWhereTheBuilderPutsEachKind()
        {
            var built = Build(new SdrfRowInput(Copied(Input(InputRow())), Assay()));
            var h = built.Header.ToList();

            Assert.That(h[0], Is.EqualTo("source name"));
            Assert.That(h.LastIndexOf("characteristics[biological replicate]"), Is.LessThan(h.IndexOf("assay name")));
            Assert.That(h.IndexOf("comment[sample preparation batch]"), Is.GreaterThan(h.IndexOf("comment[data file]")));
            Assert.That(h.IndexOf("factor value[disease]"), Is.GreaterThan(h.IndexOf("comment[sample preparation batch]")));
            Assert.That(h.IndexOf("factor value[treatment]"), Is.EqualTo(h.Count - 1));
            Assert.That(SdrfValidator.Validate(built).IsValid, Is.True,
                string.Join(" | ", SdrfValidator.Validate(built).Errors.Select(e => e.ToString())));
        }

        /// <summary>Nine corpus files repeat characteristics[organism part]; the per-field inputs keep only the first.</summary>
        [Test]
        public void ARepeatedColumnIsWrittenTwiceInOrder()
        {
            var columns = new[] { "source name", "characteristics[organism part]", "characteristics[organism part]", "comment[data file]" };
            var header = new SdrfHeader(columns);
            var input = new SdrfDocument(header, new[] { new SdrfRow(header, new[] { "S1", "NT=liver;AC=UBERON:0002107", "NT=heart;AC=UBERON:0000948", "a.raw" }) });

            var built = Build(new SdrfRowInput(Copied(input), Assay()));

            Assert.That(built.Results[0].All("characteristics[organism part]"),
                Is.EqualTo(new[] { "NT=liver;AC=UBERON:0002107", "NT=heart;AC=UBERON:0000948" }));
        }

        /// <summary>
        /// organism and biological replicate are written from dedicated fields on the per-field path, and the
        /// replicate there is a number. A copy carries both as written, so a pooled reference stays pooled (G49).
        /// </summary>
        [Test]
        public void OrganismAndAPooledReplicateAreCopiedAsWritten()
        {
            var built = Build(new SdrfRowInput(Copied(Input(InputRow(replicate: "pooled"))), Assay()));

            Assert.That(built.Results[0]["characteristics[organism]"], Is.EqualTo("NT=Mus musculus;AC=NCBITaxon:10090"));
            Assert.That(built.Results[0]["characteristics[biological replicate]"], Is.EqualTo("pooled"));
            Assert.That(built.Header.IndexesOf("characteristics[organism]"), Has.Count.EqualTo(1), "written once, not also from Organism");
        }

        /// <summary>The input's own not available is a statement; RequireSampleMetadata (the default) must not refuse it.</summary>
        [Test]
        public void ACopiedReservedWordSatisfiesRequireSampleMetadata()
        {
            var built = SdrfBuilder.Build(new[] { new SdrfRowInput(Copied(Input(InputRow())), Assay()) },
                new SdrfBuilderOptions { RequireSampleMetadata = true });

            Assert.That(built.Results[0]["characteristics[disease]"], Is.EqualTo("not available"));
        }

        /// <summary>
        /// The input's comment[cleavage agent details] says Glu-C, but this search used trypsin: an assay column is
        /// the search's to write, and is never copied.
        /// </summary>
        [Test]
        public void ABuilderOwnedCommentInTheInputIsNotCopied()
        {
            var built = Build(new SdrfRowInput(Copied(Input(InputRow())), Assay()));

            Assert.That(built.Header.IndexesOf("comment[cleavage agent details]"), Has.Count.EqualTo(1));
            var perField = SdrfBuilder.Build(new[] { new SdrfRowInput(new SdrfSample
            {
                SourceName = "S",
                Organism = new CvParam("NCBITaxon", "NCBITaxon:10090", "Mus musculus", ""),
                Label = new CvParam("", "", "label free sample", "")
            }, Assay()) }, new SdrfBuilderOptions { RequireSampleMetadata = false });
            Assert.That(built.Results[0]["comment[cleavage agent details]"],
                Is.EqualTo(perField.Results[0]["comment[cleavage agent details]"]), "the search's trypsin, as the per-field path writes it");
            Assert.That(built.Results[0]["comment[cleavage agent details]"], Does.Not.Contain("Glu-C"));
            Assert.That(SdrfSampleCopy.FromRow(Input(InputRow()), Input(InputRow()).Results[0]).Cells.Select(c => c.Key),
                Has.None.EqualTo("comment[label]").And.None.EqualTo("assay name").And.None.EqualTo("comment[data file]"));
        }

        private static IEnumerable<TestCaseData> FieldsThatConflictWithACopy()
        {
            yield return new TestCaseData(nameof(SdrfSample.Organism),
                (Func<SdrfSample, SdrfSample>)(s => s with { Organism = new CvParam("NCBITaxon", "NCBITaxon:9606", "Homo sapiens", "") }));
            yield return new TestCaseData(nameof(SdrfSample.Characteristics),
                (Func<SdrfSample, SdrfSample>)(s => s with { Characteristics = new Dictionary<string, CvParam> { ["characteristics[cell type]"] = new CvParam("CL", "CL:0000182", "hepatocyte", "") } }));
            yield return new TestCaseData(nameof(SdrfSample.RawCharacteristics),
                (Func<SdrfSample, SdrfSample>)(s => s with { RawCharacteristics = new Dictionary<string, string> { ["characteristics[sex]"] = "male" } }));
            yield return new TestCaseData(nameof(SdrfSample.BiologicalReplicate),
                (Func<SdrfSample, SdrfSample>)(s => s with { BiologicalReplicate = 2 }));
            yield return new TestCaseData(nameof(SdrfSample.BiologicalReplicate) + " (null)",
                (Func<SdrfSample, SdrfSample>)(s => s with { BiologicalReplicate = null }));
            yield return new TestCaseData(nameof(SdrfSample.FactorValue),
                (Func<SdrfSample, SdrfSample>)(s => s with { FactorValue = "treated" }));
            yield return new TestCaseData(nameof(SdrfSample.FactorValueColumn),
                (Func<SdrfSample, SdrfSample>)(s => s with { FactorValueColumn = "factor value[treatment]" }));
            yield return new TestCaseData(nameof(SdrfSample.FactorValues),
                (Func<SdrfSample, SdrfSample>)(s => s with { FactorValues = new Dictionary<string, string> { ["factor value[treatment]"] = "treated" } }));
            yield return new TestCaseData(nameof(SdrfSample.Comments),
                (Func<SdrfSample, SdrfSample>)(s => s with { Comments = new Dictionary<string, string> { ["comment[characteristics source]"] = "config" } }));
        }

        [TestCaseSource(nameof(FieldsThatConflictWithACopy))]
        public void ACopyPlusAPerFieldSampleValueThrowsNamingTheField(string field, Func<SdrfSample, SdrfSample> set)
        {
            var sample = set(Copied(Input(InputRow())));

            var e = Assert.Throws<ArgumentException>(() => Build(new SdrfRowInput(sample, Assay())));
            Assert.That(e!.Message, Does.Contain(field.Split(' ')[0]));
        }

        [Test]
        public void ASourceNameThatIsNotTheCopiedOneThrows()
        {
            var sample = Copied(Input(InputRow())) with { SourceName = "Sample 1" };

            var e = Assert.Throws<ArgumentException>(() => Build(new SdrfRowInput(sample, Assay())));
            Assert.That(e!.Message, Does.Contain("Mouse 1"));
        }

        [Test]
        public void TwoRowsCopyingDifferentColumnSequencesThrow()
        {
            var other = new[] { "source name", "characteristics[age]", "characteristics[organism]", "comment[data file]" };
            var header = new SdrfHeader(other);
            var second = new SdrfDocument(header, new[] { new SdrfRow(header, new[] { "Mouse 2", "12M", "NT=Mus musculus;AC=NCBITaxon:10090", "b.raw" }) });

            Assert.Throws<ArgumentException>(() => Build(
                new SdrfRowInput(Copied(Input(InputRow())), Assay("a.raw")),
                new SdrfRowInput(Copied(second), Assay("b.raw"))));
        }

        [Test]
        public void SomeRowsCopyingAndOthersNotThrows()
        {
            var perField = new SdrfSample
            {
                SourceName = "Mouse 2",
                Organism = new CvParam("NCBITaxon", "NCBITaxon:10090", "Mus musculus", ""),
                Label = new CvParam("", "", "label free sample", "")
            };

            Assert.Throws<ArgumentException>(() => Build(
                new SdrfRowInput(Copied(Input(InputRow())), Assay("a.raw")),
                new SdrfRowInput(perField, Assay("b.raw"))));
        }

        /// <summary>A short row (PXD059974 has 17) reaches the builder as a blank cell, which it refuses, naming the column.</summary>
        [Test]
        public void ABlankOrMissingCopiedCellThrowsNamingTheColumn()
        {
            var shortRow = InputRow().Take(4).ToArray();

            var e = Assert.Throws<ArgumentException>(() => Build(new SdrfRowInput(Copied(Input(shortRow)), Assay())));
            Assert.That(e!.Message, Does.Contain("characteristics[disease]"));
        }

        [Test]
        public void ACopiedCellWithATabThrows()
        {
            var row = InputRow();
            row[3] = "12M\textra";

            var e = Assert.Throws<ArgumentException>(() => Build(new SdrfRowInput(Copied(Input(row)), Assay())));
            Assert.That(e!.Message, Does.Contain("characteristics[age]"));
        }

        [Test]
        public void FromRowRefusesARowFromAnotherHeaderAndAHeaderWithoutOneSourceName()
        {
            var input = Input(InputRow());
            var otherHeader = new SdrfHeader(new[] { "source name", "comment[data file]" });

            Assert.Throws<ArgumentException>(() => SdrfSampleCopy.FromRow(input, new SdrfRow(otherHeader, new[] { "S", "a.raw" })));

            var noSource = new SdrfHeader(new[] { "characteristics[organism]", "comment[data file]" });
            var doc = new SdrfDocument(noSource, new[] { new SdrfRow(noSource, new[] { "x", "a.raw" }) });
            Assert.Throws<ArgumentException>(() => SdrfSampleCopy.FromRow(doc, doc.Results[0]));
            Assert.Throws<ArgumentNullException>(() => SdrfSampleCopy.FromRow(null!, input.Results[0]));
            Assert.Throws<ArgumentNullException>(() => SdrfSampleCopy.FromRow(input, null!));
        }

        /// <summary>
        /// The MetaMorpheus flow, end to end: restrict the input SDRF to a calibrated search's files with the public
        /// join, copy each kept row, and build. The copied sample half is the input's; the searched file is the search's.
        /// </summary>
        [Test]
        public void TheJoinAndTheCopyTogetherCarryACalibratedSearchsSamples()
        {
            var input = Input(InputRow("Mouse 1", "a.raw"), InputRow("Mouse 2", "b.raw", "2"), InputRow("Mouse 3", "c.raw", "3"));
            var scoped = SdrfSearchScope.Restrict(input, new[] { @"C:\out\a-calib.mzML", @"C:\out\c-calib.mzML" });
            Assert.That(scoped.SearchedWithoutRow, Is.Empty);

            var rows = scoped.Document.Results.Select(r =>
            {
                var copy = SdrfSampleCopy.FromRow(scoped.Document, r);
                var sample = new SdrfSample { SourceName = copy.SourceName, CopiedFrom = copy, Label = new CvParam("", "", "label free sample", "") };
                return new SdrfRowInput(sample, Assay(r["comment[data file]"]!) with { SearchedDataFileName = r["comment[searched data file]"] });
            }).ToArray();
            var built = Build(rows);

            Assert.That(built.Results.Select(r => r["source name"]), Is.EqualTo(new[] { "Mouse 1", "Mouse 3" }));
            Assert.That(built.Results.Select(r => r["characteristics[biological replicate]"]), Is.EqualTo(new[] { "1", "3" }));
            Assert.That(built.Results.Select(r => r["comment[searched data file]"]), Is.EqualTo(new[] { "a-calib.mzML", "c-calib.mzML" }));
            Assert.That(built.Header.IndexesOf("comment[searched data file]"), Has.Count.EqualTo(1), "the join's column is the builder's, not copied");
        }
    }
}
