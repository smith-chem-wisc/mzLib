using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
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
            "characteristics[disease]", "characteristics[biological replicate]", "material type", "assay name",
            "comment[label]", "comment[cleavage agent details]", "comment[data file]", "comment[sample preparation batch]",
            "factor value[disease]", "factor value[treatment]"
        };

        private static SdrfDocument Input(params string[][] rows)
        {
            var header = new SdrfHeader(Columns);
            return new SdrfDocument(header, rows.Select(r => new SdrfRow(header, r)));
        }

        private static string[] InputRow(string source = "Mouse 1", string file = "a.raw", string replicate = "1") => new[]
        {
            source, "NT=Mus musculus;AC=NCBITaxon:10090", OrganismPart, "12M", "not available", replicate, "tissue",
            "run " + file, "label free sample", "NT=Glu-C;AC=MS:1001917", file, "batch_20220601", "normal", "vehicle"
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

        /// <summary>
        /// source name first, then characteristics, material type, the assay block, the copied comments, and factors
        /// last.
        /// </summary>
        [Test]
        public void CopiedColumnsAreGroupedWhereTheBuilderPutsEachKind()
        {
            var built = Build(new SdrfRowInput(Copied(Input(InputRow())), Assay()));
            var h = built.Header.ToList();

            Assert.That(h[0], Is.EqualTo("source name"));
            Assert.That(h.LastIndexOf("characteristics[biological replicate]"), Is.LessThan(h.IndexOf("material type")));
            Assert.That(h.IndexOf("material type"), Is.LessThan(h.IndexOf("assay name")));
            Assert.That(built.Results[0]["material type"], Is.EqualTo("tissue"));
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

        /// <summary>
        /// A short row reaches the builder as blank cells. Cut inside the characteristics, it leaves a blank sample
        /// fact, which the builder refuses, naming the column.
        /// </summary>
        [Test]
        public void AShortRowCutInsideTheCharacteristicsThrowsNamingTheColumn()
        {
            var shortRow = InputRow().Take(4).ToArray();

            var e = Assert.Throws<ArgumentException>(() => Build(new SdrfRowInput(Copied(Input(shortRow)), Assay())));
            Assert.That(e!.Message, Does.Contain("characteristics[disease]"));
        }

        [TestCase("source name")]
        [TestCase("characteristics[age]")]
        [TestCase("material type")]
        [TestCase("factor value[treatment]")]
        public void ABlankSampleFactThrowsNamingTheColumn(string column)
        {
            var row = InputRow();
            row[Array.IndexOf(Columns, column)] = " ";

            var e = Assert.Throws<ArgumentException>(() => Build(new SdrfRowInput(Copied(Input(row)), Assay())));
            Assert.That(e!.Message, Does.Contain($"blank '{column}'"));
        }

        /// <summary>
        /// A blank comment says nothing about the sample, so it is carried blank, as the input has it. Refusing it
        /// refused 357 of 12,298 curated SDRFs, 322 of them only for comment[file uri] (pcruzparri, #1451).
        /// </summary>
        [Test]
        public void ABlankCopiedCommentIsWrittenBlank()
        {
            var fixture = new SdrfDocument(Fixture("PXD000070.sdrf.tsv"));
            int uri = fixture.Header.IndexOf("comment[file uri]");
            Assert.That(uri, Is.GreaterThanOrEqualTo(0), "the fixture has comment[file uri]");
            var rows = fixture.Results.Select((r, i) =>
                new SdrfRow(fixture.Header, r.Cells.Select((c, j) => i == 0 && j == uri ? "" : c).ToArray())).ToList();
            var input = new SdrfDocument(fixture.Header, rows);

            var built = Build(input.Results.Select((r, i) => new SdrfRowInput(Copied(input, i), Assay($"f{i}.raw"))).ToArray());

            Assert.That(built.Results[0]["comment[file uri]"], Is.EqualTo(""));
            Assert.That(built.Results[1]["comment[file uri]"], Is.EqualTo(fixture.Results[1]["comment[file uri]"]));
            Assert.That(built.Results[0]["characteristics[organism]"], Is.EqualTo(fixture.Results[0]["characteristics[organism]"]));
        }

        /// <summary>
        /// Every real SDRF fixture in this repo, copied row by row: each input column is in the output, and each copied
        /// cell comes back byte-equal. The two named exceptions are the documented ones (SdrfSampleCopy, "What is not
        /// carried"): the accession, which comes through the option, and a case variant of a sample column.
        /// </summary>
        [Test]
        public void EveryRealFixtureIsCarriedWholeExceptTheNamedColumns()
        {
            var notCarried = new Dictionary<string, string[]> { ["PXD026824.sdrf.tsv"] = new[] { "Material Type" } };
            var files = Directory.GetFiles(Fixture(), "*.sdrf.tsv", SearchOption.AllDirectories);
            Assert.That(files, Has.Length.GreaterThanOrEqualTo(11));

            Assert.Multiple(() =>
            {
                foreach (var path in files)
                {
                    var input = new SdrfDocument(path);
                    string? accession = input.Header.Contains(PxAccession) ? input.Results[0][PxAccession] : null;
                    var copies = input.Results.Select(r => SdrfSampleCopy.FromRow(input, r)).ToList();

                    var built = SdrfBuilder.Build(
                        copies.Select((c, i) => new SdrfRowInput(Copied(input, i), Assay($"f{i}.raw"))).ToArray(),
                        new SdrfBuilderOptions { ProteomeXchangeAccession = accession });

                    string name = Path.GetFileName(path);
                    Assert.That(input.Header.Distinct().Except(built.Header),
                        Is.EquivalentTo(notCarried.GetValueOrDefault(name) ?? Array.Empty<string>()), name);
                    for (int i = 0; i < copies.Count; i++)
                        foreach (var column in copies[i].Cells.Select(c => c.Key).Distinct())
                            Assert.That(built.Results[i].All(column),
                                Is.EqualTo(copies[i].Cells.Where(c => c.Key == column).Select(c => c.Value)), $"{name} row {i} {column}");
                }
            });
        }

        /// <summary>The accession is one of the builder's own columns: without the option, a copy does not carry it.</summary>
        [Test]
        public void TheAccessionIsCarriedOnlyThroughTheOption()
        {
            var input = new SdrfDocument(Fixture("PXD000070.sdrf.tsv"));
            var rows = input.Results.Select((r, i) => new SdrfRowInput(Copied(input, i), Assay($"f{i}.raw"))).ToArray();

            Assert.That(input.Header.Except(SdrfBuilder.Build(rows).Header), Is.EqualTo(new[] { PxAccession }));
            var withOption = SdrfBuilder.Build(rows, new SdrfBuilderOptions { ProteomeXchangeAccession = "PXD000070" });
            Assert.That(withOption.Results.Select(r => r[PxAccession]), Is.All.EqualTo("PXD000070"));
        }

        private const string PxAccession = "comment[proteomexchange accession number]";

        private static string Fixture(params string[] parts) => Path.Combine(
            new[] { TestContext.CurrentContext.TestDirectory, "FileReadingTests", "ExternalFileTypes" }.Concat(parts).ToArray());

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
