using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using MzLibUtil;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Tests for the SDRF-Proteomics reader/writer.
    ///
    /// Each fixture here pins a property measured across the 1,236-file curated corpus at
    /// bigbio/sdrf-annotated-datasets rather than a property of the specification, because the two
    /// disagree in ways that would break a reader written from the spec alone:
    ///
    ///   PXD000070 -- typical file; repeats comment[modification parameters] 8 times, and uses a
    ///                PRIDE CV term (not PSI-MS) for the acquisition method.
    ///   PXD026824 -- contains literal double-quote characters inside cell values, which a
    ///                CSV parser with default quoting would strip.
    ///   PXD059974 -- ragged: a 46-column header with 17 of its 23 rows carrying only 42 cells.
    ///
    /// The load-bearing test is <see cref="RoundTrip_IsByteIdentical"/>. A reader that normalizes
    /// anything -- line endings, short rows, quoting, header case -- fails it, and normalizing a
    /// curated annotation silently is the failure mode this whole component exists to avoid.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrf
    {
        private static string FixtureDir => Path.Combine(
            TestContext.CurrentContext.TestDirectory, "FileReadingTests", "ExternalFileTypes");

        private static string Fixture(string name) => Path.Combine(FixtureDir, name);

        private static IEnumerable<string> AllFixtures()
        {
            yield return "PXD000070.sdrf.tsv";
            yield return "PXD026824.sdrf.tsv";
            yield return "PXD059974.sdrf.tsv";
        }

        [Test]
        [TestCaseSource(nameof(AllFixtures))]
        public void RoundTrip_IsByteIdentical(string fixtureName)
        {
            string source = Fixture(fixtureName);
            var document = new SdrfDocument(source);

            string outputPath = Path.Combine(
                TestContext.CurrentContext.TestDirectory, $"roundtrip_{Guid.NewGuid():N}.sdrf.tsv");
            try
            {
                document.WriteResults(outputPath);

                var expected = File.ReadAllBytes(source);
                var actual = File.ReadAllBytes(outputPath);

                Assert.That(actual, Is.EqualTo(expected),
                    $"round trip of {fixtureName} was not byte-identical");
            }
            finally
            {
                if (File.Exists(outputPath))
                    File.Delete(outputPath);
            }
        }

        [Test]
        public void Read_ParsesHeaderAndRows()
        {
            var document = new SdrfDocument(Fixture("PXD000070.sdrf.tsv"));

            Assert.That(document.Header[0], Is.EqualTo("source name"));
            Assert.That(document.Header.Contains("characteristics[organism]"), Is.True);
            Assert.That(document.Results, Is.Not.Empty);
            Assert.That(document.Results.Count, Is.EqualTo(6));
        }

        [Test]
        public void RepeatedColumns_AreAllPreserved()
        {
            // The house [Name(...)] pattern cannot express this: the column name repeats, and each
            // occurrence carries a different modification.
            var document = new SdrfDocument(Fixture("PXD000070.sdrf.tsv"));

            var indexes = document.Header.IndexesOf("comment[modification parameters]");
            Assert.That(indexes.Count, Is.EqualTo(8));

            var mods = document.Results[0].All("comment[modification parameters]");
            Assert.That(mods.Count, Is.EqualTo(8));
            Assert.That(mods[0], Does.StartWith("NT=Carbamidomethyl"));
            Assert.That(mods.Any(m => m.Contains("UNIMOD:35")), Is.True);

            // The single-value indexer returns the FIRST occurrence, not a join of them.
            Assert.That(document.Results[0]["comment[modification parameters]"],
                Is.EqualTo(mods[0]));
        }

        [Test]
        public void LiteralQuotes_AreNotTreatedAsCsvQuoting()
        {
            // SDRF defines no quoting or escaping. These double quotes are data.
            var document = new SdrfDocument(Fixture("PXD026824.sdrf.tsv"));

            bool anyQuoted = document.Results
                .SelectMany(r => r.Cells)
                .Any(c => c.StartsWith("\"") && c.EndsWith("\"") && c.Length > 1);

            Assert.That(anyQuoted, Is.True,
                "expected at least one cell whose literal double quotes survived parsing");
        }

        [Test]
        public void ShortRows_ArePreservedNotPadded()
        {
            // PXD059974 is genuinely malformed. The reader preserves it so the file round-trips;
            // reporting it is the validator's job (PR-3).
            var document = new SdrfDocument(Fixture("PXD059974.sdrf.tsv"));

            Assert.That(document.Header.Count, Is.EqualTo(46));

            var widths = document.Results.Select(r => r.Cells.Count).Distinct().OrderBy(x => x).ToList();
            Assert.That(widths, Is.EqualTo(new[] { 42, 46 }));

            // Reaching past the end of a short row yields null rather than throwing.
            var shortRow = document.Results.First(r => r.Cells.Count == 42);
            Assert.That(shortRow["comment[isolation window width]"], Is.Null);
        }

        [Test]
        public void HeaderCase_IsNotNormalized()
        {
            // The spec says lowercase and case-sensitive, but the corpus contains mixed-case
            // headers. Rewriting a user's column name is not the reader's call.
            var document = new SdrfDocument(Fixture("PXD059974.sdrf.tsv"));

            Assert.That(document.Header.Any(h => h.Any(char.IsUpper)), Is.True);
            Assert.That(document.Header.Contains("comment[MS max charge]"), Is.True);
            Assert.That(document.Header.Contains("comment[ms max charge]"), Is.False);
        }

        [Test]
        public void Write_IsDeterministic()
        {
            // Same document, written twice, must produce identical bytes. This is what makes the
            // accumulated corpus joinable: a document that serialized differently between runs
            // would show as a spurious change on every re-search.
            var document = new SdrfDocument(Fixture("PXD000070.sdrf.tsv"));

            string first = Path.Combine(TestContext.CurrentContext.TestDirectory, $"det1_{Guid.NewGuid():N}.sdrf.tsv");
            string second = Path.Combine(TestContext.CurrentContext.TestDirectory, $"det2_{Guid.NewGuid():N}.sdrf.tsv");
            try
            {
                document.WriteResults(first);
                document.WriteResults(second);
                Assert.That(File.ReadAllBytes(second), Is.EqualTo(File.ReadAllBytes(first)));
            }
            finally
            {
                foreach (var p in new[] { first, second })
                    if (File.Exists(p)) File.Delete(p);
            }
        }

        [Test]
        public void InMemoryDocument_Writes()
        {
            var header = new SdrfHeader(new[] { "source name", "assay name", "comment[data file]" });
            var rows = new List<SdrfRow>
            {
                new SdrfRow(header, new[] { "Sample 1", "run 1", "a.raw" }),
                new SdrfRow(header, new[] { "Sample 2", "run 2", "b.raw" })
            };
            var document = new SdrfDocument(header, rows);

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, $"mem_{Guid.NewGuid():N}.sdrf.tsv");
            try
            {
                document.WriteResults(outputPath);

                var text = File.ReadAllText(outputPath);
                Assert.That(text, Is.EqualTo(
                    "source name\tassay name\tcomment[data file]\r\n" +
                    "Sample 1\trun 1\ta.raw\r\n" +
                    "Sample 2\trun 2\tb.raw\r\n"));

                // No byte-order mark -- not one corpus file has one.
                var bytes = File.ReadAllBytes(outputPath);
                Assert.That(bytes.Take(3).ToArray(), Is.Not.EqualTo(new byte[] { 0xEF, 0xBB, 0xBF }));
            }
            finally
            {
                if (File.Exists(outputPath)) File.Delete(outputPath);
            }
        }

        [Test]
        public void Write_AppendsExtensionWhenMissing()
        {
            var header = new SdrfHeader(new[] { "source name" });
            var document = new SdrfDocument(header, new[] { new SdrfRow(header, new[] { "Sample 1" }) });

            string stem = Path.Combine(TestContext.CurrentContext.TestDirectory, $"stem_{Guid.NewGuid():N}");
            string expected = stem + ".sdrf.tsv";
            try
            {
                document.WriteResults(stem);
                Assert.That(File.Exists(expected), Is.True);
            }
            finally
            {
                if (File.Exists(expected)) File.Delete(expected);
            }
        }

        /// <summary>
        /// The line ending must come from the FILE, not from the platform or from a survey.
        ///
        /// Regression for a real defect: the writer pinned CRLF because a corpus survey reported
        /// "all 1,236 files are CRLF". That survey measured a working tree on a machine with
        /// core.autocrlf=true; upstream the corpus is LF. The pinned writer therefore flipped every
        /// real file's endings, and the byte-identical round-trip test only passed because git was
        /// silently supplying the CRs it was asserting on.
        /// </summary>
        [Test]
        [TestCase("\n")]
        [TestCase("\r\n")]
        [TestCase("\r")]
        public void RoundTrip_PreservesLineEnding(string lineEnding)
        {
            string source = Path.Combine(TestContext.CurrentContext.TestDirectory, $"eol_{Guid.NewGuid():N}.sdrf.tsv");
            string output = Path.Combine(TestContext.CurrentContext.TestDirectory, $"eol_out_{Guid.NewGuid():N}.sdrf.tsv");
            string text = string.Join(lineEnding, "source name\tassay name", "S1\trun 1", "S2\trun 2") + lineEnding;
            try
            {
                File.WriteAllText(source, text, new UTF8Encoding(false));

                var document = new SdrfDocument(source);
                Assert.That(document.Results.Count, Is.EqualTo(2));
                document.WriteResults(output);

                Assert.That(File.ReadAllBytes(output), Is.EqualTo(File.ReadAllBytes(source)));
            }
            finally
            {
                foreach (var p in new[] { source, output })
                    if (File.Exists(p)) File.Delete(p);
            }
        }

        [Test]
        public void RoundTrip_PreservesAbsentTrailingNewline()
        {
            string source = Path.Combine(TestContext.CurrentContext.TestDirectory, $"noeol_{Guid.NewGuid():N}.sdrf.tsv");
            string output = Path.Combine(TestContext.CurrentContext.TestDirectory, $"noeol_out_{Guid.NewGuid():N}.sdrf.tsv");
            try
            {
                File.WriteAllText(source, "source name\tassay name\nS1\trun 1", new UTF8Encoding(false));

                new SdrfDocument(source).WriteResults(output);

                Assert.That(File.ReadAllBytes(output), Is.EqualTo(File.ReadAllBytes(source)));
            }
            finally
            {
                foreach (var p in new[] { source, output })
                    if (File.Exists(p)) File.Delete(p);
            }
        }

        [Test]
        public void LoneCarriageReturn_IsALineBreak()
        {
            // Treating only LF as a break collapsed an old-Mac file into one line: a header of one
            // oddly-named column and zero rows, which the validator then reported as "NoRows" -- a
            // confident, completely wrong diagnosis of a file full of data.
            string source = Path.Combine(TestContext.CurrentContext.TestDirectory, $"cr_{Guid.NewGuid():N}.sdrf.tsv");
            try
            {
                File.WriteAllText(source, "source name\tassay name\rS1\trun 1\r", new UTF8Encoding(false));

                var document = new SdrfDocument(source);
                Assert.That(document.Header.Count, Is.EqualTo(2));
                Assert.That(document.Results.Count, Is.EqualTo(1));
                Assert.That(document.Results[0]["source name"], Is.EqualTo("S1"));
            }
            finally
            {
                if (File.Exists(source)) File.Delete(source);
            }
        }

        [Test]
        public void Write_RejectsHeaderContainingTab()
        {
            // The more dangerous of the two omissions: authored column names are built from run
            // metadata, and a tab in one emits N+1 columns, leaving every row one narrower than the
            // header and the whole document silently misaligned.
            var header = new SdrfHeader(new[] { "source name", "characteristics[organism\tpart]" });
            var document = new SdrfDocument(header, new[] { new SdrfRow(header, new[] { "S1", "liver" }) });

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, $"badhdr_{Guid.NewGuid():N}.sdrf.tsv");
            try
            {
                Assert.That(() => document.WriteResults(outputPath), Throws.TypeOf<MzLibException>());
            }
            finally
            {
                if (File.Exists(outputPath)) File.Delete(outputPath);
            }
        }

        [Test]
        public void Write_RejectingABadCell_DoesNotTouchTheOutputFile()
        {
            // File.Create truncates. Validating inside the write loop meant the guard destroyed the
            // file it was refusing to write -- and writing in place over a document deleted the
            // original, then threw "SDRF file is empty" when the lazy load re-read the wreckage.
            var header = new SdrfHeader(new[] { "source name", "assay name" });
            var document = new SdrfDocument(header, new[]
            {
                new SdrfRow(header, new[] { "S1", "run 1" }),
                new SdrfRow(header, new[] { "S2\tbad", "run 2" })
            });

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, $"keep_{Guid.NewGuid():N}.sdrf.tsv");
            const string existing = "PRECIOUS EXISTING CONTENT";
            try
            {
                File.WriteAllText(outputPath, existing);

                Assert.That(() => document.WriteResults(outputPath), Throws.TypeOf<MzLibException>());
                Assert.That(File.ReadAllText(outputPath), Is.EqualTo(existing),
                    "the guard truncated the very file it refused to write");
            }
            finally
            {
                if (File.Exists(outputPath)) File.Delete(outputPath);
            }
        }

        /// <summary>
        /// Regression: the "validate before opening the stream" fix did not cover a document with
        /// ZERO rows, which is the case its own comment described.
        ///
        /// ResultFile.Results reloads whenever the list is EMPTY, not merely unloaded, so a
        /// header-only document re-reads from disk on every access. WriteResults touched Results
        /// again after File.Create had truncated the file, so an in-place write destroyed the
        /// original and then threw "SDRF file is empty" about the wreckage.
        /// </summary>
        [Test]
        public void WriteInPlace_OnAHeaderOnlyDocument_DoesNotDestroyIt()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, $"hdr_{Guid.NewGuid():N}.sdrf.tsv");
            const string original = "source name\tassay name\ttechnology type\r\n";
            try
            {
                File.WriteAllText(path, original, new UTF8Encoding(false));

                var document = new SdrfDocument(path);
                Assert.That(document.Header.Count, Is.EqualTo(3));
                Assert.That(document.Results, Is.Empty);

                document.WriteResults(path);

                Assert.That(File.ReadAllText(path), Is.EqualTo(original),
                    "an in-place write of a header-only document must round-trip, not truncate");
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [Test]
        public void WriteInPlace_WithRows_RoundTrips()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, $"inplace_{Guid.NewGuid():N}.sdrf.tsv");
            const string original = "source name\tassay name\r\nS1\trun 1\r\nS2\trun 2\r\n";
            try
            {
                File.WriteAllText(path, original, new UTF8Encoding(false));
                new SdrfDocument(path).WriteResults(path);
                Assert.That(File.ReadAllText(path), Is.EqualTo(original));
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [Test]
        public void Write_RejectsCellContainingTab()
        {
            // There is no escape mechanism, so a tab inside a cell would silently shift every
            // following column. Fail loudly rather than emit a file that reads back wrong.
            var header = new SdrfHeader(new[] { "source name", "assay name" });
            var document = new SdrfDocument(header,
                new[] { new SdrfRow(header, new[] { "Sample\t1", "run 1" }) });

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, $"bad_{Guid.NewGuid():N}.sdrf.tsv");
            try
            {
                Assert.That(() => document.WriteResults(outputPath), Throws.TypeOf<MzLibException>());
            }
            finally
            {
                if (File.Exists(outputPath)) File.Delete(outputPath);
            }
        }

        [Test]
        public void FileReader_ResolvesSdrfByExtension()
        {
            Assert.That(Fixture("PXD000070.sdrf.tsv").ParseFileType(), Is.EqualTo(SupportedFileType.Sdrf));
            Assert.That(SupportedFileType.Sdrf.GetResultFileType(), Is.EqualTo(typeof(SdrfDocument)));

            var loaded = FileReader.ReadFile<SdrfDocument>(Fixture("PXD000070.sdrf.tsv"));
            Assert.That(loaded.Results, Is.Not.Empty);
        }

        /// <summary>
        /// Round-trips EVERY file in the curated corpus. [Explicit] so it never runs in CI: it needs
        /// a local clone of bigbio/sdrf-annotated-datasets (~1,236 files), which is not committed.
        ///
        ///     git clone https://github.com/bigbio/sdrf-annotated-datasets
        ///     set MZLIB_SDRF_CORPUS=&lt;clone&gt;\datasets
        ///     dotnet test --filter "FullyQualifiedName~RoundTrip_EntireCorpus" -- NUnit.Explicit=true
        ///
        /// This is the regression that actually guards the format layer. The three committed
        /// fixtures pin the failure modes already found; this catches the ones that have not been.
        /// </summary>
        [Test]
        [Explicit("Requires a local clone of bigbio/sdrf-annotated-datasets; set MZLIB_SDRF_CORPUS.")]
        public void RoundTrip_EntireCorpus()
        {
            string? corpus = Environment.GetEnvironmentVariable("MZLIB_SDRF_CORPUS");
            if (string.IsNullOrWhiteSpace(corpus) || !Directory.Exists(corpus))
                Assert.Ignore($"MZLIB_SDRF_CORPUS not set or not found: '{corpus}'");

            var files = Directory.GetFiles(corpus, "*.sdrf.tsv", SearchOption.AllDirectories);
            Assert.That(files, Is.Not.Empty, "corpus directory contained no .sdrf.tsv files");

            var failures = new List<string>();
            string scratch = Path.Combine(Path.GetTempPath(), $"sdrf_corpus_{Guid.NewGuid():N}");
            Directory.CreateDirectory(scratch);
            try
            {
                foreach (var file in files)
                {
                    string outputPath = Path.Combine(scratch, "rt.sdrf.tsv");
                    try
                    {
                        new SdrfDocument(file).WriteResults(outputPath);
                        if (!File.ReadAllBytes(outputPath).SequenceEqual(File.ReadAllBytes(file)))
                            failures.Add($"{Path.GetFileName(file)}: bytes differ");
                    }
                    catch (Exception e)
                    {
                        failures.Add($"{Path.GetFileName(file)}: {e.GetType().Name}: {e.Message}");
                    }
                }
            }
            finally
            {
                if (Directory.Exists(scratch)) Directory.Delete(scratch, true);
            }

            TestContext.Progress.WriteLine($"round-tripped {files.Length} files, {failures.Count} failures");
            foreach (var f in failures.Take(25))
                TestContext.Progress.WriteLine("  " + f);

            Assert.That(failures, Is.Empty,
                $"{failures.Count} of {files.Length} corpus files did not round-trip byte-identically");
        }

        [Test]
        public void EmptyFile_ThrowsMzLibException()
        {
            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, $"empty_{Guid.NewGuid():N}.sdrf.tsv");
            File.WriteAllText(outputPath, string.Empty);
            try
            {
                var document = new SdrfDocument(outputPath);
                Assert.That(() => document.LoadResults(), Throws.TypeOf<MzLibException>());
            }
            finally
            {
                if (File.Exists(outputPath)) File.Delete(outputPath);
            }
        }
    }
}
