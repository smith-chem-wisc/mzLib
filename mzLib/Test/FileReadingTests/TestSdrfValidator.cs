using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Tests for <see cref="SdrfValidator"/>.
    ///
    /// Two halves. The unit tests below build a minimal synthetic document and check that each rule
    /// fires exactly when it should. <see cref="CorpusCalibration"/> then runs the validator over
    /// the whole curated corpus and pins the observed rates, because a rule that condemns most of
    /// the community reference set is a wrong rule, not 1,236 wrong files.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrfValidator
    {
        private static string FixtureDir => Path.Combine(
            TestContext.CurrentContext.TestDirectory, "FileReadingTests", "ExternalFileTypes");

        /// <summary>
        /// A header carrying every Required column (the 13 present in all 1,236 curated files) plus
        /// the Recommended ones, so a document built on it is genuinely warning-free. Shared, because
        /// several tests are about ONE rule and should not fail merely for being skeletal.
        /// </summary>
        private static SdrfHeader CompleteHeader => new(new[]
        {
            "source name", "characteristics[organism]", "characteristics[organism part]",
            "characteristics[disease]", "characteristics[cell type]",
            "characteristics[biological replicate]",
            "assay name", "technology type",
            "comment[proteomics data acquisition method]", "comment[label]", "comment[instrument]",
            "comment[cleavage agent details]", "comment[modification parameters]",
            "comment[fraction identifier]", "comment[technical replicate]", "comment[data file]",
            "factor value[disease]"
        });

        /// <summary>One valid row over <see cref="CompleteHeader"/>, plus any extra trailing cells.</summary>
        private static SdrfRow CompleteRow(SdrfHeader header, string sample, params string[] extras)
        {
            var cells = new List<string>
            {
                sample, "homo sapiens", "liver", "normal", "epithelial cell", "1",
                "run " + sample, "proteomic profiling by mass spectrometry",
                "NT=Data-dependent acquisition;AC=PRIDE:0000627", "label free sample",
                "NT=Q Exactive;AC=MS:1001911", "NT=Trypsin;AC=MS:1001251",
                "NT=Oxidation;AC=UNIMOD:35;TA=M;MT=Variable", "1", "1", sample + ".raw",
                "normal"
            };
            cells.AddRange(extras);
            return new SdrfRow(header, cells);
        }

        /// <summary>A minimal document that should validate with no errors.</summary>
        private static SdrfDocument MinimalValid()
        {
            var header = CompleteHeader;
            return new SdrfDocument(header,
                new[] { CompleteRow(header, "Sample 1"), CompleteRow(header, "Sample 2") });
        }

        [Test]
        public void MinimalDocument_IsValidWithNoWarnings()
        {
            var result = SdrfValidator.Validate(MinimalValid());
            Assert.That(result.IsValid, Is.True, result.ToString());
            Assert.That(result.Messages, Is.Empty,
                "unexpected: " + string.Join(" | ", result.Messages.Select(m => m.ToString())));
        }

        [Test]
        public void EmptyHeader_IsOneErrorRatherThanOnePerMissingColumn()
        {
            var result = SdrfValidator.Validate(
                new SdrfDocument(new SdrfHeader(Array.Empty<string>()), Array.Empty<SdrfRow>()));

            Assert.That(result.Errors.Any(e => e.Rule == "EmptyHeader"), Is.True, result.ToString());

            // The early return matters: without it a document with no columns would report one
            // RequiredColumn error for each of the 13 universal columns, burying the real problem.
            Assert.That(result.Errors.Any(e => e.Rule == "RequiredColumn"), Is.False,
                "EmptyHeader should short-circuit the required-column sweep");
        }

        [Test]
        public void MisCasedRequiredColumn_SaysSoInsteadOfJustSayingMissing()
        {
            // The severity is right -- the spec makes lowercase a MUST and names
            // "Characteristics[organism]" as an explicit failure, so a consumer joining on column
            // names really does not have this column. But "missing comment[data file]" on a file
            // that visibly contains "Comment[data file]" sends a curator hunting for the wrong bug.
            var names = CompleteHeader.ToList();
            names[names.IndexOf("comment[data file]")] = "Comment[data file]";
            var header = new SdrfHeader(names);

            var result = SdrfValidator.Validate(
                new SdrfDocument(header, new[] { CompleteRow(header, "Sample 1") }));

            var missing = result.Errors.Single(e => e.Rule == "RequiredColumn");
            Assert.That(missing.ColumnName, Is.EqualTo("comment[data file]"));
            Assert.That(missing.Message, Does.Contain("Comment[data file]"),
                "the message must name the column that IS present");
            Assert.That(missing.Message, Does.Contain("differs only in casing"));

            // ColumnNameCase still reports it independently; the hint does not replace that.
            Assert.That(result.Warnings.Any(w => w.Rule == "ColumnNameCase"), Is.True);
        }

        [Test]
        public void GenuinelyMissingRequiredColumn_GetsNoCasingHint()
        {
            var names = CompleteHeader.Where(n => n != "comment[data file]").ToList();
            var header = new SdrfHeader(names);

            var result = SdrfValidator.Validate(
                new SdrfDocument(header, new[] { new SdrfRow(header, names) }));

            var missing = result.Errors.Single(
                e => e.Rule == "RequiredColumn" && e.ColumnName == "comment[data file]");
            Assert.That(missing.Message, Does.Not.Contain("casing"));
        }

        [Test]
        public void HeaderWithNoRows_WarnsButIsStillValid()
        {
            var result = SdrfValidator.Validate(
                new SdrfDocument(CompleteHeader, Array.Empty<SdrfRow>()));

            Assert.That(result.Warnings.Any(w => w.Rule == "NoRows"), Is.True, result.ToString());
            Assert.That(result.IsValid, Is.True, "a header with no rows is incomplete, not invalid");
        }

        [Test]
        public void BlankColumnName_WarnsWithoutInvalidatingTheDocument()
        {
            // Almost always a trailing tab on the header line; one corpus file carries 23 of them.
            // It must not be an Error, or those files would all fail for a whitespace artefact.
            var header = new SdrfHeader(CompleteHeader.Concat(new[] { "" }));
            var result = SdrfValidator.Validate(
                new SdrfDocument(header, new[] { CompleteRow(header, "Sample 1", "not applicable") }));

            Assert.That(result.Warnings.Any(w => w.Rule == "EmptyColumnName"), Is.True, result.ToString());
            Assert.That(result.IsValid, Is.True);
        }

        [Test]
        public void EmptyCellWarns_AndCellsPastTheHeaderAreNamedByPosition()
        {
            // A row wider than its header is an Error, because every cell past the last header
            // entry is being read under no name at all -- reported positionally so the message
            // still says where it is.
            var header = CompleteHeader;
            var result = SdrfValidator.Validate(
                new SdrfDocument(header, new[] { CompleteRow(header, "Sample 1", "", "extra") }));

            Assert.That(result.Warnings.Any(w => w.Rule == "EmptyCell"), Is.True, result.ToString());
            Assert.That(result.Errors.Any(e => e.Rule == "RowWidth"), Is.True,
                "a row with more cells than the header is a width error");
            Assert.That(result.Messages.Any(m => m.ColumnName == "(column 17)"), Is.True,
                "cells beyond the header are named by position: " + result.ToString());
        }

        [Test]
        public void MissingRequiredColumn_IsError()
        {
            var header = new SdrfHeader(new[] { "assay name", "technology type" });
            var document = new SdrfDocument(header, new[] { new SdrfRow(header, new[] { "run 1", "x" }) });

            var result = SdrfValidator.Validate(document);
            Assert.That(result.IsValid, Is.False);
            Assert.That(result.Errors.Any(e => e.Rule == "RequiredColumn" && e.ColumnName == "source name"));
        }

        [Test]
        public void MissingRecommendedColumn_IsWarningNotError()
        {
            // Every REQUIRED column present; only a Recommended one (cell type) removed.
            var names = CompleteHeader.Where(n => n != "characteristics[cell type]").ToArray();
            var header = new SdrfHeader(names);
            var cells = CompleteRow(CompleteHeader, "S1").Cells.ToList();
            cells.RemoveAt(4);
            var document = new SdrfDocument(header, new[] { new SdrfRow(header, cells) });

            var result = SdrfValidator.Validate(document);
            Assert.That(result.IsValid, Is.True, "a missing recommended column must not invalidate");
            Assert.That(result.Warnings.Any(w => w.Rule == "RecommendedColumn"));
        }

        [Test]
        public void RaggedRow_IsError()
        {
            var header = new SdrfHeader(new[] { "source name", "assay name", "technology type" });
            var document = new SdrfDocument(header, new[] { new SdrfRow(header, new[] { "s", "a" }) });

            var result = SdrfValidator.Validate(document);
            var finding = result.Errors.FirstOrDefault(e => e.Rule == "RowWidth");
            Assert.That(finding, Is.Not.Null);
            Assert.That(finding!.RowIndex, Is.EqualTo(0));
            Assert.That(finding.LineNumber, Is.EqualTo(2), "line number is row index + 2 (header is line 1)");
        }

        [Test]
        public void DuplicateRowKey_IsError()
        {
            var header = new SdrfHeader(new[] { "source name", "assay name", "technology type", "comment[label]" });
            var document = new SdrfDocument(header, new[]
            {
                new SdrfRow(header, new[] { "Sample 1", "run 1", "t", "label free sample" }),
                new SdrfRow(header, new[] { "Sample 1", "run 1", "t", "label free sample" })
            });

            var result = SdrfValidator.Validate(document);
            Assert.That(result.Errors.Any(e => e.Rule == "RowKeyUniqueness"));
        }

        [Test]
        public void SameSourceDifferentLabel_IsNotDuplicate()
        {
            // The whole point of including comment[label] in the key: one sample, two TMT channels.
            var header = new SdrfHeader(new[] { "source name", "assay name", "technology type", "comment[label]" });
            var document = new SdrfDocument(header, new[]
            {
                new SdrfRow(header, new[] { "Sample 1", "run 1", "t", "TMT126" }),
                new SdrfRow(header, new[] { "Sample 1", "run 1", "t", "TMT127N" })
            });

            var result = SdrfValidator.Validate(document);
            Assert.That(result.Errors.Any(e => e.Rule == "RowKeyUniqueness"), Is.False);
        }

        [Test]
        public void UppercaseColumnName_IsWarning()
        {
            var names = CompleteHeader.ToList();
            names.Add("comment[MS min charge]");
            var header = new SdrfHeader(names);
            var document = new SdrfDocument(header, new[] { CompleteRow(header, "S1", "2") });

            var result = SdrfValidator.Validate(document);
            Assert.That(result.IsValid, Is.True);
            Assert.That(result.Warnings.Any(w => w.Rule == "ColumnNameCase"));
        }

        [Test]
        public void SpaceBeforeBracket_IsWarning()
        {
            var header = new SdrfHeader(new[] { "source name", "assay name", "technology type", "characteristics [organism]" });
            var document = new SdrfDocument(header, new[] { new SdrfRow(header, new[] { "s", "a", "t", "homo sapiens" }) });

            var result = SdrfValidator.Validate(document);
            Assert.That(result.Warnings.Any(w => w.Rule == "MalformedColumnName"));
        }

        [Test]
        public void CasedReservedWord_IsWarning()
        {
            var header = new SdrfHeader(new[] { "source name", "assay name", "technology type" });
            var document = new SdrfDocument(header, new[] { new SdrfRow(header, new[] { "s", "a", "Not Available" }) });

            var result = SdrfValidator.Validate(document);
            Assert.That(result.Warnings.Any(w => w.Rule == "ReservedWordCase"));
        }

        [Test]
        public void NonIntegerReplicate_IsWarning()
        {
            var header = new SdrfHeader(new[]
                { "source name", "assay name", "technology type", "comment[technical replicate]" });
            var document = new SdrfDocument(header, new[] { new SdrfRow(header, new[] { "s", "a", "t", "one" }) });

            var result = SdrfValidator.Validate(document);
            Assert.That(result.Warnings.Any(w => w.Rule == "NonIntegerValue"));
        }

        [Test]
        public void ReservedWordInIntegerColumn_IsAccepted()
        {
            var header = new SdrfHeader(new[]
                { "source name", "assay name", "technology type", "characteristics[biological replicate]" });
            var document = new SdrfDocument(header, new[] { new SdrfRow(header, new[] { "s", "a", "t", "pooled" }) });

            var result = SdrfValidator.Validate(document);
            Assert.That(result.Warnings.Any(w => w.Rule == "NonIntegerValue"), Is.False,
                "'pooled' is the spec's own way to say the replicate is not a number");
        }

        [Test]
        public void FactorValueBeforeComment_IsOrderingWarning()
        {
            var header = new SdrfHeader(new[]
                { "source name", "factor value[disease]", "assay name", "technology type", "comment[data file]" });
            var document = new SdrfDocument(header,
                new[] { new SdrfRow(header, new[] { "s", "normal", "a", "t", "x.raw" }) });

            var result = SdrfValidator.Validate(document);
            Assert.That(result.Warnings.Any(w => w.Rule == "ColumnOrdering"));
        }

        /// <summary>
        /// Regression: the original rule compared only two pairs of extremes, fired ZERO times
        /// across all 1,236 curated files, and could not detect the violation its own summary
        /// described. Both cases below produced no message.
        /// </summary>
        [Test]
        public void SampleMetadataAfterDataFileMetadata_IsOrderingWarning()
        {
            var header = new SdrfHeader(new[]
                { "source name", "assay name", "technology type", "characteristics[organism]", "comment[data file]" });
            var document = new SdrfDocument(header,
                new[] { new SdrfRow(header, new[] { "s", "a", "t", "homo sapiens", "x.raw" }) });

            var result = SdrfValidator.Validate(document);
            Assert.That(result.Warnings.Any(w => w.Rule == "ColumnOrdering"), Is.True,
                "characteristics[...] after assay name is out of block order");
        }

        [Test]
        public void SourceNameNotFirst_IsOrderingWarning()
        {
            var header = new SdrfHeader(new[]
                { "characteristics[organism]", "source name", "assay name", "technology type", "comment[data file]" });
            var document = new SdrfDocument(header,
                new[] { new SdrfRow(header, new[] { "homo sapiens", "s", "a", "t", "x.raw" }) });

            var result = SdrfValidator.Validate(document);
            Assert.That(result.Warnings.Any(w => w.Rule == "ColumnOrdering"), Is.True,
                "nothing previously checked the position of source name at all");
        }

        [Test]
        public void MixedCaseFactorValue_DoesNotDisarmOrderingOrFalselyTrigger()
        {
            // Ordinal matching meant "Factor Value[...]" (which the corpus contains) silently
            // disarmed the rule, while a mixed-case pair elsewhere produced a false positive.
            var header = new SdrfHeader(new[]
                { "source name", "assay name", "technology type", "comment[data file]", "Factor Value[organism part]" });
            var document = new SdrfDocument(header,
                new[] { new SdrfRow(header, new[] { "s", "a", "t", "x.raw", "liver" }) });

            var result = SdrfValidator.Validate(document);
            Assert.That(result.Warnings.Any(w => w.Rule == "ColumnOrdering"), Is.False,
                "this header IS correctly ordered; only its casing is wrong");
            Assert.That(result.Warnings.Any(w => w.Rule == "ColumnNameCase"), Is.True,
                "the casing is reported separately, which is where it belongs");
        }

        [Test]
        public void ShortRows_DoNotManufactureDuplicateKeyErrors()
        {
            // Two rows too short to reach assay name both keyed as "" and were reported as
            // duplicates of one another -- a second, wrong diagnosis on top of the RowWidth error
            // they already had.
            var header = new SdrfHeader(new[] { "source name", "assay name", "technology type" });
            var document = new SdrfDocument(header, new[]
            {
                new SdrfRow(header, new[] { "S1" }),
                new SdrfRow(header, new[] { "S2" })
            });

            var result = SdrfValidator.Validate(document);

            Assert.That(result.Errors.Count(e => e.Rule == "RowWidth"), Is.EqualTo(2));
            Assert.That(result.Errors.Any(e => e.Rule == "RowKeyUniqueness"), Is.False,
                "a row that cannot be keyed must not be reported as a duplicate");
        }

        [Test]
        public void RaggedCorpusFile_IsReportedAsError()
        {
            // PXD059974 is the known-malformed corpus file. The reader loads it happily; the
            // validator is what says so.
            var document = new SdrfDocument(Path.Combine(FixtureDir, "PXD059974.sdrf.tsv"));
            var result = SdrfValidator.Validate(document);

            Assert.That(result.IsValid, Is.False);
            Assert.That(result.Errors.Count(e => e.Rule == "RowWidth"), Is.EqualTo(17),
                "17 of its 23 rows are short");
        }

        [Test]
        public void TypicalCorpusFile_HasNoErrors()
        {
            var document = new SdrfDocument(Path.Combine(FixtureDir, "PXD000070.sdrf.tsv"));
            var result = SdrfValidator.Validate(document);

            Assert.That(result.IsValid, Is.True,
                "PXD000070 is a well-formed curated file: " +
                string.Join(" | ", result.Errors.Select(e => e.ToString())));
        }

        /// <summary>
        /// Runs the validator over the entire curated corpus and reports how often each rule fires.
        /// [Explicit]; needs MZLIB_SDRF_CORPUS. See TestSdrf.RoundTrip_EntireCorpus.
        ///
        /// This is rule calibration, not a pass/fail check on the corpus. The assertion is that the
        /// great majority of curated files carry no ERROR -- if that stops being true, a rule has
        /// been miscategorised and the validator has become useless against real data.
        /// </summary>
        [Test]
        [Explicit("Requires a local clone of bigbio/sdrf-annotated-datasets; set MZLIB_SDRF_CORPUS.")]
        public void CorpusCalibration()
        {
            string? corpus = Environment.GetEnvironmentVariable("MZLIB_SDRF_CORPUS");
            if (string.IsNullOrWhiteSpace(corpus) || !Directory.Exists(corpus))
                Assert.Ignore($"MZLIB_SDRF_CORPUS not set or not found: '{corpus}'");

            var files = Directory.GetFiles(corpus, "*.sdrf.tsv", SearchOption.AllDirectories);

            var filesByRule = new Dictionary<string, int>(StringComparer.Ordinal);
            var severityByRule = new Dictionary<string, SdrfValidationSeverity>(StringComparer.Ordinal);
            int filesWithErrors = 0, filesClean = 0;

            foreach (var file in files)
            {
                var result = SdrfValidator.Validate(new SdrfDocument(file));
                if (!result.IsValid) filesWithErrors++;
                if (result.Messages.Count == 0) filesClean++;

                foreach (var rule in result.Messages.Select(m => m.Rule).Distinct(StringComparer.Ordinal))
                {
                    filesByRule.TryGetValue(rule, out int n);
                    filesByRule[rule] = n + 1;
                    severityByRule[rule] = result.Messages.First(m => m.Rule == rule).Severity;
                }
            }

            TestContext.Progress.WriteLine($"corpus files          : {files.Length}");
            TestContext.Progress.WriteLine($"entirely clean        : {filesClean}");
            TestContext.Progress.WriteLine($"with >=1 ERROR        : {filesWithErrors} " +
                                           $"({100.0 * filesWithErrors / files.Length:F1}%)");
            TestContext.Progress.WriteLine("rule                     severity  files   %");
            foreach (var kv in filesByRule.OrderByDescending(k => k.Value))
            {
                TestContext.Progress.WriteLine(
                    $"  {kv.Key,-22} {severityByRule[kv.Key],-9} {kv.Value,5}  {100.0 * kv.Value / files.Length,5:F1}");
            }

            // Observed 2026-08-07 against corpus @ 4f823dcd, after RequiredColumns was corrected
            // from 3 entries to the 13 that are genuinely universal: 1,236 files, 1,089 entirely
            // clean (88.1%), and STILL exactly 1 with an error (PXD059974, genuinely ragged).
            //
            // That the error count did not move is the whole point. Ten columns were promoted from
            // Recommended to Required and not one additional curated file became invalid, which is
            // what "require only what is universal" is supposed to guarantee. The clean rate fell
            // because RecommendedColumn now fires on 69 files for the three columns that really are
            // sometimes absent (disease 4, cell type 15, modification parameters 60).
            //
            // The earlier 3-column list came from reading the survey's OCCURRENCE counts as FILE
            // counts.
            //
            // The ordering rule is the reason to keep watching this number. It previously fired on
            // ZERO files -- dead, and reading as coverage. Repairing it naively then made it fire on
            // 31, of which 28 were wrong, because ranking unrecognised columns forced legitimate
            // sample columns like "material type" after the comment block.
            //
            // The bounds are loose enough to survive corpus growth but tight enough that promoting
            // any warning to an error, or adding an over-eager rule, fails here rather than quietly
            // condemning the reference set.
            double errorRate = 100.0 * filesWithErrors / files.Length;
            double cleanRate = 100.0 * filesClean / files.Length;

            Assert.That(errorRate, Is.LessThan(1.0),
                $"{errorRate:F1}% of the curated corpus is reported invalid; a rule is miscategorised " +
                "as Error. Only genuinely unreadable files should error.");
            Assert.That(cleanRate, Is.GreaterThan(85.0),
                $"only {cleanRate:F1}% of the curated corpus is warning-free; a rule is over-eager.");
        }
    }
}
