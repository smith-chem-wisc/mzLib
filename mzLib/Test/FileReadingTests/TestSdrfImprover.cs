using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MzLibUtil;
using NUnit.Framework;
using Readers;
using UsefulProteomicsDatabases;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Tests for improving a DEPOSITED SDRF with a draft (sdrf D35): every stated value is kept, gaps are
    /// filled and marked, disagreements are reported and never applied.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrfImprover
    {
        private static readonly string[] Columns =
        {
            "source name", "characteristics[organism]", "characteristics[disease]", "characteristics[biological replicate]",
            "assay name", "comment[label]", "comment[data file]", "comment[technical replicate]", "factor value[disease]"
        };

        private static SdrfDocument Doc(string[] columns, params string[][] rows)
        {
            var header = new SdrfHeader(columns);
            return new SdrfDocument(header, rows.Select(r => new SdrfRow(header, r)));
        }

        private static SdrfDocument Deposited() => Doc(Columns,
            new[] { "p1", "homo sapiens", "not available", "1", "run 1", "label free sample", "NEG1.raw", "1", "normal" },
            new[] { "p1", "homo sapiens", "not available", "1", "run 1r", "label free sample", "NEG1rep.raw", "1", "normal" },
            new[] { "p2", "homo sapiens", "not available", "1", "run 2", "label free sample", "POS1.raw", "1", "COVID-19" });

        private static PrideProject Project() => new()
        {
            Accession = "PXD020394",
            ProjectDescription = "SARS-CoV-2 positive and negative swabs.",
            Organisms = { new CvParam("NEWT", "NEWT:9606", "Homo sapiens (human)", "") },
            Diseases = { new CvParam("DOID", "DOID:0080600", "Covid-19", "") },
            Instruments = { new CvParam("MS", "MS:1001911", "Q Exactive", "") },
        };

        private static readonly string[] Files = { "NEG1.raw", "NEG1rep.raw", "POS1.raw", "NEG2.raw", "NEG2rep.raw", "POS2.raw" };

        private static SdrfImprovement Improve() =>
            SdrfImprover.Improve(Deposited(), SdrfDrafter.Draft(Project(), Files));

        private static SdrfRow Row(SdrfDocument d, string file) => d.Results.Single(r => r["comment[data file]"] == file);

        [Test]
        public void AStatedValueIsNeverChangedAndADisagreementIsReported()
        {
            var i = Improve();

            // The depositor wrote technical replicate 1 for NEG1rep; the names say it is the second injection.
            Assert.That(Row(i.Document, "NEG1rep.raw")["comment[technical replicate]"], Is.EqualTo("1"));
            var d = i.Disagreements.Single(x => x.DataFile == "NEG1rep.raw" && x.Column == "comment[technical replicate]");
            Assert.That((d.Deposited, d.Drafted), Is.EqualTo(("1", "2")));
            Assert.That(d.Evidence, Is.Not.Empty);
        }

        [Test]
        public void ANotAvailableCellIsFilledAndMarked()
        {
            var i = Improve();

            var neg = Row(i.Document, "NEG1.raw");
            Assert.That(neg["characteristics[disease]"], Does.Contain("normal"));
            Assert.That(neg["comment[disease source]"], Is.EqualTo("inferred"));
            Assert.That(neg["comment[characteristics source]"], Is.EqualTo("deposited"));
        }

        [Test]
        public void AMissingColumnIsAddedInItsBlockAndMarked()
        {
            var i = Improve();

            var h = i.Document.Header.ToList();
            Assert.That(h, Does.Contain("comment[instrument]"));
            Assert.That(h.IndexOf("comment[instrument]"), Is.GreaterThan(h.IndexOf("assay name")));
            Assert.That(h.IndexOf("comment[instrument]"), Is.LessThan(h.IndexOf("factor value[disease]")));
            Assert.That(Row(i.Document, "POS1.raw")["comment[instrument]"], Does.Contain("MS:1001911"));
            Assert.That(Row(i.Document, "POS1.raw")["comment[instrument source]"], Is.EqualTo("pride project record"));
            Assert.That(i.Document.Results.All(r => r.Cells.Count == i.Document.Header.Count), "never ragged");
        }

        [Test]
        public void ARawFileTheSdrfDoesNotListGetsADraftedRow()
        {
            var i = Improve();

            var added = Row(i.Document, "NEG2rep.raw");
            Assert.That(added["comment[characteristics source]"], Is.EqualTo("inferred"));
            Assert.That(added["comment[label]"], Is.EqualTo("label free sample"), "a column constant across the deposit is carried");
            Assert.That(added["source name"], Is.Not.EqualTo("p1").And.Not.EqualTo("p2"), "a drafted sample never takes a deposited name");
            Assert.That(i.AddedRows, Is.EqualTo(3));
        }

        [Test]
        public void ADepositedFactorColumnIsKeptAndNoDraftedFactorIsAdded()
        {
            var i = Improve();

            Assert.That(i.Document.Header.Count(h => h.StartsWith("factor value[", StringComparison.Ordinal)), Is.EqualTo(1));
            Assert.That(Row(i.Document, "POS1.raw")["factor value[disease]"], Is.EqualTo("COVID-19"));
        }

        [Test]
        public void ADepositWithNoFactorGetsTheDraftedOne()
        {
            var noFactor = Columns.Take(Columns.Length - 1).ToArray();
            var dep = Doc(noFactor,
                new[] { "p1", "homo sapiens", "normal", "1", "run 1", "label free sample", "NEG1.raw", "1" },
                new[] { "p2", "homo sapiens", "covid", "1", "run 2", "label free sample", "POS1.raw", "1" });

            var i = SdrfImprover.Improve(dep, SdrfDrafter.Draft(Project(), Files));

            Assert.That(i.Document.Header, Does.Contain("factor value[condition]"));
            Assert.That(Row(i.Document, "POS1.raw")["factor value[condition]"], Is.EqualTo("POS").IgnoreCase);
        }

        [Test]
        public void AMultiplexedFileGetsOnlyItsPerFileFactsFilled()
        {
            var cols = new[] { "source name", "characteristics[organism]", "characteristics[biological replicate]", "assay name",
                "comment[label]", "comment[data file]", "comment[fraction identifier]" };
            var dep = Doc(cols,
                new[] { "a", "homo sapiens", "not available", "run 1", "TMT126", "NEG1.raw", "not available" },
                new[] { "b", "homo sapiens", "not available", "run 1", "TMT127N", "NEG1.raw", "not available" });

            var i = SdrfImprover.Improve(dep, SdrfDrafter.Draft(Project(), Files));

            var rows = i.Document.Results.Where(r => r["comment[data file]"] == "NEG1.raw").ToList();
            Assert.That(rows.Select(r => r["comment[instrument]"]), Is.All.Contains("MS:1001911"), "per-file: filled");
            Assert.That(rows.Select(r => r["characteristics[biological replicate]"]), Is.All.EqualTo("not available"),
                "per-channel: the drafter reads files, not channels, so it does not guess");
        }

        [Test]
        public void TheImprovedDocumentValidatesNoWorseThanTheDeposited()
        {
            var dep = Deposited();
            var before = SdrfValidator.Validate(dep).Errors.Count();

            var after = SdrfValidator.Validate(Improve().Document).Errors.ToList();

            Assert.That(after.Count, Is.LessThanOrEqualTo(before), string.Join("\n", after.Select(e => e.ToString())));
        }

        /// <summary>
        /// Found by improving the whole corpus: PRIDE's project-level instrument ("LTQ Orbitrap") disputed
        /// 24,000 per-file curated values ("Q Exactive"). A project summary is weaker evidence than a
        /// per-file statement: it may fill a gap, never dispute.
        /// </summary>
        [Test]
        public void AProjectRecordValueFillsAGapButNeverDisputesAStatedOne()
        {
            var cols = new[] { "source name", "assay name", "comment[instrument]", "comment[data file]" };
            var dep = Doc(cols, new[] { "p1", "run 1", "NT=Orbitrap Fusion;AC=MS:1002416", "NEG1.raw" });

            var i = SdrfImprover.Improve(dep, SdrfDrafter.Draft(Project(), Files));

            Assert.That(i.Disagreements.Any(x => x.Column == "comment[instrument]"), Is.False);
        }

        /// <summary>
        /// Found by improving the whole corpus: 44,000 "disagreements" were the draft's DEFAULT fraction 1
        /// against real fractions. A default is not a reading: it neither fills a gap nor disputes a value.
        /// </summary>
        [Test]
        public void ADraftDefaultNeitherFillsNorDisputes()
        {
            var cols = new[] { "source name", "assay name", "comment[data file]", "comment[fraction identifier]" };
            var dep = Doc(cols,
                new[] { "p1", "run 1", "NEG1.raw", "7" },
                new[] { "p1", "run 2", "POS1.raw", "not available" });

            var i = SdrfImprover.Improve(dep, SdrfDrafter.Draft(Project(), Files));

            Assert.That(i.Disagreements.Any(x => x.Column == "comment[fraction identifier]"), Is.False);
            Assert.That(Row(i.Document, "POS1.raw")["comment[fraction identifier]"], Is.EqualTo("not available"));
        }

        /// <summary>Found by improving the whole corpus (PXD001587): one run listed under two extensions.</summary>
        [Test]
        public void DraftedRowsNeverDuplicateARowKey()
        {
            var files = new[] { "NEG1.raw", "NEG9.raw", "NEG9.mzXML" };

            var i = SdrfImprover.Improve(Deposited(), SdrfDrafter.Draft(Project(), files));

            var keys = i.Document.Results.Select(r => (r["source name"], r["assay name"], r["comment[label]"])).ToList();
            Assert.That(keys.Distinct().Count(), Is.EqualTo(keys.Count));
            Assert.That(SdrfValidator.Validate(i.Document).Errors.Where(e => e.Rule == "RowKeyUniqueness"), Is.Empty);
        }

        /// <summary>
        /// G29, found by the draft -> improve -> restrict chain: a deposit missing columns the specification
        /// requires stayed missing them. They are added, filled "not available" -- except technology type,
        /// whose one specified value every mass-spectrometry deposit has.
        /// </summary>
        [Test]
        public void EveryRequiredColumnIsPresentAfterImprovement()
        {
            var cols = new[] { "source name", "characteristics[biological replicate]", "assay name", "comment[data file]" };
            var dep = Doc(cols, new[] { "p1", "1", "run 1", "NEG1.raw" }, new[] { "p2", "1", "run 2", "POS1.raw" });

            var i = SdrfImprover.Improve(dep, SdrfDrafter.Draft(Project(), Files));

            Assert.That(SdrfValidator.Validate(i.Document).Errors.Where(e => e.Rule == "RequiredColumn"), Is.Empty);
            var h = i.Document.Header.ToList();
            Assert.That(h.IndexOf("technology type"), Is.EqualTo(h.IndexOf("assay name") + 1));
            Assert.That(h.IndexOf("characteristics[organism part]"), Is.LessThan(h.IndexOf("assay name")));
            Assert.That(Row(i.Document, "NEG1.raw")["technology type"], Is.EqualTo("proteomic profiling by mass spectrometry"));
            Assert.That(Row(i.Document, "NEG1.raw")["comment[cleavage agent details]"], Is.EqualTo("not available"));
            Assert.That(i.Document.Results.All(r => r.Cells.Count == h.Count), "never ragged");
        }

        [Test]
        public void AColumnWrittenInOtherCasingIsNeverShadowedByALowercaseCopy()
        {
            var cols = new[] { "Source Name", "Characteristics[organism]", "characteristics[biological replicate]", "Assay Name",
                "Technology Type", "Comment[label]", "comment[data file]" };
            var dep = Doc(cols,
                new[] { "p1", "not available", "1", "run 1", "proteomic profiling by mass spectrometry", "label free sample", "NEG1.raw" },
                new[] { "p2", "homo sapiens", "1", "run 2", "proteomic profiling by mass spectrometry", "label free sample", "POS1.raw" });

            var i = SdrfImprover.Improve(dep, SdrfDrafter.Draft(Project(), Files));

            var h = i.Document.Header.ToList();
            foreach (var name in cols)
                Assert.That(h.Count(x => string.Equals(x, name, StringComparison.OrdinalIgnoreCase)), Is.EqualTo(1), name);
            var neg = Row(i.Document, "NEG1.raw");
            Assert.That((neg["Source Name"], neg["Assay Name"], neg["Comment[label]"]), Is.EqualTo(("p1", "run 1", "label free sample")));
            Assert.That(neg["Characteristics[organism]"], Is.EqualTo("homo sapiens"), "the gap is filled in the depositor's own column");
            Assert.That(SdrfValidator.Validate(i.Document).Errors.Where(e => e.Rule == "RequiredColumn").Select(e => e.Message),
                Has.Some.Contains("differs only in casing"), "the curator still hears about the casing");
        }

        [Test]
        public void ADraftedRowCarriesOnlyAssayWideColumns()
        {
            var cols = new[] { "source name", "characteristics[organism]", "characteristics[individual]", "assay name",
                "comment[label]", "comment[cleavage agent details]", "comment[data file]", "comment[file uri]" };
            var dep = Doc(cols, new[] { "p1", "homo sapiens", "patient 7", "run 1", "label free sample", "NT=Trypsin;AC=MS:1001251", "S1.raw", "ftp://x/S1.raw" });

            var i = SdrfImprover.Improve(dep, SdrfDrafter.Draft(Project(), new[] { "S1.raw", "S2.raw" }));

            var added = Row(i.Document, "S2.raw");
            Assert.That(added["comment[file uri]"], Is.EqualTo("not available"), "never another file's download link");
            Assert.That(added["characteristics[individual]"], Is.EqualTo("not available"), "a sample the deposit did not describe");
            Assert.That(added["comment[label]"], Is.EqualTo("label free sample"));
            Assert.That(added["comment[cleavage agent details]"], Is.EqualTo("NT=Trypsin;AC=MS:1001251"));
        }

        [Test]
        public void ARunListedUnderTwoExtensionsIsNotAMultiplexedFile()
        {
            var cols = new[] { "source name", "characteristics[organism]", "assay name", "comment[label]", "comment[data file]" };
            var dep = Doc(cols,
                new[] { "p1", "not available", "run 1", "label free sample", "S1.raw" },
                new[] { "p1", "not available", "run 1m", "label free sample", "S1.mzML" },
                new[] { "p2", "not available", "run 2", "label free sample", "S2.raw" });

            var i = SdrfImprover.Improve(dep, SdrfDrafter.Draft(Project(), new[] { "S1.raw", "S2.raw" }));

            Assert.That(new[] { "S1.raw", "S1.mzML", "S2.raw" }.Select(f => Row(i.Document, f)["characteristics[organism]"]),
                Is.All.Contains("homo sapiens"));
        }

        [Test]
        public void AFillFollowsTheFormItsColumnIsWrittenIn()
        {
            var cols = new[] { "source name", "characteristics[organism]", "characteristics[disease]", "assay name", "comment[label]", "comment[data file]" };
            var freeText = Doc(cols,
                new[] { "p1", "homo sapiens", "not available", "run 1", "label free sample", "NEG1.raw" },
                new[] { "p2", "not available", "COVID-19", "run 2", "label free sample", "POS1.raw" });
            var terms = Doc(cols,
                new[] { "p1", "NT=homo sapiens;AC=NCBITaxon:9606", "not available", "run 1", "label free sample", "NEG1.raw" },
                new[] { "p2", "not available", "NT=COVID-19;AC=MONDO:0100096", "run 2", "label free sample", "POS1.raw" });

            var asText = SdrfImprover.Improve(freeText, SdrfDrafter.Draft(Project(), Files)).Document;
            var asTerms = SdrfImprover.Improve(terms, SdrfDrafter.Draft(Project(), Files)).Document;

            Assert.That(Row(asText, "POS1.raw")["characteristics[organism]"], Is.EqualTo("homo sapiens"));
            Assert.That(Row(asText, "NEG1.raw")["characteristics[disease]"], Is.EqualTo("normal"));
            Assert.That(Row(asText, "NEG2.raw")["characteristics[organism]"], Is.EqualTo("homo sapiens"), "a drafted row too");
            Assert.That(Row(asTerms, "POS1.raw")["characteristics[organism]"], Does.Contain("AC=NCBITaxon:9606"));
            Assert.That(Row(asTerms, "NEG1.raw")["characteristics[disease]"], Does.Contain("AC=PATO:0000461"));
        }

        [Test]
        public void ARepeatedColumnIsCarriedByPositionNotByName()
        {
            var cols = new[] { "source name", "assay name", "comment[data file]", "comment[modification parameters]", "comment[modification parameters]" };
            var dep = Doc(cols, new[] { "p1", "run 1", "NEG1.raw", "NT=Oxidation", "NT=Carbamidomethyl" });

            var i = SdrfImprover.Improve(dep, SdrfDrafter.Draft(Project(), Files));

            var added = Row(i.Document, "POS1.raw");
            Assert.That(added.All("comment[modification parameters]"), Is.EqualTo(new[] { "NT=Oxidation", "NT=Carbamidomethyl" }));
        }

        /// <summary>
        /// Improves every curated corpus SDRF that has a cached PRIDE listing, and reports what happened.
        /// A crash is a failure; everything else is a count. Needs MZLIB_SDRF_CORPUS (the datasets/ folder)
        /// and MZLIB_PRIDE_CACHE (one folder per accession holding project.json and files.json, as the
        /// sdrf project's benchmark writes them).
        /// </summary>
        [Test]
        [Explicit("Requires the curated SDRF corpus and a PRIDE metadata cache; set MZLIB_SDRF_CORPUS and MZLIB_PRIDE_CACHE.")]
        public void CorpusImproveReport()
        {
            string? corpus = Environment.GetEnvironmentVariable("MZLIB_SDRF_CORPUS");
            string? cache = Environment.GetEnvironmentVariable("MZLIB_PRIDE_CACHE");
            if (string.IsNullOrWhiteSpace(corpus) || string.IsNullOrWhiteSpace(cache)) Assert.Ignore("MZLIB_SDRF_CORPUS / MZLIB_PRIDE_CACHE not set.");

            int files = 0, crashed = 0, filled = 0, addedRows = 0, disagreements = 0, worse = 0, joined = 0;
            var crashes = new List<string>();
            var byColumn = new Dictionary<string, int>(); var examples = new List<string>(); var samples = new List<string>(); int manyAdded = 0, requiredFixed = 0;
            foreach (var path in System.IO.Directory.GetFiles(corpus!, "*.sdrf.tsv", System.IO.SearchOption.AllDirectories))
            {
                string acc = System.IO.Path.GetFileName(System.IO.Path.GetDirectoryName(path))!;
                string projectJson = System.IO.Path.Combine(cache!, acc, "project.json");
                if (!System.IO.File.Exists(projectJson)) continue;
                files++;
                try
                {
                    var project = Newtonsoft.Json.JsonConvert.DeserializeObject<PrideProject>(System.IO.File.ReadAllText(projectJson))!;
                    var listing = Newtonsoft.Json.JsonConvert.DeserializeObject<List<PrideArchiveFile>>(
                        System.IO.File.ReadAllText(System.IO.Path.Combine(cache!, acc, "files.json")))!;
                    var raw = listing.Where(f => string.Equals(f.FileCategory?.Value, "RAW", StringComparison.OrdinalIgnoreCase)).Select(f => f.FileName).ToList();
                    var deposited = new SdrfDocument(path);
                    if (!deposited.Header.Contains("comment[data file]") || raw.Count == 0) continue;
                    joined++;
                    int before = SdrfValidator.Validate(deposited).Errors.Count();
                    int requiredBefore = SdrfValidator.Validate(deposited).Errors.Count(e => e.Rule == "RequiredColumn");
                    var i = SdrfImprover.Improve(deposited, SdrfDrafter.Draft(project, raw));
                    filled += i.FilledCells; addedRows += i.AddedRows; disagreements += i.Disagreements.Count;
                    foreach (var g in i.Disagreements.GroupBy(x => x.Column)) byColumn[g.Key] = byColumn.GetValueOrDefault(g.Key) + g.Count();
                    if (i.AddedRows > deposited.Results.Count) { manyAdded++; if (examples.Count < 6) examples.Add($"{acc}: {deposited.Results.Count} rows, {i.AddedRows} added; deposited e.g. '{deposited.Results[0]["comment[data file]"]}', PRIDE e.g. '{raw[0]}'"); }
                    foreach (var x in i.Disagreements.Where(x => x.Column == "characteristics[organism]" || x.Column == "comment[instrument]").Take(1))
                        if (samples.Count < 8) samples.Add($"{acc} {x.Column}: deposited '{x.Deposited}' vs drafted '{x.Drafted}'");
                    var afterErrors = SdrfValidator.Validate(i.Document).Errors.ToList();
                    if (requiredBefore > 0 && afterErrors.All(e => e.Rule != "RequiredColumn")) requiredFixed++;
                    if (afterErrors.Count > before)
                    {
                        worse++;
                        var beforeRules = SdrfValidator.Validate(deposited).Errors.Select(e => e.Rule).ToHashSet();
                        foreach (var e in afterErrors.Where(e => !beforeRules.Contains(e.Rule)).Take(2)) examples.Add($"WORSE {acc}: {e}");
                    }
                }
                catch (Exception e)
                {
                    crashed++;
                    if (crashes.Count < 15)
                        crashes.Add($"{acc}: {e.GetType().Name}: {e.StackTrace?.Split('\n').FirstOrDefault(l => l.Contains("Sdrf"))?.Trim()}");
                }
            }
            TestContext.Progress.WriteLine($"{files} corpus SDRFs with a cached PRIDE project; {joined} improved; {crashed} crashed");
            TestContext.Progress.WriteLine($"cells filled {filled}; rows added {addedRows}; disagreements reported {disagreements}; validator worse on {worse}");
            TestContext.Progress.WriteLine($"files that were missing a required column and now have them all: {requiredFixed}");
            foreach (var c in crashes) TestContext.Progress.WriteLine("  CRASH " + c);
            foreach (var (c, n) in byColumn.OrderByDescending(kv => kv.Value)) TestContext.Progress.WriteLine($"  disagree {n,7} {c}");
            TestContext.Progress.WriteLine($"  deposits with more rows added than they had: {manyAdded}");
            foreach (var e in examples) TestContext.Progress.WriteLine("  ADDED " + e);
            foreach (var e in samples) TestContext.Progress.WriteLine("  DIFF " + e);
            Assert.That(crashed, Is.Zero);
        }

        [Test]
        public void MalformedArgumentsThrow()
        {
            Assert.Throws<ArgumentNullException>(() => SdrfImprover.Improve(null!, SdrfDrafter.Draft(Project(), Files)));
            Assert.Throws<ArgumentNullException>(() => SdrfImprover.Improve(Deposited(), null!));
        }
    }
}
