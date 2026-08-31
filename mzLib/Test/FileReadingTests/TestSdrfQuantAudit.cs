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
    /// Covers <see cref="SdrfQuantAuditor"/>, which reads an SDRF for what it says about
    /// quantification rather than for whether it is well formed.
    ///
    /// The defects it looks for are all SILENT: a kit-only file parses, validates, and then yields
    /// one channel for a ten-channel experiment; a duplicated (data file, channel) pair makes a
    /// dictionary keyed on it drop measurements; a channel-level file declaring no isobaric
    /// modification cannot be searched from its own annotation. None is a structural error, so
    /// <see cref="SdrfValidator"/> is right not to report them.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrfQuantAudit
    {
        private static string Fixture(string name) =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, "FileReadingTests",
                "ExternalFileTypes", "TmtAuditFixture", name);

        private static string FixtureDataDirectory() =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, "FileReadingTests",
                "ExternalFileTypes", "TmtAuditFixture", "data");

        [Test]
        public void ChannelLevelDocumentReportsItsChannels()
        {
            var audit = SdrfQuantAuditor.Audit(Fixture("TmtChannelLevel.sdrf.tsv"));

            Assert.That(audit.Kind, Is.EqualTo(SdrfQuantKind.ChannelLevel));
            Assert.That(audit.Channels, Is.EquivalentTo(new[] { "TMT126", "TMT127N", "TMT128N", "TMT129C" }));
            Assert.That(audit.ImpliedPlex, Is.EqualTo(4));
            Assert.That(audit.DeclaresIsobaricModification, Is.True);
        }

        /// <summary>
        /// The distinction the whole type exists for. A kit-only file names the reagent on every row
        /// and never enumerates a channel, so there is no channel-level design in it to recover --
        /// and read as though there were, every row claims the same channel.
        /// </summary>
        [Test]
        public void KitOnlyDocumentIsNotMistakenForAChannelLevelOne()
        {
            var audit = SdrfQuantAuditor.Audit(Fixture("TmtKitOnly.sdrf.tsv"));

            Assert.That(audit.Kind, Is.EqualTo(SdrfQuantKind.KitOnly));
            Assert.That(audit.Channels, Is.Empty);
            Assert.That(audit.KitLabels, Is.EquivalentTo(new[] { "TMT10PLEX" }));
            Assert.That(audit.DeclaresIsobaricModification, Is.False,
                "and it declares no isobaric modification either, so nothing in it says the run was labelled");

            var channel = audit.Facts.Single(f => f.Number == 2);
            Assert.That(channel.Presence, Is.EqualTo(SdrfFactPresence.Unparseable));
            Assert.That(channel.Detail, Does.Contain("kit"));
        }

        /// <summary>
        /// A label-free document is a legitimate answer, not a failure. The audit is pointed at
        /// whatever a user just downloaded, so "this is not an isobaric experiment" has to be
        /// something it can say.
        /// </summary>
        [Test]
        public void LabelFreeDocumentIsReportedAsNotIsobaric()
        {
            string labelFree = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "FileReadingTests", "ExternalFileTypes", "PXD000070.sdrf.tsv");

            var audit = SdrfQuantAuditor.Audit(labelFree);

            Assert.That(audit.Kind, Is.EqualTo(SdrfQuantKind.NotIsobaric));
            Assert.That(audit.Channels, Is.Empty);
            Assert.That(audit.KitLabels, Is.Empty);
        }

        /// <summary>
        /// One curated corpus file writes both label forms, so the form is a per-cell decision. A
        /// reader that sniffed the document once would misread every cell of the other form.
        /// </summary>
        [Test]
        public void LabelFormIsCountedPerCellSoAMixedDocumentIsVisible()
        {
            var audit = SdrfQuantAuditor.Audit(Fixture("TmtChannelLevel.sdrf.tsv"));

            Assert.That(audit.AccessionedLabelCells, Is.EqualTo(2), "the two NT=TMT127N cells");
            Assert.That(audit.BareLabelCells, Is.EqualTo(7));
            Assert.That(audit.MixesLabelForms, Is.True);
        }

        /// <summary>
        /// Both key orders occur in the corpus -- NT first and AC first -- so the name is taken by
        /// parsing the cell, never by matching a prefix.
        /// </summary>
        [TestCase("NT=TMT127N;AC=PRIDE:0000519")]
        [TestCase("AC=PRIDE:0000519;NT=TMT127N")]
        public void AccessionedLabelsAreReadInEitherKeyOrder(string cell)
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, $"keyorder_{Guid.NewGuid():N}.sdrf.tsv");
            File.WriteAllText(path,
                "source name\tcomment[data file]\tcomment[label]\r\n"
                + $"s1\tfile.raw\t{cell}\r\n");
            try
            {
                var audit = SdrfQuantAuditor.Audit(path);
                Assert.That(audit.Channels, Is.EqualTo(new[] { "TMT127N" }));
                Assert.That(audit.AccessionedLabelCells, Is.EqualTo(1));
            }
            finally
            {
                File.Delete(path);
            }
        }

        /// <summary>
        /// Case varies in the corpus -- itraq114, TMT10Plex, tmt126 -- and those are the same reagents
        /// as their canonical spellings. Folding them is right for identity; reporting the variant is
        /// right for anyone joining on the raw string.
        /// </summary>
        [Test]
        public void SpellingVariantsAreFoldedForIdentityAndStillReported()
        {
            var audit = SdrfQuantAuditor.Audit(Fixture("TmtChannelLevel.sdrf.tsv"));

            Assert.That(audit.Channels, Does.Contain("TMT128N"), "tmt128N is that channel, not another one");
            Assert.That(audit.NonCanonicalChannelSpellings, Does.Contain("tmt128N"));
        }

        /// <summary>
        /// A channel is one measurement of one file, so a repeated (data file, channel) pair means the
        /// pair does not identify a row -- and any dictionary keyed on it silently keeps one. Five
        /// files in the curated corpus do this, one of them 551 times.
        /// </summary>
        [Test]
        public void RepeatedFileChannelPairsAreReported()
        {
            var audit = SdrfQuantAuditor.Audit(Fixture("TmtChannelLevel.sdrf.tsv"));

            Assert.That(audit.DuplicateFileLabelPairs, Has.Count.EqualTo(1));
            Assert.That(audit.DuplicateFileLabelPairs[0].DataFile, Is.EqualTo("plexA_fr1.raw"));
            Assert.That(audit.DuplicateFileLabelPairs[0].Label, Is.EqualTo("TMT126"));
            Assert.That(audit.DuplicateFileLabelPairs[0].Rows, Is.EqualTo(2));
        }

        /// <summary>
        /// A reserved word is the spec-correct way to say nothing, so it is absence rather than an
        /// answer -- a validator is right to accept it and a consumer still learns nothing. A value
        /// that is neither reserved nor parseable is a different problem and gets a different verdict.
        /// </summary>
        [Test]
        public void ReservedWordsReadAsAbsentAndBadValuesAsUnparseable()
        {
            var audit = SdrfQuantAuditor.Audit(Fixture("TmtKitOnly.sdrf.tsv"));

            var sampleType = audit.Facts.Single(f => f.Number == 10);
            Assert.That(sampleType.Presence, Is.EqualTo(SdrfFactPresence.Absent),
                "'not available' is spec-correct and still answers nothing");

            var fraction = audit.Facts.Single(f => f.Number == 9);
            Assert.That(fraction.Presence, Is.EqualTo(SdrfFactPresence.Unparseable));
            Assert.That(fraction.Detail, Does.Contain("19b"));
        }

        /// <summary>
        /// The plex is reported as unstated rather than guessed. The project measured the source-name
        /// partition at 3 correct in 9 where ground truth exists, failing by over-splitting, and no
        /// bulk dataset in the corpus carries a usable batch column -- so inventing one here would
        /// manufacture plexes that are not in the file.
        /// </summary>
        [Test]
        public void AnUnstatedPlexIsReportedRatherThanInferred()
        {
            var audit = SdrfQuantAuditor.Audit(Fixture("TmtChannelLevel.sdrf.tsv"));

            Assert.That(audit.PlexSource, Is.EqualTo(SdrfPlexSource.None));
            Assert.That(audit.Plexes, Is.Empty);
            Assert.That(audit.ToReport(), Does.Contain("none stated"));
        }

        /// <summary>
        /// Disk matching cannot be exercised by the curated corpus, which carries no data files at
        /// all -- hence the zero-byte stub beside the fixture. One of the two named files is
        /// deliberately absent, because "0 of 2 found" and "not looked" must not read the same.
        /// </summary>
        [Test]
        public void DataFilesAreMatchedOnDiskWhenADirectoryIsGiven()
        {
            var audit = SdrfQuantAuditor.Audit(Fixture("TmtChannelLevel.sdrf.tsv"), FixtureDataDirectory());

            Assert.That(audit.SearchedForDataFiles, Is.True);
            Assert.That(audit.MatchedDataFiles, Is.EqualTo(new[] { "plexA_fr1.raw" }));
            Assert.That(audit.UnmatchedDataFiles, Is.EqualTo(new[] { "plexA_fr2.raw" }));
        }

        [Test]
        public void NoDirectoryMeansNotLookedRatherThanNothingFound()
        {
            var audit = SdrfQuantAuditor.Audit(Fixture("TmtChannelLevel.sdrf.tsv"));

            Assert.That(audit.SearchedForDataFiles, Is.False);
            Assert.That(audit.MatchedDataFiles, Is.Empty);
            Assert.That(audit.UnmatchedDataFiles, Is.Empty);
        }

        [Test]
        public void ReportNamesTheDesignTheChannelsAndTheFacts()
        {
            string report = SdrfQuantAuditor.Audit(Fixture("TmtChannelLevel.sdrf.tsv")).ToReport();

            Assert.That(report, Does.Contain("ChannelLevel"));
            Assert.That(report, Does.Contain("TMT127N"));
            Assert.That(report, Does.Contain("reporter m/z"));
            Assert.That(report, Does.Contain("duplicate rows"));
        }

        /// <summary>
        /// All eleven facts are reported for every document, present or not. A fact that vanished
        /// from the report when it was missing would be indistinguishable from one nobody checked.
        /// </summary>
        [Test]
        public void EveryDocumentReportsAllElevenFacts()
        {
            foreach (string fixture in new[] { "TmtChannelLevel.sdrf.tsv", "TmtKitOnly.sdrf.tsv" })
            {
                var audit = SdrfQuantAuditor.Audit(Fixture(fixture));
                Assert.That(audit.Facts.Select(f => f.Number), Is.EqualTo(Enumerable.Range(1, 11)), fixture);
            }
        }

        /// <summary>
        /// Fact 3 has no SDRF column at all, which is a different answer from "the column is empty"
        /// and has to survive as one -- it is why the reporter m/z has to come from a channel table.
        /// </summary>
        [Test]
        public void ReporterMzIsReportedAsHavingNoColumnRatherThanAsUnfilled()
        {
            var fact = SdrfQuantAuditor.Audit(Fixture("TmtChannelLevel.sdrf.tsv")).Facts.Single(f => f.Number == 3);

            Assert.That(fact.Column, Is.Null);
            Assert.That(fact.Presence, Is.EqualTo(SdrfFactPresence.Absent));
            Assert.That(fact.Detail, Does.Contain("no column"));
        }

        /// <summary>
        /// Classifies the whole curated corpus and compares the counts to a checked-in expectations
        /// table, so that drift in the corpus OR in this classifier shows up as a diff rather than as
        /// a number nobody re-derived.
        ///
        /// This is the milestone's reason for existing. Three independent surveys of this corpus
        /// returned 49, 51 and 52 channel-level TMT files and 31 versus 32 kit-only; the qualitative
        /// findings agreed and the integers did not. A count that is not produced by committed code
        /// against a committed expectation is a number somebody remembered.
        ///
        /// [Explicit]; needs MZLIB_SDRF_CORPUS, following TestSdrfValidator's corpus calibration.
        /// </summary>
        [Test]
        [Explicit("Requires a local clone of bigbio/sdrf-annotated-datasets; set MZLIB_SDRF_CORPUS.")]
        public void ClassifiesTheEntireCorpusAgainstTheExpectationsTable()
        {
            string corpus = Environment.GetEnvironmentVariable("MZLIB_SDRF_CORPUS");
            if (string.IsNullOrWhiteSpace(corpus) || !Directory.Exists(corpus))
                Assert.Ignore($"MZLIB_SDRF_CORPUS not set or not found: '{corpus}'");

            var files = Directory.GetFiles(CuratedRoot(corpus), "*.sdrf.tsv", SearchOption.AllDirectories);
            Assert.That(files, Is.Not.Empty, "no .sdrf.tsv under the corpus root");

            // An unreadable file is counted, not tolerated silently and not thrown out of the run.
            // One curated file is zero bytes, and a count that quietly skipped it would be the same
            // kind of unre-derived number this test exists to stop.
            var audits = new List<SdrfQuantAudit>(files.Length);
            int unreadable = 0;
            foreach (string file in files)
            {
                try
                {
                    audits.Add(SdrfQuantAuditor.Audit(file));
                }
                catch (MzLibUtil.MzLibException)
                {
                    unreadable++;
                }
            }

            var channelLevel = audits.Where(a => a.Kind == SdrfQuantKind.ChannelLevel).ToList();
            var kitOnly = audits.Where(a => a.Kind == SdrfQuantKind.KitOnly).ToList();

            var observed = new Dictionary<string, int>
            {
                ["files"] = audits.Count,
                ["channel-level"] = channelLevel.Count,
                ["kit-only"] = kitOnly.Count,
                ["not-isobaric"] = audits.Count(a => a.Kind == SdrfQuantKind.NotIsobaric),
                ["channel-level, no isobaric modification"] = channelLevel.Count(a => !a.DeclaresIsobaricModification),
                ["channel-level, plex column present"] = channelLevel.Count(a => a.PlexSource == SdrfPlexSource.Column),
                ["channel-level, duplicate (file, channel)"] = channelLevel.Count(a => a.DuplicateFileLabelPairs.Count > 0),
                ["channel-level, ragged rows"] = channelLevel.Count(a => a.RaggedRows > 0),
                ["channel-level, any accessioned label"] = channelLevel.Count(a => a.AccessionedLabelCells > 0),
                ["documents mixing label forms"] = audits.Count(a => a.MixesLabelForms),
                ["channel-level, non-canonical spelling"] = channelLevel.Count(a => a.NonCanonicalChannelSpellings.Count > 0),
                ["unreadable"] = unreadable,

                // The family split is what makes this corpus's counts comparable across surveys.
                // 51 channel-level and 31 kit-only is right for the TMT family alone -- which is the
                // number the project's earlier TMT-scoped surveys were reaching for -- and iTRAQ adds
                // 10 and 6 on top. Pinning both means a future disagreement is a scope question with
                // an answer rather than three integers and no way to choose.
                ["channel-level, TMT family"] = channelLevel.Count(a => a.Families.Contains("TMT")),
                ["kit-only, TMT family"] = kitOnly.Count(a => a.Families.Contains("TMT")),
                ["channel-level, iTRAQ only"] = channelLevel.Count(a => !a.Families.Contains("TMT")),
                ["kit-only, iTRAQ only"] = kitOnly.Count(a => !a.Families.Contains("TMT")),
            };

            foreach (var pair in observed)
                TestContext.Progress.WriteLine($"{pair.Key,-45}{pair.Value}");

            var expected = ReadExpectations();
            var drift = new List<string>();
            foreach (var pair in expected)
            {
                if (!observed.TryGetValue(pair.Key, out int actual))
                {
                    drift.Add($"{pair.Key}: no longer measured");
                }
                else if (actual != pair.Value)
                {
                    drift.Add($"{pair.Key}: expected {pair.Value}, observed {actual}");
                }
            }

            Assert.That(drift, Is.Empty,
                "the corpus or the classifier moved. Re-derive, decide which changed, and update "
                + "SdrfQuantAuditCorpusExpectations.tsv in the same commit:\n  "
                + string.Join("\n  ", drift));
        }

        /// <summary>
        /// The named duplicate cases, kept as their own test because they are the numbers that make
        /// the defect concrete rather than statistical.
        ///
        /// The project's earlier survey recorded PXD040455 at 551 repeated (data file, channel) pairs,
        /// and that is confirmed exactly here. It also called it the worst in the corpus, and that
        /// part is wrong: PXD017291-mixed-label carries 781. Both are pinned, so neither claim has to
        /// be taken on trust again.
        /// </summary>
        [Test]
        [Explicit("Requires a local clone of bigbio/sdrf-annotated-datasets; set MZLIB_SDRF_CORPUS.")]
        public void TheNamedDuplicatedDocumentsStillCarryTheCountsWeMeasured()
        {
            string corpus = Environment.GetEnvironmentVariable("MZLIB_SDRF_CORPUS");
            if (string.IsNullOrWhiteSpace(corpus) || !Directory.Exists(corpus))
                Assert.Ignore($"MZLIB_SDRF_CORPUS not set or not found: '{corpus}'");

            var byName = Directory.GetFiles(CuratedRoot(corpus), "*.sdrf.tsv", SearchOption.AllDirectories)
                .Where(f => new FileInfo(f).Length > 0)
                .Select(f => SdrfQuantAuditor.Audit(f))
                .Where(a => a.DuplicateFileLabelPairs.Count > 0)
                .ToDictionary(a => Path.GetFileName(a.FilePath)!, a => a);

            Assert.That(byName["PXD040455.sdrf.tsv"].DuplicateFileLabelPairs, Has.Count.EqualTo(551));
            Assert.That(byName["PXD017291-mixed-label.sdrf.tsv"].DuplicateFileLabelPairs, Has.Count.EqualTo(781),
                "the worst in the corpus, which the earlier survey missed");

            // PXD007160 is the other shape: few pairs, but one of them spread over 21 separate rows.
            var pxd007160 = byName["PXD007160.sdrf.tsv"];
            Assert.That(pxd007160.DuplicateFileLabelPairs, Has.Count.EqualTo(80));
            Assert.That(pxd007160.DuplicateFileLabelPairs[0].Rows, Is.EqualTo(21));

            Assert.That(byName, Has.Count.EqualTo(5), "five documents duplicate at all");
        }

        /// <summary>
        /// The curated corpus is the <c>datasets</c> subtree. The clone also carries a <c>sandbox</c>
        /// directory of 295 further files that are explicitly not part of it, and counting those was
        /// how at least one earlier survey of this corpus came back with numbers nobody could
        /// reproduce. Accepting either the clone root or the subtree itself removes that as a way to
        /// get a different answer by pointing the variable one level out.
        /// </summary>
        private static string CuratedRoot(string corpus)
        {
            string datasets = Path.Combine(corpus, "datasets");
            return Directory.Exists(datasets) ? datasets : corpus;
        }

        private static Dictionary<string, int> ReadExpectations()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "FileReadingTests",
                "ExternalFileTypes", "TmtAuditFixture", "SdrfQuantAuditCorpusExpectations.tsv");

            var expected = new Dictionary<string, int>();
            foreach (string line in File.ReadLines(path))
            {
                if (line.Length == 0 || line[0] == '#') continue;
                string[] parts = line.Split('\t');
                if (parts.Length == 2 && int.TryParse(parts[1], out int value))
                    expected[parts[0]] = value;
            }
            Assert.That(expected, Is.Not.Empty, $"no expectations read from {path}");
            return expected;
        }
    }
}
