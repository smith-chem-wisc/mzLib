using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Tests for restricting an SDRF to the files one search read (sdrf D36): the SDRF that governs the
    /// quantification names exactly the searched files, and joins each to its acquisition by the MAP-12 rule.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrfSearchScope
    {
        private static readonly string[] Columns =
            { "source name", "characteristics[biological replicate]", "assay name", "comment[label]", "comment[data file]", "comment[technical replicate]" };

        private static SdrfDocument Doc(params string[][] rows)
        {
            var header = new SdrfHeader(Columns);
            return new SdrfDocument(header, rows.Select(r => new SdrfRow(header, r)));
        }

        private static SdrfDocument Deposit() => Doc(
            new[] { "S1", "1", "run 1", "label free sample", "WT_1.raw", "1" },
            new[] { "S2", "2", "run 2", "label free sample", "WT_2.raw", "1" },
            new[] { "S3", "3", "run 3", "label free sample", "WT_3.raw", "1" });

        private static SdrfRow Row(SdrfDocument d, string file) => d.Results.Single(r => r["comment[data file]"] == file);

        [Test]
        public void OnlyTheSearchedFilesKeepARowAndTheRestAreReported()
        {
            var s = SdrfSearchScope.Restrict(Deposit(), new[] { @"C:\data\WT_1.raw", @"C:\data\WT_3.raw" });

            Assert.That(s.Document.Results.Select(r => r["comment[data file]"]), Is.EquivalentTo(new[] { "WT_1.raw", "WT_3.raw" }));
            Assert.That(s.DroppedDataFiles, Is.EqualTo(new[] { "WT_2.raw" }));
        }

        [Test]
        public void ADerivativeTheSearchReadIsJoinedToItsAcquisitionAndNamed()
        {
            var s = SdrfSearchScope.Restrict(Deposit(), new[] { @"C:\out\WT_1-calib.mzML", @"C:\out\wt_2-calib-averaged.mzML" });

            Assert.That(Row(s.Document, "WT_1.raw")["comment[searched data file]"], Is.EqualTo("WT_1-calib.mzML"));
            Assert.That(Row(s.Document, "WT_2.raw")["comment[searched data file]"], Is.EqualTo("wt_2-calib-averaged.mzML"), "chained suffixes, any case");
            var h = s.Document.Header.ToList();
            Assert.That(h.IndexOf("comment[searched data file]"), Is.EqualTo(h.IndexOf("comment[data file]") + 1));
        }

        [Test]
        public void AKnownAcquiredNameWinsOverTheSuffixRule()
        {
            var acquired = new Dictionary<string, string> { [@"C:\out\sample_A.mzML"] = "WT_3.raw" };

            var s = SdrfSearchScope.Restrict(Deposit(), new[] { @"C:\out\sample_A.mzML" }, acquired);

            Assert.That(Row(s.Document, "WT_3.raw")["comment[searched data file]"], Is.EqualTo("sample_A.mzML"));
        }

        [Test]
        public void ASearchedFileWithNoRowIsReportedNotInvented()
        {
            var s = SdrfSearchScope.Restrict(Deposit(), new[] { "WT_1.raw", "Blank_01.raw" });

            Assert.That(s.SearchedWithoutRow, Is.EqualTo(new[] { "Blank_01.raw" }));
            Assert.That(s.Document.Results.Count, Is.EqualTo(1));
        }

        [Test]
        public void EveryChannelOfAMultiplexedFileIsKept()
        {
            var header = new SdrfHeader(Columns);
            var tmt = new SdrfDocument(header, new[]
            {
                new SdrfRow(header, new[] { "a", "1", "run 1", "TMT126", "plex1.raw", "1" }),
                new SdrfRow(header, new[] { "b", "2", "run 1", "TMT127N", "plex1.raw", "1" }),
                new SdrfRow(header, new[] { "c", "1", "run 2", "TMT126", "plex2.raw", "1" }),
            });

            var s = SdrfSearchScope.Restrict(tmt, new[] { "plex1-calib.mzML" });

            Assert.That(s.Document.Results.Count, Is.EqualTo(2));
        }

        [Test]
        public void TheGroupingIsCarriedExactlyNeverRenumbered()
        {
            var s = SdrfSearchScope.Restrict(Deposit(), new[] { "WT_1.raw", "WT_3.raw" });

            Assert.That(s.Document.Results.Select(r => r["characteristics[biological replicate]"]), Is.EqualTo(new[] { "1", "3" }),
                "ranking for quantification is the design reader's job (MAP-33), not this one's");
        }

        [Test]
        public void TwoSearchedFilesClaimingOneAcquisitionAreReported()
        {
            var s = SdrfSearchScope.Restrict(Deposit(), new[] { "WT_1.raw", "WT_1-calib.mzML" });

            Assert.That(s.Ambiguous, Is.Not.Empty);
            Assert.That(s.Document.Results.Count(r => r["comment[data file]"] == "WT_1.raw"), Is.EqualTo(1));
        }

        [Test]
        public void AnEarlierSearchsColumnIsReplaced()
        {
            var header = new SdrfHeader(Columns.Append("comment[searched data file]"));
            var old = new SdrfDocument(header, new[] { new SdrfRow(header, new[] { "S1", "1", "run 1", "label free sample", "WT_1.raw", "1", "WT_1-old.mzML" }) });

            var s = SdrfSearchScope.Restrict(old, new[] { "WT_1-calib.mzML" });

            Assert.That(s.Document.Results.Single()["comment[searched data file]"], Is.EqualTo("WT_1-calib.mzML"));
            Assert.That(s.Document.Header.Count(h => h == "comment[searched data file]"), Is.EqualTo(1));
        }

        /// <summary>
        /// The whole chain for one deposit: draft from PRIDE, improve the deposited SDRF with it, restrict to a
        /// calibrated search of part of it. The grouping the search quantifies with is exactly the improved
        /// SDRF's (the D36 invariant, on the SDRF side), and the result validates.
        /// </summary>
        [Test]
        public void DraftImproveAndRestrictCarryTheGroupingIntact()
        {
            var project = new UsefulProteomicsDatabases.PrideProject { Accession = "PXD0", ProjectDescription = "Wild type (WT) and knockout (KO) cultures." };
            var raw = new[] { "WT_1.raw", "WT_2.raw", "WT_3.raw", "KO_1.raw", "KO_2.raw", "KO_3.raw" };
            var improved = SdrfImprover.Improve(Deposit(), SdrfDrafter.Draft(project, raw)).Document;
            var searched = new[] { "WT_1-calib.mzML", "WT_3-calib.mzML", "KO_2-calib.mzML" };

            var s = SdrfSearchScope.Restrict(improved, searched);

            Assert.That(s.Document.Results.Select(r => r["comment[searched data file]"]), Is.EquivalentTo(searched));
            string[] grouping = { "source name", "characteristics[biological replicate]", "comment[technical replicate]" };
            foreach (var row in s.Document.Results)
            {
                var before = improved.Results.Single(r => r["comment[data file]"] == row["comment[data file]"]);
                Assert.That(grouping.Select(c => row[c]), Is.EqualTo(grouping.Select(c => before[c])));
            }
            // The fixture's deposit lacks required columns; restricting must add no error of its own.
            var rulesBefore = SdrfValidator.Validate(improved).Errors.Select(e => e.Rule).ToHashSet();
            Assert.That(SdrfValidator.Validate(s.Document).Errors.Where(e => !rulesBefore.Contains(e.Rule)), Is.Empty);
        }

        [Test]
        public void MalformedArgumentsThrow()
        {
            Assert.Throws<ArgumentNullException>(() => SdrfSearchScope.Restrict(null!, new[] { "a.raw" }));
            Assert.Throws<ArgumentNullException>(() => SdrfSearchScope.Restrict(Deposit(), null!));
            Assert.Throws<ArgumentException>(() => SdrfSearchScope.Restrict(Deposit(), Array.Empty<string>()));
        }
    }
}
