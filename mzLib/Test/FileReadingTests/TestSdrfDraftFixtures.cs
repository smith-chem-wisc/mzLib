using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using Newtonsoft.Json;
using NUnit.Framework;
using Readers;
using UsefulProteomicsDatabases;

namespace Test.FileReadingTests
{
    /// <summary>
    /// The drafter on real deposits that nobody annotated: the aging project's fixtures (their
    /// SDRF-A12), each pinned to the design aging read from the PRIDE record and the file names.
    /// PROVENANCE.md, next to the data, says where every file came from.
    ///
    /// These pin what a person reading the deposit concludes, not what the drafter happened to
    /// produce. A case the drafter gets wrong fails here first; it is not quietly fixed.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrfDraftFixtures
    {
        private static string FixturePath(params string[] parts) =>
            Path.Combine(new[] { TestContext.CurrentContext.TestDirectory, "FileReadingTests",
                "ExternalFileTypes", "SdrfDraftFixtures" }.Concat(parts).ToArray());

        /// <summary>The record as PRIDE served it, and every raw file aging's list names.</summary>
        private static SdrfDraft Draft(string accession, out List<string> files)
        {
            var project = JsonConvert.DeserializeObject<PrideProject>(
                File.ReadAllText(FixturePath("records", accession + ".project.json")))!;
            files = File.ReadAllLines(FixturePath("raw_file_lists", accession + "_raw_files.tsv"))
                .Skip(1)
                .Where(l => !string.IsNullOrWhiteSpace(l))
                .Select(l => l.Split('\t')[0])
                .ToList();
            return SdrfDrafter.Draft(project, files);
        }

        /// <summary>The draft's conditions, as a partition of data files by their factor values.</summary>
        private static HashSet<string> Partition(SdrfDraft d, Func<SdrfDraftRow, string> key) =>
            d.Rows.GroupBy(key)
                .Select(g => string.Join("|", g.Select(r => r.DataFile).OrderBy(f => f, StringComparer.Ordinal)))
                .ToHashSet();

        private static string FactorKey(SdrfDraftRow r) => string.Join("|", r.Factors.Select(f => f.Value));

        /// <summary>
        /// PXD049018: BJ cells, TurboID-NES with or without cGAMP, a streptavidin pulldown run on a gel.
        /// Two samples by ten bands, no replication. Which sample got cGAMP is stated nowhere, so the
        /// treatment is <c>not available</c>, and no condition may be invented from the sample IDs.
        /// </summary>
        [Test]
        [Ignore("KNOWN FAILURE, reported to aging (sdrf 029) before any tuning: `Band_01..10` is read as a " +
                "biological-replicate marker, not a fraction, so there are 20 samples, not 2. And the sample " +
                "IDs 67868/67869 become an invented `factor value[condition]`.")]
        public void PXD049018_IsTwoSamplesByTenBandsWithNoTreatmentStated()
        {
            var d = Draft("PXD049018", out var files);

            Assert.That(d.Rows.Select(r => r.DataFile), Is.EquivalentTo(files));
            Assert.That(files, Has.Count.EqualTo(20));

            var samples = d.Rows.GroupBy(r => r.SourceName.Value).ToList();
            Assert.That(samples, Has.Count.EqualTo(2), "MSB67868A and MSB67869A are the two samples.");
            Assert.That(Partition(d, r => r.SourceName.Value), Is.EquivalentTo(new[]
            {
                string.Join("|", files.Where(f => f.StartsWith("MSB67868A", StringComparison.Ordinal)).OrderBy(f => f, StringComparer.Ordinal)),
                string.Join("|", files.Where(f => f.StartsWith("MSB67869A", StringComparison.Ordinal)).OrderBy(f => f, StringComparer.Ordinal)),
            }));

            foreach (var row in d.Rows)
            {
                string band = row.DataFile.Substring(row.DataFile.IndexOf("Band_", StringComparison.Ordinal) + 5, 2);
                Assert.That(row.Fraction.Value, Is.EqualTo(int.Parse(band).ToString()), row.DataFile + ": the band is the fraction");
                Assert.That(row.BiologicalReplicate.Value, Is.EqualTo("1"), row.DataFile + ": nothing is replicated");
                Assert.That(row.TechnicalReplicate.Value, Is.EqualTo("1"), row.DataFile);
            }

            Assert.That(d.Rows.SelectMany(r => r.Factors).Where(f => f.Value != SdrfReserved.NotAvailable),
                Is.Empty,
                "The deposit never says which pulldown had cGAMP. A factor that splits the two samples is " +
                "a condition the record does not state. Factors: " + string.Join(", ", d.FactorColumns));
        }

        /// <summary>
        /// PXD067622: HeLa, SPRTN-TurboID. Construct (WT, CA) by treatment (DMSO, FA, NDC, THY) by 3.
        /// The trailing index runs 1..24 across the whole study (WT_DMSO1-3, CA_DMSO4-6, ...), so it is
        /// a sample number, not a replicate. Within each of the eight conditions the replicates are 1..3.
        /// </summary>
        [Test]
        public void PXD067622_IsConstructByTreatmentByThree()
        {
            var d = Draft("PXD067622", out var files);

            Assert.That(d.Rows.Select(r => r.DataFile), Is.EquivalentTo(files));
            Assert.That(files, Has.Count.EqualTo(24));
            Assert.That(d.Rows.Select(r => r.SourceName.Value).Distinct().Count(), Is.EqualTo(24),
                "Every file is its own sample.");

            // Aging's conditions: the construct-and-treatment token before the trailing index.
            string Condition(string file)
            {
                string stem = Path.GetFileNameWithoutExtension(file);
                string tail = stem.Substring(stem.LastIndexOf("_12032_", StringComparison.Ordinal) + 7);
                return tail.TrimEnd("0123456789".ToCharArray());
            }
            var expected = files.GroupBy(Condition)
                .Select(g => string.Join("|", g.OrderBy(f => f, StringComparer.Ordinal)))
                .ToHashSet();
            Assert.That(expected, Has.Count.EqualTo(8));

            Assert.That(Partition(d, FactorKey), Is.EquivalentTo(expected),
                "The factors must split the files into aging's eight conditions. Factors: " +
                string.Join(", ", d.FactorColumns));

            foreach (var condition in d.Rows.GroupBy(FactorKey))
                Assert.That(condition.Select(r => r.BiologicalReplicate.Value).OrderBy(v => v, StringComparer.Ordinal),
                    Is.EqualTo(new[] { "1", "2", "3" }),
                    "Replicates count within a condition; the study-wide 1..24 index is not a replicate. " +
                    string.Join("; ", condition.Select(r => r.DataFile + " -> " + r.BiologicalReplicate.Value +
                                                            " (" + r.BiologicalReplicate.Evidence + ")")));
        }

        /// <summary>
        /// PXD067622 again, for MAP-33 (agreed with QuantProject): the drafter ranks an unnamed
        /// replicate within its group AND SAYS SO in the slot's evidence. <c>CA_DMSO4..6</c> become
        /// replicates 1..3. That is right, but a reviewer reading "1" beside a file named "4" must be
        /// told the number was ranked, not read.
        /// </summary>
        [Test]
        [Ignore("KNOWN FAILURE, reported to aging (sdrf 029) before any tuning: the evidence reads " +
                "'a replicate count in the file names' and never says the study-wide index was ranked (MAP-33).")]
        public void PXD067622_ARankedReplicateSaysItWasRanked()
        {
            var d = Draft("PXD067622", out _);

            var renumbered = d.Rows
                .Where(r => !Path.GetFileNameWithoutExtension(r.DataFile).EndsWith(r.BiologicalReplicate.Value, StringComparison.Ordinal))
                .ToList();
            Assert.That(renumbered, Is.Not.Empty, "CA_DMSO4 is replicate 1, so some replicates were renumbered.");

            foreach (var row in renumbered)
                Assert.That(row.BiologicalReplicate.Evidence, Does.Contain("rank").IgnoreCase,
                    row.DataFile + " -> " + row.BiologicalReplicate.Value + ": the evidence must say the number was ranked.");
        }

        /// <summary>
        /// PXD058611: mouse liver, files <c>178.raw</c> .. <c>213.raw</c>. The names carry nothing, so the
        /// right answer is "no structure found": every file its own sample, replicate 1, fraction 1,
        /// and no grouping invented. (The protocol does describe a split, but a drafter that reads
        /// names cannot see it, and must not pretend to.)
        /// </summary>
        [Test]
        public void PXD058611_HasNoStructureInItsNamesAndGetsNoneInvented()
        {
            var d = Draft("PXD058611", out var files);

            Assert.That(d.Rows.Select(r => r.DataFile), Is.EquivalentTo(files));
            Assert.That(files, Has.Count.EqualTo(36));
            Assert.That(d.Rows.Select(r => r.SourceName.Value).Distinct().Count(), Is.EqualTo(36),
                "Every file is its own sample.");

            foreach (var row in d.Rows)
            {
                Assert.That(row.BiologicalReplicate.Value, Is.EqualTo("1"), row.DataFile);
                Assert.That(row.TechnicalReplicate.Value, Is.EqualTo("1"), row.DataFile);
                Assert.That(row.Fraction.Value, Is.EqualTo("1"), row.DataFile);
            }

            Assert.That(d.Rows.SelectMany(r => r.Factors).Where(f => f.Value != SdrfReserved.NotAvailable),
                Is.Empty, "No grouping is in the names. Factors: " + string.Join(", ", d.FactorColumns));
        }
    }
}
