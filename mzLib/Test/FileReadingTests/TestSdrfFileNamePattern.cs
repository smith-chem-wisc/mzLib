using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Tests for reading a design out of one deposit's file names. The refusals matter as much as the
    /// reads: the plex audit measured file-name partitioning at 3 correct of 9, failing by
    /// over-splitting, so every case that should NOT yield structure is pinned here too.
    ///
    /// The deposit-shaped sets follow the patterns aging reported for PXD049018 (2 lysates x 10 gel
    /// bands), PXD067622 (genotype x treatment x 3) and PXD024803 / PXD032040 (IP with IgG controls).
    /// They are written from that description; the real PRIDE listings replace them when aging sends
    /// them (SDRF-A12).
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrfFileNamePattern
    {
        private static List<string> Names(string pattern, params (string, int)[] ranges) =>
            ranges.SelectMany(r => Enumerable.Range(1, r.Item2).Select(i => string.Format(pattern, r.Item1, i))).ToList();

        private static SdrfFileNameSlot Slot(SdrfFileNameStructure s, SdrfFileNameRole role) =>
            s.Slots.Single(x => x.Role == role);

        [Test]
        public void GelBandsAreFractionsOfTheirLysateNotSamples()
        {
            var names = Names("Lysate{0}_Band{1:00}.raw", ("1", 10), ("2", 10));

            var s = SdrfFileNamePattern.Read(names);

            Assert.That(s.Found, Is.True, s.NoStructureReason);
            Assert.That(s.Files.Select(f => f.SampleKey).Distinct().Count(), Is.EqualTo(2), "20 files, 2 samples");
            Assert.That(s.Files.Select(f => f.Fraction), Is.EquivalentTo(Enumerable.Range(1, 10).Concat(Enumerable.Range(1, 10)).Select(i => (int?)i)));
            Assert.That(Slot(s, SdrfFileNameRole.Fraction).Evidence, Does.Contain("'Band'"));
            Assert.That(Slot(s, SdrfFileNameRole.Sample).Evidence, Does.Contain("'Lysate'"));
            Assert.That(s.Files.All(f => f.Replicate == null && f.BiologicalReplicate == null));
        }

        [Test]
        public void GenotypeByTreatmentByThreeGivesTwoFactorsAndAnUnnamedReplicate()
        {
            var names = new[] { "WT_DMSO", "WT_Drug", "KO_DMSO", "KO_Drug" }
                .SelectMany(arm => Enumerable.Range(1, 3).Select(i => $"{arm}{i}.raw")).ToList();

            var s = SdrfFileNamePattern.Read(names);

            Assert.That(s.Found, Is.True, s.NoStructureReason);
            Assert.That(s.Slots.Count(x => x.Role == SdrfFileNameRole.Factor), Is.EqualTo(2));
            Assert.That(Slot(s, SdrfFileNameRole.Replicate).Levels, Is.EqualTo(new[] { "1", "2", "3" }));
            Assert.That(s.Files.Select(f => f.SampleKey).Distinct().Count(), Is.EqualTo(12));
            var wtDmso2 = s.Files.Single(f => f.FileName == "WT_DMSO2.raw");
            Assert.That(wtDmso2.FactorLevels, Is.EqualTo(new[] { "WT", "DMSO" }));
            Assert.That(wtDmso2.Replicate, Is.EqualTo(2));
            Assert.That(wtDmso2.BiologicalReplicate, Is.Null, "a bare trailing number does not say biological or technical");
        }

        /// <summary>
        /// A bait gene with a digit in its name beside IgG leaves the names ragged when every word+number
        /// is split; splitting only the last part aligns them without breaking the gene name.
        /// </summary>
        [Test]
        public void AnIgGControlIsItsOwnLevelAndNeverMergesIntoTheBait()
        {
            var names = Names("{0}_IP_{1}.raw", ("CHD6", 3), ("IgG", 3));

            var s = SdrfFileNamePattern.Read(names);

            Assert.That(s.Found, Is.True, s.NoStructureReason);
            Assert.That(Slot(s, SdrfFileNameRole.Factor).Levels, Is.EqualTo(new[] { "CHD6", "IgG" }));
            Assert.That(Slot(s, SdrfFileNameRole.Factor).Evidence, Does.Contain("'IgG' read as a control"));
            Assert.That(s.Files.Where(f => f.IsControl).Select(f => f.FileName), Is.EquivalentTo(Names("{0}_IP_{1}.raw", ("IgG", 3))));
            Assert.That(s.Files.Where(f => f.IsControl).Select(f => f.SampleKey)
                .Intersect(s.Files.Where(f => !f.IsControl).Select(f => f.SampleKey)), Is.Empty);
        }

        [Test]
        public void NamedIndicesKeepTheirKind()
        {
            var names = new List<string>();
            foreach (int sample in new[] { 1, 2 })
                foreach (int br in new[] { 1, 2 })
                    foreach (int inj in new[] { 1, 2 })
                        names.Add($"S{sample}_BR{br}_inj{inj}.raw");

            var s = SdrfFileNamePattern.Read(names);

            Assert.That(s.Found, Is.True, s.NoStructureReason);
            var f = s.Files.Single(x => x.FileName == "S2_BR1_inj2.raw");
            Assert.That((f.BiologicalReplicate, f.TechnicalReplicate), Is.EqualTo(((int?)1, (int?)2)));
            Assert.That(s.Files.Select(x => x.SampleKey).Distinct().Count(), Is.EqualTo(4),
                "a re-injection is the same sample; a biological replicate is not");
        }

        [Test]
        public void NumberingThatContinuesAcrossArmsIsRenumberedWithinEachAndSaysSo()
        {
            var names = new[] { "WT_1", "WT_2", "WT_3", "KO_4", "KO_5", "KO_6" };

            var s = SdrfFileNamePattern.Read(names);

            Assert.That(s.Found, Is.True, s.NoStructureReason);
            Assert.That(s.Files.Single(f => f.FileName == "KO_4").Replicate, Is.EqualTo(1));
            Assert.That(Slot(s, SdrfFileNameRole.Replicate).Evidence, Does.Contain("renumbered"));
        }

        [Test]
        public void ADateThatSplitsTheFilesIsABatchAndNotPartOfTheSample()
        {
            var names = new List<string>();
            var dates = new[] { "20210304", "20210611" };
            for (int rep = 1; rep <= dates.Length; rep++)
                foreach (string arm in new[] { "WT", "KO" })
                    names.Add($"{dates[rep - 1]}_{arm}_{rep}.raw");

            var s = SdrfFileNamePattern.Read(names);

            Assert.That(s.Found, Is.True, s.NoStructureReason);
            Assert.That(Slot(s, SdrfFileNameRole.Batch).Levels, Is.EqualTo(dates));
            var ko2 = s.Files.Single(f => f.FileName == "20210611_KO_2.raw");
            Assert.That(ko2.Batch, Is.EqualTo("20210611"));
            Assert.That(ko2.SampleKey, Does.Not.Contain("20210611"), "the batch is not part of the sample");
        }

        /// <summary>
        /// The same arm in two batches with nothing else varying is either one sample run twice or two
        /// samples; the names cannot say which, so neither is guessed.
        /// </summary>
        [Test]
        public void TheSameArmInTwoBatchesWithNothingElseIsNotGuessed()
        {
            var names = new[] { "20210304_WT", "20210304_KO", "20210611_WT", "20210611_KO" };

            var s = SdrfFileNamePattern.Read(names);

            Assert.That(s.Found, Is.False);
            Assert.That(s.NoStructureReason, Does.Contain("same design cell"));
        }

        [Test]
        public void ADateOnEveryFileIsAcquisitionNotDesignAndIsIgnored()
        {
            var names = new[] { "20210301_WT_1", "20210302_WT_2", "20210303_KO_1", "20210304_KO_2" };

            var s = SdrfFileNamePattern.Read(names);

            Assert.That(s.Found, Is.True, s.NoStructureReason);
            Assert.That(s.Slots.Any(x => x.Role == SdrfFileNameRole.Batch), Is.False);
            Assert.That(s.Files.All(f => f.Batch == null));
        }

        // ---- The refusals: each of these is an over-split the plex audit warned about. ----

        [TestCase(new[] { "HeLa_A", "HeLa_B", "HeLa_C", "HeLa_D" }, "identifier", TestName = "A word that names every file differently is an identifier")]
        [TestCase(new[] { "QE_run001", "QE_run002", "QE_run003", "QE_run004" }, "run number", TestName = "A run number alone is not a replicate")]
        [TestCase(new[] { "WT_1", "KO_2", "WT_3", "KO_4", "WT_5", "KO_6" }, "contiguously", TestName = "Interleaved run order is not a replicate count")]
        [TestCase(new[] { "A_1_1", "A_1_2", "A_2_1", "A_2_2", "B_1_1", "B_1_2", "B_2_1", "B_2_2" }, "two numbered parts", TestName = "Two unnamed indices cannot be told apart")]
        [TestCase(new[] { "WT_rep1", "KO_replicate_1_extra" }, "shape", TestName = "Names of different shapes cannot be compared")]
        [TestCase(new[] { "only.raw" }, "one file", TestName = "One file has no siblings")]
        [TestCase(new[] { "run.raw", "run.mzML" }, "extension", TestName = "Names differing only in extension")]
        [TestCase(new[] { "20210301_x", "20210302_x" }, "acquisition date", TestName = "A date alone is not a design")]
        public void WhenTheEvidenceDoesNotDecideTheAnswerIsNoStructure(string[] names, string reasonContains)
        {
            var s = SdrfFileNamePattern.Read(names);

            Assert.That(s.Found, Is.False);
            Assert.That(s.NoStructureReason, Does.Contain(reasonContains));
            Assert.That(s.Slots, Is.Empty);
            Assert.That(s.Files.Select(f => f.SampleKey), Is.EqualTo(names), "with no structure, every file is its own sample");
        }

        [TestCase("a/b/Sample_1.raw", "Sample_1")]
        [TestCase("X.mzML.gz", "X")]
        [TestCase("X.raw.gz", "X")]
        [TestCase("X.d", "X")]
        [TestCase("X.wiff2", "X")]
        [TestCase("X.tsv", "X.tsv")]
        public void TheStemDropsTheFolderAndDataExtensionsOnly(string name, string stem)
        {
            Assert.That(SdrfFileNamePattern.Stem(name), Is.EqualTo(stem));
        }

        [Test]
        public void MalformedArgumentsThrowWeakEvidenceDoesNot()
        {
            Assert.Throws<ArgumentNullException>(() => SdrfFileNamePattern.Read(null!));
            Assert.Throws<ArgumentException>(() => SdrfFileNamePattern.Read(new[] { "a.raw", " " }));
            var e = Assert.Throws<ArgumentException>(() => SdrfFileNamePattern.Read(new[] { "a.raw", "A.RAW" }));
            Assert.That(e!.Message, Does.Contain("listed 2 times"));
        }
    }
}
