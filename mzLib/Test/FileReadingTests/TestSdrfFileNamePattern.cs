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
    /// reads: the nearest measurement of this kind of inference (a source-name partition, 3 right of 9)
    /// failed by over-splitting, so every case that should NOT yield structure is pinned here too.
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

        /// <summary>
        /// Fractions are matched ACROSS samples by index (FlashLFQ transfers between fractions at most
        /// one apart), so a lysate whose first band was never uploaded must keep its bands' own numbers.
        /// Renumbering it from 1 would line its band 2 up with the other lysate's band 1 (MAP-34).
        /// </summary>
        [Test]
        public void AMissingBandLeavesAGapAndDoesNotShiftTheOthers()
        {
            var names = Names("Lysate{0}_Band{1:00}.raw", ("2", 10));
            names.AddRange(Enumerable.Range(2, 9).Select(i => $"Lysate1_Band{i:00}.raw"));
            names.Remove("Lysate2_Band05.raw");

            var s = SdrfFileNamePattern.Read(names);

            Assert.That(s.Found, Is.True, s.NoStructureReason);
            Assert.That(s.Files.Single(f => f.FileName == "Lysate1_Band02.raw").Fraction, Is.EqualTo(2));
            Assert.That(s.Files.Single(f => f.FileName == "Lysate2_Band06.raw").Fraction, Is.EqualTo(6));
            Assert.That(Slot(s, SdrfFileNameRole.Fraction).Evidence, Does.Not.Contain("renumbered"));
        }

        [Test]
        public void BandsCountedFromZeroAreShiftedTogetherNotPerSample()
        {
            var names = Enumerable.Range(0, 3).Select(i => $"Lysate1_Band{i:00}")
                .Concat(Enumerable.Range(1, 2).Select(i => $"Lysate2_Band{i:00}")).ToList();

            var s = SdrfFileNamePattern.Read(names);

            Assert.That(s.Found, Is.True, s.NoStructureReason);
            Assert.That(s.Files.Single(f => f.FileName == "Lysate1_Band00").Fraction, Is.EqualTo(1));
            Assert.That(s.Files.Single(f => f.FileName == "Lysate2_Band01").Fraction, Is.EqualTo(2),
                "one offset for the whole deposit, so band 01 is fraction 2 in both lysates");
            Assert.That(Slot(s, SdrfFileNameRole.Fraction).Evidence, Does.Contain("shifted by 1"));
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
        /// A six-digit number that parses as yyMMdd is read as a date. That is the common case in PRIDE
        /// (12,358 such tokens in raw-file names across 167 cached deposits), and a sample ID that happens
        /// to parse would be misread, so the evidence names the format for a curator to check.
        /// </summary>
        [Test]
        public void ASixDigitDateSaysItWasReadAsYyMMdd()
        {
            var names = new[] { "210304_WT_1", "210304_KO_1", "210611_WT_2", "210611_KO_2" };

            var s = SdrfFileNamePattern.Read(names);

            Assert.That(s.Found, Is.True, s.NoStructureReason);
            Assert.That(Slot(s, SdrfFileNameRole.Batch).Evidence, Does.Contain("yyMMdd"));
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

        // ---- Pilot 1 (blind grading of 60 deposits, 2026-09-23): shapes a human reads at a glance. ----

        [Test]
        public void AReinjectionMarkerAfterTheNumberIsATechnicalReplicateOfTheSameSample()
        {
            var names = new[] { "NEG", "POS" }.SelectMany(arm => Enumerable.Range(1, 3)
                .SelectMany(i => new[] { $"{arm}{i}.raw", $"{arm}{i}rep.raw" })).ToList();

            var s = SdrfFileNamePattern.Read(names);

            Assert.That(s.Found, Is.True, s.NoStructureReason);
            var first = s.Files.Single(f => f.FileName == "NEG2.raw");
            var again = s.Files.Single(f => f.FileName == "NEG2rep.raw");
            Assert.That(again.SampleKey, Is.EqualTo(first.SampleKey), "a re-injection is the same sample");
            Assert.That((first.TechnicalReplicate, again.TechnicalReplicate), Is.EqualTo(((int?)1, (int?)2)));
            Assert.That(s.Files.Select(f => f.SampleKey).Distinct().Count(), Is.EqualTo(6));
            Assert.That(s.Slots.Single(x => x.Role == SdrfFileNameRole.Factor).Levels, Is.EqualTo(new[] { "NEG", "POS" }));
        }

        [TestCase("AF_17_Control", "AF_17_Replicate_Control")]
        [TestCase("HP_C10", "HP_C10_rr")]
        [TestCase("Cat_1_long", "Cat_1_long2")]
        [TestCase("Run1_0_13C", "Run1_0_13C_2")]
        public void AReinjectionMarkerIsReadWhereverItStandsWhenTheUnmarkedNameExists(string first, string again)
        {
            var names = new[] { first, again, first.Replace("1", "9"), again.Replace("1", "9") }.Distinct().ToList();

            var s = SdrfFileNamePattern.Read(names);

            Assert.That(s.Found, Is.True, s.NoStructureReason);
            var a = s.Files.Single(f => f.FileName == first);
            var b = s.Files.Single(f => f.FileName == again);
            Assert.That(b.SampleKey, Is.EqualTo(a.SampleKey));
            Assert.That(b.TechnicalReplicate, Is.EqualTo(2));
        }

        [Test]
        public void AMarkerWithNoUnmarkedTwinIsNotAReinjection()
        {
            var names = new[] { "WT_rep1", "WT_rep2", "KO_rep1", "KO_rep2" };

            var s = SdrfFileNamePattern.Read(names);

            Assert.That(s.Files.All(f => f.TechnicalReplicate == null), "rep1..2 with no bare twin is a replicate count, not a re-injection");
            Assert.That(s.Files.Single(f => f.FileName == "KO_rep2").Replicate, Is.EqualTo(2));
        }

        [Test]
        public void ADepositMixingFamiliesReadsEachFamilyOnItsOwn()
        {
            var names = new List<string> { "Blank.raw", "QC_standard_mix.raw" };
            names.AddRange(Names("{0}_{1}.raw", ("WT", 3), ("KO", 3)));

            var s = SdrfFileNamePattern.Read(names);

            Assert.That(s.Found, Is.True, s.NoStructureReason);
            Assert.That(s.Files.Single(f => f.FileName == "KO_2.raw").FactorLevels, Is.EqualTo(new[] { "KO" }));
            var blank = s.Files.Single(f => f.FileName == "Blank.raw");
            Assert.That(blank.FactorLevels, Is.Empty);
            Assert.That(s.Files.Count(f => f.SampleKey == blank.SampleKey), Is.EqualTo(1), "a lone file is its own sample");
        }

        [Test]
        public void EveryReadingSaysWhichFamilyItWasReadIn()
        {
            var names = new List<string> { "Blank.raw" };
            names.AddRange(Names("{0}_{1}.raw", ("WT", 2), ("KO", 2)));
            names.AddRange(Names("20180222_ZJ_{0}_{1}.raw", ("MG", 2)));

            var s = SdrfFileNamePattern.Read(names);

            var wt = s.Files.Single(f => f.FileName == "WT_1.raw");
            var mg = s.Files.Single(f => f.FileName == "20180222_ZJ_MG_1.raw");
            Assert.That(wt.Family, Is.GreaterThanOrEqualTo(0));
            Assert.That(mg.Family, Is.Not.EqualTo(wt.Family), "a differently shaped name is another family");
            Assert.That(s.Files.Single(f => f.FileName == "Blank.raw").Family, Is.EqualTo(-1), "a lone run belongs to no family");
            Assert.That(s.Files.Where(f => f.FileName.StartsWith("WT") || f.FileName.StartsWith("KO")).Select(f => f.Family).Distinct().Count(), Is.EqualTo(1));
        }

        [Test]
        public void ASidecarFileFoldsIntoItsRun()
        {
            var names = Names("{0}_{1}.wiff", ("WT", 2), ("KO", 2));
            names.AddRange(names.Select(n => n + ".scan").ToList());

            var s = SdrfFileNamePattern.Read(names);

            Assert.That(s.Found, Is.True, s.NoStructureReason);
            var run = s.Files.Single(f => f.FileName == "KO_2.wiff");
            var scan = s.Files.Single(f => f.FileName == "KO_2.wiff.scan");
            Assert.That(scan with { FileName = run.FileName }, Is.EqualTo(run));
        }

        [Test]
        public void ADashSeparatedDateIsOneDate()
        {
            var names = new[] { "2018-02-28_WT_1", "2018-02-28_WT_2", "2018-03-05_KO_1", "2018-03-05_KO_2" };

            var s = SdrfFileNamePattern.Read(names);

            Assert.That(s.Found, Is.True, s.NoStructureReason);
            Assert.That(s.Files.Single(f => f.FileName == "2018-03-05_KO_2").Replicate, Is.EqualTo(2));
            Assert.That(s.Slots.Single(x => x.Role == SdrfFileNameRole.Factor).Levels, Is.EqualTo(new[] { "KO", "WT" }));
        }

        [Test]
        public void ANumberThatDoesNotCountIsACategoryWhenEachValueIsShared()
        {
            var names = new[] { "0", "5", "25" }.SelectMany(c => new[] { $"Bsub_{c}_R1", $"Bsub_{c}_R2" }).ToList();

            var s = SdrfFileNamePattern.Read(names);

            Assert.That(s.Found, Is.True, s.NoStructureReason);
            Assert.That(s.Slots.Single(x => x.Role == SdrfFileNameRole.Factor).Levels, Is.EquivalentTo(new[] { "0", "5", "25" }));
            Assert.That(s.Files.Single(f => f.FileName == "Bsub_25_R2").Replicate, Is.EqualTo(2));
        }

        [Test]
        public void InterleavedRunNumbersIdentifySamplesButAreNotAReplicateCount()
        {
            var names = new[] { "WT_1", "KO_2", "WT_3", "KO_4", "WT_5", "KO_6" };

            var s = SdrfFileNamePattern.Read(names);

            Assert.That(s.Found, Is.True, s.NoStructureReason);
            Assert.That(s.Slots.Single(x => x.Role == SdrfFileNameRole.Factor).Levels, Is.EqualTo(new[] { "KO", "WT" }));
            Assert.That(s.Files.All(f => f.Replicate == null && f.BiologicalReplicate == null), "1,3,5 does not count replicates");
            Assert.That(s.Files.Select(f => f.SampleKey).Distinct().Count(), Is.EqualTo(6));
        }

        /// <summary>A number that is alone in every group passes "contiguous" trivially; 1..1 counts nothing.</summary>
        [Test]
        public void ANumberAloneInEveryGroupIsNotAReplicateCount()
        {
            var names = new[] { "WT_a_3", "WT_b_5", "KO_a_7", "KO_b_9" };

            var s = SdrfFileNamePattern.Read(names);

            Assert.That(s.Slots.Any(x => x.Role == SdrfFileNameRole.Replicate), Is.False);
            Assert.That(s.Files.All(f => f.Replicate == null));
        }

        [Test]
        public void OfTwoUnnamedNumbersTheLastCountsAndTheSharedOneIsACategory()
        {
            var names = new[] { "A_1_1", "A_1_2", "A_2_1", "A_2_2", "B_1_1", "B_1_2", "B_2_1", "B_2_2" };

            var s = SdrfFileNamePattern.Read(names);

            Assert.That(s.Found, Is.True, s.NoStructureReason);
            Assert.That(s.Slots.Count(x => x.Role == SdrfFileNameRole.Factor), Is.EqualTo(2));
            Assert.That(s.Files.Single(f => f.FileName == "B_2_1").Replicate, Is.EqualTo(1));
        }

        /// <summary>
        /// A number that is the only part to vary is refused when unnamed (a run number looks the same),
        /// but a word before it that states what it counts settles that: one sample injected twice.
        /// </summary>
        [Test]
        public void ANamedTechnicalReplicateAloneIsOneSampleInjectedSeveralTimes()
        {
            var names = new[] { "Sample_inj1.raw", "Sample_inj2.raw", "Sample_inj3.raw" };

            var s = SdrfFileNamePattern.Read(names);

            Assert.That(s.Found, Is.True, s.NoStructureReason);
            Assert.That(s.Files.Select(f => f.TechnicalReplicate), Is.EqualTo(new int?[] { 1, 2, 3 }));
            Assert.That(s.Files.Select(f => f.SampleKey).Distinct().Single(), Is.Empty, "the whole deposit is one sample");
            Assert.That(Slot(s, SdrfFileNameRole.TechnicalReplicate).Evidence, Does.Contain("'inj'"));
        }

        [Test]
        public void ANamedBiologicalReplicateAloneIsOneSamplePerReplicate()
        {
            var names = new[] { "HeLa_BR1.raw", "HeLa_BR2.raw", "HeLa_BR3.raw" };

            var s = SdrfFileNamePattern.Read(names);

            Assert.That(s.Found, Is.True, s.NoStructureReason);
            Assert.That(s.Files.Select(f => f.BiologicalReplicate), Is.EqualTo(new int?[] { 1, 2, 3 }));
            Assert.That(s.Files.Select(f => f.SampleKey).Distinct().Count(), Is.EqualTo(3));
        }

        [Test]
        public void NumericLevelsAreListedInNumericOrder()
        {
            var names = Enumerable.Range(1, 12).SelectMany(c => new[] { $"Bsub_{c}_R1", $"Bsub_{c}_R2" }).ToList();

            var s = SdrfFileNamePattern.Read(names);

            Assert.That(s.Found, Is.True, s.NoStructureReason);
            Assert.That(s.Slots.Single(x => x.Role == SdrfFileNameRole.Factor).Levels,
                Is.EqualTo(Enumerable.Range(1, 12).Select(i => i.ToString())));
        }

        // ---- The refusals: each of these would be an over-split. ----

        [TestCase(new[] { "HeLa_A", "HeLa_B", "HeLa_C", "HeLa_D" }, "identifier", TestName = "A word that names every file differently is an identifier")]
        [TestCase(new[] { "QE_run001", "QE_run002", "QE_run003", "QE_run004" }, "run number", TestName = "A run number alone is not a replicate")]
        [TestCase(new[] { "WT_rep1", "KO_replicate_1_extra" }, "family", TestName = "Names of different shapes cannot be compared")]
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
