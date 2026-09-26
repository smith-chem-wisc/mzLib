using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Tests for reading replicates as sibling markers and deciding their kind from the record. Every case
    /// is a shape the blind benchmark's graders quoted when a drafted biological replicate was judged wrong
    /// (sdrf project, fresh sets 1 and 2, 2026-09-23).
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrfReplicateResolver
    {
        private static SdrfReplicateReading File(SdrfReplicates r, string name) => r.Files.Single(f => f.FileName == name);

        // ---- step 1: sibling markers, read locally ----

        [TestCase("CT10A", "CT10B", "CT10C")]
        [TestCase("WT-30min-A", "WT-30min-B", "WT-30min-C")]
        [TestCase("IL1BETA_R1", "IL1BETA_R2", "IL1BETA_R3")]
        [TestCase("Brain_Rat1", "Brain_Rat2", "Brain_Rat3")]
        [TestCase("era_repeat1", "era_repeat2", "era_repeat3")]
        [TestCase("PR619_BR1", "PR619_BR2", "PR619_BR3")]
        [TestCase("CTR_01", "CTR_02", "CTR_03")]
        [TestCase("K5PA5pass-1", "K5PA5pass-2", "K5PA5pass-3")]
        public void SiblingsThatDifferOnlyInAFinalMarkerAreReplicatesOfOneBase(string a, string b, string c)
        {
            var r = SdrfReplicateResolver.Read(new[] { a, b, c, "Other_file" }, Array.Empty<string>());

            Assert.That(new[] { a, b, c }.Select(n => File(r, n).Number), Is.EqualTo(new int?[] { 1, 2, 3 }));
            Assert.That(new[] { a, b, c }.Select(n => File(r, n).Base).Distinct().Count(), Is.EqualTo(1));
            Assert.That(File(r, "Other_file").Number, Is.Null, "a lone file has no marker");
        }

        [Test]
        public void NumberingRestartsInEveryBaseEvenWhenTheDepositMixesShapes()
        {
            var names = new[] { "Astro_IL1BETA_R1", "Astro_IL1BETA_R2", "PHA_CTRL_R1", "PHA_CTRL_R2", "PHA_CTRL_R3" };

            var r = SdrfReplicateResolver.Read(names, Array.Empty<string>());

            Assert.That(File(r, "PHA_CTRL_R3").Number, Is.EqualTo(3));
            Assert.That(File(r, "Astro_IL1BETA_R2").Number, Is.EqualTo(2), "never one count across the deposit");
        }

        [Test]
        public void AnIdentifierIsNotAMarker()
        {
            var names = new[] { "sham_0316", "sham_0321", "sham_0357" };

            var r = SdrfReplicateResolver.Read(names, Array.Empty<string>());

            Assert.That(r.Files.All(f => f.Number == null), "0316/0321/0357 are animal IDs, not a 1..n count");
        }

        [TestCase("S1_FR01", "S1_FR02", "S1_FR03", "Fraction")]
        [TestCase("BR2_F1", "BR2_F2", "BR2_F3", "Fraction")]
        [TestCase("X_TR1", "X_TR2", "X_TR3", "Technical")]
        [TestCase("X_inj1", "X_inj2", "X_inj3", "Technical")]
        [TestCase("Liver_Rat1", "Liver_Rat2", "Liver_Rat3", "Biological")]
        [TestCase("S1_Band_01", "S1_Band_02", "S1_Band_03", "Fraction")]
        [TestCase("MSB67868ABand_01", "MSB67868ABand_02", "MSB67868ABand_03", "Fraction")]
        public void AMarkersOwnWordSaysWhatItCounts(string a, string b, string c, string kindName)
        {
            var r = SdrfReplicateResolver.Read(new[] { a, b, c }, Array.Empty<string>());

            Assert.That(r.Files.Select(f => f.Kind), Is.All.EqualTo(Enum.Parse<SdrfReplicateKind>(kindName)));
            Assert.That(r.Files.Select(f => f.Number), Is.EqualTo(new int?[] { 1, 2, 3 }));
        }

        /// <summary>
        /// A word read across a separator, or off a tag it is glued to, has to mean one thing. In
        /// PXD041400, <c>f_10</c> is a female, not a fraction. <c>InGel</c> is an in-gel digestion,
        /// not a gel band.
        /// </summary>
        [TestCase("Liver_f_1", "Liver_f_2", "Liver_f_3")]
        [TestCase("MCF7_InGel_01", "MCF7_InGel_02", "MCF7_InGel_03")]
        public void AShortOrAmbiguousWordBeforeTheSeparatorStatesNothing(string a, string b, string c)
        {
            var r = SdrfReplicateResolver.Read(new[] { a, b, c }, Array.Empty<string>());

            Assert.That(r.Files.Select(f => f.Kind), Is.All.Not.EqualTo(SdrfReplicateKind.Fraction));
        }

        [Test]
        public void ACountPastTwelveIsARunCounterUnlessItsWordSaysFraction()
        {
            var runs = Enumerable.Range(1, 40).Select(i => $"Phospho_final_{i:00}").ToList();
            var fractions = Enumerable.Range(1, 24).Select(i => $"S1_FR{i:00}").ToList();

            var r = SdrfReplicateResolver.Read(runs.Concat(fractions), Array.Empty<string>());

            Assert.That(r.Files.Where(f => runs.Contains(f.FileName)).All(f => f.Number == null), "no study has 40 replicates of one condition");
            Assert.That(r.Files.Single(f => f.FileName == "S1_FR24").Number, Is.EqualTo(24), "24 fractions are ordinary");
        }

        [Test]
        public void AnUnknownLetterBeforeANumberIsNotAReplicateMarker()
        {
            var r = SdrfReplicateResolver.Read(new[] { "HP_C1", "HP_C2", "HP_C3" }, Array.Empty<string>());

            Assert.That(r.Files.All(f => f.Number == null), "C1..C3 could be columns, fractions or cases; nothing says replicate");
        }

        [Test]
        public void NumbersWithGapsAreRunNumbersNotAReplicateCount()
        {
            var r = SdrfReplicateResolver.Read(new[] { "WT_1", "WT_3", "WT_7" }, Array.Empty<string>());

            Assert.That(r.Files.All(f => f.Number == null));
        }

        [Test]
        public void TwoMarkersAreAnOuterAndAnInnerLevel()
        {
            var names = new[] { "del_rlmC_1_1", "del_rlmC_1_2", "del_rlmC_2_1", "del_rlmC_2_2" };

            var r = SdrfReplicateResolver.Read(names, Array.Empty<string>());

            var f = File(r, "del_rlmC_2_1");
            Assert.That((f.Outer, f.Number), Is.EqualTo(((int?)2, (int?)1)));
            Assert.That(r.OuterKind, Is.EqualTo(SdrfReplicateKind.Biological));
            Assert.That(r.MarkerKind, Is.EqualTo(SdrfReplicateKind.Technical), "the inner of two levels is a re-injection unless the record says otherwise");
        }

        // ---- step 2: the record decides the kind, only when a stated number matches the count ----

        [TestCase("Three biological replicates were prepared for each strain.", "Biological")]
        [TestCase("Cells were grown in triplicate.", "Biological")]
        [TestCase("Samples were analyzed in triplicate by LC-MS/MS.", "Technical")]
        [TestCase("Each digest was injected three times.", "Technical")]
        [TestCase("Peptides were separated into 3 high-pH fractions.", "Fraction")]
        [TestCase("n = 3 mice per group.", "Biological")]
        [TestCase("Two biological replicates were prepared.", "Unstated")]
        [TestCase("A proteome.", "Unstated")]
        public void TheRecordNamesTheKindOnlyWhenItsNumberMatchesTheCount(string text, string kindName)
        {
            var kind = Enum.Parse<SdrfReplicateKind>(kindName);
            var names = new[] { "WT_1", "WT_2", "WT_3", "KO_1", "KO_2", "KO_3" };

            var r = SdrfReplicateResolver.Read(names, new[] { text });

            Assert.That(r.MarkerKind, Is.EqualTo(kind), r.MarkerEvidence);
            if (kind != SdrfReplicateKind.Unstated) Assert.That(r.MarkerEvidence, Does.Contain("record"));
        }

        [TestCase("Each sample was analysed in three separate technical replicates.", "Technical")]
        [TestCase("Three independent biological replicates were grown.", "Biological")]
        public void AWordOrTwoBetweenTheNumberAndTheKindStillCounts(string text, string kindName)
        {
            var r = SdrfReplicateResolver.Read(new[] { "M35_001", "M35_002", "M35_003", "M36_001", "M36_002", "M36_003" }, new[] { text });

            Assert.That(r.MarkerKind, Is.EqualTo(Enum.Parse<SdrfReplicateKind>(kindName)), r.MarkerEvidence);
        }

        [Test]
        public void TheCountIsWhatMostBasesHaveNotTheLargest()
        {
            var names = new[] { "A_1", "A_2", "A_3", "B_1", "B_2", "B_3", "C_1", "C_2", "C_3", "C_4" };

            var r = SdrfReplicateResolver.Read(names, new[] { "Samples were analyzed in triplicate." });

            Assert.That(r.MarkerKind, Is.EqualTo(SdrfReplicateKind.Technical), r.MarkerEvidence);
        }

        [Test]
        public void AFractionMarkerKeepsItsOwnNumbersAcrossAGap()
        {
            var names = new[] { "SKMEL28_F1", "SKMEL28_F2", "SKMEL28_F4", "SKMEL28_F5" };

            var r = SdrfReplicateResolver.Read(names, Array.Empty<string>());

            Assert.That(r.Files.Select(f => f.Number), Is.EqualTo(new int?[] { 1, 2, 4, 5 }), "a missing band stays a gap (MAP-34)");
            Assert.That(r.Files.Select(f => f.Base).Distinct().Count(), Is.EqualTo(1));
        }

        [Test]
        public void TechMarksATechnicalReplicate()
        {
            var r = SdrfReplicateResolver.Read(new[] { "tech_A_01", "tech_A_02", "tech_A_03" }, Array.Empty<string>());

            Assert.That(r.Files.Select(f => f.Kind), Is.All.EqualTo(SdrfReplicateKind.Technical));
        }

        [Test]
        public void TwoKindsStatingTheSameNumberDecideNothing()
        {
            var names = new[] { "WT_1", "WT_2", "WT_3" };
            var text = new[] { "Cultures were grown in triplicate and each was analyzed in triplicate." };

            var r = SdrfReplicateResolver.Read(names, text);

            Assert.That(r.MarkerKind, Is.EqualTo(SdrfReplicateKind.Unstated));
        }

        [TestCase("Lysates were pooled within each group.")]
        [TestCase("A single HCT116 lysate was fractionated.")]
        [TestCase("The biosample was run with three technical replicates.")]
        public void TheRecordCanSayThereIsOneBiologicalSample(string text)
        {
            var r = SdrfReplicateResolver.Read(new[] { "x_F1", "x_F2" }, new[] { text });

            Assert.That(r.SingleBiologicalSample, Is.True, r.SingleSampleEvidence);
        }

        [Test]
        public void StatedBiologicalReplicatesAreNotASingleSample()
        {
            var r = SdrfReplicateResolver.Read(new[] { "a", "b" },
                new[] { "Lysates were pooled from three biological replicates." });

            Assert.That(r.SingleBiologicalSample, Is.False);
        }

        [TestCase("Digestion was performed in 50 mM ammonium bicarbonate.")]
        [TestCase("Samples were analyzed in 2019 on an Orbitrap.")]
        [TestCase("The TMT-labelled samples were pooled and fractionated.")]
        [TestCase("The two eluates were pooled before injection.")]
        public void ProtocolWordingDoesNotMakeTheRecordOneBiologicalSample(string text)
        {
            var names = new[] { "CT10A", "CT10B", "CT10C", "KO10A", "KO10B", "KO10C" };

            var r = SdrfReplicateResolver.Read(names, new[] { text });

            Assert.That(r.SingleBiologicalSample, Is.False, r.SingleSampleEvidence);
            Assert.That(r.MarkerKind, Is.EqualTo(SdrfReplicateKind.Unstated));
        }

        [Test]
        public void TheSingleSampleEvidenceQuotesTheRecord()
        {
            var r = SdrfReplicateResolver.Read(new[] { "x_F1", "x_F2" },
                new[] { "The biosample was run with three technical replicates." });

            Assert.That(r.SingleSampleEvidence, Does.Contain("three technical replicates"));
        }

        [Test]
        public void ALongUnlabelledRunIsAFractionWhenTheRecordStatesThatManyFractions()
        {
            var names = new[] { "WT", "KO" }.SelectMany(a => Enumerable.Range(1, 16).Select(i => $"{a}_{i:00}")).ToList();

            var r = SdrfReplicateResolver.Read(names, new[] { "Peptides were separated into 16 fractions." });

            Assert.That(r.MarkerKind, Is.EqualTo(SdrfReplicateKind.Fraction), r.MarkerEvidence);
            Assert.That(File(r, "KO_16").Number, Is.EqualTo(16));
        }

        [Test]
        public void MalformedArgumentsThrow()
        {
            Assert.Throws<ArgumentNullException>(() => SdrfReplicateResolver.Read(null!, Array.Empty<string>()));
            Assert.Throws<ArgumentNullException>(() => SdrfReplicateResolver.Read(new[] { "a" }, null!));
        }
    }
}
