using System;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class DigestionParameterTests
    {
        private static Stopwatch Stopwatch { get; set; }

        [SetUp]
        public static void Setup()
        {
            Stopwatch = new Stopwatch();
            Stopwatch.Start();
        }

        [TearDown]
        public static void TearDown()
        {
            Console.WriteLine($"Analysis time: {Stopwatch.Elapsed.Hours}h {Stopwatch.Elapsed.Minutes}m {Stopwatch.Elapsed.Seconds}s");
        }

        [Test]
        public static void TestDigestionParamsClone()
        {
            DigestionParams digestionParams = new DigestionParams(
                protease: "top-down",
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain,
                keepNGlycopeptide: true,
                keepOGlycopeptide: true);

            DigestionParams digestionParamsClone = (DigestionParams)digestionParams.Clone();
            Assert.AreEqual(digestionParams, digestionParamsClone);
            Assert.AreEqual(digestionParams.InitiatorMethionineBehavior, digestionParamsClone.InitiatorMethionineBehavior);
            Assert.AreEqual(digestionParams.MaxMissedCleavages, digestionParamsClone.MaxMissedCleavages);
            Assert.AreEqual(digestionParams.MaxModificationIsoforms, digestionParamsClone.MaxModificationIsoforms);
            Assert.AreEqual(digestionParams.MinLength, digestionParamsClone.MinLength);
            Assert.AreEqual(digestionParams.MaxLength, digestionParamsClone.MaxLength);
            Assert.AreEqual(digestionParams.MaxMods, digestionParamsClone.MaxMods);
            Assert.AreEqual(digestionParams.Protease, digestionParamsClone.Protease);
            Assert.AreEqual(digestionParams.SearchModeType, digestionParamsClone.SearchModeType);
            Assert.AreEqual(digestionParams.FragmentationTerminus, digestionParamsClone.FragmentationTerminus);
            Assert.AreEqual(digestionParams.GeneratehUnlabeledProteinsForSilac, digestionParamsClone.GeneratehUnlabeledProteinsForSilac);
            Assert.AreEqual(digestionParams.KeepNGlycopeptide, digestionParamsClone.KeepNGlycopeptide);
            Assert.AreEqual(digestionParams.KeepOGlycopeptide, digestionParamsClone.KeepOGlycopeptide);
            Assert.AreEqual(digestionParams.SpecificProtease, digestionParamsClone.SpecificProtease);
            NUnit.Framework.Assert.That(!ReferenceEquals(digestionParams, digestionParamsClone));

            digestionParams = new DigestionParams(
                protease: "top-down",
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain,
                keepNGlycopeptide: true,
                keepOGlycopeptide: true,
                maxModificationIsoforms: 5,
                maxModsForPeptides: 6,
                maxPeptideLength: 7,
                searchModeType: CleavageSpecificity.None,
                fragmentationTerminus: FragmentationTerminus.C,
                generateUnlabeledProteinsForSilac: false);

            digestionParamsClone = (DigestionParams)digestionParams.Clone();
            Assert.AreEqual(digestionParams, digestionParamsClone);
            Assert.AreEqual(digestionParams.InitiatorMethionineBehavior, digestionParamsClone.InitiatorMethionineBehavior);
            Assert.AreEqual(digestionParams.MaxMissedCleavages, digestionParamsClone.MaxMissedCleavages);
            Assert.AreEqual(digestionParams.MaxModificationIsoforms, digestionParamsClone.MaxModificationIsoforms);
            Assert.AreEqual(digestionParams.MinLength, digestionParamsClone.MinLength);
            Assert.AreEqual(digestionParams.MaxLength, digestionParamsClone.MaxLength);
            Assert.AreEqual(digestionParams.MaxMods, digestionParamsClone.MaxMods);
            Assert.AreEqual(digestionParams.Protease, digestionParamsClone.Protease);
            Assert.AreEqual(digestionParams.SearchModeType, digestionParamsClone.SearchModeType);
            Assert.AreEqual(digestionParams.FragmentationTerminus, digestionParamsClone.FragmentationTerminus);
            Assert.AreEqual(digestionParams.GeneratehUnlabeledProteinsForSilac, digestionParamsClone.GeneratehUnlabeledProteinsForSilac);
            Assert.AreEqual(digestionParams.KeepNGlycopeptide, digestionParamsClone.KeepNGlycopeptide);
            Assert.AreEqual(digestionParams.KeepOGlycopeptide, digestionParamsClone.KeepOGlycopeptide);
            Assert.AreEqual(digestionParams.SpecificProtease, digestionParamsClone.SpecificProtease);
            NUnit.Framework.Assert.That(!ReferenceEquals(digestionParams, digestionParamsClone));
        }

        [Test]
        public static void TestDigestionParamsCloneWithNewTerminus()
        {
            DigestionParams digestionParams = new DigestionParams(
                protease: "top-down",
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain,
                keepNGlycopeptide: true,
                keepOGlycopeptide: true);

            DigestionParams digestionParamsClone = (DigestionParams)digestionParams.Clone(FragmentationTerminus.N);
            Assert.AreNotEqual(digestionParams, digestionParamsClone);
            Assert.AreEqual(digestionParams.InitiatorMethionineBehavior, digestionParamsClone.InitiatorMethionineBehavior);
            Assert.AreEqual(digestionParams.MaxMissedCleavages, digestionParamsClone.MaxMissedCleavages);
            Assert.AreEqual(digestionParams.MaxModificationIsoforms, digestionParamsClone.MaxModificationIsoforms);
            Assert.AreEqual(digestionParams.MinLength, digestionParamsClone.MinLength);
            Assert.AreEqual(digestionParams.MaxLength, digestionParamsClone.MaxLength);
            Assert.AreEqual(digestionParams.MaxMods, digestionParamsClone.MaxMods);
            Assert.AreEqual(digestionParams.Protease, digestionParamsClone.Protease);
            Assert.AreEqual(digestionParams.SearchModeType, digestionParamsClone.SearchModeType);
            Assert.AreEqual(FragmentationTerminus.N, digestionParamsClone.FragmentationTerminus);
            Assert.AreEqual(digestionParams.GeneratehUnlabeledProteinsForSilac, digestionParamsClone.GeneratehUnlabeledProteinsForSilac);
            Assert.AreEqual(digestionParams.KeepNGlycopeptide, digestionParamsClone.KeepNGlycopeptide);
            Assert.AreEqual(digestionParams.KeepOGlycopeptide, digestionParamsClone.KeepOGlycopeptide);
            Assert.AreEqual(digestionParams.SpecificProtease, digestionParamsClone.SpecificProtease);
            NUnit.Framework.Assert.That(!ReferenceEquals(digestionParams, digestionParamsClone));

            digestionParams = new DigestionParams(
                protease: "top-down",
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain,
                keepNGlycopeptide: true,
                keepOGlycopeptide: true,
                maxModificationIsoforms: 5,
                maxModsForPeptides: 6,
                maxPeptideLength: 7,
                searchModeType: CleavageSpecificity.None,
                fragmentationTerminus: FragmentationTerminus.None,
                generateUnlabeledProteinsForSilac: false);

            digestionParamsClone = (DigestionParams)digestionParams.Clone(FragmentationTerminus.N);
            Assert.AreNotEqual(digestionParams, digestionParamsClone);
            Assert.AreEqual(digestionParams.InitiatorMethionineBehavior, digestionParamsClone.InitiatorMethionineBehavior);
            Assert.AreEqual(digestionParams.MaxMissedCleavages, digestionParamsClone.MaxMissedCleavages);
            Assert.AreEqual(digestionParams.MaxModificationIsoforms, digestionParamsClone.MaxModificationIsoforms);
            Assert.AreEqual(digestionParams.MinLength, digestionParamsClone.MinLength);
            Assert.AreEqual(digestionParams.MaxLength, digestionParamsClone.MaxLength);
            Assert.AreEqual(digestionParams.MaxMods, digestionParamsClone.MaxMods);
            Assert.AreEqual(ProteaseDictionary.Dictionary["singleN"], digestionParamsClone.Protease);
            Assert.AreEqual(digestionParams.SearchModeType, digestionParamsClone.SearchModeType);
            Assert.AreEqual(FragmentationTerminus.N, digestionParamsClone.FragmentationTerminus);
            Assert.AreEqual(digestionParams.GeneratehUnlabeledProteinsForSilac, digestionParamsClone.GeneratehUnlabeledProteinsForSilac);
            Assert.AreEqual(digestionParams.KeepNGlycopeptide, digestionParamsClone.KeepNGlycopeptide);
            Assert.AreEqual(digestionParams.KeepOGlycopeptide, digestionParamsClone.KeepOGlycopeptide);
            Assert.AreEqual(digestionParams.SpecificProtease, digestionParamsClone.SpecificProtease);
            NUnit.Framework.Assert.That(!ReferenceEquals(digestionParams, digestionParamsClone));

            digestionParamsClone = (DigestionParams)digestionParams.Clone(FragmentationTerminus.C);
            Assert.AreNotEqual(digestionParams, digestionParamsClone);
            Assert.AreEqual(digestionParams.InitiatorMethionineBehavior, digestionParamsClone.InitiatorMethionineBehavior);
            Assert.AreEqual(digestionParams.MaxMissedCleavages, digestionParamsClone.MaxMissedCleavages);
            Assert.AreEqual(digestionParams.MaxModificationIsoforms, digestionParamsClone.MaxModificationIsoforms);
            Assert.AreEqual(digestionParams.MinLength, digestionParamsClone.MinLength);
            Assert.AreEqual(digestionParams.MaxLength, digestionParamsClone.MaxLength);
            Assert.AreEqual(digestionParams.MaxMods, digestionParamsClone.MaxMods);
            Assert.AreEqual(ProteaseDictionary.Dictionary["singleC"], digestionParamsClone.Protease);
            Assert.AreEqual(digestionParams.SearchModeType, digestionParamsClone.SearchModeType);
            Assert.AreEqual(FragmentationTerminus.C, digestionParamsClone.FragmentationTerminus);
            Assert.AreEqual(digestionParams.GeneratehUnlabeledProteinsForSilac, digestionParamsClone.GeneratehUnlabeledProteinsForSilac);
            Assert.AreEqual(digestionParams.KeepNGlycopeptide, digestionParamsClone.KeepNGlycopeptide);
            Assert.AreEqual(digestionParams.KeepOGlycopeptide, digestionParamsClone.KeepOGlycopeptide);
            Assert.AreEqual(digestionParams.SpecificProtease, digestionParamsClone.SpecificProtease);
            NUnit.Framework.Assert.That(!ReferenceEquals(digestionParams, digestionParamsClone));
        }
    }
}
