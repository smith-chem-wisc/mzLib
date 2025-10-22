using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Omics.BioPolymer;
using Omics.Modifications;

namespace Test.DatabaseTests.VariantTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class SequenceVariationTryAddTests
    {
        private static SequenceVariation MakeVariant(string refSeq, string altSeq, string vcf, Dictionary<int, List<Modification>> mods = null) =>
            new SequenceVariation(
                oneBasedPosition: 25,
                originalSequence: refSeq,
                variantSequence: altSeq,
                description: "TryAddBase",
                variantCallFormatDataString: vcf,
                oneBasedModifications: mods);

        private static Modification Mod(string id) =>
            new Modification(_originalId: id, _modificationType: "TestType");
            
        [Test]
        public void TryAdd_ReferenceNoOpCaughtThenAltAdded()
        {
            string vcf = "1\t300\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/1:5,7:12";
            var baseVar = MakeVariant("A", "T", vcf);
            var results = baseVar.SplitPerGenotype(includeReferenceForHeterozygous: true);

            Assert.That(results.Count, Is.EqualTo(1));
            Assert.That(results[0].Description.Contains("HeterozygousAlt"), Is.True);
            Assert.That(results[0].Description.Contains("HeterozygousRef"), Is.False,
                "No-op reference variant should be suppressed.");
        }

        [Test]
        public void TryAdd_HomozygousReference_NoVariantAdded()
        {
            string vcf = "1\t300\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/0:10,0:10";
            var baseVar = MakeVariant("A", "T", vcf);
            var results = baseVar.SplitPerGenotype(emitReferenceForHomozygousRef: true);
            Assert.That(results, Is.Empty);
        }
        [Test]
        public void TryAdd_InvalidTerminationModifications_Caught_ExceptionPath()
        {
            // This scenario throws during construction of the BASE SequenceVariation (before SplitPerGenotype)
            // because a termination ('*') variant forbids modifications at or after the begin position.
            // So the failure happens prior to TryAdd; we assert that here explicitly.
            var mods = new Dictionary<int, List<Modification>>
            {
                { 25, new List<Modification>{ Mod("StopMod") } }
            };
            string vcf = "1\t300\t.\tA\t*\t.\tPASS\tANN=*\tGT:AD:DP\t0/1:3,9:12";

            Assert.Throws<ArgumentException>(() =>
                    MakeVariant("A", "*", vcf, mods),
                "Expected constructor to reject termination variant with in?span modification site.");
        }
        [Test]
        public void TryAdd_MixedAltIndex_VariantAdded_WhenSkipDisabled()
        {
            string vcf = "1\t300\t.\tA\tT,G\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/2:4,0,6:10";
            var baseVar = MakeVariant("A", "T", vcf);
            var results = baseVar.SplitPerGenotype(skipIfAltIndexMismatch: false);

            Assert.That(results.Count, Is.EqualTo(1));
            Assert.That(results[0].Description.Contains("MixedAltIndex(StoredAltOnly)"), Is.True);
        }
        [Test]
        public void TryAdd_NoOpBaseVariant_RejectedByConstructor()
        {
            // A variant with identical original and variant sequences and no modifications is invalid by design.
            // The constructor should throw before any SplitPerGenotype logic executes.
            string vcf = "1\t300\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/1:5,5:10";
            Assert.Throws<ArgumentException>(
                () => MakeVariant("A", "A", vcf),
                "Expected rejection of no?op variant (OriginalSequence == VariantSequence with no modifications).");
        }

        [Test]
        public void TryAdd_ClonesModDictionary_OnSuccessfulAdd()
        {
            var mods = new Dictionary<int, List<Modification>>
            {
                { 10, new List<Modification>{ Mod("ModA"), Mod("ModB") } }
            };
            string vcf = "1\t300\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/1:6,6:12";
            var baseVar = MakeVariant("A", "T", vcf, mods);
            var results = baseVar.SplitPerGenotype();

            Assert.That(results.Count, Is.EqualTo(1));
            Assert.That(results[0].OneBasedModifications, Is.Not.Null);
            Assert.That(results[0].OneBasedModifications.ContainsKey(10), Is.True);
            Assert.That(ReferenceEquals(results[0].OneBasedModifications, baseVar.OneBasedModifications), Is.False,
                "Expected cloned modification map, not original reference.");
        }

        [Test]
        public void TryAdd_HomozygousAlt_SingleAdd()
        {
            string vcf = "1\t300\t.\tA\tT\t.\tPASS\tANN=T|.\tGT:AD:DP\t1/1:0,10:10";
            var baseVar = MakeVariant("A", "T", vcf);
            var results = baseVar.SplitPerGenotype();

            Assert.That(results.Count, Is.EqualTo(1));
            Assert.That(results[0].Description.Contains("HomozygousAlt"), Is.True);
        }

        [Test]
        public void TryAdd_ContainsDifferentAlt_SkippedWhenFlagTrue()
        {
            string vcf = "1\t300\t.\tA\tT,G\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/2:3,0,7:10";
            var baseVar = MakeVariant("A", "T", vcf);
            var results = baseVar.SplitPerGenotype();
            Assert.That(results, Is.Empty);
        }

        [Test]
        public void TryAdd_ContainsDifferentAlt_AddsWhenFlagFalse()
        {
            string vcf = "1\t300\t.\tA\tT,G\t.\tPASS\tANN=T|.\tGT:AD:DP\t0/2:3,0,7:10";
            var baseVar = MakeVariant("A", "T", vcf);
            var results = baseVar.SplitPerGenotype(skipIfAltIndexMismatch: false);
            Assert.That(results.Count, Is.EqualTo(1));
            Assert.That(results[0].Description.Contains("MixedAltIndex(StoredAltOnly)"), Is.True);
        }
    }
}