using System;
using System.Collections.Generic;
using NUnit.Framework;
using Omics.BioPolymer;
using Omics.Modifications;

namespace Test.DatabaseTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class SequenceVariationInvalidModificationTests
    {
        private static Modification CreateTestModification(char residue)
        {
            ModificationMotif.TryGetMotif(residue.ToString(), out var motif);
            return new Modification("testMod", null, "testType", null, motif, "Anywhere.", null, 0.0,
                null, null, null, null, null, null);
        }

        private static VariantCallFormat CreateTestVcf()
        {
            // Minimal valid VCF-like line (tab-delimited) for constructing VariantCallFormat
            return new VariantCallFormat("1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30");
        }

        [Test]
        public void Constructor_DeletionWithModificationInsideRemovedRegion_Throws()
        {
            // Original single residue at position 4 is deleted (variant sequence empty)
            // A modification is (incorrectly) specified at that deleted position (4)
            var mod = CreateTestModification('P');
            var vcf = CreateTestVcf();
            var mods = new Dictionary<int, List<Modification>>
            {
                { 4, new List<Modification> { mod } } // position 4 no longer exists after deletion
            };

            Assert.Throws<ArgumentException>(() =>
                new SequenceVariation(4, 4, "P", "", "deletion invalid mod", vcf, mods));
        }

        [Test]
        public void Constructor_StopGainedWithDownstreamModification_Throws()
        {
            // Variant introduces termination (*) at position 4; any modification at or after 4 is invalid.
            var mod = CreateTestModification('P');
            var vcf = CreateTestVcf();
            var mods = new Dictionary<int, List<Modification>>
            {
                { 5, new List<Modification> { mod } } // downstream of premature stop
            };

            Assert.Throws<ArgumentException>(() =>
                new SequenceVariation(4, 4, "P", "*", "stop gained invalid mod", vcf, mods));
        }

        [Test]
        public void Constructor_InsertionWithValidInternalModification_DoesNotThrow()
        {
            // Insertion: original 'P' (len 1) replaced by 'PPP' (len 3) at position 4; new span 4..6
            // Modification at 5 is valid (inside new inserted span)
            var mod = CreateTestModification('P');
            var vcf = CreateTestVcf();
            var mods = new Dictionary<int, List<Modification>>
            {
                { 5, new List<Modification> { mod } } // valid within expanded span
            };

            Assert.DoesNotThrow(() =>
                new SequenceVariation(4, 4, "P", "PPP", "insertion with valid mod", vcf, mods));
        }
    }
}