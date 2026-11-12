using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Omics.BioPolymer;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Assert = NUnit.Framework.Legacy.ClassicAssert;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class VariantCreatedCleavageTests
    {
        private const string Vcf = @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30";

        [Test]
        [TestCase("MAEAAPTIDER", 5, 5, "A", "K", "trypsin|P", 6, 11, false)] // K|P is not cut
        [TestCase("MAEKPTIDER", 5, 5, "P", "A", "trypsin|P", 1, 4, true)] // P->A unblocks K|A
        [TestCase("MSAPAGPALE", 3, 5, "APA", "GPA", "collagenase", 6, 10, true)] // GPX|GPX completed before the peptide
        [TestCase("MSGPAAPALE", 6, 8, "APA", "GPA", "collagenase", 1, 5, true)] // GPX|GPX completed after the peptide
        [TestCase("MAEKAPTIDER", 5, 5, "A", "G", "trypsin|P", 1, 4, false)] // K|A was already cut before the variant
        public static void AVariantIdentifiesAPeptideWhenItCreatesTheCleavage(string sequence, int variantBegin,
            int variantEnd, string original, string variant, string protease, int peptideStart, int peptideEnd,
            bool identifies)
        {
            var protein = new Protein(sequence, "P", sequenceVariations: new List<SequenceVariation>
                { new SequenceVariation(variantBegin, variantEnd, original, variant, Vcf, null) });
            var variantProtein = protein.GetVariantBioPolymers().Single(p => p.AppliedSequenceVariations.Count == 1);
            var peptide = new PeptideWithSetModifications(variantProtein, new DigestionParams(protease: protease),
                peptideStart, peptideEnd, CleavageSpecificity.Full, "full", 0, new Dictionary<int, Modification>(), 0);

            Assert.AreEqual((false, identifies),
                peptide.IntersectsAndIdentifiesVariation(variantProtein.AppliedSequenceVariations.Single()));
        }
    }
}
