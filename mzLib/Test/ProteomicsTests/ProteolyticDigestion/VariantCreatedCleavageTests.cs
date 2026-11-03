using System;
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
    /// <summary>
    /// A variant beside a peptide identifies it when the variant creates the cleavage that makes the peptide.
    /// Whether a bond is a cleavage is the protease's whole rule -- a preventing residue, a recognition sequence
    /// longer than one residue, a cut inside the motif -- not a comparison of the single residue next to the bond.
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class VariantCreatedCleavageTests
    {
        private const string Vcf = @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30";

        private static (bool intersects, bool identifies) Identify(string sequence, int variantBegin, int variantEnd,
            string original, string variant, string protease, int peptideStart, int peptideEnd)
        {
            var protein = new Protein(sequence, "P", sequenceVariations: new List<SequenceVariation>
                { new SequenceVariation(variantBegin, variantEnd, original, variant, Vcf, null) });
            var variantProtein = protein.GetVariantBioPolymers().Single(p => p.AppliedSequenceVariations.Count == 1);
            var peptide = new PeptideWithSetModifications(variantProtein, new DigestionParams(protease: protease),
                peptideStart, peptideEnd, CleavageSpecificity.Full, "full", 0, new Dictionary<int, Modification>(), 0);
            return peptide.IntersectsAndIdentifiesVariation(variantProtein.AppliedSequenceVariations.Single());
        }

        [Test]
        public static void AVariantThatMakesAKBeforeAProlineDoesNotIdentifyThePeptideAfterIt()
        {
            // MAEA[A->K]PTIDER: trypsin|P does not cut K|P, so the variant does not make PTIDER
            Assert.AreEqual((false, false), Identify("MAEAAPTIDER", 5, 5, "A", "K", "trypsin|P", 6, 11));
        }

        [Test]
        public static void AVariantThatRemovesAProlineAfterKIdentifiesThePeptideBeforeIt()
        {
            // MAEK[P->A]TIDER: trypsin|P now cuts K|A, so the variant makes MAEK
            Assert.AreEqual((false, true), Identify("MAEKPTIDER", 5, 5, "P", "A", "trypsin|P", 1, 4));
        }

        [Test]
        public static void AVariantThatCompletesACollagenaseMotifBeforeThePeptideIdentifiesIt()
        {
            // MS[APA->GPA]GPALE: GPX|GPX now cuts after residue 5, making GPALE
            Assert.AreEqual((false, true), Identify("MSAPAGPALE", 3, 5, "APA", "GPA", "collagenase", 6, 10));
        }

        [Test]
        public static void AVariantThatCompletesACollagenaseMotifAfterThePeptideIdentifiesIt()
        {
            // MSGPA[APA->GPA]LE: GPX|GPX now cuts after residue 5, making MSGPA
            Assert.AreEqual((false, true), Identify("MSGPAAPALE", 6, 8, "APA", "GPA", "collagenase", 1, 5));
        }

        [Test]
        public static void AVariantBesideAnUnchangedCollagenaseCutDoesNotIdentifyThePeptide()
        {
            // MSGPA[A->V]PALE... the cut after residue 5 needs GPX on both sides; residue 6 becomes V and kills it,
            // so the variant cannot have made MSGPA. The control: the same bond, a variant that breaks the motif.
            Assert.AreEqual((false, false), Identify("MSGPAGPALE", 6, 6, "G", "V", "collagenase", 1, 5));
        }

        /// <summary>
        /// <see cref="DigestionAgent.IsDigestionSite"/> is <see cref="DigestionAgent.GetDigestionSiteIndices"/>
        /// asked at one bond, for every shipped protease, at every bond including both termini.
        /// </summary>
        [Test]
        public static void IsDigestionSiteAgreesWithGetDigestionSiteIndicesForEveryProtease()
        {
            string[] sequences =
            {
                "PEPTIDEK", "KPEPTIDER", "RPAAAAK", "MSGPAGPALEGPXGPAK", "TATSASTKPRPDEWYFLMNQ", "K", "P",
                "ACDEFGHIKLMNPQRSTVWYKRPKPRDNGSTX", "GPAGPAGPAGPA",
            };
            foreach (var protease in ProteaseDictionary.Dictionary.Values)
            {
                foreach (var sequence in sequences)
                {
                    var sites = new HashSet<int>(protease.GetDigestionSiteIndices(sequence));
                    for (int i = 0; i <= sequence.Length; i++)
                    {
                        Assert.AreEqual(sites.Contains(i), protease.IsDigestionSite(sequence, i),
                            $"{protease.Name} on {sequence} after residue {i}");
                    }
                }
            }
        }
    }
}
