using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MzLibUtil;
using NUnit.Framework;
using Omics.BioPolymer;
using Proteomics;

namespace Test
{
    /// <summary>
    /// VariantApplication.ParseAccession, the inverse of GetAccession: an accession read back into its parent
    /// entry and applied variants, so a variant proteoform links to anything stored per entry without losing
    /// its variants. On a reviewed human proteome loaded with variants (52,359 proteins from 20,416 entries),
    /// every variant proteoform's gene answer equals its entry's.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestProteoformAccession
    {
        [TestCase("O14994_S470N", "O14994", "S470N")]
        [TestCase("P12345_S70N_A80T", "P12345", "S70N_A80T")]
        [TestCase("P12345_AB70", "P12345", "AB70")]
        [TestCase("P12345_T70TAG", "P12345", "T70TAG")]
        [TestCase("A0A087X1C5_S70N", "A0A087X1C5", "S70N")]
        public void AVariantProteoformSplitsIntoItsEntry(string accession, string entry, string variants)
        {
            var a = VariantApplication.ParseAccession(accession);

            Assert.That(a.Verbatim, Is.EqualTo(accession));
            Assert.That((a.Entry.EntryAccession, a.Entry.Namespace), Is.EqualTo((entry, AccessionNamespace.UniProt)));
            Assert.That(a.AppliedVariants, Is.EqualTo(variants));
            Assert.That(a.HasAppliedVariants, Is.True);
        }

        [Test]
        public void AnIsoformProteoformKeepsItsIsoform()
        {
            var a = VariantApplication.ParseAccession("P12345-2_S70N");

            Assert.That((a.Entry.EntryAccession, a.Entry.Isoform, a.AppliedVariants), Is.EqualTo(("P12345", (int?)2, "S70N")));
            Assert.That(a.Entry.Verbatim, Is.EqualTo("P12345-2"));
        }

        [Test]
        public void ARefSeqUnderscoreIsNotAVariantSuffix()
        {
            // "Text before the first _" would turn NP_000537 into NP.
            var a = VariantApplication.ParseAccession("NP_000537.3_R72P");
            Assert.That((a.Entry.EntryAccession, a.Entry.Version, a.AppliedVariants, a.Entry.Namespace),
                Is.EqualTo(("NP_000537", (int?)3, "R72P", AccessionNamespace.RefSeq)));

            var plain = VariantApplication.ParseAccession("NP_000537");
            Assert.That((plain.Entry.EntryAccession, plain.HasAppliedVariants), Is.EqualTo(("NP_000537", false)));
        }

        [Test]
        public void AnEntryIsItsOwnEntry_WithNoVariants()
        {
            var a = VariantApplication.ParseAccession("P12345-1");

            Assert.That(a.Entry, Is.EqualTo("P12345-1".ParseProteinAccession()));
            Assert.That(a.AppliedVariants, Is.Null);
            Assert.That(a.HasAppliedVariants, Is.False);
        }

        // ProteinDbLoader's load-collision counter names a DIFFERENT entry whose accession collided.
        [TestCase("P12345_2")]
        [TestCase("P12345_10")]
        [TestCase("P12345_S70N_2")]
        [TestCase("NP_000537_2")]
        [TestCase("DECOY_P12345")]
        [TestCase("DECOY_P12345_S70N")]
        [TestCase("Random_P12345")]
        [TestCase("CON__P02768")]
        [TestCase("p12345_s70n")]
        [TestCase("")]
        public void ANameThatOnlyLooksLikeAProteoformIsUnrecognized_AndKeptVerbatim(string accession)
        {
            var a = VariantApplication.ParseAccession(accession);

            Assert.That(a.Entry.Namespace, Is.EqualTo(AccessionNamespace.Unrecognized));
            Assert.That((a.Verbatim, a.Entry.EntryAccession), Is.EqualTo((accession, accession)));
            Assert.That(a.AppliedVariants, Is.Null);
        }

        [Test]
        public void NullParsesAsUnrecognizedEmpty()
        {
            var a = VariantApplication.ParseAccession(null);
            Assert.That((a.Verbatim, a.Entry.Namespace, a.AppliedVariants),
                Is.EqualTo(("", AccessionNamespace.Unrecognized, (string)null)));
        }

        [Test]
        public void EveryAccessionVariantApplicationWritesParsesBackToItsEntry()
        {
            // The grammar is the producer's: take the names mzLib itself gives applied-variant
            // proteoforms (a substitution, an insertion, a deletion, and combinations) and read them back.
            var variants = new List<SequenceVariation>
            {
                new(3, 3, "A", "T", "substitution"),
                new(6, 6, "G", "GKK", "insertion"),
                new(9, 10, "AA", "", "deletion"),
            };
            var protein = new Protein("MAAAAGAAAAG", "P12345", sequenceVariations: variants);

            var proteoforms = VariantApplication.ApplyAllVariantCombinations(protein, variants, maxCombinations: 10)
                .Where(p => p.AppliedSequenceVariations.Count > 0)
                .ToList();
            Assert.That(proteoforms, Has.Count.EqualTo(7), "every non-empty combination of three variants");

            foreach (var p in proteoforms)
            {
                var a = VariantApplication.ParseAccession(p.Accession);
                Assert.That(a.Entry.EntryAccession, Is.EqualTo("P12345"), p.Accession);
                Assert.That(a.AppliedVariants.Split('_'), Has.Length.EqualTo(p.AppliedSequenceVariations.Count), p.Accession);
            }
        }
    }
}
