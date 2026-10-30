using NUnit.Framework;
using Omics.Modifications;
using Assert = NUnit.Framework.Legacy.ClassicAssert;

namespace Test.Omics.Modifications
{
    /// <summary>
    /// The glycoproteases invert trypsin's relationship with modifications: their cleavage does not
    /// exist unless a glycan is present at the position the rule governs, so an unglycosylated Ser/Thr
    /// is not a site at all. Acting on that needs the modification's glycosylation CLASS, because the
    /// requirement is chemical rather than nominal -- no single modification id can express "O-GalNAc
    /// with core 1 or core 2, but not O-GlcNAc".
    ///
    /// These tests pin the classification itself: which fields it reads, which it refuses to guess
    /// from, and that an N-linked glycan never satisfies an O-linked requirement. They are the mirror
    /// of CleavageBlockingModifications_ClassifyOnlyChargeNeutralizingAcylations, and like it they test
    /// the chemistry half alone -- no digestion is involved yet.
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class CleavagePromotingModificationTests
    {
        private static Modification Mod(string originalId, string motif, string modificationType = "Test",
            string featureType = null, string locationRestriction = "Anywhere.")
        {
            Assert.IsTrue(ModificationMotif.TryGetMotif(motif, out ModificationMotif parsedMotif),
                $"test setup: '{motif}' is not a legal modification motif");
            return new Modification(_originalId: originalId, _modificationType: modificationType,
                _featureType: featureType, _target: parsedMotif, _locationRestriction: locationRestriction,
                _monoisotopicMass: 203.079373);
        }

        // --- Source 1: ModificationType, which is what MetaMorpheus's Glycan type sets ---

        [Test]
        [TestCase("O-linked glycosylation", GlycosylationClass.OLinked)]
        [TestCase("N-linked glycosylation", GlycosylationClass.NLinked)]
        [TestCase("Other glycosylation", GlycosylationClass.Other)]
        public static void ClassifyGlycosylation_ReadsTheModificationTypeMetaMorpheusSets(
            string modificationType, GlycosylationClass expected)
        {
            // MetaMorpheus's EngineLayer.Glycan (a Modification subclass) sets exactly these strings,
            // and they are also Unimod's own classification vocabulary. A glyco search's glycans are
            // the ones that have to classify correctly here.
            Modification glycan = Mod("H2N2A2F1", "T", modificationType);
            Assert.AreEqual(expected, CleavagePromotingModifications.ClassifyGlycosylation(glycan));
        }

        [Test]
        public static void ClassifyGlycosylation_ModificationTypeWins_WhenItDisagreesWithTheTargetResidue()
        {
            // An O-linked glycan declared on Asn is chemically odd, but if the modification type says
            // O-linked we believe it rather than re-deriving the linkage from the residue. The type is
            // the authoritative statement; the residue is only consulted when it is absent.
            Modification oddlyTargeted = Mod("H1N1", "N", "O-linked glycosylation");
            Assert.AreEqual(GlycosylationClass.OLinked,
                CleavagePromotingModifications.ClassifyGlycosylation(oddlyTargeted));
        }

        // --- Source 2: UniProt's CARBOHYD feature key, disambiguated by the target residue ---

        [Test]
        [TestCase("N", GlycosylationClass.NLinked)]
        [TestCase("S", GlycosylationClass.OLinked)]
        [TestCase("T", GlycosylationClass.OLinked)]
        public static void ClassifyGlycosylation_CarbohydFeature_TakesLinkageFromTheTargetResidue(
            string motif, GlycosylationClass expected)
        {
            // UniProt's ptmlist says "this is a glycosylation site" without saying which linkage, so
            // the residue decides: Asn carries N-linked glycans, Ser/Thr O-linked ones.
            Modification uniProtGlycan = Mod("Glycosylation site", motif, "UniProt", featureType: "CARBOHYD");
            Assert.AreEqual(expected, CleavagePromotingModifications.ClassifyGlycosylation(uniProtGlycan));
        }

        [Test]
        public static void ClassifyGlycosylation_CarbohydOnAnUnexpectedResidue_IsOtherNotGuessed()
        {
            // C-mannosylated tryptophan is a real UniProt CARBOHYD entry. It is a glycan, but it is not
            // a linkage any protease rule keys on, so it must classify as Other rather than being
            // forced into N or O.
            Modification cMannosylTryptophan = Mod("C-linked (Man) tryptophan", "W", "UniProt", featureType: "CARBOHYD");
            Assert.AreEqual(GlycosylationClass.Other,
                CleavagePromotingModifications.ClassifyGlycosylation(cMannosylTryptophan));
        }

        // --- What must NOT classify as a glycan ---

        [Test]
        public static void ClassifyGlycosylation_PhosphoSerine_IsNotAGlycan()
        {
            // The guard that matters most. A rule keyed on "sits on Ser/Thr" would call phospho an
            // O-glycan and let it satisfy an O-glycoprotease, inventing cleavage sites wholesale.
            // Ser/Thr is necessary but nowhere near sufficient.
            Modification phospho = Mod("Phosphoserine", "S", "UniProt");
            Assert.AreEqual(GlycosylationClass.None,
                CleavagePromotingModifications.ClassifyGlycosylation(phospho));
        }

        [Test]
        [TestCase("N6-acetyllysine", "K")]
        [TestCase("Oxidation", "M")]
        [TestCase("Carbamidomethyl", "C")]
        public static void ClassifyGlycosylation_OrdinaryModifications_AreNone(string id, string motif)
        {
            Assert.AreEqual(GlycosylationClass.None,
                CleavagePromotingModifications.ClassifyGlycosylation(Mod(id, motif, "UniProt")));
        }

        [Test]
        public static void ClassifyGlycosylation_CompositionalNameAlone_IsNotEnough()
        {
            // Deliberately NOT classified: a glycan-shaped id with no modification type and no feature
            // type. Glycan names are compositional rather than descriptive, so substring matching on
            // them would be guesswork -- unlike the blocking set, where "acetyl"/"succinyl" are
            // genuinely descriptive stems. This test pins the absence of an id fallback, so that adding
            // one later is a deliberate decision rather than an accident.
            Modification bareComposition = Mod("Hex(1)HexNAc(1)", "T", "Unimod");
            Assert.AreEqual(GlycosylationClass.None,
                CleavagePromotingModifications.ClassifyGlycosylation(bareComposition));
        }

        [Test]
        public static void ClassifyGlycosylation_GuardsAgainstMalformedModifications()
        {
            Assert.AreEqual(GlycosylationClass.None, CleavagePromotingModifications.ClassifyGlycosylation(null));

            // A CARBOHYD entry with no target motif has no residue to read the linkage from.
            Modification noTarget = new Modification(_originalId: "Glycosylation site",
                _modificationType: "UniProt", _featureType: "CARBOHYD", _target: null,
                _locationRestriction: "Anywhere.", _monoisotopicMass: 203.079373);
            Assert.AreEqual(GlycosylationClass.Other,
                CleavagePromotingModifications.ClassifyGlycosylation(noTarget));
        }

        // --- PromotedResidue ---

        [Test]
        public static void PromotedResidue_ContextBearingMotif_ReturnsTheModifiedResidue()
        {
            // The N-glycosylation sequon "Nxs" carries lower-case context around the modified residue,
            // so the upper-case letter is the one the glycan sits on.
            Modification sequon = Mod("Glycosylation site", "Nxs", "UniProt", featureType: "CARBOHYD");
            Assert.AreEqual('N', CleavagePromotingModifications.PromotedResidue(sequon));
            Assert.AreEqual(GlycosylationClass.NLinked,
                CleavagePromotingModifications.ClassifyGlycosylation(sequon));
        }

        [Test]
        public static void PromotedResidue_NoTarget_ReturnsTheNullCharacter()
        {
            Assert.AreEqual('\0', CleavagePromotingModifications.PromotedResidue(null));
        }

        // --- Satisfies: the inversion guard ---

        [Test]
        public static void Satisfies_AnNLinkedGlycan_DoesNotSatisfyAnOLinkedRequirement()
        {
            // The whole reason the class has two members rather than one "is a glycan" flag. An
            // O-glycoprotease must not be satisfied by an N-glycan, and flavastacin must not be
            // satisfied by an O-glycan.
            Modification nLinked = Mod("H5N2", "N", "N-linked glycosylation");
            Modification oLinked = Mod("H1N1", "T", "O-linked glycosylation");

            Assert.IsTrue(CleavagePromotingModifications.Satisfies(nLinked, GlycosylationClass.NLinked));
            Assert.IsFalse(CleavagePromotingModifications.Satisfies(nLinked, GlycosylationClass.OLinked));
            Assert.IsTrue(CleavagePromotingModifications.Satisfies(oLinked, GlycosylationClass.OLinked));
            Assert.IsFalse(CleavagePromotingModifications.Satisfies(oLinked, GlycosylationClass.NLinked));
        }

        [Test]
        public static void Satisfies_RequiringNone_IsNeverSatisfied()
        {
            // A rule that requires nothing must not consult Satisfies at all. Making None unsatisfiable
            // means a miswired rule fails closed (no cleavage) rather than admitting every modification.
            Modification oLinked = Mod("H1N1", "T", "O-linked glycosylation");
            Modification phospho = Mod("Phosphoserine", "S", "UniProt");

            Assert.IsFalse(CleavagePromotingModifications.Satisfies(oLinked, GlycosylationClass.None));
            Assert.IsFalse(CleavagePromotingModifications.Satisfies(phospho, GlycosylationClass.None));
            Assert.IsFalse(CleavagePromotingModifications.Satisfies(null, GlycosylationClass.None));
        }

        [Test]
        public static void Satisfies_NullModification_IsNeverSatisfied()
        {
            Assert.IsFalse(CleavagePromotingModifications.Satisfies(null, GlycosylationClass.OLinked));
            Assert.IsFalse(CleavagePromotingModifications.Satisfies(null, GlycosylationClass.NLinked));
        }
    }
}
