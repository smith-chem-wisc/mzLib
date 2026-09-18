using NUnit.Framework;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Linq;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    /// <summary>
    /// <see cref="DigestionProduct.GetCleavageObligatedSites"/>: where a glycan MUST be, given that the
    /// peptide exists at all.
    /// </summary>
    /// <remarks>
    /// <para>This is the half of the cleavage requirement that a glyco search can use, and the only half.
    /// A glyco search identifies the backbone naked and resolves the glycan afterwards from the precursor
    /// mass, so the requirement can never filter its peptides -- but it can say where the glycan has to
    /// have been. If OpeRATOR made this peptide, residue 1 carries a glycan; that is a consequence of the
    /// peptide existing, not a hypothesis to score.</para>
    ///
    /// <para>The tests below pin the conservative cases as hard as the positive ones, because this
    /// obligation is consumed as a CONSTRAINT: a site wrongly declared obligated forbids the correct
    /// localization, which is far worse than declaring nothing.</para>
    /// </remarks>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class CleavageObligationTests
    {
        private static Protease StcELike(string name)
        {
            var requirement = CleavageRequirement.NonPrime(2, GlycosylationClass.OLinked);
            var motifs = new List<DigestionMotif>
            {
                new("TXT", null, 2, null, requirement),
                new("TXS", null, 2, null, requirement),
                new("SXT", null, 2, null, requirement),
                new("SXS", null, 2, null, requirement),
            };
            var protease = new Protease(name, CleavageSpecificity.Full, null, null, motifs);
            ProteaseDictionary.Dictionary[name] = protease;
            return protease;
        }

        private static Protease OgpALike(string name)
        {
            var requirement = CleavageRequirement.Prime(1, GlycosylationClass.OLinked);
            var motifs = new List<DigestionMotif>
            {
                new("T", null, 0, null, requirement),
                new("S", null, 0, null, requirement),
            };
            var protease = new Protease(name, CleavageSpecificity.Full, null, null, motifs);
            ProteaseDictionary.Dictionary[name] = protease;
            return protease;
        }

        /// <summary>StcE motifs plus the two proline-blind tryptic ones, as the shipped composite has.</summary>
        private static Protease StcETrypsinLike(string name)
        {
            var requirement = CleavageRequirement.NonPrime(2, GlycosylationClass.OLinked);
            var motifs = new List<DigestionMotif>
            {
                new("TXT", null, 2, null, requirement),
                new("TXS", null, 2, null, requirement),
                new("SXT", null, 2, null, requirement),
                new("SXS", null, 2, null, requirement),
                new("K", null, 1, null),
                new("R", null, 1, null),
            };
            var protease = new Protease(name, CleavageSpecificity.Full, null, null, motifs);
            ProteaseDictionary.Dictionary[name] = protease;
            return protease;
        }

        /// <summary>
        /// A product spanning the given residues. Built through the public PeptideWithSetModifications
        /// constructor because ProteolyticPeptide's is internal, and because the obligation is defined on
        /// DigestionProduct, which both share.
        /// </summary>
        private static DigestionProduct Product(Protein parent, int start, int end) =>
            new PeptideWithSetModifications(parent, null, start, end, CleavageSpecificity.Full, "test", 0,
                new Dictionary<int, Modification>(), 0);

        [Test]
        public static void ANonPrimeRequirementObligatesTheResidueBeforeTheCTerminalCut()
        {
            // StcE: T5-Q6-S7 in RPPITQSSL, cut after Q6, glycan required at P2 = Thr5. For the peptide
            // RPPITQ (1..6) the C-terminal cut is the protease's, and Thr5 is its second-to-last residue.
            // Two-based: Thr5 is key 5 - 1 + 2 = 6.
            Protease stcE = StcELike("obligation-StcE");
            var protein = new Protein("RPPITQSSL", "OBL");

            List<int> obligated = Product(protein, 1, 6).GetCleavageObligatedSites(stcE);

            CollectionAssert.AreEqual(new[] { 6 }, obligated,
                "the residue at P2 of the C-terminal cut must be obligated");
        }

        [Test]
        public static void APrimeSideRequirementObligatesTheFirstResidue()
        {
            // OpeRATOR family: cut N-terminal to a glycosylated Ser/Thr, so for the peptide starting at
            // that residue the obligation is its OWN first residue -- two-based key 2.
            Protease ogpA = OgpALike("obligation-OgpA");
            var protein = new Protein("AHGVTSAPDTRK", "OBL");

            // TSAPDTRK starts at residue 5 (Thr), whose N-terminal cut OpeRATOR made.
            List<int> obligated = Product(protein, 5, 12).GetCleavageObligatedSites(ogpA);

            CollectionAssert.AreEqual(new[] { 2 }, obligated,
                "P1' of the N-terminal cut is the peptide's own first residue");
        }

        [Test]
        public static void ASubsiteLyingInTheNeighbouringPeptideIsNotObligated()
        {
            // The mirror of the two cases above, and the reason the API returns positions rather than a
            // verdict. For StcE the N-terminal cut's P2 lies in the PREVIOUS peptide, so this peptide
            // cannot localize it and must not claim to. SSL (7..9) ends at the protein C-terminus, which
            // is not a cut at all, so it obligates nothing whatsoever.
            Protease stcE = StcELike("obligation-StcE-neighbour");
            var protein = new Protein("RPPITQSSL", "OBL");

            List<int> obligated = Product(protein, 7, 9).GetCleavageObligatedSites(stcE);

            CollectionAssert.IsEmpty(obligated,
                "P2 of this peptide's N-terminal cut is Thr5, which is not in this peptide");
        }

        [Test]
        public static void ACutAnOrdinaryMotifExplainsObligatesNothing()
        {
            // StcE-trypsin. A cut after Arg is explained by the tryptic motif, which requires nothing, so
            // no glycan is implied anywhere -- even though the composite also carries glycan-requiring
            // motifs. Getting this wrong would forbid localizations on every tryptic peptide in the digest.
            Protease composite = StcETrypsinLike("obligation-StcE-trypsin");
            var protein = new Protein("AAARTTTAAA", "OBL");

            List<int> obligated = Product(protein, 1, 4).GetCleavageObligatedSites(composite);

            CollectionAssert.IsEmpty(obligated, "a tryptic cut needs no glycan to explain it");
        }

        [Test]
        public static void ProteinTerminiObligateNothing()
        {
            Protease ogpA = OgpALike("obligation-OgpA-termini");
            var protein = new Protein("TSAPDTRK", "OBL");

            // The whole protein: neither terminus was produced by the protease.
            List<int> obligated = Product(protein, 1, 8).GetCleavageObligatedSites(ogpA);

            CollectionAssert.IsEmpty(obligated,
                "a peptide that spans the whole sequence was not cut out of anything");
        }

        [Test]
        public static void AnInitiatorMethionineTerminusObligatesNothing()
        {
            // Residue 2 of a Met-initiated sequence: that N-terminus came from initiator-methionine
            // removal, not from the protease, so nothing has to explain it.
            Protease ogpA = OgpALike("obligation-OgpA-initmet");
            var protein = new Protein("MTSAPDTRK", "OBL");

            List<int> obligated = Product(protein, 2, 9).GetCleavageObligatedSites(ogpA);

            CollectionAssert.IsEmpty(obligated, "initiator methionine removal is not a protease cut");
        }

        [Test]
        public static void AnOrdinaryProteaseObligatesNothingAndPaysNothing()
        {
            Protease trypsin = ProteaseDictionary.Dictionary["trypsin"];
            var protein = new Protein("AAAKBBBKCCC", "OBL");

            Assert.IsFalse(trypsin.HasCleavageRequirement);
            CollectionAssert.IsEmpty(Product(protein, 5, 8).GetCleavageObligatedSites(trypsin));
        }

        [Test]
        public static void AnObligationCarriesTheConditionsThatWouldSatisfyIt()
        {
            // The obligation says WHERE a modification must be; the conditions say WHICH ones would do.
            // A glyco search needs both, because OpeRATOR does not merely require "a glycan" at P1' -- it
            // requires at least core 1, and is blocked again by core 2 -- and the localizer holds real
            // glycans it can put through IsSatisfiedBy directly.
            Protease opeRator = ProteaseDictionary.Dictionary["OpeRATOR"];
            var protein = new Protein("AHGVTSAPDTRK", "COND");

            var obligations = Product(protein, 5, 12).GetCleavageObligations(opeRator);

            Assert.AreEqual(1, obligations.Count, "one obligated site");
            Assert.IsTrue(obligations.ContainsKey(2), "the peptide's first residue");

            var conditions = obligations[2];
            Assert.AreEqual(1, conditions.Count(c => !c.IsForbidden), "one required condition");
            Assert.AreEqual(1, conditions.Count(c => c.IsForbidden),
                "and the forbidden one, which still says what may NOT sit on an obligated site");

            CleavageRequirement floor = conditions.First(c => !c.IsForbidden);
            Assert.IsNotNull(floor.MinimumComposition, "OpeRATOR's required condition carries a floor");

            // The whole point, expressed the way the localizer will use it.
            Assert.IsTrue(floor.IsSatisfiedBy(Glycan("Hex1HexNAc1")), "core 1 satisfies OpeRATOR");
            Assert.IsFalse(floor.IsSatisfiedBy(Glycan("HexNAc1")),
                "the Tn antigen is a lone HexNAc and does not reach core 1");

            CleavageRequirement blocked = conditions.First(c => c.IsForbidden);
            Assert.IsTrue(blocked.IsSatisfiedBy(Glycan("Hex1HexNAc2")),
                "core 2 matches the forbidding condition, so a localizer must refuse it here");
            Assert.IsFalse(blocked.IsSatisfiedBy(Glycan("Hex1HexNAc1")), "core 1 does not");
        }

        private static Modification Glycan(string composition)
        {
            ModificationMotif.TryGetMotif("T", out ModificationMotif motif);
            return new Modification(_originalId: composition, _modificationType: "O-linked glycosylation",
                _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 203.079373,
                _monosaccharideComposition: MonosaccharideComposition.Parse(composition));
        }

        [Test]
        public static void InAPrimeSideCoDigestATrypticCutCancelsTheObligation()
        {
            // The case real data pointed at. Every published protocol for these enzymes is a co-digest
            // with trypsin, so a peptide can begin with Ser/Thr and owe that N-terminus to TRYPSIN rather
            // than to the glycoprotease -- and then nothing says its first residue carries a glycan.
            // Obligating it anyway would force a glycan onto a residue that need not have one, which is
            // the over-constraining direction this API must never take.
            Protease composite = ProteaseDictionary.Dictionary["OpeRATOR-trypsin|P"];

            // AAAK | STGLDAAA -- the cut after Lys4 is explained by trypsin, which requires nothing.
            var trypticOrigin = new Protein("AAAKSTGLDAAA", "CODIGEST");
            CollectionAssert.IsEmpty(Product(trypticOrigin, 5, 12).GetCleavageObligatedSites(composite),
                "a tryptic cut needs no glycan, even when the peptide happens to start with Ser");

            // AAAG | STGLDAAA -- nothing tryptic fits, so only the OgpA motif explains the cut.
            var glycoOrigin = new Protein("AAAGSTGLDAAA", "CODIGEST");
            CollectionAssert.AreEqual(new[] { 2 }, Product(glycoOrigin, 5, 12).GetCleavageObligatedSites(composite),
                "with no tryptic explanation the glycoprotease rule stands and residue 1 is obligated");
        }

        [Test]
        public static void EveryShippedCoDigestKeepsItsTrypticMotifsRuleFree()
        {
            // If a tryptic motif inherited the glycan rule, trypsin would need a glycan and would be
            // switched off inside the composite. Checked for the PRIME-side composites specifically,
            // because their subsite arithmetic differs from StcE's non-prime one.
            foreach (string name in new[] { "OpeRATOR-trypsin|P", "IMPa-trypsin|P", "SmE-trypsin|P" })
            {
                Assert.IsTrue(ProteaseDictionary.Dictionary.ContainsKey(name), name + " is missing");
                Protease composite = ProteaseDictionary.Dictionary[name];

                var ruleFree = composite.DigestionMotifs.Where(m => !m.HasCleavageRequirement).ToList();
                CollectionAssert.AreEquivalent(new[] { "K", "R" },
                    ruleFree.Select(m => m.InducingCleavage).ToList(),
                    name + ": the two tryptic motifs, and only those, must be rule-free");

                Assert.AreEqual(2, composite.DigestionMotifs.Count(m => m.HasCleavageRequirement),
                    name + ": both glycoprotease motifs must carry the rule");
            }
        }

        [Test]
        public static void TheObligationSurvivesARealDigestAndAgreesWithThePeptideItCameFrom()
        {
            // End to end: digest with OpeRATOR, then ask each product where its glycan had to be. Every
            // peptide whose N-terminus the protease made must obligate its own first residue, and the
            // N-terminal peptide of the protein must obligate nothing.
            Protease ogpA = OgpALike("obligation-OgpA-endtoend");
            ModificationMotif.TryGetMotif("T", out ModificationMotif motifT);
            var glycan = new Modification(_originalId: "H1N1", _modificationType: "O-linked glycosylation",
                _target: motifT, _locationRestriction: "Anywhere.", _monoisotopicMass: 203.079373);

            var protein = new Protein("AHGVTSAPDTRK", "OBL",
                oneBasedModifications: new Dictionary<int, List<Modification>> { { 5, new List<Modification> { glycan } } });

            var parameters = new DigestionParams(protease: "obligation-OgpA-endtoend", maxMissedCleavages: 0,
                minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain,
                respectCleavagePromotingModifications: true);

            var products = protein.Digest(parameters, new List<Modification>(), new List<Modification>()).ToList();

            var cutOut = products.FirstOrDefault(x => x.OneBasedStartResidue == 5);
            Assert.IsNotNull(cutOut, "expected the peptide beginning at the glycosylated Thr5; got "
                + string.Join(" | ", products.Select(x => x.BaseSequence)));
            CollectionAssert.AreEqual(new[] { 2 }, cutOut.GetCleavageObligatedSites(ogpA),
                "this peptide exists only because Thr5 carried a glycan, so its first residue is obligated");

            var nTerminal = products.FirstOrDefault(x => x.OneBasedStartResidue == 1);
            Assert.IsNotNull(nTerminal);
            CollectionAssert.IsEmpty(nTerminal.GetCleavageObligatedSites(ogpA),
                "the protein's own N-terminus obligates nothing");
        }
    }
}
