using System.Collections.Generic;
using System.IO;
using System.Linq;
using MzLibUtil;
using NUnit.Framework;
using Omics.BioPolymer;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class CleavagePromotingDigestionTests
    {
        // StcE cuts RPPITQ|SSL only when Thr5, at P2, carries an O-glycan.
        private const string StcESubstrate = "RPPITQSSL";

        private static Modification Mod(string residue, string modificationType = "O-linked glycosylation", string featureType = null)
        {
            ModificationMotif.TryGetMotif(residue ?? "X", out ModificationMotif target);
            return new Modification(_originalId: "HexNAc", _modificationType: modificationType, _featureType: featureType,
                _target: residue is null ? null : target, _locationRestriction: "Anywhere.", _monoisotopicMass: 203.079373);
        }

        private static Protein WithAnnotation(string sequence, int residue, Modification modification) =>
            new(sequence, "P", oneBasedModifications: new Dictionary<int, List<Modification>> { { residue, new List<Modification> { modification } } });

        private static List<PeptideWithSetModifications> Digest(Protein protein, string protease, bool promoting = true,
            List<Modification> variableMods = null, int maxMissedCleavages = 0, bool blocking = false, int maxLength = int.MaxValue,
            InitiatorMethionineBehavior initiatorMethionine = InitiatorMethionineBehavior.Retain, int maxMods = 2,
            List<Modification> fixedMods = null) =>
            protein.Digest(new DigestionParams(protease, maxMissedCleavages, minPeptideLength: 1, maxPeptideLength: maxLength,
                        initiatorMethionineBehavior: initiatorMethionine, maxModsForPeptides: maxMods, respectCleavageBlockingModifications: blocking,
                        respectCleavagePromotingModifications: promoting),
                    fixedMods ?? new List<Modification>(), variableMods ?? new List<Modification>())
                .Cast<PeptideWithSetModifications>().ToList();

        private static List<string> BaseSequences(IEnumerable<PeptideWithSetModifications> peptides) =>
            peptides.Select(p => p.BaseSequence).Distinct().OrderBy(s => s, System.StringComparer.Ordinal).ToList();

        [Test]
        public static void ParseAll_ReadsSubsitesAsTheLiteratureWritesThem()
        {
            List<CleavageRequirement> requirements = CleavageRequirement.ParseAll("P2:O-glycan; p1' : o-glycan");
            Assert.AreEqual((2, false), (requirements[0].Subsite, requirements[0].IsPrimeSide));
            Assert.AreEqual((1, true), (requirements[1].Subsite, requirements[1].IsPrimeSide));
            Assert.IsFalse(requirements[0].IsForbidden);
            Assert.IsTrue(CleavageRequirement.ParseAll("P1:!O-glycan")[0].IsForbidden);
            CollectionAssert.IsEmpty(CleavageRequirement.ParseAll(""));
        }

        [Test]
        [TestCase("P1:!N-glycan")]
        [TestCase("P2")]
        [TestCase("X2:O-glycan")]
        [TestCase("Pa:O-glycan")]
        [TestCase("P0:O-glycan")]
        [TestCase("P2:N-glycan")]
        public static void ParseAll_RefusesWhatItCannotRead(string text)
        {
            Assert.Throws<MzLibException>(() => CleavageRequirement.ParseAll(text));
        }

        [Test]
        [TestCase("T", "O-linked glycosylation", null, true)]
        [TestCase("S", "UniProt", "CARBOHYD", true)]
        [TestCase("T", "UniProt", "CARBOHYD", true)]
        [TestCase("N", "UniProt", "CARBOHYD", false)]
        [TestCase(null, "UniProt", "CARBOHYD", false)]
        [TestCase("T", "Common Biological", null, false)]
        public static void OnlyAnOGlycanMeetsTheRequirement(string residue, string modificationType, string featureType, bool satisfies)
        {
            Assert.AreEqual(satisfies, new CleavageRequirement(2, false).IsSatisfiedBy(Mod(residue, modificationType, featureType)));
        }

        [Test]
        public static void InACoDigestTheRequirementBelongsOnlyToTheGlycoproteaseMotifs()
        {
            List<DigestionMotif> motifs = ProteaseDictionary.Dictionary["StcE-trypsin"].DigestionMotifs;
            CollectionAssert.AreEqual(new[] { 1, 1, 1, 1, 0, 0 }, motifs.Select(m => m.CleavageRequirements.Count));
            // IMPa's P1:!O-glycan names trypsin's own P1, and still must not give K| and R| the glycan rule.
            CollectionAssert.AreEqual(new[] { 2, 2, 0, 0 },
                ProteaseDictionary.Dictionary["IMPa-trypsin|P"].DigestionMotifs.Select(m => m.CleavageRequirements.Count));
            Assert.AreEqual(2, motifs[0].CleavageRequirements[0].Subsite);
        }

        [Test]
        public static void ARequirementNoMotifReachesBelongsToEveryMotif()
        {
            string path = Path.Combine(TestContext.CurrentContext.WorkDirectory, "requirement-on-every-motif.tsv");
            File.WriteAllText(path, "Name\tMotif\tSpecificity\tCleavage Requirement\nK-needs-glycan\tK|,R|\tfull\tP3:O-glycan\n");
            List<DigestionMotif> motifs = ProteaseDictionary.LoadProteaseDictionary(path)["K-needs-glycan"].DigestionMotifs;
            File.Delete(path);
            Assert.IsTrue(motifs.All(m => m.CleavageRequirements.Count == 1));
        }

        [Test]
        public static void TheFlagDefaultsOffAndSurvivesCloneEqualityAndToString()
        {
            Assert.IsFalse(new DigestionParams().RespectCleavagePromotingModifications);
            var on = new DigestionParams("StcE", respectCleavagePromotingModifications: true);
            var nonSpecific = new DigestionParams("StcE", searchModeType: CleavageSpecificity.None, respectCleavagePromotingModifications: true);
            Assert.IsTrue(((DigestionParams)on.Clone()).RespectCleavagePromotingModifications);
            Assert.IsTrue(((DigestionParams)nonSpecific.Clone()).RespectCleavagePromotingModifications);
            Assert.AreNotEqual(new DigestionParams("StcE"), on);
            Assert.AreNotEqual(new DigestionParams("StcE").GetHashCode(), on.GetHashCode());
            Assert.IsTrue(on.ToString().EndsWith(",True"));
        }

        [Test]
        public static void ThePolicyRemovesTheSitesNoGlycanCanReach()
        {
            var protein = new Protein(StcESubstrate, "P");
            var parameters = new DigestionParams("StcE", respectCleavagePromotingModifications: true);
            CollectionAssert.AreEqual(new[] { 0, 9 }, CleavageBlockingPolicy.For(parameters, null, protein).CleavageSites);
            CollectionAssert.AreEqual(new[] { 0, 6, 9 }, CleavageBlockingPolicy.For(parameters, new List<Modification> { Mod("T") }, protein).CleavageSites);
            Assert.IsNull(CleavageBlockingPolicy.For(new DigestionParams("StcE"), null, protein).CleavageSites);
        }

        [Test]
        public static void WithNothingToCarryTheGlycan_StcECutsOnlyWithTheFlagOff()
        {
            var protein = new Protein(StcESubstrate, "P");
            CollectionAssert.AreEqual(new[] { "RPPITQ", "SSL" }, BaseSequences(Digest(protein, "StcE", promoting: false)));
            CollectionAssert.AreEqual(new[] { StcESubstrate }, BaseSequences(Digest(protein, "StcE")));
        }

        [Test]
        public static void AGlycanAtP2JustifiesTheCut_AndTheUnglycosylatedFormReadsThrough()
        {
            List<PeptideWithSetModifications> peptides = Digest(WithAnnotation(StcESubstrate, 5, Mod("T")), "StcE");
            CollectionAssert.AreEqual(new[] { "RPPITQ", "RPPITQSSL", "SSL" }, BaseSequences(peptides));
            Assert.IsTrue(peptides.Where(p => p.BaseSequence == "RPPITQ").All(p => p.AllModsOneIsNterminus.ContainsKey(6)));
            PeptideWithSetModifications readThrough = peptides.Single(p => p.BaseSequence == StcESubstrate);
            Assert.IsFalse(readThrough.AllModsOneIsNterminus.ContainsKey(6));
            Assert.AreEqual(0, readThrough.MissedCleavages);
        }

        [Test]
        public static void AVariableGlycanMakesTheSiteFeasible()
        {
            CollectionAssert.AreEqual(new[] { "RPPITQ", "RPPITQSSL", "SSL" },
                BaseSequences(Digest(new Protein(StcESubstrate, "P"), "StcE", variableMods: new List<Modification> { Mod("T") })));
        }

        [Test]
        public static void ASiteNoGlycanCanReachIsNotDiscountedAgain()
        {
            // The StcE site after Gln6 has left the site list, so the read-through's one missed cleavage is Lys10.
            PeptideWithSetModifications readThrough = Digest(new Protein("GPPITQSSLKAAR", "P"), "StcE-trypsin", maxMissedCleavages: 1)
                .Single(p => p.BaseSequence == "GPPITQSSLKAAR");
            Assert.AreEqual(1, readThrough.MissedCleavages);
        }

        [Test]
        public static void APhosphateAtP2DoesNotJustifyTheCut()
        {
            CollectionAssert.AreEqual(new[] { StcESubstrate },
                BaseSequences(Digest(WithAnnotation(StcESubstrate, 5, Mod("T", "Common Biological")), "StcE")));
        }

        [Test]
        public static void APrimeSideRequirementIsMetByAPermanentOrAnAnnotatedGlycan()
        {
            var requirement = new List<CleavageRequirement> { new(1, isPrimeSide: true) };
            ProteaseDictionary.Dictionary["OgpA-like"] = new Protease("OgpA-like", CleavageSpecificity.Full, null, null,
                new List<DigestionMotif> { new("T", null, 0, null, requirement), new("S", null, 0, null, requirement) });

            var permanent = new Protein("ASAAA", "P", oneBasedFixedModifications: new Dictionary<int, Modification> { { 2, Mod("S") } });
            CollectionAssert.AreEqual(new[] { "A", "SAAA" }, BaseSequences(Digest(permanent, "OgpA-like")));

            List<PeptideWithSetModifications> annotated = Digest(WithAnnotation("ASAAA", 2, Mod("S")), "OgpA-like");
            CollectionAssert.AreEqual(new[] { "A", "ASAAA", "SAAA" }, BaseSequences(annotated));
            Assert.IsTrue(annotated.Where(p => p.BaseSequence == "SAAA").All(p => p.AllModsOneIsNterminus.ContainsKey(2)));
        }

        [Test]
        [TestCase(true)]
        [TestCase(false)]
        public static void IMPaDoesNotCutBetweenTwoAdjacentGlycosites_SmEDoes(bool annotatedOnSer3)
        {
            // GA|S|TAAA with glycans possible on Ser3 and Thr4. Cutting S|T needs Thr4 glycosylated and, for IMPa, Ser3 not.
            var annotations = new Dictionary<int, List<Modification>> { { 4, new List<Modification> { Mod("T") } } };
            if (annotatedOnSer3)
            {
                annotations[3] = new List<Modification> { Mod("S") };
            }
            var protein = new Protein("GASTAAA", "P", oneBasedModifications: annotations);
            List<Modification> variableMods = annotatedOnSer3 ? null : new List<Modification> { Mod("S") };
            CollectionAssert.AreEqual(new[] { "GA", "S", "TAAA" }, BaseSequences(Digest(protein, "IMPa", promoting: false)));
            CollectionAssert.AreEqual(new[] { "GA", "GAS", "GASTAAA", "S", "STAAA", "TAAA" },
                BaseSequences(Digest(protein, "SmE", variableMods: variableMods)));

            List<PeptideWithSetModifications> impa = Digest(protein, "IMPa", variableMods: variableMods);
            CollectionAssert.AreEqual(new[] { "GA", "GAS", "GASTAAA", "STAAA", "TAAA" }, BaseSequences(impa));
            Assert.IsFalse(impa.Single(p => p.BaseSequence == "GAS").AllModsOneIsNterminus.ContainsKey(4));
        }

        [Test]
        [TestCase(true)]
        [TestCase(false)]
        public static void AGlycanThatMustSitAtP1_StopsIMPaOnBothSidesOfTheBond(bool permanent)
        {
            // The product starting at the bond cannot see P1, so only a glycan certain to be there decides its side.
            var protein = new Protein("GASTAAA", "P",
                oneBasedModifications: new Dictionary<int, List<Modification>> { { 4, new List<Modification> { Mod("T") } } },
                oneBasedFixedModifications: permanent ? new Dictionary<int, Modification> { { 3, Mod("S") } } : null);
            List<Modification> fixedMods = permanent ? null : new List<Modification> { Mod("S") };
            CollectionAssert.AreEqual(new[] { "GA", "STAAA" }, BaseSequences(Digest(protein, "IMPa", fixedMods: fixedMods)));
        }

        [Test]
        public static void RemovingTheInitiatorMethionineIsNotACut()
        {
            List<string> products = BaseSequences(Digest(WithAnnotation("M" + StcESubstrate, 6, Mod("T")), "StcE",
                initiatorMethionine: InitiatorMethionineBehavior.Variable));
            CollectionAssert.IsSubsetOf(new[] { "MRPPITQ", "RPPITQ" }, products);
        }

        [Test]
        public static void AChainIsNotACut_AndReadsThroughTheSiteItsGlycanIsMissingFrom()
        {
            var protein = new Protein("AAAAK" + StcESubstrate + "GGG", "P", proteolysisProducts: new List<TruncationProduct> { new(6, 14, "chain") });
            PeptideWithSetModifications chain = Digest(protein, "StcE", variableMods: new List<Modification> { Mod("T") }, maxMissedCleavages: 1)
                .Single(p => p.OneBasedStartResidue == 6 && p.OneBasedEndResidue == 14);
            Assert.AreEqual(0, chain.MissedCleavages);
        }

        [Test]
        public static void WithBothCorrections_BlockedAndUnoccupiedSitesAreDiscountedTogether()
        {
            // TA|S A K6 TA|S: two StcE sites either side of an acetylated, so blocked, K6.
            ModificationMotif.TryGetMotif("K", out ModificationMotif lysine);
            var acetylK = new Modification(_originalId: "N6-acetyllysine", _modificationType: "Test", _target: lysine,
                _locationRestriction: "Anywhere.", _monoisotopicMass: 42.010565);
            PeptideWithSetModifications readThrough = Digest(new Protein("GTASAKGTASGR", "P"), "StcE-trypsin", blocking: true,
                    variableMods: new List<Modification> { Mod("T"), acetylK }, maxLength: 50, maxMods: 1)
                .Single(p => p.BaseSequence == "GTASAKGTASGR" && p.AllModsOneIsNterminus.Keys.SequenceEqual(new[] { 7 }));
            Assert.AreEqual(0, readThrough.MissedCleavages);
        }

        [Test]
        public static void WithoutALengthLimit_ReadThroughsAreCappedSoTheyGrowLinearlyWithTheMucin()
        {
            int Spans(int repeats) => Digest(new Protein(string.Concat(Enumerable.Repeat("PTTTPITTTKVTPTPTPTGTQT", repeats)), "P"), "StcE",
                    variableMods: new List<Modification> { Mod("T"), Mod("S") }, maxMissedCleavages: 2)
                .Select(p => (p.OneBasedStartResidue, p.OneBasedEndResidue)).Distinct().Count();
            Assert.Less(Spans(8), 3 * Spans(4));
        }
    }
}
