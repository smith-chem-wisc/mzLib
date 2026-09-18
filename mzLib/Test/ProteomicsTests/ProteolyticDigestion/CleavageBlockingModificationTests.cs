using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    /// <summary>
    /// Digestion places modifications AFTER cleavage and historically never checked whether a
    /// modification abolished the site it was cut at, so it reported peptidoforms ending in an
    /// acylated lysine -- a cleavage trypsin cannot perform -- often at zero missed cleavages.
    /// These tests pin the opt-in correction and, importantly, that the REAL peptide (which reads
    /// through the blocked residue) survives even when no missed cleavages are allowed.
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class CleavageBlockingModificationTests
    {
        // Trypsin cuts after K8 and after the C-terminal R16. So at 0 missed cleavages the ordinary
        // products are PEPTIDEK and AAAAAAAR; PEPTIDEKAAAAAAAR is the 1-missed-cleavage read-through.
        private const string Sequence = "PEPTIDEKAAAAAAAR";

        private static Modification MakeKModification(string originalId, double mass = 100.01604)
        {
            ModificationMotif.TryGetMotif("K", out ModificationMotif motif);
            return new Modification(_originalId: originalId, _modificationType: "Test", _target: motif,
                _locationRestriction: "Anywhere.", _monoisotopicMass: mass);
        }

        private static List<PeptideWithSetModifications> Digest(bool respectBlockingMods, int maxMissedCleavages,
            params Modification[] variableMods)
        {
            var protein = new Protein(Sequence, "accession");
            var digestionParams = new DigestionParams(protease: "trypsin", maxMissedCleavages: maxMissedCleavages,
                minPeptideLength: 7, respectCleavageBlockingModifications: respectBlockingMods);

            return protein.Digest(digestionParams, new List<Modification>(), variableMods.ToList())
                .Cast<PeptideWithSetModifications>()
                .ToList();
        }

        /// <summary>True if the peptidoform's C-terminal residue carries a modification.</summary>
        private static bool EndsInModifiedResidue(PeptideWithSetModifications peptide) =>
            peptide.AllModsOneIsNterminus.ContainsKey(peptide.BaseSequence.Length + 1);

        [Test]
        public static void CleavageBlockingModifications_ClassifyOnlyChargeNeutralizingAcylations()
        {
            // Acylations remove the epsilon-amine's positive charge, so trypsin cannot cleave.
            Assert.IsTrue(MakeKModification("N6-succinyllysine").BlocksCleavage);
            Assert.IsTrue(MakeKModification("N6-acetyllysine").BlocksCleavage);
            Assert.IsTrue(MakeKModification("Acetyl").BlocksCleavage, "the Unimod short name must classify too");
            Assert.IsTrue(MakeKModification("GG").BlocksCleavage, "the ubiquitin remnant blocks cleavage");

            // Methylation retains the charge; it impairs rather than abolishes cleavage, so it is
            // deliberately excluded (it shows up as a missed cleavage, not an impossible peptide).
            Assert.IsFalse(MakeKModification("N6-methyllysine").BlocksCleavage);
            Assert.IsFalse(MakeKModification("Trimethyl").BlocksCleavage);

            // The residue gate: a name that DOES match a blocking acyl stem, but on a non-cleavage
            // residue, must still be false -- so this fails if the classifier ever stopped checking the
            // target motif (which the phospho case below cannot catch, since "Phospho" matches no stem).
            ModificationMotif.TryGetMotif("S", out ModificationMotif serine);
            var succinylSerine = new Modification(_originalId: "N-succinylserine", _modificationType: "Test",
                _target: serine, _locationRestriction: "Anywhere.", _monoisotopicMass: 100.01604);
            Assert.IsFalse(succinylSerine.BlocksCleavage, "succinyl on serine is not on a cleavage residue");

            var phosphoSerine = new Modification(_originalId: "Phospho", _modificationType: "Test",
                _target: serine, _locationRestriction: "Anywhere.", _monoisotopicMass: 79.96633);
            Assert.IsFalse(phosphoSerine.BlocksCleavage);

            // The location gate: the same acyl group as an N-terminal (alpha-amine) modification on K
            // leaves the side chain charged, so it must NOT classify as blocking.
            ModificationMotif.TryGetMotif("K", out ModificationMotif lysine);
            var nTermAcetylK = new Modification(_originalId: "N-acetyllysine", _modificationType: "Test",
                _target: lysine, _locationRestriction: "N-terminal.", _monoisotopicMass: 42.01057);
            Assert.IsFalse(nTermAcetylK.BlocksCleavage, "alpha-amine N-terminal acetylation does not block cleavage");
        }

        [Test]
        public static void Digestion_PeptidoformEndingInBlockedLysine_IsDroppedWhenRespected()
        {
            var succinyl = MakeKModification("N6-succinyllysine");

            var withoutCorrection = Digest(respectBlockingMods: false, maxMissedCleavages: 2, succinyl);
            var withCorrection = Digest(respectBlockingMods: true, maxMissedCleavages: 2, succinyl);

            // The historical behaviour: PEPTIDEK with a succinylated C-terminal K, reported as an
            // ordinary peptide. This is the chemically impossible entry the issue is about.
            Assert.IsTrue(withoutCorrection.Any(p => p.BaseSequence == "PEPTIDEK" && EndsInModifiedResidue(p)));

            // With the correction it is gone -- no peptidoform ends in a blocked residue at an
            // internal cut anywhere in the digest.
            Assert.IsFalse(withCorrection.Any(p => p.BaseSequence == "PEPTIDEK" && EndsInModifiedResidue(p)));
            Assert.IsFalse(withCorrection.Any(p => EndsInModifiedResidue(p) && p.OneBasedEndResidueInProtein < Sequence.Length));

            // The unmodified form of the same peptide is untouched -- an unmodified K is still a site.
            Assert.IsTrue(withCorrection.Any(p => p.BaseSequence == "PEPTIDEK" && !EndsInModifiedResidue(p)));
        }

        /// <summary>
        /// The contingency that makes dropping safe. The real peptide carrying a succinylated K reads
        /// THROUGH that residue, which costs a missed cleavage under the ordinary modification-blind
        /// count -- so at MaxMissedCleavages = 0 it would never be generated and the peptide would be
        /// lost entirely rather than merely mis-reported. Discounting the blocked site restores it.
        /// </summary>
        [Test]
        public static void Digestion_ZeroMissedCleavages_StillYieldsReadThroughOfABlockedSite()
        {
            var succinyl = MakeKModification("N6-succinyllysine");

            var withCorrection = Digest(respectBlockingMods: true, maxMissedCleavages: 0, succinyl);

            // The read-through, with the blocking modification now internal, is present at 0 missed
            // cleavages because the blocked K no longer counts as a cleavage site for this peptidoform.
            var readThrough = withCorrection.Where(p => p.BaseSequence == "PEPTIDEKAAAAAAAR").ToList();
            Assert.IsTrue(readThrough.Any(p => p.AllModsOneIsNterminus.ContainsKey(9)),
                "the read-through carrying the blocked K must survive at zero missed cleavages");

            // But a genuine missed cleavage is still not allowed: the UNMODIFIED read-through has a
            // real open site inside it and must remain excluded at 0 missed cleavages.
            Assert.IsFalse(readThrough.Any(p => p.AllModsOneIsNterminus.Count == 0),
                "an unmodified read-through is still a real missed cleavage and must not appear");

            // And the impossible form is still gone.
            Assert.IsFalse(withCorrection.Any(p => p.BaseSequence == "PEPTIDEK" && EndsInModifiedResidue(p)));
        }

        /// <summary>
        /// Review concern (nbollis, seconded by Alexander-Sol): with the flag on, does the missed
        /// cleavage count coming OUT of digestion still match the digestion parameters going in?
        ///
        /// It must. A blocked residue is not a cleavage site for the peptidoform carrying it, so it is
        /// not a missed cleavage either -- the count reports cleavages that could have happened and did
        /// not, not Lys/Arg residues. Digestion enumerates a wider span internally to reach the
        /// read-through forms; this pins that the slack is discounted again before a peptide is emitted
        /// and so never reaches the caller.
        /// </summary>
        [Test]
        public static void ReadThroughOfABlockedSite_ReportsZeroMissedCleavages_NotOne()
        {
            var succinyl = MakeKModification("N6-succinyllysine");

            var withCorrection = Digest(respectBlockingMods: true, maxMissedCleavages: 0, succinyl);

            PeptideWithSetModifications readThrough = withCorrection
                .Single(p => p.BaseSequence == "PEPTIDEKAAAAAAAR" && p.AllModsOneIsNterminus.ContainsKey(9));

            // The span physically crosses K8, so the modification-blind count is 1. But K8 is blocked
            // for THIS peptidoform, so no cleavage was missed there and the reported count is 0.
            Assert.That(readThrough.MissedCleavages, Is.EqualTo(0),
                "a blocked site is not a cleavage site, so it is not a missed cleavage either");
        }

        /// <summary>
        /// The general form of the same invariant, across every peptidoform and every budget: the
        /// generation slack must never leak out as a peptide claiming more missed cleavages than the
        /// caller allowed. This is the assertion that fails if the slack is ever widened without
        /// discounting it again on the way out.
        /// </summary>
        [Test]
        public static void NoEmittedPeptide_EverReportsMoreMissedCleavagesThanRequested(
            [Values(0, 1, 2)] int maxMissedCleavages)
        {
            var succinyl = MakeKModification("N6-succinyllysine");

            var withCorrection = Digest(respectBlockingMods: true, maxMissedCleavages, succinyl);

            Assert.That(withCorrection, Is.Not.Empty, "the digestion must actually produce peptides");
            Assert.That(withCorrection.Select(p => p.MissedCleavages), Has.All.LessThanOrEqualTo(maxMissedCleavages),
                $"no peptide may report more than the requested {maxMissedCleavages} missed cleavages");
        }

        /// <summary>
        /// The counterpart: discounting applies only to residues that are actually blocked. An ordinary
        /// missed cleavage still counts, flag on or off, so the fix above cannot be a blanket zeroing.
        /// </summary>
        [Test]
        public static void AnUnblockedMissedCleavage_StillCountsWithTheFlagOn()
        {
            var succinyl = MakeKModification("N6-succinyllysine");

            var withCorrection = Digest(respectBlockingMods: true, maxMissedCleavages: 1, succinyl);

            PeptideWithSetModifications unmodifiedReadThrough = withCorrection
                .Single(p => p.BaseSequence == "PEPTIDEKAAAAAAAR" && p.AllModsOneIsNterminus.Count == 0);

            Assert.That(unmodifiedReadThrough.MissedCleavages, Is.EqualTo(1),
                "an open K is a real missed cleavage and must still be counted");
        }

        [Test]
        public static void Digestion_WithoutTheFlag_ReproducesHistoricalBehaviourExactly()
        {
            var succinyl = MakeKModification("N6-succinyllysine");

            var defaultParams = new DigestionParams(protease: "trypsin", maxMissedCleavages: 2, minPeptideLength: 7);
            Assert.IsFalse(defaultParams.RespectCleavageBlockingModifications, "the flag must default to off");

            var protein = new Protein(Sequence, "accession");
            var flagOff = protein.Digest(defaultParams, new List<Modification>(), new List<Modification> { succinyl })
                .Cast<PeptideWithSetModifications>().ToList();

            // Pin the concrete historical (modification-blind) digest of PEPTIDEKAAAAAAAR with variable
            // succinyl on K, rather than comparing one flag-off run to another (which only re-checks the
            // default). Trypsin gives PEPTIDEK, AAAAAAAR and the 1-missed-cleavage read-through
            // PEPTIDEKAAAAAAAR; the two K-bearing peptides each appear unmodified and succinylated.
            Assert.That(flagOff, Is.Not.Empty);
            CollectionAssert.AreEquivalent(
                new[] { "PEPTIDEK", "PEPTIDEK", "AAAAAAAR", "PEPTIDEKAAAAAAAR", "PEPTIDEKAAAAAAAR" },
                flagOff.Select(p => p.BaseSequence));

            // The decisive check the old assertion missed: flag-off KEEPS the succinyl-on-C-terminal-K
            // PEPTIDEK form -- the exact peptidoform the correction drops -- so this genuinely reproduces
            // pre-PR output instead of matching a second run of the same new code path.
            Assert.That(flagOff.Any(p => p.BaseSequence == "PEPTIDEK" && EndsInModifiedResidue(p)), Is.True,
                "modification-blind digestion keeps the impossible peptidoform");
        }

        [Test]
        public static void Citrulline_BlocksCleavage_OnArginineOnly()
        {
            ModificationMotif.TryGetMotif("R", out ModificationMotif arg);
            var citrullineR = new Modification(_originalId: "Citrulline", _modificationType: "Test",
                _target: arg, _locationRestriction: "Anywhere.", _monoisotopicMass: 0.98402);
            Assert.IsTrue(citrullineR.BlocksCleavage, "citrulline removes arginine's charge, so trypsin does not cleave");

            // The same name on K is not the arginine chemistry and must not classify.
            ModificationMotif.TryGetMotif("K", out ModificationMotif lys);
            var citrullineK = new Modification(_originalId: "Citrulline", _modificationType: "Test",
                _target: lys, _locationRestriction: "Anywhere.", _monoisotopicMass: 0.98402);
            Assert.IsFalse(citrullineK.BlocksCleavage);
        }

        [Test]
        public static void Clone_PreservesRespectCleavageBlockingModifications()
        {
            // Clone() has two return paths (None-specificity and everything else); both must carry the
            // flag or a cloned params silently reverts to modification-blind digestion.
            foreach (CleavageSpecificity mode in new[] { CleavageSpecificity.Full, CleavageSpecificity.None })
            {
                var dp = new DigestionParams(protease: "trypsin", searchModeType: mode,
                    respectCleavageBlockingModifications: true);
                var clone = (DigestionParams)dp.Clone();
                Assert.IsTrue(clone.RespectCleavageBlockingModifications, $"clone lost the flag for {mode}");
                Assert.AreEqual(dp, clone, $"clone compares unequal for {mode}");
            }
        }

        [Test]
        public static void BlockedResidueAtProteinCTerminus_Survives_BecauseNoCleavageHappensThere()
        {
            // The protein's own C-terminal K is not a cleavage site, so an acylation there ends a real
            // peptide and must not be dropped -- the position distinction the C-terminal drop respects.
            var succinyl = MakeKModification("N6-succinyllysine");
            var protein = new Protein("AAAAAAAPEPTIDEK", "accession");   // single K, at the protein C-terminus
            var dp = new DigestionParams(protease: "trypsin", maxMissedCleavages: 2, minPeptideLength: 7,
                respectCleavageBlockingModifications: true);

            var digest = protein.Digest(dp, new List<Modification>(), new List<Modification> { succinyl })
                .Cast<PeptideWithSetModifications>().ToList();

            Assert.That(digest.Any(p => p.BaseSequence == "AAAAAAAPEPTIDEK" && EndsInModifiedResidue(p)), Is.True,
                "an acylation on the protein's terminal K ends a real peptide and must survive");
        }

        /// <summary>
        /// A semi or nonspecific search is left ENTIRELY alone by the flag -- not just its non-Full
        /// subset. The drop and the generation slack are two halves of one exchange (impossible
        /// peptidoform out, read-through form in), and ProteinDigestion grants the slack only when
        /// SearchModeType is Full. Applying the drop without it performed half the trade: at a
        /// missed-cleavage budget too small to reach the read-through, the peptide was dropped and never
        /// replaced, so it became unidentifiable.
        ///
        /// This asserts the WHOLE digest, both specificity categories. The earlier version of this test
        /// compared only the non-Full subset, which is exactly why it passed while the Full-tagged
        /// peptides inside a semi search were being dropped unreplaced.
        /// </summary>
        [Test]
        public static void SemiSpecificSearch_IsEntirelyUnaffectedByTheFlag(
            [Values(0, 1, 2)] int maxMissedCleavages)
        {
            var succinyl = MakeKModification("N6-succinyllysine");
            var protein = new Protein("PEPTIDEKAAAAAAAR", "accession");

            List<string> WholeDigest(bool respect) =>
                protein.Digest(new DigestionParams(protease: "trypsin", maxMissedCleavages: maxMissedCleavages,
                        minPeptideLength: 7, searchModeType: CleavageSpecificity.Semi,
                        respectCleavageBlockingModifications: respect),
                        new List<Modification>(), new List<Modification> { succinyl })
                    .Cast<PeptideWithSetModifications>()
                    .Select(p => p.FullSequence).OrderBy(s => s).ToList();

            CollectionAssert.AreEqual(WholeDigest(respect: false), WholeDigest(respect: true),
                "a semi search must be untouched by the flag, not half-corrected");
        }

        /// <summary>
        /// The counterpart in Full mode, stated as the exchange it is: the impossible peptidoform leaves
        /// and the read-through that replaces it arrives, so the peptide remains identifiable. This is
        /// what a semi search does NOT get, and the pair of tests marks the boundary.
        /// </summary>
        [Test]
        public static void FullSearch_TradesTheImpossibleFormForTheReadThrough()
        {
            var succinyl = MakeKModification("N6-succinyllysine");

            var without = Digest(respectBlockingMods: false, maxMissedCleavages: 0, succinyl);
            var with = Digest(respectBlockingMods: true, maxMissedCleavages: 0, succinyl);

            Assert.Multiple(() =>
            {
                Assert.That(without.Any(p => p.BaseSequence == "PEPTIDEK" && EndsInModifiedResidue(p)), Is.True,
                    "flag off keeps the impossible form");
                Assert.That(with.Any(p => p.BaseSequence == "PEPTIDEK" && EndsInModifiedResidue(p)), Is.False,
                    "flag on drops the impossible form");
                Assert.That(with.Any(p => p.BaseSequence == "PEPTIDEKAAAAAAAR" && p.AllModsOneIsNterminus.ContainsKey(9)), Is.True,
                    "and replaces it with the read-through, so the peptide is not simply lost");
            });
        }

        [Test]
        public static void MultipleBlockedSites_ReadThroughSurvives_WithSlackTiedToMaxMods()
        {
            // Three succinylated K's in a row: the read-through spanning all of them costs three missed
            // cleavages, so a fixed slack of 2 would lose it. Slack = MaxMods (here 3) reaches it.
            var succinyl = MakeKModification("N6-succinyllysine");
            var protein = new Protein("AAAAKAAAAKAAAAKAAAAR", "accession");   // K at 5, 10, 15; R at 20
            var dp = new DigestionParams(protease: "trypsin", maxMissedCleavages: 0, minPeptideLength: 7,
                maxModsForPeptides: 3, respectCleavageBlockingModifications: true);

            var digest = protein.Digest(dp, new List<Modification>(), new List<Modification> { succinyl })
                .Cast<PeptideWithSetModifications>().ToList();

            // The full-length read-through carrying all three succinyl-K's must survive at 0 missed cleavages.
            Assert.That(
                digest.Any(p => p.BaseSequence == "AAAAKAAAAKAAAAKAAAAR" && p.AllModsOneIsNterminus.Count == 3),
                Is.True, "the read-through past three blocked sites must survive when slack = MaxMods");
        }
        /// <summary>
        /// The classifier's three defensive guards. A modification reaches this code straight from a
        /// database read, so a null instance, a null Target and a nameless modification are all shapes
        /// it has to survive rather than throw on -- and each returns false, since nothing that cannot
        /// be identified can be asserted to block a cleavage.
        /// </summary>
        [Test]
        public static void NeutralizesCleavageResidue_GuardsAgainstMalformedModifications()
        {
            Assert.IsFalse(CleavageBlockingModifications.NeutralizesCleavageResidue(null),
                "a null modification cannot block a cleavage");

            var noTarget = new Modification(_originalId: "N6-succinyllysine", _modificationType: "Test",
                _target: null, _locationRestriction: "Anywhere.", _monoisotopicMass: 100.01604);
            Assert.IsFalse(noTarget.BlocksCleavage,
                "a blocking NAME with no target motif must not classify -- the residue is unknown");

            // OriginalId null leaves IdWithMotif null too, so the classifier has no name to match on.
            ModificationMotif.TryGetMotif("K", out ModificationMotif lysine);
            var nameless = new Modification(_originalId: null, _modificationType: "Test",
                _target: lysine, _locationRestriction: "Anywhere.", _monoisotopicMass: 100.01604);
            Assert.IsFalse(nameless.BlocksCleavage,
                "a modification with no id cannot be classified as blocking");
        }

        /// <summary>
        /// Citrulline is the common name, but the reaction is deimination and databases ship both
        /// spellings; the arginine branch accepts either. An arginine modification matching neither
        /// must still be false, which is what stops the branch degrading into "any Arg mod blocks".
        /// </summary>
        [Test]
        public static void Deimination_BlocksCleavage_AndAnUnrelatedArginineModDoesNot()
        {
            ModificationMotif.TryGetMotif("R", out ModificationMotif arg);

            var deiminated = new Modification(_originalId: "Deimination", _modificationType: "Test",
                _target: arg, _locationRestriction: "Anywhere.", _monoisotopicMass: 0.98402);
            Assert.IsTrue(deiminated.BlocksCleavage, "deimination is citrullination under its other name");

            var methylArg = new Modification(_originalId: "Methyl", _modificationType: "Test",
                _target: arg, _locationRestriction: "Anywhere.", _monoisotopicMass: 14.01565);
            Assert.IsFalse(methylArg.BlocksCleavage,
                "methylarginine keeps the guanidinium charge, so trypsin still cleaves");
        }

        /// <summary>
        /// With the flag on, a variable modification that does NOT block cleavage has to be ignored by
        /// the blocking rule rather than discounting a missed cleavage. Methylation is the case that
        /// matters: it sits on the same residue trypsin cuts at, so a rule keyed on position rather
        /// than on <see cref="Modification.BlocksCleavage"/> would silently treat it as blocking.
        /// </summary>
        [Test]
        public static void ANonBlockingModOnACleavageResidue_DoesNotDiscountAMissedCleavage()
        {
            var methyl = MakeKModification("N6-methyllysine", mass: 14.01565);

            var digest = Digest(respectBlockingMods: true, maxMissedCleavages: 1, methyl);

            // The methylated read-through still spans the K, and that K is still a real missed cleavage.
            PeptideWithSetModifications methylatedReadThrough = digest
                .Single(p => p.BaseSequence == "PEPTIDEKAAAAAAAR" && p.AllModsOneIsNterminus.Count == 1);
            Assert.That(methylatedReadThrough.MissedCleavages, Is.EqualTo(1),
                "a non-blocking modification must not discount the missed cleavage at that residue");

            // And the peptidoform ending in methyl-K survives, because methylation does not abolish the site.
            Assert.That(digest.Any(p => p.BaseSequence == "PEPTIDEK" && EndsInModifiedResidue(p)), Is.True,
                "a methylated C-terminal K is a cleavage trypsin can still perform");
        }

        /// <summary>
        /// Review round 4 (acesnik). The acyl test has to run BEFORE any methyl exclusion, because
        /// Unimod ships names carrying both chemistries on the same lysine. Testing "methyl" first
        /// classified "Methyl+Acetyl:2H(3)" and "Methyl:2H(3)+Acetyl:2H(3)" -- both site K, position
        /// Anywhere in unimod.xml -- as non-blocking, when the acetyl on the epsilon-amine is exactly
        /// the chemistry this classifier exists for.
        /// </summary>
        [Test]
        public static void AnAcylationIsBlocking_EvenWhenItsNameAlsoCarriesAMethyl()
        {
            Assert.IsTrue(MakeKModification("Methyl+Acetyl:2H(3)").BlocksCleavage,
                "the acetyl on the epsilon-amine decides; the methyl in the name must not veto it");
            Assert.IsTrue(MakeKModification("Methyl:2H(3)+Acetyl:2H(3)").BlocksCleavage);

            // And the exclusion the reorder must not break: methylation with no acyl group is still
            // non-blocking, because the epsilon-amine keeps its charge.
            Assert.IsFalse(MakeKModification("Methyl").BlocksCleavage);
            Assert.IsFalse(MakeKModification("Dimethyl").BlocksCleavage);
            Assert.IsFalse(MakeKModification("Trimethyl").BlocksCleavage);
        }

        /// <summary>
        /// Review round 4 (acesnik). Of the four members that had to carry the flag, ToString was the
        /// one no test asserted -- Equals was covered only indirectly, through Clone's AreEqual. It
        /// matters downstream: MetaMorpheus keys its fragment-index cache on IndexingEngine.ToString()
        /// and SameSettings compares that string verbatim, so a flag missing from this stringification
        /// would let a flag-off index be silently reused for a flag-on run.
        /// </summary>
        [Test]
        public static void DigestionParams_DifferingOnlyInTheFlag_AreDistinguishableByEveryIdentityMember()
        {
            var off = new DigestionParams(protease: "trypsin", respectCleavageBlockingModifications: false);
            var on = new DigestionParams(protease: "trypsin", respectCleavageBlockingModifications: true);

            Assert.AreNotEqual(off.ToString(), on.ToString(),
                "the fragment-index cache key must distinguish the two settings");
            Assert.IsFalse(off.Equals(on), "params differing in the flag must not compare equal");
            Assert.AreNotEqual(off.GetHashCode(), on.GetHashCode());

            // The other direction, so the assertions above cannot pass by everything being unequal.
            var offAgain = new DigestionParams(protease: "trypsin", respectCleavageBlockingModifications: false);
            Assert.AreEqual(off.ToString(), offAgain.ToString());
            Assert.IsTrue(off.Equals(offAgain));
            Assert.AreEqual(off.GetHashCode(), offAgain.GetHashCode());
        }

        /// <summary>
        /// Review round 4 (acesnik). The generation slack widens enumeration on every protein whenever
        /// the flag and Full mode are set, and it was being paid in searches where no configured
        /// modification could block anything -- at MetaMorpheus defaults, enumerating at 4 missed
        /// cleavages instead of 2, roughly 1.7x the unmodified peptides before modification
        /// combinatorics, and carrying into the fragment index.
        ///
        /// The cost is invisible in the finished digest, because the open-site filter trims the surplus
        /// again on the way out: an output comparison passes either way and pins nothing. So this reads
        /// the enumeration itself, through ProteinDigestion.Digestion, which is where the slack is spent.
        /// </summary>
        [Test]
        public static void WithNoBlockingModificationConfigured_NoGenerationSlackIsPaid()
        {
            ModificationMotif.TryGetMotif("S", out ModificationMotif serine);
            var phospho = new Modification(_originalId: "Phospho", _modificationType: "Test", _target: serine,
                _locationRestriction: "Anywhere.", _monoisotopicMass: 79.96633);
            var methyl = MakeKModification("N6-methyllysine", mass: 14.01565);
            var succinyl = MakeKModification("N6-succinyllysine");

            var protein = new Protein("AAAAKAAAAKAAAASAAAAKAAAAKAAAAR", "accession");

            int Enumerated(bool respect, params Modification[] variableMods)
            {
                var dp = new DigestionParams(protease: "trypsin", maxMissedCleavages: 2, minPeptideLength: 5,
                    respectCleavageBlockingModifications: respect);
                return new ProteinDigestion(dp, new List<Modification>(), variableMods.ToList())
                    .Digestion(protein).Count();
            }

            int baseline = Enumerated(respect: false, phospho, methyl);

            Assert.That(Enumerated(respect: true, phospho, methyl), Is.EqualTo(baseline),
                "nothing configured can block a trypsin cleavage, so no slack may be bought");

            // The counterpart, so the gate cannot pass by disabling the slack outright: once a blocking
            // modification IS configured, the slack is bought and enumeration widens.
            Assert.That(Enumerated(respect: true, phospho, succinyl), Is.GreaterThan(baseline),
                "a configured blocking modification must still buy the slack the read-through needs");
        }

        /// <summary>
        /// And the finished digest is unchanged too -- the invariant a consumer sees, as distinct from
        /// the enumeration cost pinned above.
        /// </summary>
        [Test]
        public static void WithNoBlockingModificationConfigured_TheDigestIsUnchanged()
        {
            ModificationMotif.TryGetMotif("S", out ModificationMotif serine);
            var phospho = new Modification(_originalId: "Phospho", _modificationType: "Test", _target: serine,
                _locationRestriction: "Anywhere.", _monoisotopicMass: 79.96633);
            var methyl = MakeKModification("N6-methyllysine", mass: 14.01565);

            var protein = new Protein("PEPTIDEKAAASAAAR", "accession");

            List<string> WholeDigest(bool respect) =>
                protein.Digest(new DigestionParams(protease: "trypsin", maxMissedCleavages: 2, minPeptideLength: 7,
                        respectCleavageBlockingModifications: respect),
                        new List<Modification>(), new List<Modification> { phospho, methyl })
                    .Cast<PeptideWithSetModifications>()
                    .Select(p => p.FullSequence + "|" + p.MissedCleavages).OrderBy(x => x).ToList();

            CollectionAssert.AreEqual(WholeDigest(respect: false), WholeDigest(respect: true),
                "no configured modification can block a trypsin cleavage, so the flag must change nothing");
        }

        /// <summary>
        /// Whether a modification blocks a cleavage is a question about the PAIR, not about the
        /// modification on its own. <see cref="DigestionAgent.CleavesCTerminalTo"/> is the protease half.
        /// </summary>
        [Test]
        public static void CleavesCTerminalTo_ReadsTheProteasesOwnMotifs()
        {
            Assert.IsTrue(ProteaseDictionary.Dictionary["trypsin"].CleavesCTerminalTo('K'));
            Assert.IsTrue(ProteaseDictionary.Dictionary["trypsin"].CleavesCTerminalTo('R'));
            Assert.IsFalse(ProteaseDictionary.Dictionary["trypsin"].CleavesCTerminalTo('E'));

            // Lys-C cuts after K only, Arg-C after R only -- each is blind to the other's blocking mod.
            Assert.IsTrue(ProteaseDictionary.Dictionary["Lys-C|P"].CleavesCTerminalTo('K'));
            Assert.IsFalse(ProteaseDictionary.Dictionary["Lys-C|P"].CleavesCTerminalTo('R'));
            Assert.IsTrue(ProteaseDictionary.Dictionary["Arg-C"].CleavesCTerminalTo('R'));
            Assert.IsFalse(ProteaseDictionary.Dictionary["Arg-C"].CleavesCTerminalTo('K'));

            // Glu-C cuts after E and nowhere near a lysine.
            Assert.IsTrue(ProteaseDictionary.Dictionary["Glu-C"].CleavesCTerminalTo('E'));
            Assert.IsFalse(ProteaseDictionary.Dictionary["Glu-C"].CleavesCTerminalTo('K'));

            // Direction: Asp-N and Lys-N cut N-terminal to their recognition residue, so they sever
            // nothing C-terminal to it and report false for every residue -- which is what leaves the
            // C-terminal-side correction inert for them rather than silently mirror-imaged.
            Assert.IsFalse(ProteaseDictionary.Dictionary["Asp-N"].CleavesCTerminalTo('D'));
            Assert.IsFalse(ProteaseDictionary.Dictionary["Lys-N"].CleavesCTerminalTo('K'));

            // The wildcard is honoured through the same matcher digestion itself uses.
            Assert.IsTrue(ProteaseDictionary.Dictionary["non-specific"].CleavesCTerminalTo('K'));
        }

        /// <summary>
        /// The consequence for digestion, and the reason the protease half is not optional. Glu-C never
        /// cleaves after a lysine, so an acylated lysine abolishes nothing in a Glu-C digest. Counting
        /// it as blocked discounted a missed cleavage the peptide genuinely has -- both under-reporting
        /// the count and letting an over-budget peptide through on the generation slack.
        /// </summary>
        [Test]
        public static void AProteaseThatDoesNotCleaveAfterTheModifiedResidue_IsUnaffected(
            [Values("Glu-C", "Arg-C")] string protease)
        {
            var succinyl = MakeKModification("N6-succinyllysine");
            var protein = new Protein("AAAAAAAEAAKAAAAAEAAAAAAA", "accession");   // E8, K11, E17

            List<string> WholeDigest(bool respect) =>
                protein.Digest(new DigestionParams(protease: protease, maxMissedCleavages: 1, minPeptideLength: 7,
                        respectCleavageBlockingModifications: respect),
                        new List<Modification>(), new List<Modification> { succinyl })
                    .Cast<PeptideWithSetModifications>()
                    .Select(p => p.FullSequence + "|" + p.MissedCleavages).OrderBy(x => x).ToList();

            CollectionAssert.AreEqual(WholeDigest(respect: false), WholeDigest(respect: true),
                protease + " does not cleave after K, so a blocked K must not discount its missed cleavages");
        }

        /// <summary>
        /// The positive control for the same gate: Lys-C does cleave after K, so the correction applies
        /// there exactly as it does for trypsin. Without this, the protease gate could pass by simply
        /// disabling the feature everywhere.
        /// </summary>
        [Test]
        public static void AProteaseThatDoesCleaveAfterTheModifiedResidue_StillGetsTheCorrection()
        {
            var succinyl = MakeKModification("N6-succinyllysine");
            var protein = new Protein("PEPTIDEKAAAAAAAK", "accession");   // Lys-C cuts after K8; K16 is the protein C-terminus
            var dp = new DigestionParams(protease: "Lys-C|P", maxMissedCleavages: 0, minPeptideLength: 7,
                respectCleavageBlockingModifications: true);

            var digest = protein.Digest(dp, new List<Modification>(), new List<Modification> { succinyl })
                .Cast<PeptideWithSetModifications>().ToList();

            Assert.That(digest.Any(p => p.BaseSequence == "PEPTIDEK" && EndsInModifiedResidue(p)), Is.False,
                "Lys-C cannot cleave after a succinylated K either, so the impossible form must be dropped");
            Assert.That(digest.Any(p => p.BaseSequence == "PEPTIDEKAAAAAAAK" && p.AllModsOneIsNterminus.ContainsKey(9)), Is.True,
                "and the read-through must replace it at zero missed cleavages");
        }

        /// <summary>
        /// The same gate at residue granularity rather than protease granularity: Lys-C cleaves after K
        /// but not after R, so a citrullinated arginine inside a Lys-C peptide blocks nothing and must
        /// not discount that peptide's missed cleavage.
        /// </summary>
        [Test]
        public static void ABlockingModOnAResidueThisProteaseIgnores_DoesNotDiscountAMissedCleavage()
        {
            ModificationMotif.TryGetMotif("R", out ModificationMotif arg);
            var citrulline = new Modification(_originalId: "Citrulline", _modificationType: "Test", _target: arg,
                _locationRestriction: "Anywhere.", _monoisotopicMass: 0.98402);

            var protein = new Protein("AAAAAAAKAARAAAAAKAAAAAAA", "accession");   // K8, R11, K17
            var dp = new DigestionParams(protease: "Lys-C|P", maxMissedCleavages: 1, minPeptideLength: 7,
                respectCleavageBlockingModifications: true);

            var digest = protein.Digest(dp, new List<Modification>(), new List<Modification> { citrulline })
                .Cast<PeptideWithSetModifications>().ToList();

            // The span crossing K8 with a citrulline on R11 really does miss the cleavage at K8: Lys-C
            // never had a site at R11 for the citrulline to abolish.
            PeptideWithSetModifications readThrough = digest
                .Single(p => p.BaseSequence == "AAAAAAAKAARAAAAAK" && p.AllModsOneIsNterminus.Count == 1);
            Assert.That(readThrough.MissedCleavages, Is.EqualTo(1),
                "citrulline on R blocks nothing under Lys-C, so the missed cleavage at K stands");
        }
    }
}
