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
    }
}
