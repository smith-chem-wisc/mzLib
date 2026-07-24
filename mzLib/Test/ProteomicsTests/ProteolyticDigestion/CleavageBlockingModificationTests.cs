using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
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

            // A modification that is not on a cleavage residue cannot block a cleavage.
            ModificationMotif.TryGetMotif("S", out ModificationMotif serine);
            var phosphoSerine = new Modification(_originalId: "Phospho", _modificationType: "Test",
                _target: serine, _locationRestriction: "Anywhere.", _monoisotopicMass: 79.96633);
            Assert.IsFalse(phosphoSerine.BlocksCleavage);
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

        [Test]
        public static void Digestion_WithoutTheFlag_ReproducesHistoricalBehaviourExactly()
        {
            var succinyl = MakeKModification("N6-succinyllysine");

            // Default is off, and the explicit false must agree with it, so existing callers and every
            // downstream search are unaffected until they opt in.
            var defaultParams = new DigestionParams(protease: "trypsin", maxMissedCleavages: 2, minPeptideLength: 7);
            Assert.IsFalse(defaultParams.RespectCleavageBlockingModifications);

            var protein = new Protein(Sequence, "accession");
            var byDefault = protein.Digest(defaultParams, new List<Modification>(), new List<Modification> { succinyl })
                .Cast<PeptideWithSetModifications>().Select(p => p.FullSequence).OrderBy(s => s).ToList();
            var explicitlyOff = Digest(respectBlockingMods: false, maxMissedCleavages: 2, succinyl)
                .Select(p => p.FullSequence).OrderBy(s => s).ToList();

            CollectionAssert.AreEqual(byDefault, explicitlyOff);
        }
    }
}
