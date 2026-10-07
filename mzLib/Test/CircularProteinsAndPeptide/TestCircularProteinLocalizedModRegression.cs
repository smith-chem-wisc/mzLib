using NUnit.Framework;
using Omics;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test.CircularSearch
{
    /// <summary>
    /// Regression tests for localized-modification handling in
    /// <see cref="CircularProtein.Digest"/>. Each test is constructed to FAIL on the
    /// pre-fix code and pass after the corresponding fix:
    ///
    ///   • Wrap-around fragments must carry XML/database-localized mods on the wrapped
    ///     (second-copy) residues — the proxy now replicates localized-mod keys onto
    ///     positions N+1..2N-1.
    ///   • Intact-ring forms must include annotated localized mods even when they are
    ///     not supplied via <c>variableModifications</c>, matching the linear path.
    ///   • The full-ring variable-mod enumeration must be bounded by
    ///     <see cref="DigestionParams.MaxModificationIsoforms"/>.
    ///   • <c>NumFixedMods</c> must not count a fixed mod that a variable mod overrides
    ///     at the same residue.
    /// </summary>
    [TestFixture]
    public static class TestCircularProteinLocalizedModRegression
    {
        private const double PhosphoMass = 79.966331;
        private const double FixedSerineMass = 100.0;
        private const double Tol = 1e-3;

        [OneTimeSetUp]
        public static void OneTimeSetup()
        {
            Loaders.LoadElements();
        }

        // ── Helpers ─────────────────────────────────────────────────────────────

        private static Modification BuildPhospho()
        {
            ModificationMotif.TryGetMotif("S", out ModificationMotif motifS);
            return new Modification(
                _originalId: "Phospho",
                _modificationType: "Test",
                _target: motifS,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: PhosphoMass);
        }

        private static Modification BuildFixedSerineMod()
        {
            ModificationMotif.TryGetMotif("S", out ModificationMotif motifS);
            return new Modification(
                _originalId: "FixedSerineMod",
                _modificationType: "Test",
                _target: motifS,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: FixedSerineMass);
        }

        /// <summary>
        /// Builds a <see cref="CircularProtein"/> from an already-canonical sequence
        /// (so positions are preserved through <see cref="CircularProtein.FromProtein"/>)
        /// carrying a single localized modification at <paramref name="oneBasedPos"/>.
        /// </summary>
        private static CircularProtein BuildCircularWithLocalizedMod(
            string canonicalSequence, int oneBasedPos, Modification mod)
        {
            var mods = new Dictionary<int, List<Modification>>
            {
                { oneBasedPos, new List<Modification> { mod } }
            };
            var source = new Protein(canonicalSequence, "TESTACC", oneBasedModifications: mods);
            return CircularProtein.FromProtein(source);
        }

        private static bool CarriesPhospho(IBioPolymerWithSetMods peptide) =>
            peptide.AllModsOneIsNterminus.Values.Any(m =>
                m.MonoisotopicMass.HasValue &&
                Math.Abs(m.MonoisotopicMass.Value - PhosphoMass) < Tol);

        // ── Fix A (Critical / Codex #1): wrap-around localized mods ───────────────

        /// <summary>
        /// Ring "ASDKEFGHR" (N=9): trypsin cuts after K(4) and R(9). With one missed
        /// cleavage the proxy yields the wrap-around product "EFGHRASDK" (ring 5..9 + 1..4),
        /// whose only serine is the wrapped residue at ring position 2 (proxy position 11).
        /// The phospho is annotated as a LOCALIZED mod (not passed as a variable mod), so it
        /// is applied by position; the wrapped S can only receive it if the proxy carries the
        /// localized-mod key at N+2 = 11. Pre-fix, that key was absent and the mod was dropped.
        /// </summary>
        [Test]
        public static void WrapAroundFragment_CarriesLocalizedModOnWrappedResidue()
        {
            var phospho = BuildPhospho();
            var circular = BuildCircularWithLocalizedMod("ASDKEFGHR", oneBasedPos: 2, phospho);

            var dp = new DigestionParams(protease: "trypsin", maxMissedCleavages: 2, minPeptideLength: 1);

            var wrapPeptides = circular
                .Digest(dp, new List<Modification>(), new List<Modification>())
                .OfType<PeptideWithSetModifications>()
                .Where(p => p.BaseSequence == "EFGHRASDK")
                .ToList();

            Assert.That(wrapPeptides, Is.Not.Empty,
                "Expected the wrap-around product \"EFGHRASDK\" to be produced.");
            Assert.That(wrapPeptides.Any(CarriesPhospho), Is.True,
                "The wrap-around peptide must carry the annotated phosphoserine on its wrapped S residue " +
                "(proxy position 11). Pre-fix the localized mod was silently dropped on the wrapped portion.");
        }

        // ── Fix B (Codex #2): full-ring includes annotated localized mods ─────────

        /// <summary>
        /// Ring "ALLLLALPGS" (N=10) has no K/R, so trypsin yields only the intact-ring
        /// product. The phospho is annotated on the S (position 10) but NOT passed as a
        /// variable mod. The linear path applies OneBasedPossibleLocalizedModifications as
        /// variable mods; the full-ring branch must too, or it drops the modified isoform.
        /// </summary>
        [Test]
        public static void FullRing_IncludesAnnotatedLocalizedMod_WithoutVariableModsArg()
        {
            var phospho = BuildPhospho();
            var circular = BuildCircularWithLocalizedMod("ALLLLALPGS", oneBasedPos: 10, phospho);

            var dp = new DigestionParams(protease: "trypsin", maxMissedCleavages: 2, minPeptideLength: 1);

            var rings = circular
                .Digest(dp, new List<Modification>(), new List<Modification>())
                .OfType<CircularPeptideWithSetModifications>()
                .ToList();

            Assert.That(rings, Is.Not.Empty, "Expected the intact-ring product.");
            Assert.That(rings.Any(CarriesPhospho), Is.True,
                "The intact-ring form must include the annotated localized phospho isoform " +
                "even though it was not supplied via variableModifications.");
        }

        // ── Fix C (Major): full-ring variable-mod isoform cap ─────────────────────

        /// <summary>
        /// Ring "ASASASLLPG" (N=10) has three serines (positions 2,4,6) and no K/R. With
        /// phospho as a variable mod and MaxMods=2 there are six distinct modified ring
        /// forms; MaxModificationIsoforms=2 must bound how many are emitted.
        /// </summary>
        [Test]
        public static void FullRing_VariableModForms_BoundedByMaxModificationIsoforms()
        {
            var phospho = BuildPhospho();
            var circular = CircularProtein.FromProtein(new Protein("ASASASLLPG", "TESTACC"));

            var dp = new DigestionParams(protease: "trypsin", maxMissedCleavages: 2,
                minPeptideLength: 1, maxModificationIsoforms: 2, maxModsForPeptides: 2);

            var phosphoRingForms = circular
                .Digest(dp, new List<Modification>(), new List<Modification> { phospho })
                .OfType<CircularPeptideWithSetModifications>()
                .Count(CarriesPhospho);

            Assert.That(phosphoRingForms, Is.GreaterThan(0),
                "At least one modified ring form should be produced.");
            Assert.That(phosphoRingForms, Is.LessThanOrEqualTo(2),
                "Modified ring forms must be capped at MaxModificationIsoforms (2); " +
                "pre-fix the enumeration was uncapped and produced six.");
        }

        // ── Fix D (Major): NumFixedMods not overcounted on override ───────────────

        /// <summary>
        /// Ring "ASLLLLLLPG" (N=10) has a single serine (position 2) and no K/R. A fixed
        /// mod and a variable phospho both target that S. In the ring form where the variable
        /// phospho occupies the S, the fixed mod is overridden, so NumFixedMods must be 0.
        /// </summary>
        [Test]
        public static void FullRing_NumFixedMods_NotOvercounted_WhenVariableOverridesFixed()
        {
            var phospho = BuildPhospho();
            var fixedSerine = BuildFixedSerineMod();
            var circular = CircularProtein.FromProtein(new Protein("ASLLLLLLPG", "TESTACC"));

            var dp = new DigestionParams(protease: "trypsin", maxMissedCleavages: 2, minPeptideLength: 1);

            var rings = circular
                .Digest(dp, new List<Modification> { fixedSerine }, new List<Modification> { phospho })
                .OfType<CircularPeptideWithSetModifications>()
                .ToList();

            var phosphoForm = rings.FirstOrDefault(CarriesPhospho);

            Assert.That(phosphoForm, Is.Not.Null,
                "Expected a ring form where the serine carries the variable phospho.");
            Assert.That(phosphoForm.NumFixedMods, Is.EqualTo(0),
                "When the variable phospho overrides the fixed mod at the same serine, " +
                "NumFixedMods must not count the overridden fixed mod (pre-fix it reported 1).");
        }
    }
}
