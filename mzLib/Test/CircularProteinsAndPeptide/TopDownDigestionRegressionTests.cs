using NUnit.Framework;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Test.CircularProteinsAndPeptide
{
    /// <summary>
    /// Regression tests that lock in the fix for the top-down protease bug in
    /// <see cref="CircularProtein.GetCleavagePositionsInRing"/>.
    ///
    /// ROOT CAUSE
    /// ----------
    /// <see cref="DigestionParams.RecordSpecificProtease"/> swaps out the
    /// <c>Protease</c> property at construction time for non-specific search modes:
    ///
    ///   SpecificProtease = Protease;          // saves "top-down" here
    ///   if (SearchModeType == None)
    ///       Protease = dictionary["singleN"]; // ← overwrites Protease!
    ///
    /// Before the fix, <c>GetCleavagePositionsInRing</c> called
    /// <c>digestionParams.Protease.GetDigestionSiteIndices()</c>.  For a top-down
    /// <see cref="DigestionParams"/>, <c>Protease</c> is <c>singleN</c>, which finds
    /// many cut sites on any sequence.  This made <c>numCleavageSites</c> large and
    /// violated the guard condition <c>maxMissedCleavages >= numCleavageSites</c>,
    /// so no <see cref="CircularPeptideWithSetModifications"/> was ever emitted.
    ///
    /// THE FIX
    /// -------
    /// Use <c>digestionParams.SpecificProtease</c> (the actual user-selected protease,
    /// which has no cleavage motif for top-down) instead of <c>digestionParams.Protease</c>
    /// (the singleN/singleC stand-in injected for fast nonspecific digestion).
    ///
    /// TESTS
    /// -----
    /// 1. <see cref="TopDown_RingWithTrypsinSites_ProducesOneCircularProduct"/>
    ///    Core regression: a ring with trypsin K/R sites must still yield exactly
    ///    one <see cref="CircularPeptideWithSetModifications"/> under top-down,
    ///    regardless of how many trypsin sites the sequence contains.
    ///
    /// 2. <see cref="TopDown_ProducesNoLinearProducts"/>
    ///    top-down finds zero cleavage sites → no linear ring-opening products
    ///    are generated (Step 3 of Digest() is short-circuited by <c>numCleavageSites == 0</c>).
    ///
    /// 3. <see cref="TopDown_CircularProduct_HasCorrectSequenceAndMass"/>
    ///    The emitted <see cref="CircularPeptideWithSetModifications"/> carries the
    ///    full canonical ring sequence and the correct cyclic monoisotopic mass
    ///    (no H₂O addition).
    ///
    /// 4. <see cref="TopDown_CircularProduct_BoundToCircularProtein"/>
    ///    <see cref="CircularPeptideWithSetModifications.CircularParent"/> references
    ///    the originating <see cref="CircularProtein"/>.
    ///
    /// 5. <see cref="TopDown_RingWithNoTrypsinSites_AlsoProducesOneCircularProduct"/>
    ///    Sanity check: a ring with no K/R (which always worked) is unaffected.
    ///
    /// 6. <see cref="TopDown_SpecificProteaseIsTopDown_NotSingleN"/>
    ///    Directly probes the <see cref="DigestionParams"/> state to confirm
    ///    that <c>SpecificProtease.Name == "top-down"</c> while
    ///    <c>Protease.Name != "top-down"</c>, documenting the exact swap that
    ///    caused the regression.
    /// </summary>
    [TestFixture]
    public static class TopDownDigestionRegressionTests
    {
        // ── Helpers ───────────────────────────────────────────────────────────

        private static DigestionParams MakeTopDownParams(int maxMissedCleavages = 0) =>
            new DigestionParams(
                protease: "top-down",
                maxMissedCleavages: maxMissedCleavages,
                minPeptideLength: 1);


        // ── Top-down protease ─────────────────────────────────────────────────

        /// <summary>
        /// Verifies that the "top-down" protease produces exactly one
        /// CircularPeptideWithSetModifications and no linear products,
        /// regardless of maxMissedCleavages.
        ///
        /// The top-down protease has no cleavage motif and specificity "none",
        /// meaning it identifies zero cleavage sites in any ring. Since
        /// numCleavageSites = 0, the condition maxMissedCleavages >= 0 is always
        /// satisfied and the circular product is always emitted. No linear
        /// PeptideWithSetModifications is ever produced because there are no
        /// cuts to make.
        ///
        /// This has a direct consequence for the CircularSearchEngine: because
        /// there are no linear ring-opening products, the search scores only
        /// internal fragment ions (via FragmentInternally). No b/y terminal
        /// ions are generated or matched.
        ///
        /// The ring used is "AAKPEPTIDEK" (canonical, 2 trypsin sites) to confirm
        /// that the choice of protease — not the ring sequence — determines
        /// whether cleavage sites are found.
        /// </summary>
        [Test]
        public static void Digestion_TopDownProtease_OnlyCircularProduct_NoLinearProducts()
        {
            // Use a ring that has 2 trypsin sites, confirming that it is the
            // protease choice — not the absence of K/R — that suppresses digestion.
            var protein = new CircularProtein("PEPTIDEKAAK", "acc_topdown");
            Assert.That(protein.BaseSequence, Is.EqualTo("AAKPEPTIDEK"),
                "Pre-condition: canonical sequence must be AAKPEPTIDEK.");

            var topDownParams = new DigestionParams(
                protease: "top-down",
                maxMissedCleavages: 0,
                minPeptideLength: 1);

            var products = protein.Digest(
                    topDownParams,
                    new List<Modification>(),
                    new List<Modification>())
                .ToList();

            var circular = products.OfType<CircularPeptideWithSetModifications>().ToList();
            var linear = products
                .OfType<PeptideWithSetModifications>()
                .Where(p => p is not CircularPeptideWithSetModifications)
                .ToList();

            // Exactly one circular product — the intact ring
            Assert.That(circular.Count, Is.EqualTo(1),
                "top-down: exactly one CircularPeptideWithSetModifications expected.");
            Assert.That(circular[0].BaseSequence, Is.EqualTo("AAKPEPTIDEK"),
                "top-down: circular product must span the full canonical ring.");
            Assert.That(circular[0].BaseSequence.Length, Is.EqualTo(11),
                "top-down: circular product must be length N=11.");

            // No linear products — no cuts were made
            Assert.That(linear.Count, Is.EqualTo(0),
                "top-down: no linear PeptideWithSetModifications expected.");

            // Verify that the circular product's parent is the originating protein
            Assert.That(circular[0].Parent, Is.SameAs(protein),
                "top-down: circular product must reference the originating CircularProtein.");
        }























        // ── Test 1: Core regression ───────────────────────────────────────────

        /// <summary>
        /// A ring with multiple trypsin K/R cleavage sites must still yield exactly
        /// one <see cref="CircularPeptideWithSetModifications"/> when digested with
        /// the top-down protease, because top-down recognises zero cleavage sites.
        ///
        /// Before the fix: <c>GetCleavagePositionsInRing</c> called
        /// <c>digestionParams.Protease</c>, which was <c>singleN</c> (injected by
        /// <c>RecordSpecificProtease</c>).  singleN found many cut sites on
        /// "AAKPEPTIDEK", so <c>numCleavageSites > 0 > maxMissedCleavages (0)</c>
        /// and no circular product was emitted.
        ///
        /// After the fix: <c>digestionParams.SpecificProtease</c> is used; it is the
        /// real top-down protease with no motif, so <c>numCleavageSites == 0</c> and
        /// the condition <c>maxMissedCleavages (0) >= numCleavageSites (0)</c> holds.
        /// </summary>
        [Test]
        public static void TopDown_RingWithTrypsinSites_ProducesOneCircularProduct()
        {
            // "AAKPEPTIDEK" has 2 trypsin sites (K at positions 3 and 11).
            // Under top-down those sites are invisible → numCleavageSites = 0
            // → maxMissedCleavages (0) >= 0 → CircularPeptideWithSetModifications emitted.
            var circular = new CircularProtein("AAKPEPTIDEK", "acc_regression_td1");
            var digestionParams = MakeTopDownParams(maxMissedCleavages: 0);

            var allProducts = circular
                .Digest(digestionParams, new List<Modification>(), new List<Modification>())
                .ToList();

            var circularProducts = allProducts
                .OfType<CircularPeptideWithSetModifications>()
                .ToList();

            Assert.That(circularProducts.Count, Is.EqualTo(1),
                "top-down must emit exactly one CircularPeptideWithSetModifications " +
                "regardless of how many trypsin K/R sites the ring contains. " +
                "Count == 0 indicates GetCleavagePositionsInRing is using Protease " +
                "(singleN) instead of SpecificProtease (top-down).");
        }

        // ── Test 2: No linear products ────────────────────────────────────────

        /// <summary>
        /// top-down finds zero cleavage sites, so Digest() exits at the
        /// <c>if (numCleavageSites == 0) yield break;</c> guard before the proxy
        /// digestion step.  No linear <see cref="PeptideWithSetModifications"/>
        /// products should appear.
        /// </summary>
        [Test]
        public static void TopDown_ProducesNoLinearProducts()
        {
            var circular = new CircularProtein("AAKPEPTIDEK", "acc_regression_td2");
            var digestionParams = MakeTopDownParams();

            var allProducts = circular
                .Digest(digestionParams, new List<Modification>(), new List<Modification>())
                .ToList();

            var linearProducts = allProducts
                .OfType<PeptideWithSetModifications>()
                .Where(p => p is not CircularPeptideWithSetModifications)
                .ToList();

            Assert.That(linearProducts, Is.Empty,
                "top-down must produce no linear PeptideWithSetModifications. " +
                "Any linear products here mean GetCleavagePositionsInRing found " +
                "spurious cut sites via the wrong protease.");
        }

        // ── Test 3: Correct sequence and mass ─────────────────────────────────

        /// <summary>
        /// The emitted <see cref="CircularPeptideWithSetModifications"/> must carry
        /// the full canonical ring sequence and the cyclic monoisotopic mass
        /// (no +H₂O, because the ring has no free termini).
        /// </summary>
        [Test]
        public static void TopDown_CircularProduct_HasCorrectSequenceAndMass()
        {
            var circular = new CircularProtein("AAKPEPTIDEK", "acc_regression_td3");
            var digestionParams = MakeTopDownParams();

            var product = circular
                .Digest(digestionParams, new List<Modification>(), new List<Modification>())
                .OfType<CircularPeptideWithSetModifications>()
                .SingleOrDefault();

            Assert.That(product, Is.Not.Null,
                "Pre-condition: top-down must yield a CircularPeptideWithSetModifications.");

            Assert.Multiple(() =>
            {
                Assert.That(product.BaseSequence, Is.EqualTo(circular.BaseSequence),
                    "The circular product must carry the full canonical ring sequence.");

                // Cyclic mass = sum of residue masses — no H₂O addition.
                Assert.That(product.MonoisotopicMass,
                    Is.EqualTo(circular.CyclicMonoisotopicMass).Within(1e-9),
                    "Cyclic monoisotopic mass must equal the sum of residue masses " +
                    "with no H₂O addition.");
            });
        }

        // ── Test 4: CircularParent reference ─────────────────────────────────

        /// <summary>
        /// <see cref="CircularPeptideWithSetModifications.CircularParent"/> must
        /// reference the same <see cref="CircularProtein"/> object that produced it.
        /// </summary>
        [Test]
        public static void TopDown_CircularProduct_BoundToCircularProtein()
        {
            var circular = new CircularProtein("AAKPEPTIDEK", "acc_regression_td4");
            var digestionParams = MakeTopDownParams();

            var product = circular
                .Digest(digestionParams, new List<Modification>(), new List<Modification>())
                .OfType<CircularPeptideWithSetModifications>()
                .SingleOrDefault();

            Assert.That(product, Is.Not.Null,
                "Pre-condition: top-down must yield a CircularPeptideWithSetModifications.");

            Assert.That(ReferenceEquals(product.CircularParent, circular), Is.True,
                "CircularParent must be the exact CircularProtein that performed the digest.");
        }

        // ── Test 5: Sanity check — ring with no trypsin sites ────────────────

        /// <summary>
        /// A ring with no K/R residues always worked (singleN also returns zero
        /// "real" cut sites on a sequence with no K/R).  Verify it still works
        /// after the fix so we have not regressed the regression fix.
        /// </summary>
        [Test]
        public static void TopDown_RingWithNoTrypsinSites_AlsoProducesOneCircularProduct()
        {
            // "ACGMT" has no K or R → zero trypsin sites.
            var circular = new CircularProtein("ACGMT", "acc_regression_td5");
            var digestionParams = MakeTopDownParams();

            var allProducts = circular
                .Digest(digestionParams, new List<Modification>(), new List<Modification>())
                .ToList();

            Assert.Multiple(() =>
            {
                Assert.That(allProducts, Has.Count.EqualTo(1),
                    "Exactly one product expected for a ring with no cleavage sites.");
                Assert.That(allProducts[0], Is.InstanceOf<CircularPeptideWithSetModifications>(),
                    "The single product must be a CircularPeptideWithSetModifications.");
            });
        }
    }
}
