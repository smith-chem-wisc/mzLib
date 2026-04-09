using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Linq;

namespace Test.CircularSearch
{
    [TestFixture]
    public static class CircularSearchEngineTests
    {
        /// <summary>
        /// Verifies the exact digestion products of CircularProtein "PEPTIDEKAAK"
        /// using trypsin, 2 missed cleavages, min peptide length 1.
        ///
        /// The circular protein canonicalises to "AAKPEPTIDEK" (alphabetically
        /// earliest rotation).  Trypsin then produces:
        ///   Linear:   AAK  (pos 1-3)
        ///   Linear:   PEPTIDEK  (pos 4-11)
        ///   Linear:   AAKPEPTIDEK  (pos 1-11, 1 missed cleavage)
        ///   Linear:   PEPTIDEKAAK  (pos 4-14, 1 missed cleavage, wraps)
        ///   Circular: AAKPEPTIDEK  (two-cut cyclic product)
        /// </summary>
        [Test]
        public static void Digest_PEPTIDEKAAK_ProducesExactExpectedProducts()
        {
            var circularProtein = CircularProtein.FromProtein(
                new Protein("PEPTIDEKAAK", accession: "TEST01"));

            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 2,
                minPeptideLength: 1);

            var allProducts = circularProtein
                .Digest(digestionParams,
                        new List<Modification>(),
                        new List<Modification>())
                .ToList();

            // ── Linear products ───────────────────────────────────────────────

            var linearProducts = allProducts
                .OfType<PeptideWithSetModifications>()
                .Where(p => p is not CircularPeptideWithSetModifications)
                .ToList();

            Assert.That(linearProducts.Count, Is.EqualTo(4),
                "Expected exactly 4 linear (PeptideWithSetModifications) products.");

            var aak = linearProducts.SingleOrDefault(p => p.BaseSequence == "AAK");
            Assert.That(aak, Is.Not.Null, "Expected linear product AAK.");
            Assert.That(aak.OneBasedStartResidue, Is.EqualTo(1));
            Assert.That(aak.OneBasedEndResidue, Is.EqualTo(3));

            var peptidek = linearProducts.SingleOrDefault(p => p.BaseSequence == "PEPTIDEK");
            Assert.That(peptidek, Is.Not.Null, "Expected linear product PEPTIDEK.");
            Assert.That(peptidek.OneBasedStartResidue, Is.EqualTo(4));
            Assert.That(peptidek.OneBasedEndResidue, Is.EqualTo(11));

            var aakpeptidek = linearProducts.SingleOrDefault(p => p.BaseSequence == "AAKPEPTIDEK");
            Assert.That(aakpeptidek, Is.Not.Null, "Expected linear product AAKPEPTIDEK.");
            Assert.That(aakpeptidek.OneBasedStartResidue, Is.EqualTo(1));
            Assert.That(aakpeptidek.OneBasedEndResidue, Is.EqualTo(11));

            var peptidekaak = linearProducts.SingleOrDefault(p => p.BaseSequence == "PEPTIDEKAAK");
            Assert.That(peptidekaak, Is.Not.Null, "Expected linear product PEPTIDEKAAK.");
            Assert.That(peptidekaak.OneBasedStartResidue, Is.EqualTo(4));
            Assert.That(peptidekaak.OneBasedEndResidue, Is.EqualTo(14));

            // ── Circular product ──────────────────────────────────────────────

            var circularProducts = allProducts
                .OfType<CircularPeptideWithSetModifications>()
                .ToList();

            Assert.That(circularProducts.Count, Is.EqualTo(1),
                "Expected exactly 1 CircularPeptideWithSetModifications.");

            var circPeptide = circularProducts.Single();
            Assert.That(circPeptide.BaseSequence, Is.EqualTo("AAKPEPTIDEK"));

            Assert.That(allProducts.Count, Is.EqualTo(5),
                "Expected exactly 5 digestion products in total.");
        }

        /// <summary>
        /// Same CircularProtein and DigestionParams as above, but with variable
        /// phosphorylation on T.
        ///
        /// Canonical sequence: A-A-K-P-E-P-T-I-D-E-K  (N=11)
        ///                     1-2-3-4-5-6-7-8-9-10-11
        /// T is at canonical position 7.
        ///
        /// Key convention: AllModsOneIsNterminus key = localIndex + 2,
        /// where localIndex = (canonical position of T) - (OneBasedStartResidue), 0-based.
        ///
        /// Products and phosphorylation keys:
        ///
        ///   AAK            [1-3]   no T → 1 form, no phospho
        ///   PEPTIDEK       [4-11]  T at canonical 7: localIndex=7-4=3, key=5  → 2 forms
        ///   AAKPEPTIDEK    [1-11]  T at canonical 7: localIndex=7-1=6, key=8  → 2 forms
        ///   PEPTIDEKAAK    [4-14]  T at canonical 7: localIndex=7-4=3, key=5  → 2 forms
        ///   Circular AAKPEPTIDEK [1-11]  T at canonical 7: localIndex=7-1=6, key=8  → 2 forms
        ///
        /// Total: 1 + 2 + 2 + 2 + 2 = 9 products.
        /// </summary>
        [Test]
        public static void Digest_PEPTIDEKAAK_WithPhosphoT_ProducesCorrectProducts()
        {
            ModificationMotif.TryGetMotif("T", out ModificationMotif tMotif);
            var phospho = new Modification(
                _originalId: "Phosphorylation on T",
                _modificationType: "Common Variable",
                _target: tMotif,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: 79.966331);

            var circularProtein = CircularProtein.FromProtein(
                new Protein("PEPTIDEKAAK", accession: "TEST01"));

            Assert.That(circularProtein.BaseSequence, Is.EqualTo("AAKPEPTIDEK"));
            Assert.That(circularProtein.BaseSequence[6], Is.EqualTo('T'),
                "T is at 1-based canonical position 7 (0-based index 6).");

            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 2,
                minPeptideLength: 1);

            var allProducts = circularProtein
                .Digest(digestionParams,
                        new List<Modification>(),
                        new List<Modification> { phospho })
                .ToList();

            Assert.That(allProducts.Count, Is.EqualTo(9),
                "Expected 9 total products: 5 unmodified + 4 phosphorylated (one per T-containing product).");

            // ── Linear products ───────────────────────────────────────────────

            var linearProducts = allProducts
                .OfType<PeptideWithSetModifications>()
                .Where(p => p is not CircularPeptideWithSetModifications)
                .ToList();

            Assert.That(linearProducts.Count, Is.EqualTo(7),
                "Expected 7 linear products: AAK(1) + PEPTIDEK(2) + AAKPEPTIDEK(2) + PEPTIDEKAAK(2).");

            // AAK [1-3] — no T → exactly 1 form, no phosphorylation
            var aakForms = linearProducts.Where(p => p.BaseSequence == "AAK").ToList();
            Assert.That(aakForms.Count, Is.EqualTo(1), "AAK has no T: 1 form only.");
            Assert.That(aakForms[0].AllModsOneIsNterminus, Is.Empty,
                "AAK: no modification expected.");

            // PEPTIDEK [4-11] — T at canonical pos 7, localIndex=3, key=5
            var peptidekForms = linearProducts.Where(p => p.BaseSequence == "PEPTIDEK").ToList();
            Assert.That(peptidekForms.Count, Is.EqualTo(2), "PEPTIDEK: 2 forms (±phospho).");
            var peptidekUnmod = peptidekForms.Single(p => p.AllModsOneIsNterminus.Count == 0);
            var peptidekPhospho = peptidekForms.Single(p => p.AllModsOneIsNterminus.Count == 1);
            Assert.That(peptidekPhospho.AllModsOneIsNterminus.ContainsKey(5), Is.True,
                "PEPTIDEK phospho key: T at canonical 7, start 4 → localIndex=3, key=5.");
            Assert.That(peptidekPhospho.AllModsOneIsNterminus[5].OriginalId,
                Does.Contain("Phosphorylation"));
            Assert.That(
                peptidekPhospho.MonoisotopicMass - peptidekUnmod.MonoisotopicMass,
                Is.EqualTo(79.966331).Within(1e-5),
                "PEPTIDEK: phospho mass shift = +79.966331 Da.");

            // AAKPEPTIDEK [1-11] — T at canonical pos 7, localIndex=6, key=8
            var aakpeptidekForms = linearProducts.Where(p => p.BaseSequence == "AAKPEPTIDEK").ToList();
            Assert.That(aakpeptidekForms.Count, Is.EqualTo(2), "AAKPEPTIDEK: 2 forms (±phospho).");
            var aakpeptidekUnmod = aakpeptidekForms.Single(p => p.AllModsOneIsNterminus.Count == 0);
            var aakpeptidekPhospho = aakpeptidekForms.Single(p => p.AllModsOneIsNterminus.Count == 1);
            Assert.That(aakpeptidekPhospho.AllModsOneIsNterminus.ContainsKey(8), Is.True,
                "AAKPEPTIDEK phospho key: T at canonical 7, start 1 → localIndex=6, key=8.");
            Assert.That(aakpeptidekPhospho.AllModsOneIsNterminus[8].OriginalId,
                Does.Contain("Phosphorylation"));
            Assert.That(
                aakpeptidekPhospho.MonoisotopicMass - aakpeptidekUnmod.MonoisotopicMass,
                Is.EqualTo(79.966331).Within(1e-5),
                "AAKPEPTIDEK: phospho mass shift = +79.966331 Da.");

            // PEPTIDEKAAK [4-14] — T at canonical pos 7, localIndex=3, key=5
            var peptidekaakForms = linearProducts.Where(p => p.BaseSequence == "PEPTIDEKAAK").ToList();
            Assert.That(peptidekaakForms.Count, Is.EqualTo(2), "PEPTIDEKAAK: 2 forms (±phospho).");
            var peptidekaakUnmod = peptidekaakForms.Single(p => p.AllModsOneIsNterminus.Count == 0);
            var peptidekaakPhospho = peptidekaakForms.Single(p => p.AllModsOneIsNterminus.Count == 1);
            Assert.That(peptidekaakPhospho.AllModsOneIsNterminus.ContainsKey(5), Is.True,
                "PEPTIDEKAAK phospho key: T at canonical 7, start 4 → localIndex=3, key=5.");
            Assert.That(peptidekaakPhospho.AllModsOneIsNterminus[5].OriginalId,
                Does.Contain("Phosphorylation"));
            Assert.That(
                peptidekaakPhospho.MonoisotopicMass - peptidekaakUnmod.MonoisotopicMass,
                Is.EqualTo(79.966331).Within(1e-5),
                "PEPTIDEKAAK: phospho mass shift = +79.966331 Da.");

            // ── Circular products ─────────────────────────────────────────────
            // AAKPEPTIDEK [1-11] — T at canonical pos 7, localIndex=6, key=8

            var circularProducts = allProducts
                .OfType<CircularPeptideWithSetModifications>()
                .ToList();

            Assert.That(circularProducts.Count, Is.EqualTo(2),
                "Expected 2 circular products: unmodified and phosphorylated AAKPEPTIDEK.");

            var circUnmod = circularProducts.Single(p => p.AllModsOneIsNterminus.Count == 0);
            var circPhospho = circularProducts.Single(p => p.AllModsOneIsNterminus.Count == 1);

            Assert.That(circUnmod.BaseSequence, Is.EqualTo("AAKPEPTIDEK"));
            Assert.That(circPhospho.BaseSequence, Is.EqualTo("AAKPEPTIDEK"));

            Assert.That(circPhospho.AllModsOneIsNterminus.ContainsKey(8), Is.True,
                "Circular AAKPEPTIDEK phospho key: T at canonical 7, start 1 → localIndex=6, key=8.");
            Assert.That(circPhospho.AllModsOneIsNterminus[8].OriginalId,
                Does.Contain("Phosphorylation"));
            Assert.That(
                circPhospho.MonoisotopicMass - circUnmod.MonoisotopicMass,
                Is.EqualTo(79.966331).Within(1e-5),
                "Circular AAKPEPTIDEK: phospho mass shift = +79.966331 Da.");
        }
    }
}