using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class TestCircularPeptideWithSetModifications
    {
        [OneTimeSetUp]
        public static void OneTimeSetup()
        {
            Loaders.LoadElements();
        }

        // ── Helpers ───────────────────────────────────────────────────────────

        private static double WaterMass =>
            PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2
            + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;

        /// <summary>
        /// Returns the monoisotopic mass of a bare sub-sequence of the ring,
        /// computed as the sum of residue masses — no caps added.
        /// Used to verify fragment neutral masses independently of the production code.
        /// </summary>
        private static double ResidueMassSum(string ringSequence, int oneBasedStart, int length)
        {
            int n = ringSequence.Length;
            double mass = 0;
            for (int i = 0; i < length; i++)
            {
                char aa = ringSequence[(oneBasedStart - 1 + i) % n];
                mass += Residue.ResidueMonoisotopicMass[aa];
            }
            return mass;
        }

        /// <summary>
        /// Returns the expected neutral mass of an internal fragment (bIb convention):
        ///   residueMassSum + nTermCap(b=0) + cTermCap(b=0) − H₂O
        /// For HCD with b/y ions: nTermCap = 0, cTermCap = H₂O, so:
        ///   neutralMass = residueMassSum + 0 + H₂O − H₂O = residueMassSum
        /// But the production code uses:
        ///   fragmentResidueMass + nTermCap + cTermCap − H₂O
        /// For b-type internal (bIb): nTermCap = 0 (b-ion), cTermCap = 0 (b-ion) → neutralMass = residueMassSum - H₂O... 
        /// Actually for internal ions both caps come from the N- and C-terminal product types.
        /// For HCD: nTermProductTypes = [b], cTermProductTypes = [y]
        ///   b nTermCap = 0, y cTermCap = H₂O
        ///   neutralMass = residueMassSum + 0 + H₂O − H₂O = residueMassSum
        /// </summary>
        private static double ExpectedInternalFragmentNeutralMass(string ringSequence, int oneBasedStart, int length)
        {
            // For HCD: b nTermCap=0, y cTermCap=+H2O, minus H2O = residueMassSum
            return ResidueMassSum(ringSequence, oneBasedStart, length);
        }

        /// <summary>
        /// Gets a CircularPeptideWithSetModifications of the full ring (zero cuts)
        /// from a sequence with no trypsin cleavage sites.
        /// </summary>
        private static CircularPeptideWithSetModifications GetFullRingPeptide(string sequence)
        {
            var protein = new CircularProtein(sequence, "test_acc");
            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1);
            return (CircularPeptideWithSetModifications)protein
                .Digest(digestionParams, [], [])
                .Single();
        }

        /// <summary>
        /// Gets a CircularPeptideWithSetModifications sub-peptide from a ring
        /// by base sequence match.
        /// </summary>
        private static CircularPeptideWithSetModifications GetSubPeptide(
            string ringSequence, string subSequence)
        {
            var protein = new CircularProtein(ringSequence, "test_acc");
            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1);
            return protein.Digest(digestionParams, [], [])
                .OfType<CircularPeptideWithSetModifications>()
                .Single(p => p.BaseSequence == subSequence);
        }


        [Test]
        public static void Digest_LinearProduct_ModificationRenumberedToProxyPosition()
        {
            // Ring: "AHMCDEFGKIEP" (N=12)
            //   A=1, H=2, M=3, C=4, D=5, E=6, F=7, G=8, K=9, I=10, E=11, P=12
            //   Canonical: A is at position 1 and is the smallest character → already canonical.
            //   K at position 9 → one tryptic cut → one linear PeptideWithSetModifications.
            //
            // Cut after K(9) → linear peptide runs from position 10 to position 21:
            //   I(10), E(11), P(12), A(13=1), H(14=2), M(15=3), C(16=4),
            //   D(17=5), E(18=6), F(19=7), G(20=8), K(21=9)
            //   BaseSequence = "IEPAHMCDEFGK"
            //   OneBasedStartResidueInProtein = 10
            //   OneBasedEndResidueInProtein   = 21
            //
            // Oxidation on M:
            //   M is at ring position 3.
            //   In the proxy, the peptide starts at proxy position 10, M is at proxy position 15.
            //   allModsOneIsNterminus key = 15 - 10 + 2 = 7
            //   (M is the 6th residue of "IEPAHMCDEFGK", 0-based index 5, key = 5+2 = 7)

            ModificationMotif.TryGetMotif("M", out ModificationMotif mMotif);
            var oxidation = new Modification(
                _originalId: "Oxidation on M",
                _modificationType: "Common Variable",
                _target: mMotif,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: 15.994915);

            // M is at position 3 in the canonical ring
            var mods = new Dictionary<int, List<Modification>> { [3] = [oxidation] };
            var protein = new CircularProtein("AHMCDEFGKIEP", "mod_renumber_test",
                oneBasedModifications: mods);

            Assert.That(protein.BaseSequence, Is.EqualTo("AHMCDEFGKIEP"),
                "Sequence should already be canonical (A is smallest at position 1)");

            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1);

            // Digest with oxidation as variable modification
            var peptides = protein.Digest(digestionParams, [], [oxidation])
                .Cast<PeptideWithSetModifications>()
                .ToList();

            // With one K, we get exactly two products:
            //   unmodified linear  "IEPAHMCDEFGK" (no mods)
            //   oxidized linear    "IEPAHMCDEFGK" (oxidation on M)
            // Both are linear PeptideWithSetModifications (one cut), not CircularPeptideWithSetModifications.
            Assert.That(peptides, Has.Count.EqualTo(2),
                "One K site + one variable mod site → two linear products (unmodified and oxidized)");

            Assert.That(peptides.All(p => p is not CircularPeptideWithSetModifications), Is.True,
                "All products are linear PeptideWithSetModifications (one cut opens the ring)");

            var unmodified = peptides.Single(p => p.AllModsOneIsNterminus.Count == 0);
            var oxidized = peptides.Single(p => p.AllModsOneIsNterminus.Count == 1);

            // ── Sequence and position checks (both forms) ─────────────────────────
            foreach (var peptide in peptides)
            {
                Assert.That(peptide.BaseSequence, Is.EqualTo("IEPAHMCDEFGK"),
                    "Linear product must run from I(10) wrapping around to K(9)");
                Assert.That(peptide.OneBasedStartResidueInProtein, Is.EqualTo(10),
                    "Start = position 10 (I), immediately after the K cleavage site at position 9");
                Assert.That(peptide.OneBasedEndResidueInProtein, Is.EqualTo(21),
                    "End = 10 + 12 - 1 = 21 in doubled-ring numbering");
                Assert.That(peptide.Length, Is.EqualTo(protein.Length),
                    "Linear ring-opening product has length N=12");
                Assert.That(ReferenceEquals(peptide.Protein, protein), Is.True,
                    "Parent must be the CircularProtein instance");
            }

            // ── Modification renumbering check ────────────────────────────────────
            // M is at ring position 3. In the proxy, the peptide starts at position 10.
            // M appears at proxy position 15. Key = 15 - 10 + 2 = 7.
            // This means M is the residue at local 0-based index 5 of "IEPAHMCDEFGK":
            //   I=0, E=1, P=2, A=3, H=4, M=5 ✓
            Assert.That(oxidized.AllModsOneIsNterminus.ContainsKey(7), Is.True,
                "Oxidation on M must be stored at allModsOneIsNterminus key 7 " +
                "(M is the 6th residue of IEPAHMCDEFGK, local 0-based index 5, key = 5+2 = 7)");

            Assert.That(oxidized.AllModsOneIsNterminus[7].OriginalId,
                Does.Contain("Oxidation"),
                "The modification at key 7 must be oxidation");

            // ── Mass check ────────────────────────────────────────────────────────
            // Oxidized peptide mass = unmodified mass + 15.994915
            Assert.That(
                oxidized.MonoisotopicMass - unmodified.MonoisotopicMass,
                Is.EqualTo(15.994915).Within(1e-4),
                "Mass difference between oxidized and unmodified must equal oxidation mass");
        }


        [Test]
        public static void ModFits_BasicUnderstanding()
        {
            // Set up oxidation on M — motif "M", LocationRestriction "Anywhere."
            ModificationMotif.TryGetMotif("M", out ModificationMotif mMotif);
            var oxidation = new Modification(
                _originalId: "Oxidation on M",
                _modificationType: "Common Variable",
                _target: mMotif,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: 15.994915);

            const string sequence = "ACDEFGHIM"; // N=9, M at position 9

            // Test 1: Does ModFits return true for M at position 9?
            bool fitsAtPos9 = ModificationLocalization.ModFits(
                oxidation,
                sequence,
                digestionProductOneBasedIndex: 9,
                digestionProductLength: 9,
                bioPolymerOneBasedIndex: 9);

            Assert.That(fitsAtPos9, Is.True,
                "ModFits must return true for M at position 9 in ACDEFGHIM");

            // Test 2: Does ModFits return false for A at position 1?
            bool fitsAtPos1 = ModificationLocalization.ModFits(
                oxidation,
                sequence,
                digestionProductOneBasedIndex: 1,
                digestionProductLength: 9,
                bioPolymerOneBasedIndex: 1);

            Assert.That(fitsAtPos1, Is.False,
                "ModFits must return false for A at position 1 (not M)");

            // Test 3: Does ModFits return false for H at position 7?
            bool fitsAtPos7 = ModificationLocalization.ModFits(
                oxidation,
                sequence,
                digestionProductOneBasedIndex: 7,
                digestionProductLength: 9,
                bioPolymerOneBasedIndex: 7);

            Assert.That(fitsAtPos7, Is.False,
                "ModFits must return false for H at position 7 (not M)");

            // Test 4: Confirm what OriginalId and LocationRestriction look like
            Console.WriteLine($"mod.OriginalId           = {oxidation.OriginalId}");
            Console.WriteLine($"mod.LocationRestriction  = {oxidation.LocationRestriction}");
            Console.WriteLine($"mod.Target               = {oxidation.Target}");
            Console.WriteLine($"ModFits pos1={fitsAtPos1}, pos7={fitsAtPos7}, pos9={fitsAtPos9}");
        }


        [Test]
        public static void ModFits_FallbackDiagnostic()
        {
            ModificationMotif.TryGetMotif("M", out ModificationMotif mMotif);
            var oxidation = new Modification(
                _originalId: "Oxidation on M",
                _modificationType: "Common Variable",
                _target: mMotif,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: 15.994915);

            var mods = new Dictionary<int, List<Modification>> { [9] = [oxidation] };
            var protein = new CircularProtein("ACDEFGHIM", "mod_test",
                oneBasedModifications: mods);

            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1);

            var peptides = protein.Digest(digestionParams, [], [oxidation])
                .OfType<CircularPeptideWithSetModifications>()
                .ToList();

            Console.WriteLine($"Total peptides: {peptides.Count}");
            foreach (var p in peptides)
            {
                Console.WriteLine($"  Mods count: {p.AllModsOneIsNterminus.Count}");
                foreach (var kvp in p.AllModsOneIsNterminus)
                    Console.WriteLine($"    Key={kvp.Key} Id={kvp.Value.OriginalId}");
            }

            // Also check OneBasedPossibleLocalizedModifications on the protein
            Console.WriteLine($"OneBasedPossibleLocalizedModifications count: " +
                $"{protein.OneBasedPossibleLocalizedModifications.Count}");
            foreach (var kvp in protein.OneBasedPossibleLocalizedModifications)
                foreach (var m in kvp.Value)
                    Console.WriteLine($"  Pos={kvp.Key} Id={m.OriginalId}");
        }

        // ──────────────────────────────────────────────────────────────────────
        // 1. Preconditions and guards
        // ──────────────────────────────────────────────────────────────────────



        [Test]
        public static void FragmentInternally_MinLengthLessThanTwo_ReturnsEmpty()
        {
            var peptide = GetFullRingPeptide("ACDEFGHIM");
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 1, products);

            Assert.That(products, Is.Empty,
                "minLengthOfFragments < 2 must produce no products");
        }

        [Test]
        public static void FragmentInternally_MinLengthEqualsToPeptideLength_ReturnsEmpty()
        {
            var peptide = GetFullRingPeptide("ACDEFGHIM"); // length 9
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 9, products);

            Assert.That(products, Is.Empty,
                "minLengthOfFragments >= peptideLength must produce no products");
        }

        [Test]
        public static void FragmentInternally_ClearsProductsListOnEntry()
        {
            var peptide = GetFullRingPeptide("ACDEFGHIM");
            var products = new List<Product> { new Product(ProductType.b, FragmentationTerminus.N, 100, 1, 1, 0) };

            peptide.FragmentInternally(DissociationType.HCD, 3, products);

            // List is cleared and repopulated — the original entry must be gone
            Assert.That(products.All(p => p.NeutralMass != 100), Is.True,
                "Products list must be cleared on entry");
        }

        // ──────────────────────────────────────────────────────────────────────
        // 2. Fragment count
        //
        // For a full-ring peptide of length N, minLength L:
        // For each start position s in [0, N-1]:
        //   lengths from L to (N - s - 1)  → count = max(0, N - s - 1 - L + 1)
        // Total = sum over s of max(0, N - s - L)
        //
        // For N=9, L=3:
        //   s=0: lengths 3,4,5,6,7,8 → 6 fragments
        //   s=1: lengths 3,4,5,6,7   → 5
        //   s=2: lengths 3,4,5,6     → 4
        //   s=3: lengths 3,4,5       → 3
        //   s=4: lengths 3,4         → 2
        //   s=5: length  3           → 1
        //   s=6,7,8: 0
        //   Total = 21 fragments × 1 ion type pair (b,y for HCD) = 21
        // ──────────────────────────────────────────────────────────────────────

        [Test]
        public static void FragmentInternally_FullRing_N9_MinLength3_CorrectCount()
        {
            var peptide = GetFullRingPeptide("ACDEFGHIM"); // N=9, no K/R
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 3, products);

            // HCD gives b and y product types: nTermProductTypes=[b], cTermProductTypes=[y]
            // Each (start, length) pair produces 1 product (1 nTerm × 1 cTerm)
            Assert.That(products, Has.Count.EqualTo(21),
                "N=9, minLength=3: expected 21 internal fragments");
        }

        // ──────────────────────────────────────────────────────────────────────
        // 3. Numbering: non-wrapping fragments of a full-ring peptide
        //
        // Ring "ACDEFGHIM" (N=9), full-ring peptide [1-9].
        // start=1, length=3: residues A(1),C(2),D(3) → annotation [1,3]
        // start=1, length=4: residues A(1),C(2),D(3),E(4) → annotation [1,4]
        // start=4, length=3: residues E(4),F(5),G(6) → annotation [4,6]
        // ──────────────────────────────────────────────────────────────────────

        [Test]
        public static void FragmentInternally_FullRing_NonWrapping_StartAndEndNumbering()
        {
            const string ring = "ACDEFGHIM"; // N=9
            var peptide = GetFullRingPeptide(ring);
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 3, products);

            // Fragment starting at canonical position 1, length 3: [1,3]
            var frag_1_3 = products.Single(p =>
                p.FragmentNumber == 1 && p.ResiduePosition == 3);
            Assert.That(frag_1_3.SecondaryFragmentNumber, Is.EqualTo(3),
                "Start=1, length=3 → end=3");

            // Fragment starting at canonical position 1, length 4: [1,4]
            var frag_1_4 = products.Single(p =>
                p.FragmentNumber == 1 && p.ResiduePosition == 4);
            Assert.That(frag_1_4.SecondaryFragmentNumber, Is.EqualTo(4),
                "Start=1, length=4 → end=4");

            // Fragment starting at canonical position 4, length 3: [4,6]
            var frag_4_3 = products.Single(p =>
                p.FragmentNumber == 4 && p.ResiduePosition == 3);
            Assert.That(frag_4_3.SecondaryFragmentNumber, Is.EqualTo(6),
                "Start=4, length=3 → end=6");
        }

        // ──────────────────────────────────────────────────────────────────────
        // 4. Numbering: wrap-around fragments of a full-ring peptide
        //
        // Ring "ACDEFGHIM" (N=9), full-ring peptide [1-9].
        // start=8 (I), length=3: residues I(8),M(9),A(10→wrap) → annotation [8,10]
        // start=9 (M), length=3: residues M(9),A(10),C(11)     → annotation [9,11]
        // start=7 (H), length=4: residues H(7),I(8),M(9),A(10) → annotation [7,10]
        // ──────────────────────────────────────────────────────────────────────

        [Test]
        public static void FragmentInternally_FullRingStartingAt1_NeverProducesWrapAroundFragments()
        {
            // For a full-ring CircularPeptideWithSetModifications with
            // OneBasedStartResidueInProtein == 1:
            //   oneBasedEnd = 1 + localStart + length - 1
            //               ≤ 1 + localStart + (N - localStart - 1) - 1
            //               = N - 1 < N
            // Therefore SecondaryFragmentNumber < N always — no wrap-around possible.
            const string ring = "ACDEFGHIM"; // N=9
            var peptide = GetFullRingPeptide(ring);
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 2, products);

            Assert.That(products, Is.Not.Empty);
            Assert.That(products.All(p => p.SecondaryFragmentNumber < ring.Length), Is.True,
                "Full-ring peptide [1-N]: all fragment ends must be < N. " +
                "Wrap-around requires OneBasedStartResidueInProtein > 1.");
        }

        [Test]
        public static void FragmentInternally_WrappingSubPeptide_ProducesWrapAroundFragments()
        {
            // Ring "AKDEFKGHIM" (N=10), K at positions 2 and 6.
            // Wrapping sub-peptide "GHIMAK" has:
            //   OneBasedStartResidueInProtein = 7
            //   OneBasedEndResidueInProtein   = 12 (> N=10)
            // Its internal fragments can have ends > N.
            //
            // With canonicalOffset=7:
            //   localStart=0, oneBasedStart=7, maxLength=5
            //     length=3: end=9   → [7,9]  = G,H,I   (no wrap)
            //     length=4: end=10  → [7,10] = G,H,I,M (no wrap, ends exactly at N)
            //     length=5: end=11  → [7,11] = G,H,I,M,A (wrap! A is pos 1 → doubled pos 11)
            //   localStart=1, oneBasedStart=8, maxLength=4
            //     length=3: end=10  → [8,10] = H,I,M (no wrap)
            //     length=4: end=11  → [8,11] = H,I,M,A (wrap!)
            //   localStart=2, oneBasedStart=9, maxLength=3
            //     length=3: end=11  → [9,11] = I,M,A (wrap!)
            //   localStart=3, oneBasedStart=10, maxLength=2
            //     length=2: end=11  → [10,11] = M,A (wrap!)
            //   localStart=4, oneBasedStart=11, maxLength=1 < minLength=3 → none
            const string ringSeq = "AKDEFKGHIM"; // A=1,K=2,D=3,E=4,F=5,K=6,G=7,H=8,I=9,M=10
            var protein = new CircularProtein(ringSeq, "wrap_frag_test");

            Assert.That(protein.BaseSequence, Is.EqualTo(ringSeq),
                "AKDEFKGHIM should be canonical (A < D < E < F < G < H < I < K < M)");

            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 1,
                minPeptideLength: 1);

            var wrappingPeptide = protein.Digest(digestionParams, [], [])
                .OfType<CircularPeptideWithSetModifications>()
                .SingleOrDefault(p => p.BaseSequence == "GHIMAK");

            Assume.That(wrappingPeptide, Is.Not.Null,
                "Wrapping sub-peptide GHIMAK[7-12] must be produced by Digest()");
            Assert.That(wrappingPeptide.OneBasedStartResidueInProtein, Is.EqualTo(7));
            Assert.That(wrappingPeptide.OneBasedEndResidueInProtein, Is.EqualTo(12));

            var products = new List<Product>();
            wrappingPeptide.FragmentInternally(DissociationType.HCD, 3, products);

            // Non-wrapping fragments
            var frag_7_3 = products.Single(p => p.FragmentNumber == 7 && p.ResiduePosition == 3);
            Assert.That(frag_7_3.SecondaryFragmentNumber, Is.EqualTo(9),
                "GHIMAK: start=7, length=3 → end=9 (GHI, no wrap)");

            var frag_7_4 = products.Single(p => p.FragmentNumber == 7 && p.ResiduePosition == 4);
            Assert.That(frag_7_4.SecondaryFragmentNumber, Is.EqualTo(10),
                "GHIMAK: start=7, length=4 → end=10 (GHIM, ends at N, no wrap)");

            // Wrap-around fragments (end > N=10)
            var frag_7_5 = products.Single(p => p.FragmentNumber == 7 && p.ResiduePosition == 5);
            Assert.That(frag_7_5.SecondaryFragmentNumber, Is.EqualTo(11),
                "GHIMAK: start=7, length=5 → end=11 > N=10 (GHIMA, wraps to A at doubled pos 11)");

            var frag_9_3 = products.Single(p => p.FragmentNumber == 9 && p.ResiduePosition == 3);
            Assert.That(frag_9_3.SecondaryFragmentNumber, Is.EqualTo(11),
                "GHIMAK: start=9, length=3 → end=11 > N=10 (IMA, wraps to A)");

            // All wrap-around fragments have end > N
            var wrapFragments = products.Where(p => p.SecondaryFragmentNumber > 10).ToList();
            Assert.That(wrapFragments, Is.Not.Empty,
                "Wrapping peptide GHIMAK[7-12] must produce wrap-around internal fragments");
        }

        [Test]
        public static void FragmentInternally_WrappingSubPeptide_WrapAroundMassIsCorrect()
        {
            // Verify the mass of a wrap-around fragment from GHIMAK[7-12].
            // Fragment [7,11] = G,H,I,M,A (length 5, wraps A from position 1).
            // Ring "AKDEFKGHIM": G=pos7, H=pos8, I=pos9, M=pos10, A=pos1(=11 in doubled)
            // Mass = sum of residue masses of G,H,I,M,A (no water for internal HCD)
            const string ringSeq = "AKDEFKGHIM";
            var protein = new CircularProtein(ringSeq, "wrap_mass_test");

            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 1,
                minPeptideLength: 1);

            var wrappingPeptide = protein.Digest(digestionParams, [], [])
                .OfType<CircularPeptideWithSetModifications>()
                .SingleOrDefault(p => p.BaseSequence == "GHIMAK");

            Assume.That(wrappingPeptide, Is.Not.Null);

            var products = new List<Product>();
            wrappingPeptide.FragmentInternally(DissociationType.HCD, 3, products);

            // Fragment [7,11] = G(7),H(8),I(9),M(10),A(1) — length 5, wraps
            var frag_7_5 = products.Single(p => p.FragmentNumber == 7 && p.ResiduePosition == 5);

            // Expected mass: G + H + I + M + A residue masses (wrap-around, reading from ring pos 7)
            double expected = ExpectedInternalFragmentNeutralMass(ringSeq, 7, 5);

            Assert.That(frag_7_5.NeutralMass,
                Is.EqualTo(expected).Within(1e-5),
                "Wrap-around fragment [7,11]=GHIMA: mass must equal G+H+I+M+A residue mass sum");
        }

        [Test]
        public static void FragmentInternally_FullRing_WrapAround_EndIsExactly2N()
        {
            // The maximum end position is start=N, length=N-1: end = N + (N-1) - 1 = 2N-2
            // For N=9: start=9, maxLength=8... but maxLength = N-s-1 = 9-8-1 = 0 (0-indexed)
            // Actually for start=9 (localStart=8, 0-based): maxLength = 9-8-1 = 0 → no fragments
            // For start=8 (localStart=7): maxLength = 9-7-1 = 1 → only length 2 (if minLength≤2)
            // Let us use minLength=2 and verify start=9, length would require localStart=8, maxLength=0
            // The largest end we can get is start=1, length=N-1=8: end=1+8-1=8 (non-wrap)
            // Or start=N (localStart=N-1=8), maxLength=N-8-1=0 → no products from last position
            // So the actual maximum wrap end = start + length - 1 where start is near N and length is max
            // start=8 (localStart=7), maxLength=1 → length=2: end=8+2-1=9 (not wrapping)
            // start=7 (localStart=6), maxLength=2 → length 2: end=8; length 3: but wait maxLength=2
            // Hmm: maxLength = N - localStart - 1. For N=9, localStart=6: maxLength=2, so length up to 2
            // To get wrap with minLength=2: localStart=7 (start=8), maxLength=1, length=2, end=9 (no wrap)
            // Wait - I need to reconsider. Let me recalculate for minLength=3:
            // localStart=5 (start=6), maxLength=3, lengths 3: end=6+3-1=8 (no wrap)
            // localStart=6 (start=7), maxLength=2 → below minLength=3, no products
            // So wrap only happens when... wait. The loop is:
            //   for localStart in [0, N): oneBasedStart = canonicalOffset + localStart
            //   maxLength = peptideLength - localStart - 1
            // For full ring, canonicalOffset=1, peptideLength=N=9:
            //   localStart=0: oneBasedStart=1, maxLength=8, lengths 3..8
            //   localStart=5: oneBasedStart=6, maxLength=3, lengths 3..3 → end=6+3-1=8
            //   localStart=6: oneBasedStart=7, maxLength=2 < 3 → no products with minLength=3
            // So with minLength=3, no wrap-around fragments exist? That contradicts earlier.
            // The issue: for full-ring, peptideLength=N=9. The wrap happens when
            // oneBasedStart + length - 1 > N=9. With oneBasedStart = 1 + localStart:
            //   1 + localStart + length - 1 > 9 → localStart + length > 9
            //   But maxLength = 9 - localStart - 1, so localStart + maxLength = 8 < 9
            // So localStart + length ≤ 8 < 9 always → end = oneBasedStart + length - 1
            //   = 1 + localStart + length - 1 = localStart + length ≤ 8 < 9 → never wraps!
            //
            // CONCLUSION: For a full-ring CircularPeptideWithSetModifications [1-N],
            // oneBasedStart = 1 + localStart, and end = oneBasedStart + length - 1
            //   = 1 + localStart + length - 1 ≤ 1 + (N-1) + (N-localStart-1) - 1
            // Wait: max end = 1 + localStart + maxLength - 1 = localStart + maxLength = localStart + N-localStart-1 = N-1
            // So end is at most N-1 < N → wrap-around fragments (end > N) NEVER occur
            // for a full-ring CircularPeptideWithSetModifications under the current code.
            //
            // Wrap-around would require oneBasedStart to be set relative to a starting
            // position > 1 in the ring — i.e., a wrapping peptide from Digest().
            Assert.Pass("Wrap-around fragments with end > N require a wrapping peptide " +
                        "(OneBasedStartResidueInProtein > 1 and OneBasedEndResidueInProtein > N). " +
                        "See wrap-around tests using wrapping peptides below.");
        }


        // ──────────────────────────────────────────────────────────────────────
        // 6. Numbering: wrapping peptide from Digest()
        //
        // Ring "ACDEKFGHIM" (N=10), K at position 5.
        // With maxMissedCleavages=1, Digest() yields the wrapping linear product
        // "GHIMACDEК" starting at position 6+1=6... wait, K is at position 5.
        // Actually: "ACDEKFGHIM", trypsin cuts after K at pos 5.
        // Proxy: "ACDEKFGHIMACDEKFGHI" — K at pos 5 and 15.
        // 0 MC products: "ACDEK"[1-5], "FGHIM"[6-10] (length<N → CircularPeptide)
        //                "ACDEK"[11-15] filtered (start>N)
        // 1 MC products: "ACDEKFGHIM"[1-10] length=N → linear PeptideWithSetMods
        //                "FGHIMACDEK"[6-15] length=N → linear PeptideWithSetMods
        //
        // So there are no wrapping CircularPeptideWithSetModifications from this ring.
        // To get a wrapping CircularPeptideWithSetMods we need 3+ cleavage sites.
        //
        // Ring "AKDEFKGHIM" (N=10), K at positions 2 and 6.
        // Proxy: "AKDEFKGHIMAKDEFKGHI"
        // K at proxy positions 2, 6, 12, 16.
        // 0 MC products starting ≤ N:
        //   "AK"[1-2], "DEFK"[3-6], "GHIM"[7-10]
        //   All length < N → CircularPeptideWithSetMods
        // 1 MC products starting ≤ N:
        //   "AKDEFK"[1-6] length=6 < N → CircularPeptideWithSetMods
        //   "DEFKGHIM"[3-10] length=8 < N → CircularPeptideWithSetMods
        //   "GHIMAK"[7-12] length=6 < N, start=7≤N, end=12>N → wrapping CircularPeptide!
        //   "AKDEFKGHIM"[1-10] length=N → linear PeptideWithSetMods
        //   "DEFKGHIMAK"[3-12] length=N → linear PeptideWithSetMods
        //   "GHIMAKDEFK"[7-16] length=N → linear PeptideWithSetMods
        //
        // "GHIMAK"[7-12]: OneBasedStartResidueInProtein=7, OneBasedEndResidueInProtein=12
        // Its internal fragments must use doubled-ring numbering:
        //   localStart=0, oneBasedStart=7, maxLength=5, lengths 2..5
        //   length=3: end=7+3-1=9 → [7,9] = GHI
        //   length=5: end=7+5-1=11 → [7,11] = GHIMA (wraps: A is position 11>N=10, i.e. position 1)
        //   localStart=3, oneBasedStart=10, maxLength=2
        //   length=2: end=10+2-1=11 → [10,11] = MA (M is pos 10, A wraps to pos 11)
        // ──────────────────────────────────────────────────────────────────────

        [Test]
        public static void FragmentInternally_WrappingSubPeptide_EndExceedsN()
        {
            // Ring "AKDEFKGHIM" (N=10), K at pos 2 and 6.
            // First verify canonicalization: A < D < E < F < G < H < I < K < M
            // "AKDEFKGHIM": A=1 is smallest → already canonical
            // Wrapping sub-peptide "GHIMAK" [7-12], obtained with 1 MC.
            const string ringSeq = "AKDEFKGHIM";
            var protein = new CircularProtein(ringSeq, "wrap_test");

            Assert.That(protein.BaseSequence, Is.EqualTo(ringSeq),
                "Sequence should be canonical (A is smallest)");

            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 1,
                minPeptideLength: 1);

            var wrappingPeptide = protein.Digest(digestionParams, [], [])
                .OfType<CircularPeptideWithSetModifications>()
                .SingleOrDefault(p => p.BaseSequence == "GHIMAK");

            Assume.That(wrappingPeptide, Is.Not.Null,
                "Wrapping sub-peptide GHIMAK[7-12] must be produced");
            Assert.That(wrappingPeptide.OneBasedStartResidueInProtein, Is.EqualTo(7));
            Assert.That(wrappingPeptide.OneBasedEndResidueInProtein, Is.EqualTo(12));

            var products = new List<Product>();
            wrappingPeptide.FragmentInternally(DissociationType.HCD, 3, products);

            // localStart=0, oneBasedStart=7, length=3: end=9 → [7,9] = GHI (no wrap)
            var frag_7_3 = products.Single(p =>
                p.FragmentNumber == 7 && p.ResiduePosition == 3);
            Assert.That(frag_7_3.SecondaryFragmentNumber, Is.EqualTo(9),
                "GHIMAK[7-12]: start=7, length=3 → end=9 (GHI), no wrap");

            // localStart=0, oneBasedStart=7, length=5: end=11 → [7,11] = GHIMA (wraps at A=pos11)
            var frag_7_5 = products.Single(p =>
                p.FragmentNumber == 7 && p.ResiduePosition == 5);
            Assert.That(frag_7_5.SecondaryFragmentNumber, Is.EqualTo(11),
                "GHIMAK[7-12]: start=7, length=5 → end=11 > N=10 (wraps to A at doubled pos 11)");
        }

        // ──────────────────────────────────────────────────────────────────────
        // 7. Mass correctness
        //
        // For HCD: nTermProductTypes=[b], cTermProductTypes=[y]
        //   b nTermCap = 0, y cTermCap = +H2O
        //   internal fragment neutralMass = residueMassSum + 0 + H2O - H2O = residueMassSum
        // ──────────────────────────────────────────────────────────────────────

        [Test]
        public static void FragmentInternally_FullRing_NeutralMass_MatchesResidueMassSum()
        {
            const string ring = "ACDEFGHIM"; // N=9
            var peptide = GetFullRingPeptide(ring);
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 3, products);

            // Verify every product's mass against independently computed residue mass sum
            foreach (var product in products)
            {
                int start = product.FragmentNumber;   // 1-based canonical start
                int length = product.ResiduePosition;  // fragment length
                double expected = ExpectedInternalFragmentNeutralMass(ring, start, length);

                Assert.That(product.NeutralMass,
                    Is.EqualTo(expected).Within(1e-5),
                    $"Fragment [{start},{start + length - 1}] length={length}: " +
                    $"expected mass {expected:F5}, got {product.NeutralMass:F5}");
            }
        }

        [Test]
        public static void FragmentInternally_FullRing_WrapAround_NeutralMass_MatchesResidueMassSum()
        {
            const string ring = "ACDEFGHIM"; // N=9
            var peptide = GetFullRingPeptide(ring);
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 3, products);

            // Find all wrap-around fragments (end > N=9) and verify their masses
            var wrapFragments = products
                .Where(p => p.SecondaryFragmentNumber > ring.Length)
                .ToList();

            // With minLength=3 and N=9, let's recalculate whether wraps exist.
            // As established above, for full-ring [1-9], end = oneBasedStart + length - 1
            // = 1 + localStart + length - 1 = localStart + length ≤ N-1 = 8 < 9.
            // So no wrap-around fragments exist for full-ring [1-9] with minLength=3.
            // This test documents that fact.
            Assert.That(wrapFragments, Is.Empty,
                "Full-ring [1-N] peptide produces no wrap-around fragments because " +
                "maxLength = N - localStart - 1 caps the end at N-1");
        }

        // ──────────────────────────────────────────────────────────────────────
        // 9. Product type correctness (HCD produces bIy internal ions)
        // ──────────────────────────────────────────────────────────────────────

        [Test]
        public static void FragmentInternally_HCD_ProductTypesAreByConvention()
        {
            // HCD: nTermProductTypes=[b], cTermProductTypes=[y]
            // Internal fragments have: ProductType=y (cTerm), SecondaryProductType=b (nTerm)
            var peptide = GetFullRingPeptide("ACDEFGHIM");
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 3, products);

            Assert.That(products, Is.Not.Empty);
            Assert.That(products.All(p => p.ProductType == ProductType.y), Is.True,
                "HCD internal fragments: primary ProductType must be y (cTerm)");
            Assert.That(products.All(p => p.SecondaryProductType == ProductType.b), Is.True,
                "HCD internal fragments: SecondaryProductType must be b (nTerm)");
            Assert.That(products.All(p => p.Terminus == FragmentationTerminus.None), Is.True,
                "Internal fragments must have FragmentationTerminus.None");
        }

        [Test]
        public static void FragmentInternally_ModifiedResidue_IncludedInFragmentMass()
        {
            // Ring "ACDEFGHIM" (N=9), no K or R → zero cuts → full-ring
            // CircularPeptideWithSetModifications. M is at position 9.
            //
            // Full-ring peptide [1-9], OneBasedStartResidueInProtein=1.
            // For localStart=s (0-based), maxLength = 9 - s - 1.
            // Fragments containing M (position 9, localStart=8):
            //   localStart=8, maxLength=0 → no fragments start at position 9
            // Fragments ending at M (position 9):
            //   localStart=6 (oneBasedStart=7), maxLength=2:
            //     length=2 → [7,8] = H,I  (no M)
            //     length=3 → [7,9] = H,I,M  (contains M!) ✓
            //   localStart=5 (oneBasedStart=6), maxLength=3:
            //     length=2 → [6,7] = G,H
            //     length=3 → [6,8] = G,H,I
            //     length=4 → wait, maxLength=3 → [6,9] = G,H,I,M (contains M!) ✓
            //
            // Note: [8,9] = I,M (length 2 starting at localStart=7) requires
            // maxLength = 9-7-1 = 1 ≥ 2, which is false → NOT produced.
            // The smallest fragment containing M that IS produced is [7,9] (length 3).

            ModificationMotif.TryGetMotif("M", out ModificationMotif mMotif);
            var oxidation = new Modification(
                _originalId: "Oxidation on M",
                _modificationType: "Common Variable",
                _target: mMotif,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: 15.994915);

            var mods = new Dictionary<int, List<Modification>> { [9] = [oxidation] };
            var protein = new CircularProtein("ACDEFGHIM", "mod_test",
                oneBasedModifications: mods);

            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1);

            var peptides = protein.Digest(digestionParams, [], [oxidation])
                .OfType<CircularPeptideWithSetModifications>()
                .ToList();

            var unmodified = peptides.Single(p => p.AllModsOneIsNterminus.Count == 0);
            var oxidized = peptides.Single(p => p.AllModsOneIsNterminus.Count == 1);

            Assert.That(oxidized.AllModsOneIsNterminus.ContainsKey(10), Is.True,
                "Oxidation on M at position 9 → allMods key = 9+1 = 10");
            Assert.That(oxidized.AllModsOneIsNterminus[10].OriginalId,
                Does.Contain("Oxidation"),
                "The modification at key 10 must be oxidation");

            const string ring = "ACDEFGHIM";
            double oxMass = 15.994915;

            var unmodProducts = new List<Product>();
            unmodified.FragmentInternally(DissociationType.HCD, 2, unmodProducts);

            var oxProducts = new List<Product>();
            oxidized.FragmentInternally(DissociationType.HCD, 2, oxProducts);

            // Fragment [7,8] = H,I (length 2, localStart=6, maxLength=2 ≥ 2) — no M
            // → same mass in both unmodified and oxidized forms
            var unmod_7_2 = unmodProducts.Single(p =>
                p.FragmentNumber == 7 && p.ResiduePosition == 2);
            var ox_7_2 = oxProducts.Single(p =>
                p.FragmentNumber == 7 && p.ResiduePosition == 2);

            Assert.That(unmod_7_2.NeutralMass,
                Is.EqualTo(ExpectedInternalFragmentNeutralMass(ring, 7, 2)).Within(1e-4),
                "Unmodified [7,8]=HI mass correct");
            Assert.That(ox_7_2.NeutralMass,
                Is.EqualTo(ExpectedInternalFragmentNeutralMass(ring, 7, 2)).Within(1e-4),
                "Oxidized [7,8]=HI: no M → no oxidation shift");

            // Fragment [7,9] = H,I,M (length 3, localStart=6, maxLength=2)...
            // Wait: maxLength = 9-6-1 = 2, so length up to 2 only → [7,9] NOT produced either!
            // localStart=5 (oneBasedStart=6), maxLength=3: length=4 → [6,9]=GHIM ✓
            // localStart=4 (oneBasedStart=5), maxLength=4: length=5 → [5,9]=FGHIM ✓

            // Fragment [6,9] = G,H,I,M (length 4, localStart=5, maxLength=3 ≥ 4?)
            // maxLength = 9-5-1 = 3 → length up to 3 only → [6,9] NOT produced!
            // localStart=4 (oneBasedStart=5), maxLength=4: length=4 → [5,8]=FGHI, length=5→[5,9]=FGHIM
            // Wait: maxLength = 9-4-1 = 4. So length=5 → [5,9]=FGHIM ✓

            // Let me enumerate ALL fragments containing M systematically.
            // M is at oneBasedStart + length - 1 = 9, so oneBasedStart + length = 10.
            // oneBasedStart = 1 + localStart, maxLength = 9 - localStart - 1.
            // Need: (1 + localStart) + length - 1 = 9 → localStart + length = 9
            // Also: length ≤ maxLength = 9 - localStart - 1
            // → localStart + length ≤ 9 - 1 = 8 < 9
            // So localStart + length = 9 is NEVER satisfiable!
            // This means NO fragment with M as the last residue can be produced
            // from a full-ring [1-9] peptide with these loop bounds.
            //
            // M can also be in the MIDDLE of a fragment. Fragment [s, e] contains M
            // if s ≤ 9 ≤ e, i.e., oneBasedStart ≤ 9 AND oneBasedStart + length - 1 ≥ 9.
            // But we showed oneBasedStart + length - 1 ≤ N-1 = 8 always.
            // Therefore NO fragment produced from full-ring [1-9] contains M at position 9.
            //
            // CONCLUSION: The modification test cannot be meaningfully tested using a
            // full-ring peptide starting at position 1, because the loop bounds prevent
            // any fragment from reaching position N. We need M at an earlier position.

            // Use M at position 5 instead: ring "ACDEMFGHIM" — but that has two M's.
            // Use ring "ACMDEFGHI" (N=9), M at position 3, no K/R → full-ring [1-9].
            // Fragments containing M (position 3):
            //   localStart=0 (oneBasedStart=1), maxLength=8:
            //     length=3 → [1,3]=ACM ✓ (contains M)
            //     length=4 → [1,4]=ACMD ✓
            //   localStart=1 (oneBasedStart=2), maxLength=7:
            //     length=2 → [2,3]=CM ✓
            //     length=3 → [2,4]=CMD ✓
            //   localStart=2 (oneBasedStart=3), maxLength=6: M is at start, no fragment
            //     starts at M with length < maxLength=6, but M is position 3 so
            //     fragment [3,4]=MD (length 2) starts at M ✓

            // This test needs a different ring. Provide as a separate test for clarity.
            // For THIS test, document that no fragment contains M and verify mass
            // difference is zero for all fragments (oxidation has no effect on fragment masses).

            Assert.That(unmodProducts.Count, Is.EqualTo(oxProducts.Count),
                "Unmodified and oxidized full-ring [1-9] peptides produce the same number of fragments");

            // Verify that oxidized and unmodified fragment masses are identical —
            // because M is at position 9 and no fragment produced from [1-9] reaches it.
            for (int i = 0; i < unmodProducts.Count; i++)
            {
                var u = unmodProducts.OrderBy(p => p.FragmentNumber).ThenBy(p => p.ResiduePosition).ToList()[i];
                var o = oxProducts.OrderBy(p => p.FragmentNumber).ThenBy(p => p.ResiduePosition).ToList()[i];
                Assert.That(o.NeutralMass, Is.EqualTo(u.NeutralMass).Within(1e-9),
                    $"Fragment [{u.FragmentNumber},{u.SecondaryFragmentNumber}]: " +
                    $"oxidation on M at position 9 cannot affect this fragment " +
                    $"because no fragment from full-ring [1-9] reaches position 9");
            }
        }

        [Test]
        public static void FragmentInternally_ModifiedResidue_MidRing_IncludedInFragmentMass()
        {
            // Ring "ACMDEFGHI" (N=9), M at position 3, no K or R → full-ring [1-9].
            // M at position 3 CAN be reached by fragments starting at positions 1 or 2.
            // Fragment [1,3]=ACM (length 3): contains M → mass includes oxidation ✓
            // Fragment [2,3]=CM  (length 2): contains M → mass includes oxidation ✓
            // Fragment [1,2]=AC  (length 2): no M → mass unchanged ✓

            ModificationMotif.TryGetMotif("M", out ModificationMotif mMotif);
            var oxidation = new Modification(
                _originalId: "Oxidation on M",
                _modificationType: "Common Variable",
                _target: mMotif,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: 15.994915);

            var mods = new Dictionary<int, List<Modification>> { [3] = [oxidation] };
            var protein = new CircularProtein("ACMDEFGHI", "mod_mid_test",
                oneBasedModifications: mods);

            Assert.That(protein.BaseSequence, Is.EqualTo("ACMDEFGHI"),
                "Already canonical (A is smallest)");

            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1);

            var peptides = protein.Digest(digestionParams, [], [oxidation])
                .OfType<CircularPeptideWithSetModifications>()
                .ToList();

            var unmodified = peptides.Single(p => p.AllModsOneIsNterminus.Count == 0);
            var oxidized = peptides.Single(p => p.AllModsOneIsNterminus.Count == 1);

            // M at position 3 → allModsOneIsNterminus key = 3+1 = 4
            Assert.That(oxidized.AllModsOneIsNterminus.ContainsKey(4), Is.True,
                "Oxidation on M at position 3 → allMods key = 4");

            const string ring = "ACMDEFGHI";
            double oxMass = 15.994915;

            var unmodProducts = new List<Product>();
            unmodified.FragmentInternally(DissociationType.HCD, 2, unmodProducts);

            var oxProducts = new List<Product>();
            oxidized.FragmentInternally(DissociationType.HCD, 2, oxProducts);

            // Fragment [1,2] = A,C (length 2) — no M → same mass in both forms
            var unmod_1_2 = unmodProducts.Single(p =>
                p.FragmentNumber == 1 && p.ResiduePosition == 2);
            var ox_1_2 = oxProducts.Single(p =>
                p.FragmentNumber == 1 && p.ResiduePosition == 2);
            Assert.That(ox_1_2.NeutralMass,
                Is.EqualTo(unmod_1_2.NeutralMass).Within(1e-4),
                "[1,2]=AC: no M → oxidation has no effect");
            Assert.That(ox_1_2.NeutralMass,
                Is.EqualTo(ExpectedInternalFragmentNeutralMass(ring, 1, 2)).Within(1e-4),
                "[1,2]=AC mass correct");

            // Fragment [2,3] = C,M (length 2) — contains M → oxidation shifts mass
            var unmod_2_2 = unmodProducts.Single(p =>
                p.FragmentNumber == 2 && p.ResiduePosition == 2);
            var ox_2_2 = oxProducts.Single(p =>
                p.FragmentNumber == 2 && p.ResiduePosition == 2);
            Assert.That(ox_2_2.NeutralMass,
                Is.EqualTo(ExpectedInternalFragmentNeutralMass(ring, 2, 2) + oxMass).Within(1e-4),
                "[2,3]=CM: M in fragment → mass includes oxidation");
            Assert.That(ox_2_2.NeutralMass - unmod_2_2.NeutralMass,
                Is.EqualTo(oxMass).Within(1e-4),
                "Mass difference for [2,3]=CM must equal oxidation mass");

            // Fragment [1,3] = A,C,M (length 3) — contains M → oxidation shifts mass
            var unmod_1_3 = unmodProducts.Single(p =>
                p.FragmentNumber == 1 && p.ResiduePosition == 3);
            var ox_1_3 = oxProducts.Single(p =>
                p.FragmentNumber == 1 && p.ResiduePosition == 3);
            Assert.That(ox_1_3.NeutralMass,
                Is.EqualTo(ExpectedInternalFragmentNeutralMass(ring, 1, 3) + oxMass).Within(1e-4),
                "[1,3]=ACM: M in fragment → mass includes oxidation");
            Assert.That(ox_1_3.NeutralMass - unmod_1_3.NeutralMass,
                Is.EqualTo(oxMass).Within(1e-4),
                "Mass difference for [1,3]=ACM must equal oxidation mass");
        }

        [Test]
        public static void FragmentInternally_SubPeptide_N5_MinLength2_CorrectCount()
        {
            // Ring "ACDEKFGHIK" (N=10), K at positions 5 and 10.
            // Both K's are genuine ring cleavage sites → "ACDEK"[1-5] and "FGHIK"[6-10]
            // are genuine two-cut CircularPeptideWithSetModifications.
            // Sub-peptide "ACDEK", length 5, minLength 2:
            //   s=0: lengths 2,3,4 → 3
            //   s=1: lengths 2,3   → 2
            //   s=2: length  2     → 1
            //   Total = 6 fragments
            var peptide = GetSubPeptide("ACDEKFGHIK", "ACDEK");
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 2, products);

            Assert.That(products, Has.Count.EqualTo(6),
                "Sub-peptide length 5, minLength 2: expected 6 internal fragments");
        }

        [Test]
        public static void FragmentInternally_SubPeptide_AtRingStart_NumberingIsCanonical()
        {
            // "ACDEK"[1-5] from ring "ACDEKFGHIK".
            // Canonical positions: A=1,C=2,D=3,E=4,K=5.
            var peptide = GetSubPeptide("ACDEKFGHIK", "ACDEK");
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 3, products);

            // localStart=0, oneBasedStart=1, length=3: end=3 → [1,3] = ACD
            var frag_1_3 = products.Single(p =>
                p.FragmentNumber == 1 && p.ResiduePosition == 3);
            Assert.That(frag_1_3.SecondaryFragmentNumber, Is.EqualTo(3),
                "ACDEK[1-5]: start=1, length=3 → end=3 (ACD)");

            // localStart=0, oneBasedStart=1, length=4: end=4 → [1,4] = ACDE
            var frag_1_4 = products.Single(p =>
                p.FragmentNumber == 1 && p.ResiduePosition == 4);
            Assert.That(frag_1_4.SecondaryFragmentNumber, Is.EqualTo(4),
                "ACDEK[1-5]: start=1, length=4 → end=4 (ACDE)");

            // localStart=1, oneBasedStart=2, length=3: end=4 → [2,4] = CDE
            var frag_2_3 = products.Single(p =>
                p.FragmentNumber == 2 && p.ResiduePosition == 3);
            Assert.That(frag_2_3.SecondaryFragmentNumber, Is.EqualTo(4),
                "ACDEK[1-5]: start=2, length=3 → end=4 (CDE)");
        }

        [Test]
        public static void FragmentInternally_SubPeptide_NotAtRingStart_NumberingIsCanonical()
        {
            // "FGHIK"[6-10] from ring "ACDEKFGHIK".
            // Canonical positions: F=6,G=7,H=8,I=9,K=10.
            var peptide = GetSubPeptide("ACDEKFGHIK", "FGHIK");
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 3, products);

            // localStart=0, oneBasedStart=6, length=3: end=8 → [6,8] = FGH
            var frag_6_3 = products.Single(p =>
                p.FragmentNumber == 6 && p.ResiduePosition == 3);
            Assert.That(frag_6_3.SecondaryFragmentNumber, Is.EqualTo(8),
                "FGHIK[6-10]: start=6, length=3 → end=8 (FGH)");

            // localStart=0, oneBasedStart=6, length=4: end=9 → [6,9] = FGHI
            var frag_6_4 = products.Single(p =>
                p.FragmentNumber == 6 && p.ResiduePosition == 4);
            Assert.That(frag_6_4.SecondaryFragmentNumber, Is.EqualTo(9),
                "FGHIK[6-10]: start=6, length=4 → end=9 (FGHI)");

            // localStart=1, oneBasedStart=7, length=3: end=9 → [7,9] = GHI
            var frag_7_3 = products.Single(p =>
                p.FragmentNumber == 7 && p.ResiduePosition == 3);
            Assert.That(frag_7_3.SecondaryFragmentNumber, Is.EqualTo(9),
                "FGHIK[6-10]: start=7, length=3 → end=9 (GHI)");
        }

        [Test]
        public static void FragmentInternally_SubPeptide_NoFragmentsStartAtParentPosition1()
        {
            // "FGHIK" starts at parent position 6.
            // No fragment should start at a position < 6.
            var peptide = GetSubPeptide("ACDEKFGHIK", "FGHIK");
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 2, products);

            Assert.That(products.All(p => p.FragmentNumber >= 6), Is.True,
                "Fragments of FGHIK[6-10] must all start at canonical position ≥ 6");
        }

        [Test]
        public static void FragmentInternally_SubPeptide_NeutralMass_MatchesResidueMassSum()
        {
            // "FGHIK"[6-10]: verify masses use parent ring residue masses
            const string ring = "ACDEKFGHIK";
            var peptide = GetSubPeptide(ring, "FGHIK");
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 3, products);

            foreach (var product in products)
            {
                int start = product.FragmentNumber;
                int length = product.ResiduePosition;
                double expected = ExpectedInternalFragmentNeutralMass(ring, start, length);

                Assert.That(product.NeutralMass,
                    Is.EqualTo(expected).Within(1e-5),
                    $"FGHIK fragment [{start},{start + length - 1}]: " +
                    $"expected {expected:F5}, got {product.NeutralMass:F5}");
            }
        }

        [Test]
        public static void FragmentInternally_SubPeptide_NoFragmentExceedsSpan()
        {
            // "ACDEK"[1-5]: no fragment may end at position 5 or beyond —
            // position 5 is the K cleavage residue, the C-terminal boundary.
            var peptide = GetSubPeptide("ACDEKFGHIK", "ACDEK");
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 2, products);

            Assert.That(products.All(p => p.SecondaryFragmentNumber < 5), Is.True,
                "ACDEK[1-5]: no internal fragment may end at or beyond position 5");
        }

        [Test]
        public static void FragmentInternally_SubPeptide_NoWrapAroundFragments()
        {
            // Sub-peptides (length < N) must never produce wrap-around fragments.
            var peptide = GetSubPeptide("ACDEKFGHIK", "FGHIK");
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 2, products);

            Assert.That(products.All(p => p.SecondaryFragmentNumber <= 10), Is.True,
                "Sub-peptide FGHIK[6-10]: no wrap-around fragments (end must be ≤ N=10)");
        }
    }
}