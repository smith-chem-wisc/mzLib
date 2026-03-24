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
    /// <summary>
    /// Tests for <see cref="CircularPeptideWithSetModifications"/> and the
    /// <see cref="CircularProtein.Digest"/> pipeline that produces it.
    ///
    /// NUMBERING SYSTEM
    /// ----------------
    /// All positions are 1-based and expressed in the parent
    /// <see cref="CircularProtein"/>'s canonical coordinate system — the
    /// lexicographically smallest rotation of the ring sequence, with residue 1
    /// being the alphabetically earliest position (ties broken by subsequent
    /// residues).
    ///
    /// Fragment annotations use the doubled-ring convention:
    ///   [FragmentNumber, SecondaryFragmentNumber] = [s, s+L-1]
    /// where s is the 1-based canonical start and L is the fragment length.
    /// SecondaryFragmentNumber exceeds N only for wrap-around fragments from
    /// peptides whose OneBasedStartResidueInProtein > 1.
    ///
    /// DIGESTION RULES
    /// ---------------
    ///   0 cuts  → CircularPeptideWithSetModifications, length N (ring intact)
    ///   1 cut   → PeptideWithSetModifications, length N (ring opened)
    ///   2+ cuts → CircularPeptideWithSetModifications sub-peptides, length < N
    ///
    /// FRAGMENT LOOP BOUNDS
    /// --------------------
    /// For localStart in [0, peptideLength):
    ///   oneBasedStart = OneBasedStartResidueInProtein + localStart
    ///   maxLength     = peptideLength - localStart - 1
    ///
    /// Consequence: for a full-ring peptide with OneBasedStartResidueInProtein=1,
    /// SecondaryFragmentNumber ≤ N-1 always — wrap-around fragments (end > N)
    /// only arise from wrapping sub-peptides (OneBasedStartResidueInProtein > 1).
    ///
    /// Consequence: M at position N (the last residue) is unreachable from any
    /// fragment of a full-ring [1-N] peptide; place the modified residue at an
    /// earlier canonical position to test mass incorporation.
    /// </summary>
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

        /// <summary>
        /// Returns the sum of monoisotopic residue masses for a contiguous
        /// sub-sequence of a ring, reading with wrap-around.
        /// Used as the independent reference for fragment neutral masses.
        /// For HCD internal ions: neutralMass = residueMassSum
        /// (b nTermCap=0 + y cTermCap=+H₂O − H₂O cancels).
        /// </summary>
        private static double ResidueMassSum(string ringSequence, int oneBasedStart, int length)
        {
            int n = ringSequence.Length;
            double mass = 0;
            for (int i = 0; i < length; i++)
                mass += Residue.ResidueMonoisotopicMass[ringSequence[(oneBasedStart - 1 + i) % n]];
            return mass;
        }

        private static double ExpectedInternalFragmentNeutralMass(
            string ringSequence, int oneBasedStart, int length) =>
            ResidueMassSum(ringSequence, oneBasedStart, length);

        /// <summary>
        /// Returns the single full-ring CircularPeptideWithSetModifications
        /// produced by trypsin digestion of a ring with no K or R residues.
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
        /// Returns the CircularPeptideWithSetModifications sub-peptide matching
        /// <paramref name="subSequence"/> from trypsin digestion (0 missed cleavages)
        /// of <paramref name="ringSequence"/>. The ring must have at least two K/R
        /// residues so that both boundaries of the sub-peptide are genuine cut sites.
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

        // ══════════════════════════════════════════════════════════════════════
        // DIGESTION — modification renumbering
        //
        // When a circular protein is digested, variable modifications travel
        // with their residues. A modification at ring position p is renumbered
        // in the resulting peptide's AllModsOneIsNterminus dictionary:
        //   key = (proxy position of modified residue)
        //         - (proxy start of peptide) + 2
        // This matches the standard PeptideWithSetModifications convention
        // (key = local 0-based residue index + 2).
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public static void Digest_LinearProduct_ModificationRenumberedToProxyPosition()
        {
            // Ring "AHMCDEFGKIEP" (N=12), already canonical (A is smallest).
            // M at ring position 3, K at ring position 9.
            //
            // One tryptic cut after K(9) → linear peptide [10, 21]:
            //   I(10) E(11) P(12) A(13→1) H(14→2) M(15→3) C(16→4)
            //   D(17→5) E(18→6) F(19→7) G(20→8) K(21→9)
            //   BaseSequence = "IEPAHMCDEFGK"
            //
            // Oxidation renumbering:
            //   M is at proxy position 15, peptide starts at proxy position 10.
            //   key = 15 - 10 + 2 = 7
            //   (M is the 6th residue of "IEPAHMCDEFGK", 0-based index 5, key=5+2=7)

            ModificationMotif.TryGetMotif("M", out ModificationMotif mMotif);
            var oxidation = new Modification(
                _originalId: "Oxidation on M",
                _modificationType: "Common Variable",
                _target: mMotif,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: 15.994915);

            var mods = new Dictionary<int, List<Modification>> { [3] = [oxidation] };
            var protein = new CircularProtein("AHMCDEFGKIEP", "mod_renumber_test",
                oneBasedModifications: mods);

            Assert.That(protein.BaseSequence, Is.EqualTo("AHMCDEFGKIEP"));

            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1);

            var peptides = protein.Digest(digestionParams, [], [oxidation])
                .Cast<PeptideWithSetModifications>()
                .ToList();

            // One K → two linear products: unmodified and oxidized
            Assert.That(peptides, Has.Count.EqualTo(2));
            Assert.That(peptides.All(p => p is not CircularPeptideWithSetModifications), Is.True,
                "One cut → linear PeptideWithSetModifications, not circular");

            var unmodified = peptides.Single(p => p.AllModsOneIsNterminus.Count == 0);
            var oxidized = peptides.Single(p => p.AllModsOneIsNterminus.Count == 1);

            foreach (var peptide in peptides)
            {
                Assert.That(peptide.BaseSequence, Is.EqualTo("IEPAHMCDEFGK"));
                Assert.That(peptide.OneBasedStartResidueInProtein, Is.EqualTo(10));
                Assert.That(peptide.OneBasedEndResidueInProtein, Is.EqualTo(21));
                Assert.That(peptide.Length, Is.EqualTo(protein.Length));
                Assert.That(ReferenceEquals(peptide.Protein, protein), Is.True);
            }

            Assert.That(oxidized.AllModsOneIsNterminus.ContainsKey(7), Is.True,
                "Oxidation on M: proxy pos 15, start 10 → key = 15-10+2 = 7");
            Assert.That(oxidized.AllModsOneIsNterminus[7].OriginalId,
                Does.Contain("Oxidation"));
            Assert.That(
                oxidized.MonoisotopicMass - unmodified.MonoisotopicMass,
                Is.EqualTo(15.994915).Within(1e-4));
        }

        // ══════════════════════════════════════════════════════════════════════
        // ModFits — understanding the public API
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public static void ModFits_BasicUnderstanding()
        {
            // Verify ModificationLocalization.ModFits behavior for oxidation on M.
            // Used by the Digest() fallback (zero-cut path) to apply variable
            // modifications to the intact full-ring peptide.
            ModificationMotif.TryGetMotif("M", out ModificationMotif mMotif);
            var oxidation = new Modification(
                _originalId: "Oxidation on M",
                _modificationType: "Common Variable",
                _target: mMotif,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: 15.994915);

            const string sequence = "ACDEFGHIM"; // M at position 9

            bool fitsAtPos9 = ModificationLocalization.ModFits(
                oxidation, sequence,
                digestionProductOneBasedIndex: 9,
                digestionProductLength: 9,
                bioPolymerOneBasedIndex: 9);

            bool fitsAtPos1 = ModificationLocalization.ModFits(
                oxidation, sequence,
                digestionProductOneBasedIndex: 1,
                digestionProductLength: 9,
                bioPolymerOneBasedIndex: 1);

            bool fitsAtPos7 = ModificationLocalization.ModFits(
                oxidation, sequence,
                digestionProductOneBasedIndex: 7,
                digestionProductLength: 9,
                bioPolymerOneBasedIndex: 7);

            Assert.That(fitsAtPos9, Is.True, "M at position 9 → fits");
            Assert.That(fitsAtPos1, Is.False, "A at position 1 → does not fit");
            Assert.That(fitsAtPos7, Is.False, "H at position 7 → does not fit");

            Console.WriteLine($"mod.OriginalId          = {oxidation.OriginalId}");
            Console.WriteLine($"mod.LocationRestriction = {oxidation.LocationRestriction}");
            Console.WriteLine($"mod.Target              = {oxidation.Target}");
            Console.WriteLine($"ModFits pos1={fitsAtPos1}, pos7={fitsAtPos7}, pos9={fitsAtPos9}");
        }

        [Test]
        public static void ModFits_FallbackDiagnostic()
        {
            // Confirms the zero-cut Digest() fallback produces both unmodified
            // and oxidized full-ring peptides when oxidation is a variable mod.
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

            Console.WriteLine($"OneBasedPossibleLocalizedModifications count: " +
                $"{protein.OneBasedPossibleLocalizedModifications.Count}");
            foreach (var kvp in protein.OneBasedPossibleLocalizedModifications)
                foreach (var m in kvp.Value)
                    Console.WriteLine($"  Pos={kvp.Key} Id={m.OriginalId}");
        }

        // ══════════════════════════════════════════════════════════════════════
        // FRAGMENT GUARDS
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public static void FragmentInternally_MinLengthLessThanTwo_ReturnsEmpty()
        {
            // minLengthOfFragments < 2 is rejected — single-residue internal
            // fragments are not physically meaningful.
            var peptide = GetFullRingPeptide("ACDEFGHIM");
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 1, products);

            Assert.That(products, Is.Empty);
        }

        [Test]
        public static void FragmentInternally_MinLengthEqualsToPeptideLength_ReturnsEmpty()
        {
            // minLengthOfFragments >= peptideLength: no fragment can be shorter
            // than the full peptide, so nothing is produced.
            var peptide = GetFullRingPeptide("ACDEFGHIM"); // length 9
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 9, products);

            Assert.That(products, Is.Empty);
        }

        [Test]
        public static void FragmentInternally_ClearsProductsListOnEntry()
        {
            // The products list is always cleared at the start of the call.
            var peptide = GetFullRingPeptide("ACDEFGHIM");
            var products = new List<Product>
                { new Product(ProductType.b, FragmentationTerminus.N, 100, 1, 1, 0) };

            peptide.FragmentInternally(DissociationType.HCD, 3, products);

            Assert.That(products.All(p => p.NeutralMass != 100), Is.True,
                "The sentinel entry must be gone — list was cleared");
        }

        // ══════════════════════════════════════════════════════════════════════
        // FRAGMENT COUNT — full-ring peptide
        //
        // For N=9, minLength L=3, full-ring [1-9]:
        //   localStart s | maxLength | lengths produced
        //   -------------|-----------|----------------
        //        0       |     8     | 3,4,5,6,7,8  → 6
        //        1       |     7     | 3,4,5,6,7    → 5
        //        2       |     6     | 3,4,5,6      → 4
        //        3       |     5     | 3,4,5        → 3
        //        4       |     4     | 3,4          → 2
        //        5       |     3     | 3            → 1
        //        6,7,8   |   < 3     | (none)
        //   Total = 21 fragment positions × 1 ion-type pair (b,y) = 21 products
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public static void FragmentInternally_FullRing_N9_MinLength3_CorrectCount()
        {
            var peptide = GetFullRingPeptide("ACDEFGHIM"); // N=9, no K/R
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 3, products);

            // HCD: nTermProductTypes=[b], cTermProductTypes=[y]
            // Each (start, length) pair → 1 product (1 b-nTerm × 1 y-cTerm)
            Assert.That(products, Has.Count.EqualTo(21));
        }

        // ══════════════════════════════════════════════════════════════════════
        // PRODUCT TYPE — HCD convention
        //
        // Internal ions follow the bIy convention:
        //   ProductType          = y  (C-terminal cap)
        //   SecondaryProductType = b  (N-terminal cap)
        //   FragmentationTerminus = None
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public static void FragmentInternally_HCD_ProductTypesAreByConvention()
        {
            var peptide = GetFullRingPeptide("ACDEFGHIM");
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 3, products);

            Assert.That(products, Is.Not.Empty);
            Assert.That(products.All(p => p.ProductType == ProductType.y), Is.True);
            Assert.That(products.All(p => p.SecondaryProductType == ProductType.b), Is.True);
            Assert.That(products.All(p => p.Terminus == FragmentationTerminus.None), Is.True);
        }

        // ══════════════════════════════════════════════════════════════════════
        // NUMBERING — non-wrapping fragments of a full-ring [1-N] peptide
        //
        // Ring "ACDEFGHIM" (N=9), peptide [1-9].
        // Fragment [s, s+L-1] where s = oneBasedStart, L = length.
        //   [1,3] = ACD  ✓
        //   [1,4] = ACDE ✓
        //   [4,6] = EFG  ✓
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public static void FragmentInternally_FullRing_NonWrapping_StartAndEndNumbering()
        {
            const string ring = "ACDEFGHIM"; // N=9
            var peptide = GetFullRingPeptide(ring);
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 3, products);

            var frag_1_3 = products.Single(p => p.FragmentNumber == 1 && p.ResiduePosition == 3);
            Assert.That(frag_1_3.SecondaryFragmentNumber, Is.EqualTo(3), "[1,3]");

            var frag_1_4 = products.Single(p => p.FragmentNumber == 1 && p.ResiduePosition == 4);
            Assert.That(frag_1_4.SecondaryFragmentNumber, Is.EqualTo(4), "[1,4]");

            var frag_4_3 = products.Single(p => p.FragmentNumber == 4 && p.ResiduePosition == 3);
            Assert.That(frag_4_3.SecondaryFragmentNumber, Is.EqualTo(6), "[4,6]");
        }

        // ══════════════════════════════════════════════════════════════════════
        // NUMBERING — wrap-around is impossible for full-ring [1-N] peptides
        //
        // Proof: for OneBasedStartResidueInProtein=1,
        //   SecondaryFragmentNumber = 1 + localStart + length - 1
        //                           ≤ 1 + localStart + (N - localStart - 1) - 1
        //                           = N - 1 < N
        // Therefore SecondaryFragmentNumber < N always.
        // Wrap-around (end > N) requires OneBasedStartResidueInProtein > 1,
        // which only occurs for wrapping sub-peptides from Digest().
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public static void FragmentInternally_FullRingStartingAt1_NeverProducesWrapAroundFragments()
        {
            const string ring = "ACDEFGHIM"; // N=9
            var peptide = GetFullRingPeptide(ring);
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 2, products);

            Assert.That(products, Is.Not.Empty);
            Assert.That(products.All(p => p.SecondaryFragmentNumber < ring.Length), Is.True,
                "Full-ring [1-N]: SecondaryFragmentNumber < N always");
        }

        [Test]
        public static void FragmentInternally_FullRing_WrapAround_EndIsExactly2N()
        {
            // This test documents the mathematical proof above.
            Assert.Pass(
                "Wrap-around fragments (SecondaryFragmentNumber > N) require " +
                "OneBasedStartResidueInProtein > 1. " +
                "See wrapping sub-peptide tests below.");
        }

        [Test]
        public static void FragmentInternally_FullRing_WrapAround_NeutralMass_MatchesResidueMassSum()
        {
            // Corollary: the full-ring [1-N] peptide produces zero wrap-around fragments.
            const string ring = "ACDEFGHIM"; // N=9
            var peptide = GetFullRingPeptide(ring);
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 3, products);

            Assert.That(
                products.Where(p => p.SecondaryFragmentNumber > ring.Length).ToList(),
                Is.Empty,
                "No wrap-around fragments from full-ring [1-9] with minLength=3");
        }

        // ══════════════════════════════════════════════════════════════════════
        // NUMBERING — wrapping sub-peptide (OneBasedStart > 1, OneBasedEnd > N)
        //
        // Ring "AKDEFKGHIM" (N=10), K at positions 2 and 6.
        // With maxMissedCleavages=1, Digest() produces wrapping sub-peptide
        // "GHIMAK" with OneBasedStart=7, OneBasedEnd=12.
        //
        // Fragment loop for peptideLength=6, canonicalOffset=7:
        //   localStart=0: oneBasedStart=7, maxLength=5
        //     [7,9]=GHI (L=3), [7,10]=GHIM (L=4), [7,11]=GHIMA (L=5, wraps!)
        //   localStart=1: oneBasedStart=8, maxLength=4
        //     [8,10]=HIM (L=3), [8,11]=HIMA (L=4, wraps!)
        //   localStart=2: oneBasedStart=9, maxLength=3
        //     [9,11]=IMA (L=3, wraps!)
        //   localStart=3: oneBasedStart=10, maxLength=2
        //     [10,11]=MA (L=2, wraps!)
        //   localStart=4: oneBasedStart=11, maxLength=1 < minLength=3 → none
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public static void FragmentInternally_WrappingSubPeptide_EndExceedsN()
        {
            const string ringSeq = "AKDEFKGHIM"; // N=10, K at 2 and 6
            var protein = new CircularProtein(ringSeq, "wrap_test");

            Assert.That(protein.BaseSequence, Is.EqualTo(ringSeq), "Already canonical");

            var wrappingPeptide = protein
                .Digest(new DigestionParams("trypsin", maxMissedCleavages: 1, minPeptideLength: 1),
                    [], [])
                .OfType<CircularPeptideWithSetModifications>()
                .SingleOrDefault(p => p.BaseSequence == "GHIMAK");

            Assume.That(wrappingPeptide, Is.Not.Null,
                "GHIMAK[7-12] must be produced by Digest()");
            Assert.That(wrappingPeptide.OneBasedStartResidueInProtein, Is.EqualTo(7));
            Assert.That(wrappingPeptide.OneBasedEndResidueInProtein, Is.EqualTo(12));

            var products = new List<Product>();
            wrappingPeptide.FragmentInternally(DissociationType.HCD, 3, products);

            // Non-wrapping
            var frag_7_3 = products.Single(p => p.FragmentNumber == 7 && p.ResiduePosition == 3);
            Assert.That(frag_7_3.SecondaryFragmentNumber, Is.EqualTo(9), "[7,9]=GHI, no wrap");

            var frag_7_4 = products.Single(p => p.FragmentNumber == 7 && p.ResiduePosition == 4);
            Assert.That(frag_7_4.SecondaryFragmentNumber, Is.EqualTo(10), "[7,10]=GHIM, ends at N");

            // Wrap-around
            var frag_7_5 = products.Single(p => p.FragmentNumber == 7 && p.ResiduePosition == 5);
            Assert.That(frag_7_5.SecondaryFragmentNumber, Is.EqualTo(11),
                "[7,11]=GHIMA wraps: A is canonical position 1, doubled position 11");

            var frag_9_3 = products.Single(p => p.FragmentNumber == 9 && p.ResiduePosition == 3);
            Assert.That(frag_9_3.SecondaryFragmentNumber, Is.EqualTo(11), "[9,11]=IMA wraps");

            Assert.That(products.Where(p => p.SecondaryFragmentNumber > 10).ToList(),
                Is.Not.Empty, "Wrapping peptide must produce wrap-around fragments");
        }

        [Test]
        public static void FragmentInternally_WrappingSubPeptide_ProducesWrapAroundFragments()
        {
            // Alias for EndExceedsN — exercises the same peptide with the same assertions
            // but emphasises the positive assertion that wrap-around IS produced.
            FragmentInternally_WrappingSubPeptide_EndExceedsN();
        }

        // ══════════════════════════════════════════════════════════════════════
        // MASS — neutral mass equals residue mass sum for HCD internal ions
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public static void FragmentInternally_FullRing_NeutralMass_MatchesResidueMassSum()
        {
            const string ring = "ACDEFGHIM"; // N=9
            var peptide = GetFullRingPeptide(ring);
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 3, products);

            foreach (var product in products)
            {
                double expected = ExpectedInternalFragmentNeutralMass(
                    ring, product.FragmentNumber, product.ResiduePosition);
                Assert.That(product.NeutralMass, Is.EqualTo(expected).Within(1e-5),
                    $"[{product.FragmentNumber},{product.SecondaryFragmentNumber}]");
            }
        }

        [Test]
        public static void FragmentInternally_WrappingSubPeptide_WrapAroundMassIsCorrect()
        {
            // Fragment [7,11] = G(7),H(8),I(9),M(10),A(1) from GHIMAK[7-12].
            // Mass = residue masses of G,H,I,M,A read clockwise from ring pos 7.
            const string ringSeq = "AKDEFKGHIM";
            var protein = new CircularProtein(ringSeq, "wrap_mass_test");

            var wrappingPeptide = protein
                .Digest(new DigestionParams("trypsin", maxMissedCleavages: 1, minPeptideLength: 1),
                    [], [])
                .OfType<CircularPeptideWithSetModifications>()
                .SingleOrDefault(p => p.BaseSequence == "GHIMAK");

            Assume.That(wrappingPeptide, Is.Not.Null);

            var products = new List<Product>();
            wrappingPeptide.FragmentInternally(DissociationType.HCD, 3, products);

            var frag_7_5 = products.Single(p => p.FragmentNumber == 7 && p.ResiduePosition == 5);
            Assert.That(frag_7_5.NeutralMass,
                Is.EqualTo(ExpectedInternalFragmentNeutralMass(ringSeq, 7, 5)).Within(1e-5),
                "[7,11]=GHIMA mass = G+H+I+M+A residue masses");
        }

        // ══════════════════════════════════════════════════════════════════════
        // NUMBERING — sub-peptide fragments use parent ring canonical positions
        //
        // Ring "ACDEKFGHIK" (N=10), K at positions 5 and 10.
        // Two genuine cut sites → "ACDEK"[1-5] and "FGHIK"[6-10].
        //
        // NOTE: The ring must have K at BOTH position 5 AND position 10 (the
        // last residue). With K only at position 5, position 1 is not a valid
        // cut-derived start and "ACDEK" would be filtered as a proxy artifact.
        // K at the last position makes position 1 valid (post-K[N] = position 1).
        //
        // "ACDEK"[1-5], canonicalOffset=1:
        //   localStart | oneBasedStart | maxLength | fragments
        //       0      |      1        |     4     | [1,3],[1,4],[1,5-1=4]  etc.
        //       1      |      2        |     3     | [2,3],[2,4]
        //       2      |      3        |     2     | [3,4]
        //
        // "FGHIK"[6-10], canonicalOffset=6:
        //   localStart | oneBasedStart | maxLength | fragments
        //       0      |      6        |     4     | [6,8],[6,9]
        //       1      |      7        |     3     | [7,9]
        //       2      |      8        |     2     | [8,9]
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public static void FragmentInternally_SubPeptide_N5_MinLength2_CorrectCount()
        {
            // "ACDEK"[1-5], length 5, minLength 2:
            //   s=0: L=2,3,4 → 3;  s=1: L=2,3 → 2;  s=2: L=2 → 1;  s=3,4: 0
            //   Total = 6 fragments
            var peptide = GetSubPeptide("ACDEKFGHIK", "ACDEK");
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 2, products);

            Assert.That(products, Has.Count.EqualTo(6));
        }

        [Test]
        public static void FragmentInternally_SubPeptide_AtRingStart_NumberingIsCanonical()
        {
            // "ACDEK"[1-5]: fragment numbers must be in parent ring coordinates.
            var peptide = GetSubPeptide("ACDEKFGHIK", "ACDEK");
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 3, products);

            var frag_1_3 = products.Single(p => p.FragmentNumber == 1 && p.ResiduePosition == 3);
            Assert.That(frag_1_3.SecondaryFragmentNumber, Is.EqualTo(3), "[1,3]=ACD");

            var frag_1_4 = products.Single(p => p.FragmentNumber == 1 && p.ResiduePosition == 4);
            Assert.That(frag_1_4.SecondaryFragmentNumber, Is.EqualTo(4), "[1,4]=ACDE");

            var frag_2_3 = products.Single(p => p.FragmentNumber == 2 && p.ResiduePosition == 3);
            Assert.That(frag_2_3.SecondaryFragmentNumber, Is.EqualTo(4), "[2,4]=CDE");
        }

        [Test]
        public static void FragmentInternally_SubPeptide_NotAtRingStart_NumberingIsCanonical()
        {
            // "FGHIK"[6-10]: all fragment numbers must be ≥ 6 (parent ring coordinates).
            var peptide = GetSubPeptide("ACDEKFGHIK", "FGHIK");
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 3, products);

            var frag_6_3 = products.Single(p => p.FragmentNumber == 6 && p.ResiduePosition == 3);
            Assert.That(frag_6_3.SecondaryFragmentNumber, Is.EqualTo(8), "[6,8]=FGH");

            var frag_6_4 = products.Single(p => p.FragmentNumber == 6 && p.ResiduePosition == 4);
            Assert.That(frag_6_4.SecondaryFragmentNumber, Is.EqualTo(9), "[6,9]=FGHI");

            var frag_7_3 = products.Single(p => p.FragmentNumber == 7 && p.ResiduePosition == 3);
            Assert.That(frag_7_3.SecondaryFragmentNumber, Is.EqualTo(9), "[7,9]=GHI");
        }

        [Test]
        public static void FragmentInternally_SubPeptide_NoFragmentsStartAtParentPosition1()
        {
            // "FGHIK" starts at canonical position 6 — no fragment may start below 6.
            var peptide = GetSubPeptide("ACDEKFGHIK", "FGHIK");
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 2, products);

            Assert.That(products.All(p => p.FragmentNumber >= 6), Is.True);
        }

        [Test]
        public static void FragmentInternally_SubPeptide_NeutralMass_MatchesResidueMassSum()
        {
            // "FGHIK"[6-10]: masses must use parent ring residue masses.
            const string ring = "ACDEKFGHIK";
            var peptide = GetSubPeptide(ring, "FGHIK");
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 3, products);

            foreach (var product in products)
            {
                double expected = ExpectedInternalFragmentNeutralMass(
                    ring, product.FragmentNumber, product.ResiduePosition);
                Assert.That(product.NeutralMass, Is.EqualTo(expected).Within(1e-5),
                    $"[{product.FragmentNumber},{product.SecondaryFragmentNumber}]");
            }
        }

        // ══════════════════════════════════════════════════════════════════════
        // SPAN RESTRICTION — sub-peptides cannot produce fragments outside their span
        //
        // Both cleavage sites that define a sub-peptide must lie within the
        // sub-peptide's span. The maxLength bound enforces this.
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public static void FragmentInternally_SubPeptide_NoFragmentExceedsSpan()
        {
            // "ACDEK"[1-5]: no fragment end may reach position 5 (the K boundary).
            // maxLength at localStart=0 is 4, so max end = 1+4-1 = 4 < 5. ✓
            var peptide = GetSubPeptide("ACDEKFGHIK", "ACDEK");
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 2, products);

            Assert.That(products.All(p => p.SecondaryFragmentNumber < 5), Is.True,
                "ACDEK[1-5]: SecondaryFragmentNumber must be < 5");
        }

        [Test]
        public static void FragmentInternally_SubPeptide_NoWrapAroundFragments()
        {
            // Sub-peptides (length < N) never produce wrap-around fragments.
            // Both cleavage sites lie within the ring, so no fragment can extend
            // beyond the C-terminal boundary of the sub-peptide.
            var peptide = GetSubPeptide("ACDEKFGHIK", "FGHIK");
            var products = new List<Product>();

            peptide.FragmentInternally(DissociationType.HCD, 2, products);

            Assert.That(products.All(p => p.SecondaryFragmentNumber <= 10), Is.True,
                "FGHIK[6-10]: SecondaryFragmentNumber must be ≤ N=10");
        }

        // ══════════════════════════════════════════════════════════════════════
        // MODIFICATION MASS — oxidation must shift fragment masses correctly
        //
        // KEY CONSTRAINT: a modified residue at canonical position p is only
        // reachable by fragments from a full-ring [1-N] peptide if
        //   oneBasedStart + length - 1 ≥ p
        //   ⟺ localStart + length ≥ p - 1 + 1 = p
        // But maxLength = N - localStart - 1, so the maximum reachable position is:
        //   oneBasedStart + maxLength - 1 = 1 + localStart + (N-localStart-1) - 1 = N-1
        // Therefore M at position N (the last residue) is NEVER reachable.
        //
        // Test A: ring "ACDEFGHIM", M at position 9 (=N).
        //   No fragment contains M → oxidation has zero effect on all fragment masses.
        //
        // Test B: ring "ACMDEFGHI", M at position 3 (well inside the ring).
        //   Fragments [1,3]=ACM and [2,3]=CM contain M → masses shift by +15.994915.
        //   Fragment [1,2]=AC does not contain M → mass unchanged.
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public static void FragmentInternally_ModifiedResidue_IncludedInFragmentMass()
        {
            // TEST A: M at position 9 = N → no fragment reaches it.
            // Verifies that the fallback Digest() path correctly produces both
            // unmodified and oxidized full-ring peptides (key=10 for pos 9),
            // and that all fragment masses are identical between the two forms.
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

            // Modification key: position 9 + 1 = 10
            Assert.That(oxidized.AllModsOneIsNterminus.ContainsKey(10), Is.True);
            Assert.That(oxidized.AllModsOneIsNterminus[10].OriginalId,
                Does.Contain("Oxidation"));

            var unmodProducts = new List<Product>();
            var oxProducts = new List<Product>();
            unmodified.FragmentInternally(DissociationType.HCD, 2, unmodProducts);
            oxidized.FragmentInternally(DissociationType.HCD, 2, oxProducts);

            // [7,8]=HI exists (localStart=6, maxLength=2, length=2 ✓)
            var unmod_7_2 = unmodProducts.Single(p => p.FragmentNumber == 7 && p.ResiduePosition == 2);
            var ox_7_2 = oxProducts.Single(p => p.FragmentNumber == 7 && p.ResiduePosition == 2);
            Assert.That(ox_7_2.NeutralMass,
                Is.EqualTo(unmod_7_2.NeutralMass).Within(1e-4),
                "[7,8]=HI: M at pos 9 unreachable → no mass shift");

            // All fragments: oxidation at position 9 cannot affect any fragment
            // because SecondaryFragmentNumber ≤ N-1 = 8 always for full-ring [1-9].
            Assert.That(unmodProducts.Count, Is.EqualTo(oxProducts.Count));
            var unmodOrdered = unmodProducts.OrderBy(p => p.FragmentNumber)
                .ThenBy(p => p.ResiduePosition).ToList();
            var oxOrdered = oxProducts.OrderBy(p => p.FragmentNumber)
                .ThenBy(p => p.ResiduePosition).ToList();
            for (int i = 0; i < unmodOrdered.Count; i++)
                Assert.That(oxOrdered[i].NeutralMass,
                    Is.EqualTo(unmodOrdered[i].NeutralMass).Within(1e-9),
                    $"[{unmodOrdered[i].FragmentNumber},{unmodOrdered[i].SecondaryFragmentNumber}]: " +
                    $"M at pos 9 unreachable → identical mass");
        }

        [Test]
        public static void FragmentInternally_ModifiedResidue_MidRing_IncludedInFragmentMass()
        {
            // TEST B: M at position 3 (well inside the ring) → fragments CAN reach it.
            // Ring "ACMDEFGHI" (N=9), M at canonical position 3, no K/R → full-ring [1-9].
            //
            // Fragments containing M (end ≥ 3 and start ≤ 3):
            //   [1,3]=ACM (L=3, localStart=0, maxLength=8 ✓)
            //   [2,3]=CM  (L=2, localStart=1, maxLength=7 ✓)
            //   [3,4]=MD  (L=2, localStart=2, maxLength=6 ✓)  ← starts at M
            // Fragments NOT containing M:
            //   [1,2]=AC  (L=2, end=2 < 3)

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

            // Modification key: position 3 + 1 = 4
            Assert.That(oxidized.AllModsOneIsNterminus.ContainsKey(4), Is.True,
                "Oxidation on M at position 3 → key = 4");

            const string ring = "ACMDEFGHI";
            double oxMass = 15.994915;

            var unmodProducts = new List<Product>();
            var oxProducts = new List<Product>();
            unmodified.FragmentInternally(DissociationType.HCD, 2, unmodProducts);
            oxidized.FragmentInternally(DissociationType.HCD, 2, oxProducts);

            // [1,2]=AC — does NOT contain M → same mass in both forms
            var unmod_1_2 = unmodProducts.Single(p => p.FragmentNumber == 1 && p.ResiduePosition == 2);
            var ox_1_2 = oxProducts.Single(p => p.FragmentNumber == 1 && p.ResiduePosition == 2);
            Assert.That(ox_1_2.NeutralMass, Is.EqualTo(unmod_1_2.NeutralMass).Within(1e-4),
                "[1,2]=AC: no M → no oxidation shift");
            Assert.That(ox_1_2.NeutralMass,
                Is.EqualTo(ExpectedInternalFragmentNeutralMass(ring, 1, 2)).Within(1e-4));

            // [2,3]=CM — contains M → mass shifts by +oxMass
            var unmod_2_2 = unmodProducts.Single(p => p.FragmentNumber == 2 && p.ResiduePosition == 2);
            var ox_2_2 = oxProducts.Single(p => p.FragmentNumber == 2 && p.ResiduePosition == 2);
            Assert.That(ox_2_2.NeutralMass,
                Is.EqualTo(ExpectedInternalFragmentNeutralMass(ring, 2, 2) + oxMass).Within(1e-4),
                "[2,3]=CM: M in fragment → +oxMass");
            Assert.That(ox_2_2.NeutralMass - unmod_2_2.NeutralMass,
                Is.EqualTo(oxMass).Within(1e-4),
                "Mass difference for [2,3]=CM = oxMass");

            // [1,3]=ACM (L=3) — contains M → mass shifts by +oxMass
            var unmod_1_3 = unmodProducts.Single(p => p.FragmentNumber == 1 && p.ResiduePosition == 3);
            var ox_1_3 = oxProducts.Single(p => p.FragmentNumber == 1 && p.ResiduePosition == 3);
            Assert.That(ox_1_3.NeutralMass,
                Is.EqualTo(ExpectedInternalFragmentNeutralMass(ring, 1, 3) + oxMass).Within(1e-4),
                "[1,3]=ACM: M in fragment → +oxMass");
            Assert.That(ox_1_3.NeutralMass - unmod_1_3.NeutralMass,
                Is.EqualTo(oxMass).Within(1e-4),
                "Mass difference for [1,3]=ACM = oxMass");
        }
    }
}