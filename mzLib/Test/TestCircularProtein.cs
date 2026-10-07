using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class TestCircularProtein
    {
        [OneTimeSetUp]
        public static void OneTimeSetup()
        {
            Loaders.LoadElements();
        }

        // ──────────────────────────────────────────────────────────────────────
        // GetCanonicalRotation
        // ──────────────────────────────────────────────────────────────────────

        [TestCase("ACGT", "ACGT", TestName = "AlreadyCanonical")]
        [TestCase("CGTA", "ACGT", TestName = "RotationByOne")]
        [TestCase("GTAC", "ACGT", TestName = "RotationByTwo")]
        [TestCase("TACG", "ACGT", TestName = "RotationByThree")]
        [TestCase("A",    "A",    TestName = "SingleCharacter")]
        [TestCase("",     "",     TestName = "EmptyString")]
        public static void GetCanonicalRotation_ReturnsSmallestRotation(string input, string expected)
        {
            Assert.That(CircularProtein.GetCanonicalRotation(input), Is.EqualTo(expected));
        }

        // ──────────────────────────────────────────────────────────────────────
        // Constructor
        // ──────────────────────────────────────────────────────────────────────

        [Test]
        public static void Constructor_CanonicalizesInputSequence()
        {
            // "TACG" is not the smallest rotation; "ACGT" is
            var circular = new CircularProtein("TACG", "acc1");

            Assert.Multiple(() =>
            {
                Assert.That(circular.BaseSequence, Is.EqualTo("ACGT"));
                Assert.That(circular.CircularSequence, Is.EqualTo(circular.BaseSequence));
                Assert.That(circular.Length, Is.EqualTo(4));
            });
        }

        [Test]
        public static void Constructor_PreservesMetadata()
        {
            var geneNames = new List<Tuple<string, string>> { Tuple.Create("primary", "CYCA") };
            var circular = new CircularProtein(
                sequence: "ACGT",
                accession: "ACC001",
                organism: "Homo sapiens",
                geneNames: geneNames,
                name: "cyclic test",
                fullName: "cyclic test protein",
                isDecoy: false,
                isContaminant: true);

            Assert.Multiple(() =>
            {
                Assert.That(circular.Accession, Is.EqualTo("ACC001"));
                Assert.That(circular.Organism, Is.EqualTo("Homo sapiens"));
                Assert.That(circular.Name, Is.EqualTo("cyclic test"));
                Assert.That(circular.FullName, Is.EqualTo("cyclic test protein"));
                Assert.That(circular.IsDecoy, Is.False);
                Assert.That(circular.IsContaminant, Is.True);
                Assert.That(circular.GeneNames, Has.Count.EqualTo(1));
            });
        }

        // ──────────────────────────────────────────────────────────────────────
        // CyclicMonoisotopicMass
        // ──────────────────────────────────────────────────────────────────────

        [Test]
        public static void CyclicMonoisotopicMass_IsLinearMassMinusWater()
        {
            // A=Ala C=Cys G=Gly T=Thr — all standard residues, already canonical
            const string sequence = "ACGT";
            var circular = new CircularProtein(sequence, "acc1");

            var linearPeptide = new PeptideWithSetModifications(sequence, new Dictionary<string, Modification>());
            double waterMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2
                             + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;

            Assert.That(circular.CyclicMonoisotopicMass,
                Is.EqualTo(linearPeptide.MonoisotopicMass - waterMass).Within(1e-5));
        }

        // ──────────────────────────────────────────────────────────────────────
        // FromProtein
        // ──────────────────────────────────────────────────────────────────────

        [Test]
        public static void FromProtein_CanonicalizesAndPreservesMetadata()
        {
            var source = new Protein(
                "TACG",
                "ACC002",
                organism: "Mus musculus",
                name: "source protein",
                isDecoy: true);

            var circular = CircularProtein.FromProtein(source);

            Assert.Multiple(() =>
            {
                Assert.That(circular, Is.InstanceOf<CircularProtein>());
                Assert.That(circular.BaseSequence, Is.EqualTo("ACGT"));
                Assert.That(circular.Accession, Is.EqualTo("ACC002"));
                Assert.That(circular.Organism, Is.EqualTo("Mus musculus"));
                Assert.That(circular.Name, Is.EqualTo("source protein"));
                Assert.That(circular.IsDecoy, Is.True);
            });
        }

        // ──────────────────────────────────────────────────────────────────────
        // Digest
        //
        // Sequences are designed so the reasoning is auditable by inspection:
        //
        //   "ACDEFK"    N=6,  one K at the last position.
        //               Proxy = "ACDEFKACDEF" → K only at proxy pos 6.
        //               One cut → one linearised full-circle peptide "ACDEFK".
        //
        //   "ACDEKFGHIK" N=10, K at positions 5 and 10 (last).
        //               Proxy = "ACDEKFGHIKACDEKFGHI" → K at proxy pos 5, 10, 15.
        //               0 MC: arcs [1-5]="ACDEK" and [6-10]="FGHIK".
        //               1 MC: adds [1-10]="ACDEKFGHIK" (end=N, non-wrapping)
        //                     and   [6-15]="FGHIKACDEK" (end>N, wrapping).
        // ──────────────────────────────────────────────────────────────────────

        [Test]
        public static void Digest_OneLysine_ProducesOnePeptideOfFullLength()
        {
            var circular = new CircularProtein("ACDEFK", "acc_digest1");
            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1);

            var peptides = circular.Digest(digestionParams, [], [])
                .Cast<PeptideWithSetModifications>().ToList();

            Assert.Multiple(() =>
            {
                Assert.That(peptides, Has.Count.EqualTo(1));
                Assert.That(peptides[0].BaseSequence, Is.EqualTo("ACDEFK"));
                Assert.That(peptides[0].Length, Is.EqualTo(circular.Length));
            });
        }

        [Test]
        public static void Digest_TwoLysines_ZeroMissedCleavages_ReturnsTwoShortPeptides()
        {
            var circular = new CircularProtein("ACDEKFGHIK", "acc_digest2");
            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1);

            var peptides = circular.Digest(digestionParams, [], [])
                .Cast<PeptideWithSetModifications>()
                .OrderBy(p => p.BaseSequence)
                .ToList();

            Assert.Multiple(() =>
            {
                Assert.That(peptides, Has.Count.EqualTo(2));
                Assert.That(peptides[0].BaseSequence, Is.EqualTo("ACDEK"));
                Assert.That(peptides[1].BaseSequence, Is.EqualTo("FGHIK"));
                Assert.That(peptides.All(p => p.Length < circular.Length), Is.True);
            });
        }

        [Test]
        public static void Digest_TwoLysines_OneMissedCleavage_ReturnsFourPeptides()
        {
            // 0 MC: "ACDEK", "FGHIK"
            // 1 MC: "ACDEKFGHIK" (full-length, non-wrapping)
            //       "FGHIKACDEK" (full-length, wrapping — starts at position 6 of circle)
            var circular = new CircularProtein("ACDEKFGHIK", "acc_digest3");
            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 1,
                minPeptideLength: 1);

            var sequences = circular.Digest(digestionParams, [], [])
                .Cast<PeptideWithSetModifications>()
                .Select(p => p.BaseSequence)
                .ToHashSet();

            Assert.Multiple(() =>
            {
                Assert.That(sequences, Has.Count.EqualTo(4));
                Assert.That(sequences, Does.Contain("ACDEK"));
                Assert.That(sequences, Does.Contain("FGHIK"));
                Assert.That(sequences, Does.Contain("ACDEKFGHIK"));
                Assert.That(sequences, Does.Contain("FGHIKACDEK"));
            });
        }

        [Test]
        public static void Digest_WrappingPeptide_HasEndResidueGreaterThanN()
        {
            // "FGHIKACDEK" spans the circular junction:
            // proxy positions [6, 15], end=15 > N=10.
            var circular = new CircularProtein("ACDEKFGHIK", "acc_digest4");
            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 1,
                minPeptideLength: 1);

            var peptides = circular.Digest(digestionParams, [], [])
                .Cast<PeptideWithSetModifications>().ToList();

            var wrapping = peptides.Single(p => p.BaseSequence == "FGHIKACDEK");
            Assert.Multiple(() =>
            {
                Assert.That(wrapping.OneBasedStartResidueInProtein, Is.LessThanOrEqualTo(circular.Length));
                Assert.That(wrapping.OneBasedEndResidueInProtein, Is.GreaterThan(circular.Length));
            });
        }

        [Test]
        public static void Digest_AllPeptidesAreBoundToCircularProtein()
        {
            // Rebinding replaces the transient proxy Protein with `this`.
            var circular = new CircularProtein("ACDEKFGHIK", "acc_digest5");
            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 1,
                minPeptideLength: 1);

            var peptides = circular.Digest(digestionParams, [], [])
                .Cast<PeptideWithSetModifications>().ToList();

            Assert.That(peptides, Is.Not.Empty);
            Assert.That(peptides.All(p => ReferenceEquals(p.Protein, circular)), Is.True);
        }

        // ──────────────────────────────────────────────────────────────────────
        // Modification localization
        //
        // Sequence: "ACDKEPGHIK" (N=10, K at 4 and 10, P at 6)
        // Proxy:    "ACDKEPGHIKACDKEPGHI" (2N-1=19 residues)
        //           P is at proxy position 6 (first copy, ≤ N → in mod dict).
        //
        // Trypsin 1 MC peptides that span the proline:
        //   "EPGHIK"     proxy [5-10]   → allModsKey = 6-5+2  = 3
        //   "ACDKEPGHIK" proxy [1-10]   → allModsKey = 6-1+2  = 7   ← same P, different key!
        //   "EPGHIKACDK" proxy [5-14]   → allModsKey = 6-5+2  = 3   (wrapping, same start)
        // ──────────────────────────────────────────────────────────────────────

        [Test]
        public static void Digest_HydroxyprolineLocalization_DiffersByCleavageSite()
        {
            ModificationMotif.TryGetMotif("P", out ModificationMotif pMotif);
            var hydroxyproline = new Modification(
                _originalId: "Oxidation on P",
                _modificationType: "Common Variable",
                _target: pMotif,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: 15.994915);

            // P is at position 6 in "ACDKEPGHIK" (K at 4 and 10)
            var mods = new Dictionary<int, List<Modification>> { [6] = [hydroxyproline] };
            var circular = new CircularProtein("ACDKEPGHIK", "acc_mod", oneBasedModifications: mods);
            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 1,
                minPeptideLength: 1);

            // Only keep peptides that actually carry the modification
            var modifiedPeptides = circular.Digest(digestionParams, [], [])
                .Cast<PeptideWithSetModifications>()
                .Where(p => p.AllModsOneIsNterminus.Count > 0)
                .ToList();

            var modEpghik     = modifiedPeptides.Single(p => p.BaseSequence == "EPGHIK");
            var modFullLength = modifiedPeptides.Single(p => p.BaseSequence == "ACDKEPGHIK");
            var modWrapping   = modifiedPeptides.Single(p => p.BaseSequence == "EPGHIKACDK");

            Assert.Multiple(() =>
            {
                // "EPGHIK" starts at proxy pos 5 → key = 6-5+2 = 3
                Assert.That(modEpghik.AllModsOneIsNterminus.Keys.Single(), Is.EqualTo(3));

                // "ACDKEPGHIK" starts at proxy pos 1 → key = 6-1+2 = 7  (same physical P, different key!)
                Assert.That(modFullLength.AllModsOneIsNterminus.Keys.Single(), Is.EqualTo(7));

                // "EPGHIKACDK" starts at proxy pos 5 → key = 6-5+2 = 3  (wrapping, same start as "EPGHIK")
                Assert.That(modWrapping.AllModsOneIsNterminus.Keys.Single(), Is.EqualTo(3));

                // All three carry the same modification identity
                Assert.That(modEpghik.AllModsOneIsNterminus.Values.Single().OriginalId, Is.EqualTo("Oxidation"));
                Assert.That(modFullLength.AllModsOneIsNterminus.Values.Single().OriginalId, Is.EqualTo("Oxidation"));
                Assert.That(modWrapping.AllModsOneIsNterminus.Values.Single().OriginalId, Is.EqualTo("Oxidation"));
            });
        }

        [Test]
        public static void Digest_NoCleavageSites_YieldsCircularPeptideWithSetModifications()
        {
            // "ACGMT" contains no K or R — trypsin finds no cleavage sites.
            // Digest() must still yield the full-length ring as a
            // CircularPeptideWithSetModifications.
            var circular = new CircularProtein("ACGMT", "acc_type1");
            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1);

            var peptides = circular.Digest(digestionParams, [], []).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(peptides, Has.Count.EqualTo(1));
                Assert.That(peptides[0], Is.InstanceOf<CircularPeptideWithSetModifications>(),
                    "Digest() must yield CircularPeptideWithSetModifications when no cleavage sites exist");
                Assert.That(((PeptideWithSetModifications)peptides[0]).BaseSequence,
                    Is.EqualTo(circular.BaseSequence),
                    "The yielded peptide must be the full-length canonical ring sequence");
            });
        }

        [Test]
        public static void Digest_NoCleavageSites_CircularParentReferenceIsCorrect()
        {
            var circular = new CircularProtein("ACGMT", "acc_type5");
            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1);

            var peptide = circular.Digest(digestionParams, [], [])
                .Cast<CircularPeptideWithSetModifications>()
                .Single();

            Assert.Multiple(() =>
            {
                Assert.That(peptide, Is.InstanceOf<CircularPeptideWithSetModifications>());
                Assert.That(ReferenceEquals(peptide.CircularParent, circular), Is.True,
                    "CircularParent must reference the originating CircularProtein");
                Assert.That(peptide.CircularParent.CyclicMonoisotopicMass,
                    Is.EqualTo(circular.CyclicMonoisotopicMass).Within(1e-9));
            });
        }
        [Test]
        public static void Digest_TwoKs_MaxMissedCleavagesTwo_ProducesCorrectMixOfTypes()
        {
            // Ring "ACDEKFGHIK" (N=10): K at positions 5 and 10 → 2 cleavage sites.
            // maxMissedCleavages: 2 >= numCleavageSites: 2
            //   → CircularPeptideWithSetModifications IS emitted (Step 2, full ring intact).
            //
            // Step 3 (proxy digestion, maxMissedCleavages=2) produces linear products:
            //   0 missed cleavages (two cuts, length < N):
            //     "ACDEK" [1-5]         → PeptideWithSetModifications
            //     "FGHIK" [6-10]        → PeptideWithSetModifications
            //   1 missed cleavage (one cut, length == N):
            //     "ACDEKFGHIK" [1-10]   → PeptideWithSetModifications
            //     "FGHIKACDEK" [6-15]   → PeptideWithSetModifications
            //   2 missed cleavages: length 15 > N → discarded.
            //
            // Total: 1 circular + 4 linear = 5 products.

            var circular = new CircularProtein("ACDEKFGHIK", "acc_type2");

            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 2,
                minPeptideLength: 1);

            var peptides = circular.Digest(digestionParams, [], [])
                .Cast<PeptideWithSetModifications>()
                .ToList();

            var circularProducts = peptides.OfType<CircularPeptideWithSetModifications>().ToList();
            var linearProducts = peptides.Where(p => p is not CircularPeptideWithSetModifications).ToList();
            var subPeptides = linearProducts.Where(p => p.Length < circular.Length).ToList();
            var fullLengthLinear = linearProducts.Where(p => p.Length == circular.Length).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(peptides, Has.Count.EqualTo(5));

                // Circular product: full ring intact, length == N
                Assert.That(circularProducts, Has.Count.EqualTo(1),
                    "maxMissedCleavages (2) >= numCleavageSites (2): exactly one circular product.");
                Assert.That(circularProducts[0].BaseSequence, Is.EqualTo("ACDEKFGHIK"));
                Assert.That(circularProducts[0].Length, Is.EqualTo(circular.Length));

                // Sub-peptides: two cuts, length < N, linear
                Assert.That(subPeptides, Has.Count.EqualTo(2));
                Assert.That(subPeptides.All(p => p is not CircularPeptideWithSetModifications), Is.True,
                    "Sub-peptides (length < N) are linear PeptideWithSetModifications.");
                Assert.That(subPeptides.Select(p => p.BaseSequence),
                    Is.EquivalentTo(new[] { "ACDEK", "FGHIK" }));

                // Full-length linear products: one cut, length == N
                Assert.That(fullLengthLinear, Has.Count.EqualTo(2));
                Assert.That(fullLengthLinear.All(p => p is not CircularPeptideWithSetModifications), Is.True,
                    "Full-length linear products (one cut) are PeptideWithSetModifications.");
                Assert.That(fullLengthLinear.Select(p => p.BaseSequence),
                    Is.EquivalentTo(new[] { "ACDEKFGHIK", "FGHIKACDEK" }));
                Assert.That(fullLengthLinear.All(p => p.Length == circular.Length), Is.True);

                // All bound to the CircularProtein parent
                Assert.That(peptides.All(p => ReferenceEquals(p.Protein, circular)), Is.True);
            });
        }

        [Test]
        public static void Digest_AllPeptidesAreBoundToCircularProteinParent()
        {
            // Every peptide yielded by CircularProtein.Digest() — regardless of type —
            // must be bound to the originating CircularProtein as its parent.
            //
            // Ring "ACDEKFGHIK" (N=10): K at positions 5 and 10 → 2 cleavage sites.
            // maxMissedCleavages: 1 < numCleavageSites: 2 → NO CircularPeptideWithSetModifications.
            // All products are linear PeptideWithSetModifications.
            //
            // Expected products (trypsin, maxMissedCleavages=1, minPeptideLength=1):
            //   0 missed cleavages: ACDEK [1-5], FGHIK [6-10]           (length < N)
            //   1 missed cleavage:  ACDEKFGHIK [1-10], FGHIKACDEK [6-15] (length == N)

            var circular = new CircularProtein("ACDEKFGHIK", "acc_type3");

            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 1,
                minPeptideLength: 1);

            var peptides = circular.Digest(digestionParams, [], [])
                .Cast<PeptideWithSetModifications>()
                .ToList();

            var subPeptides = peptides.Where(p => p.Length < circular.Length).ToList();
            var fullLengthPeptides = peptides.Where(p => p.Length == circular.Length).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(peptides, Is.Not.Empty);

                // No CircularPeptideWithSetModifications: maxMissedCleavages (1) < numCleavageSites (2)
                Assert.That(peptides.All(p => p is not CircularPeptideWithSetModifications), Is.True,
                    "maxMissedCleavages (1) < numCleavageSites (2): no circular product expected.");

                // Sub-peptides (length < N) are linear PeptideWithSetModifications — two cuts opened the ring
                Assert.That(subPeptides, Is.Not.Empty);
                Assert.That(subPeptides.All(p => p is not CircularPeptideWithSetModifications), Is.True,
                    "Sub-peptides (length < N) are linear PeptideWithSetModifications, never circular.");

                // Full-length products (length == N) are single-cut ring-opening linear products
                Assert.That(fullLengthPeptides, Is.Not.Empty);
                Assert.That(fullLengthPeptides.All(p => p is not CircularPeptideWithSetModifications), Is.True,
                    "Full-length products (length == N, one cut) are linear PeptideWithSetModifications.");

                // All products are bound to the CircularProtein parent
                Assert.That(peptides.All(p => ReferenceEquals(p.Protein, circular)), Is.True,
                    "Every product from CircularProtein.Digest() must reference the CircularProtein parent.");
            });
        }

        [Test]
        public static void Digest_CircularPeptide_CircularParentIsTheSameInstance()
        {
            // The CircularParent property of every CircularPeptideWithSetModifications
            // must reference the exact same CircularProtein instance that was digested.
            // Linear products (PeptideWithSetModifications) must also reference the same
            // CircularProtein instance via their Protein property.
            //
            // Ring "ACDEKFGHIK" (N=10): K at positions 5 and 10 → 2 cleavage sites.
            // maxMissedCleavages: 2 >= numCleavageSites: 2 → CircularPeptideWithSetModifications IS emitted.
            //
            // Expected products:
            //   Circular: ACDEKFGHIK [1-10]          (maxMissedCleavages absorbs both sites)
            //   Linear:   ACDEK [1-5]                (0 missed cleavages)
            //   Linear:   FGHIK [6-10]               (0 missed cleavages)
            //   Linear:   ACDEKFGHIK [1-10]          (1 missed cleavage, length == N)
            //   Linear:   FGHIKACDEK [6-15]          (1 missed cleavage, wraps, length == N)

            var circular = new CircularProtein("ACDEKFGHIK", "acc_type4");

            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 2,
                minPeptideLength: 1);

            var peptides = circular.Digest(digestionParams, [], [])
                .Cast<PeptideWithSetModifications>()
                .ToList();

            var circularPeptides = peptides.OfType<CircularPeptideWithSetModifications>().ToList();
            var linearPeptides = peptides.Where(p => p is not CircularPeptideWithSetModifications).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(peptides, Is.Not.Empty);

                // CircularPeptideWithSetModifications: CircularParent must be this instance
                Assert.That(circularPeptides, Is.Not.Empty,
                    "maxMissedCleavages (2) >= numCleavageSites (2): circular product must be emitted.");
                Assert.That(circularPeptides.All(p => ReferenceEquals(p.CircularParent, circular)), Is.True,
                    "CircularParent must be the same instance as the digested CircularProtein.");

                // Linear products: Protein must be this instance
                Assert.That(linearPeptides, Is.Not.Empty);
                Assert.That(linearPeptides.All(p => ReferenceEquals(p.Protein, circular)), Is.True,
                    "All linear products must also reference the CircularProtein as their Protein.");
            });
        }

        [Test]
        public static void Digest_CircularPeptideMonoisotopicMass_IsLinearMassMinusWater()
        {
            // "ACDEFGHIM" has no K or R — trypsin makes zero cuts — the ring is
            // intact and yields a CircularPeptideWithSetModifications of length N.
            // Its mass must equal the linear peptide mass minus one water molecule.
            const string sequence = "ACDEFGHIM";
            var circular = new CircularProtein(sequence, "acc_mass");
            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1);

            var peptide = circular.Digest(digestionParams, [], [])
                .OfType<CircularPeptideWithSetModifications>()
                .Single();

            var linearPeptide = new PeptideWithSetModifications(
                sequence, new Dictionary<string, Modification>());

            double waterMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2
                             + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;

            Assert.Multiple(() =>
            {
                Assert.That(peptide.MonoisotopicMass,
                    Is.EqualTo(linearPeptide.MonoisotopicMass - waterMass).Within(1e-5),
                    "Circular peptide mass must equal linear peptide mass minus H₂O");

                Assert.That(peptide.MonoisotopicMass,
                    Is.EqualTo(circular.CyclicMonoisotopicMass).Within(1e-5),
                    "Circular peptide mass must equal parent CyclicMonoisotopicMass");
            });
        }
        // ──────────────────────────────────────────────────────────────────────
        // Single cleavage site — linear ring-opening product
        // ──────────────────────────────────────────────────────────────────────

        [Test]
        public static void Digest_OneLysine_ZeroMissedCleavages_ProducesLinearPeptide()
        {
            // One K at position 6: a single tryptic cut opens the ring, producing
            // a linear PeptideWithSetModifications, NOT a CircularPeptideWithSetModifications.
            // The ring-opening ion is mass-equivalent to the precursor and carries
            // no additional circular identity.
            var circular = new CircularProtein("ACDEFK", "acc_linear1");
            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1);

            var peptides = circular.Digest(digestionParams, [], [])
                .Cast<PeptideWithSetModifications>().ToList();

            Assert.Multiple(() =>
            {
                Assert.That(peptides, Has.Count.EqualTo(1));
                Assert.That(peptides[0].BaseSequence, Is.EqualTo("ACDEFK"));
                Assert.That(peptides[0].Length, Is.EqualTo(circular.Length));

                // Must be a plain linear peptide, not circular
                Assert.That(peptides[0], Is.Not.InstanceOf<CircularPeptideWithSetModifications>(),
                    "A single-cut ring-opening product must be a linear PeptideWithSetModifications");

                // But it must still be bound to the CircularProtein parent
                // so that canonical numbering is preserved
                Assert.That(ReferenceEquals(peptides[0].Protein, circular), Is.True,
                    "Linear ring-opening product must reference the CircularProtein parent");
            });
        }

        [Test]
        public static void Digest_OneLysine_ZeroMissedCleavages_NumberingIsCanonical()
        {
            // The single K is at position 6 (last residue) in the canonical sequence "ACDEFK".
            // The ring-opening product spans the entire canonical sequence [1, 6].
            var circular = new CircularProtein("ACDEFK", "acc_linear2");
            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1);

            var peptide = circular.Digest(digestionParams, [], [])
                .Cast<PeptideWithSetModifications>()
                .Single();

            Assert.Multiple(() =>
            {
                Assert.That(peptide.OneBasedStartResidueInProtein, Is.EqualTo(1),
                    "Start residue must be 1 in the canonical ring numbering");
                Assert.That(peptide.OneBasedEndResidueInProtein, Is.EqualTo(circular.Length),
                    "End residue must equal N in the canonical ring numbering");
            });
        }

        [Test]
        public static void Digest_OneLysine_ZeroMissedCleavages_NumberingReflectsCleavageSitePosition()
        {
            // "KACDEF" canonicalizes to "ACDEFK" (A < K).
            // The K moves from position 1 (input) to position 6 (canonical).
            // The ring-opening product must start at canonical position 1.
            var circular = new CircularProtein("KACDEF", "acc_linear3");
            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1);

            var peptide = circular.Digest(digestionParams, [], [])
                .Cast<PeptideWithSetModifications>()
                .Single();

            Assert.Multiple(() =>
            {
                // After canonicalization, the sequence is "ACDEFK"
                Assert.That(circular.BaseSequence, Is.EqualTo("ACDEFK"));
                // The product spans the full canonical sequence
                Assert.That(peptide.BaseSequence, Is.EqualTo("ACDEFK"));
                Assert.That(peptide.OneBasedStartResidueInProtein, Is.EqualTo(1));
                Assert.That(peptide.OneBasedEndResidueInProtein, Is.EqualTo(6));
                // Not circular
                Assert.That(peptide, Is.Not.InstanceOf<CircularPeptideWithSetModifications>());
            });
        }

        [Test]
        public static void Digest_TwoLysines_ZeroMissedCleavages_ProducesOnlyLinearSubPeptides()
        {
            // Ring "ACDEKFGHIK" (N=10): K at positions 5 and 10 → 2 cleavage sites.
            // maxMissedCleavages: 0 < numCleavageSites: 2
            //   → no CircularPeptideWithSetModifications emitted.
            //   → two cuts produce two linear sub-peptides, both length < N.
            //
            // Sub-peptides are PeptideWithSetModifications — the ring has been opened
            // and they are linear fragments with free termini. They are NOT circular.

            var circular = new CircularProtein("ACDEKFGHIK", "acc_linear4");

            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1);

            var peptides = circular.Digest(digestionParams, [], [])
                .Cast<PeptideWithSetModifications>()
                .ToList();

            Assert.Multiple(() =>
            {
                Assert.That(peptides, Has.Count.EqualTo(2));

                Assert.That(peptides.All(p => p is not CircularPeptideWithSetModifications), Is.True,
                    "Sub-peptides from two cuts are linear PeptideWithSetModifications, not circular.");

                Assert.That(peptides.All(p => p.Length < circular.Length), Is.True,
                    "Both sub-peptides must be shorter than the full ring length N.");
            });
        }

        [Test]
        public static void Digest_ZeroCuts_NoCleavageSites_YieldsCircularPeptideOfLengthN()
        {
            // "ACDEFG" has no K or R — trypsin makes zero cuts.
            // The ring is intact → CircularPeptideWithSetModifications of length N.
            // This is the correct way to obtain a full-length CircularPeptideWithSetModifications:
            // not via missed cleavages, but via absence of cleavage sites entirely.
            var circular = new CircularProtein("ACDEFGM", "acc_zero_cuts");
            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1);

            var peptides = circular.Digest(digestionParams, [], []).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(peptides, Has.Count.EqualTo(1));
                Assert.That(peptides[0], Is.InstanceOf<CircularPeptideWithSetModifications>(),
                    "Zero cuts → ring intact → CircularPeptideWithSetModifications");
                Assert.That(((PeptideWithSetModifications)peptides[0]).Length,
                    Is.EqualTo(circular.Length));
                Assert.That(((PeptideWithSetModifications)peptides[0]).BaseSequence,
                    Is.EqualTo(circular.BaseSequence));
            });
        }
        // ──────────────────────────────────────────────────────────────────────
        // Zero cuts: no cleavage sites → CircularPeptideWithSetModifications, length N
        // ──────────────────────────────────────────────────────────────────────

        [Test]
        public static void Digest_ZeroCuts_YieldsCircularPeptideOfLengthN()
        {
            // No K or R → trypsin makes zero cuts → ring is intact.
            var circular = new CircularProtein("ACDEFGHIM", "acc_cuts0");
            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1);

            var peptides = circular.Digest(digestionParams, [], []).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(peptides, Has.Count.EqualTo(1));
                Assert.That(peptides[0], Is.InstanceOf<CircularPeptideWithSetModifications>(),
                    "Zero cuts → ring intact → CircularPeptideWithSetModifications");
                Assert.That(((PeptideWithSetModifications)peptides[0]).Length,
                    Is.EqualTo(circular.Length));
                Assert.That(((PeptideWithSetModifications)peptides[0]).OneBasedStartResidueInProtein,
                    Is.EqualTo(1));
                Assert.That(((PeptideWithSetModifications)peptides[0]).OneBasedEndResidueInProtein,
                    Is.EqualTo(circular.Length));
            });
        }

        // ──────────────────────────────────────────────────────────────────────
        // One cut: exactly one cleavage site → PeptideWithSetModifications, length N
        // ──────────────────────────────────────────────────────────────────────
        [Test]
        public static void Digest_OneCut_YieldsLinearPeptideOfLengthN_KAtEnd()
        {
            // Ring "ACDEFGHIK" (N=9): K at position 9 (last residue) → 1 cleavage site.
            //
            // maxMissedCleavages: 0 < numCleavageSites: 1
            //   → condition maxMissedCleavages >= numCleavageSites is NOT satisfied
            //   → no CircularPeptideWithSetModifications emitted.
            //
            // The single cut opens the ring at K(9), producing one linear product:
            //   "ACDEFGHIK" [1-9], length N, PeptideWithSetModifications.
            //
            // NOTE: maxMissedCleavages: 1 would satisfy 1 >= 1 and additionally emit
            // a CircularPeptideWithSetModifications, giving 2 products total — not 1.
            // maxMissedCleavages: 0 is required to isolate the single-cut linear product.

            var circular = new CircularProtein("ACDEFGHIK", "acc_cuts1a");

            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1);

            var peptides = circular.Digest(digestionParams, [], [])
                .Cast<PeptideWithSetModifications>()
                .ToList();

            Assert.Multiple(() =>
            {
                Assert.That(peptides, Has.Count.EqualTo(1),
                    "One cleavage site with maxMissedCleavages=0 produces exactly one linear product.");

                Assert.That(peptides[0], Is.Not.InstanceOf<CircularPeptideWithSetModifications>(),
                    "One cut → linear PeptideWithSetModifications, never circular.");

                Assert.That(peptides[0].Length, Is.EqualTo(circular.Length),
                    "The single-cut product spans the full ring: length == N.");

                Assert.That(peptides[0].BaseSequence, Is.EqualTo("ACDEFGHIK"));
                Assert.That(peptides[0].OneBasedStartResidueInProtein, Is.EqualTo(1));
                Assert.That(peptides[0].OneBasedEndResidueInProtein, Is.EqualTo(9));

                Assert.That(ReferenceEquals(peptides[0].Protein, circular), Is.True,
                    "The linear product must reference the CircularProtein as its parent.");
            });
        }

        [Test]
        public static void Digest_OneCut_YieldsLinearPeptideOfLengthN_KNotAtEnd()
        {
            // "ACDEFKGHIM": K at position 6, N=10.
            // One cut after K at pos 6 → single linear product "GHIMACDEFK",
            // start=7, end=16. No sub-peptides: "ACDEFK" [1-6] is excluded
            // because its end (6) ≤ N but it requires the proxy endpoint as its
            // second boundary, not a genuine second cut within the ring.
            var circular = new CircularProtein("ACDEFKGHIM", "acc_cuts1b");
            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1);

            var peptides = circular.Digest(digestionParams, [], [])
                .Cast<PeptideWithSetModifications>()
                .ToList();

            Assert.Multiple(() =>
            {
                Assert.That(peptides, Has.Count.EqualTo(1),
                    "One cleavage site produces exactly one product");
                Assert.That(peptides[0], Is.Not.InstanceOf<CircularPeptideWithSetModifications>(),
                    "One cut → linear PeptideWithSetModifications");
                Assert.That(peptides[0].BaseSequence, Is.EqualTo("GHIMACDEFK"),
                    "Linear product runs from residue after K at pos 6, wrapping to K");
                Assert.That(peptides[0].OneBasedStartResidueInProtein, Is.EqualTo(7),
                    "Starts at position 7 — immediately after the K cleavage site");
                Assert.That(peptides[0].OneBasedEndResidueInProtein, Is.EqualTo(16),
                    "Ends at position 16 = 7 + 10 - 1 in doubled-ring numbering");
                Assert.That(ReferenceEquals(peptides[0].Protein, circular), Is.True);
            });
        }
        // ──────────────────────────────────────────────────────────────────────
        // Two cuts (maxMissedCleavages=0): all products are linear
        // PeptideWithSetModifications, length < N.
        // ──────────────────────────────────────────────────────────────────────
        [Test]
        public static void Digest_TwoCuts_YieldsOnlySubPeptides_AllLinearType()
        {
            // Ring "ACDEFKGHIK" (N=10): K at positions 6 and 10 → 2 cleavage sites.
            // maxMissedCleavages: 0 < numCleavageSites: 2
            //   → no CircularPeptideWithSetModifications emitted.
            //   → two cuts produce two linear sub-peptides, both length < N.

            var circular = new CircularProtein("ACDEFKGHIK", "acc_cuts2");

            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1);

            var peptides = circular.Digest(digestionParams, [], [])
                .Cast<PeptideWithSetModifications>()
                .ToList();

            Assert.Multiple(() =>
            {
                Assert.That(peptides, Has.Count.EqualTo(2));

                Assert.That(peptides.All(p => p.Length < circular.Length), Is.True,
                    "Two cuts → all products shorter than N.");

                Assert.That(peptides.All(p => p is not CircularPeptideWithSetModifications), Is.True,
                    "Two cuts → all products are linear PeptideWithSetModifications, never circular.");

                Assert.That(peptides.Select(p => p.BaseSequence),
                    Is.EquivalentTo(new[] { "ACDEFK", "GHIK" }));
            });
        }

        [Test]
        public static void Digest_TwoCuts_WithMissedCleavages_LinearProductsAreStillLengthN()
        {
            // Ring "ACDEFKGHIK" (N=10): K at positions 6 and 10 → 2 cleavage sites.
            // maxMissedCleavages: 1 < numCleavageSites: 2
            //   → no CircularPeptideWithSetModifications emitted.
            //   → all products are linear PeptideWithSetModifications.
            //
            // Expected products:
            //   0 missed cleavages (two cuts):
            //     "ACDEFK"     [1-6]   length < N
            //     "GHIK"       [7-10]  length < N
            //   1 missed cleavage (one cut):
            //     "ACDEFKGHIK" [1-10]  length == N  (cut after K@10, K@6 is the missed cleavage)
            //     "GHIKACDEFK" [7-16]  length == N  (cut after K@6, K@10 is the missed cleavage)

            var circular = new CircularProtein("ACDEFKGHIK", "acc_cuts2mc");

            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 1,
                minPeptideLength: 1);

            var peptides = circular.Digest(digestionParams, [], [])
                .Cast<PeptideWithSetModifications>()
                .ToList();

            var subPeptides = peptides.Where(p => p.Length < circular.Length).ToList();
            var fullLengthLinear = peptides.Where(p => p.Length == circular.Length).ToList();

            Assert.Multiple(() =>
            {
                // No circular product: maxMissedCleavages (1) < numCleavageSites (2)
                Assert.That(peptides.All(p => p is not CircularPeptideWithSetModifications), Is.True,
                    "maxMissedCleavages (1) < numCleavageSites (2): no circular product.");

                // Sub-peptides from two cuts — linear, length < N
                Assert.That(subPeptides, Has.Count.EqualTo(2));
                Assert.That(subPeptides.All(p => p is not CircularPeptideWithSetModifications), Is.True,
                    "Sub-peptides (length < N) are linear PeptideWithSetModifications.");
                Assert.That(subPeptides.Select(p => p.BaseSequence),
                    Is.EquivalentTo(new[] { "ACDEFK", "GHIK" }));

                // Full-length products from one cut — linear, length == N
                Assert.That(fullLengthLinear, Has.Count.EqualTo(2));
                Assert.That(fullLengthLinear.All(p => p is not CircularPeptideWithSetModifications), Is.True,
                    "Full-length products (length == N, one cut) are linear PeptideWithSetModifications.");
                Assert.That(fullLengthLinear.Select(p => p.BaseSequence),
                    Is.EquivalentTo(new[] { "ACDEFKGHIK", "GHIKACDEFK" }));

                // Canonical numbering for full-length linear products
                var fromK10 = fullLengthLinear.Single(p => p.BaseSequence == "ACDEFKGHIK");
                Assert.That(fromK10.OneBasedStartResidueInProtein, Is.EqualTo(1));
                Assert.That(fromK10.OneBasedEndResidueInProtein, Is.EqualTo(10));

                var fromK6 = fullLengthLinear.Single(p => p.BaseSequence == "GHIKACDEFK");
                Assert.That(fromK6.OneBasedStartResidueInProtein, Is.EqualTo(7));
                Assert.That(fromK6.OneBasedEndResidueInProtein, Is.EqualTo(16));

                // All products are bound to the CircularProtein parent
                Assert.That(peptides.All(p => ReferenceEquals(p.Protein, circular)), Is.True);
            });
        }
        // ──────────────────────────────────────────────────────────────────────
        // Modification remapping (fix_001)
        // ──────────────────────────────────────────────────────────────────────

        private static Modification CreateTestMod(char targetResidue, string id = "TestMod", double mass = 42.01057)
        {
            ModificationMotif.TryGetMotif(targetResidue.ToString(), out var motif);
            return new Modification(
                _originalId: id,
                _modificationType: "Common",
                _target: motif,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: mass);
        }

        [Test]
        public static void Constructor_NonCanonicalInputWithModification_ModificationRemappedToCorrectResidue()
        {
            // "GKAEC" → canonical "AECGK", offset = 2
            // Original position 3 (A) should become canonical position 1 (A).
            var testMod = CreateTestMod('A');
            var mods = new Dictionary<int, List<Modification>>
            {
                { 3, new List<Modification> { testMod } }
            };

            var cp = new CircularProtein("GKAEC", "acc_remap1", oneBasedModifications: mods);

            Assert.That(cp.BaseSequence, Is.EqualTo("AECGK"));
            Assert.That(cp.OneBasedPossibleLocalizedModifications.ContainsKey(1), Is.True,
                "Modification should be remapped to position 1 (residue A)");
            Assert.That(cp.OneBasedPossibleLocalizedModifications.ContainsKey(3), Is.False,
                "Position 3 in canonical is residue C, not the original target");
        }

        [Test]
        public static void Constructor_AlreadyCanonicalInput_ModificationsUnchanged()
        {
            // "AECGK" position 2 is 'E'
            var testMod = CreateTestMod('E');
            var mods = new Dictionary<int, List<Modification>>
            {
                { 2, new List<Modification> { testMod } }
            };

            var cp = new CircularProtein("AECGK", "acc_remap2", oneBasedModifications: mods);

            Assert.That(cp.BaseSequence, Is.EqualTo("AECGK"));
            Assert.That(cp.OneBasedPossibleLocalizedModifications.ContainsKey(2), Is.True);
            Assert.That(cp.OneBasedPossibleLocalizedModifications[2], Contains.Item(testMod));
        }

        [Test]
        public static void Constructor_NullModifications_DoesNotThrow()
        {
            Assert.DoesNotThrow(() =>
                new CircularProtein("GKAEC", "acc_remap3", oneBasedModifications: null));
        }

        [Test]
        public static void Constructor_EmptyModifications_DoesNotThrow()
        {
            var mods = new Dictionary<int, List<Modification>>();
            var cp = new CircularProtein("GKAEC", "acc_remap4", oneBasedModifications: mods);

            Assert.That(cp.OneBasedPossibleLocalizedModifications, Is.Empty);
        }

        [Test]
        public static void Constructor_MultipleModificationsRemapped_AllPositionsCorrect()
        {
            // "GKAEC" → "AECGK", offset = 2
            // pos 1 (G) → ((1-1-2+5)%5)+1 = 4
            // pos 4 (E) → ((4-1-2+5)%5)+1 = 2
            var testMod = CreateTestMod('G');
            var mod2 = CreateTestMod('E', "OtherMod", 15.99491);

            var mods = new Dictionary<int, List<Modification>>
            {
                { 1, new List<Modification> { testMod } },
                { 4, new List<Modification> { mod2 } }
            };

            var cp = new CircularProtein("GKAEC", "acc_remap5", oneBasedModifications: mods);

            Assert.That(cp.BaseSequence, Is.EqualTo("AECGK"));
            Assert.That(cp.OneBasedPossibleLocalizedModifications.ContainsKey(4), Is.True,
                "G was at original pos 1, should be at canonical pos 4");
            Assert.That(cp.OneBasedPossibleLocalizedModifications[4], Contains.Item(testMod));
            Assert.That(cp.OneBasedPossibleLocalizedModifications.ContainsKey(2), Is.True,
                "E was at original pos 4, should be at canonical pos 2");
            Assert.That(cp.OneBasedPossibleLocalizedModifications[2], Contains.Item(mod2));
        }

        [Test]
        public static void Constructor_ModificationOnLastResidue_WrapsCorrectly()
        {
            // "GKAEC" → "AECGK", offset = 2
            // pos 5 (C) → ((5-1-2+5)%5)+1 = 3
            var testMod = CreateTestMod('C');
            var mods = new Dictionary<int, List<Modification>>
            {
                { 5, new List<Modification> { testMod } }
            };

            var cp = new CircularProtein("GKAEC", "acc_remap6", oneBasedModifications: mods);

            Assert.That(cp.BaseSequence, Is.EqualTo("AECGK"));
            Assert.That(cp.OneBasedPossibleLocalizedModifications.ContainsKey(3), Is.True,
                "C was at original pos 5, should be at canonical pos 3");
        }

        [Test]
        public static void FromProtein_SourceHasModificationsAndNonCanonicalSequence_ModificationsRemapped()
        {
            // "GKAEC" position 3 is 'A'
            var testMod = CreateTestMod('A');
            var mods = new Dictionary<int, List<Modification>>
            {
                { 3, new List<Modification> { testMod } }
            };

            var source = new Protein("GKAEC", "acc_remap7", oneBasedModifications: mods);
            var cp = CircularProtein.FromProtein(source);

            Assert.That(cp.BaseSequence, Is.EqualTo("AECGK"));
            Assert.That(cp.OneBasedPossibleLocalizedModifications.ContainsKey(1), Is.True,
                "Modification at original pos 3 (A) should map to canonical pos 1 (A)");
        }

        [Test]
        public static void Constructor_SingleResidue_ModificationStaysAtPositionOne()
        {
            var testMod = CreateTestMod('A');
            var mods = new Dictionary<int, List<Modification>>
            {
                { 1, new List<Modification> { testMod } }
            };

            var cp = new CircularProtein("A", "acc_remap8", oneBasedModifications: mods);

            Assert.That(cp.BaseSequence, Is.EqualTo("A"));
            Assert.That(cp.OneBasedPossibleLocalizedModifications.ContainsKey(1), Is.True);
        }

        [Test]
        public static void GetCanonicalRotationWithOffset_ReturnsCorrectOffset()
        {
            var (canonical, offset) = CircularProtein.GetCanonicalRotationWithOffset("GKAEC");

            Assert.That(canonical, Is.EqualTo("AECGK"));
            Assert.That(offset, Is.EqualTo(2));
        }

        [Test]
        public static void GetCanonicalRotationWithOffset_AlreadyCanonical_OffsetIsZero()
        {
            var (canonical, offset) = CircularProtein.GetCanonicalRotationWithOffset("AECGK");

            Assert.That(canonical, Is.EqualTo("AECGK"));
            Assert.That(offset, Is.EqualTo(0));
        }

        [TestCase("KAEC", 1, 4, 'K', TestName = "ModRemap_K_at_pos1_to_pos4")]
        [TestCase("KAEC", 2, 1, 'A', TestName = "ModRemap_A_at_pos2_to_pos1")]
        [TestCase("KAEC", 3, 2, 'E', TestName = "ModRemap_E_at_pos3_to_pos2")]
        [TestCase("KAEC", 4, 3, 'C', TestName = "ModRemap_C_at_pos4_to_pos3")]
        public static void Constructor_ModAtPosition_LandsOnSameAminoAcid(
            string sequence, int originalPos, int expectedCanonicalPos, char expectedResidue)
        {
            var testMod = CreateTestMod(sequence[originalPos - 1]);
            var mods = new Dictionary<int, List<Modification>>
            {
                { originalPos, new List<Modification> { testMod } }
            };

            var cp = new CircularProtein(sequence, "acc_tc", oneBasedModifications: mods);

            Assert.That(cp.OneBasedPossibleLocalizedModifications.ContainsKey(expectedCanonicalPos), Is.True,
                $"Mod should be at canonical position {expectedCanonicalPos}");
            Assert.That(cp.BaseSequence[expectedCanonicalPos - 1], Is.EqualTo(expectedResidue),
                $"Canonical position {expectedCanonicalPos} should be residue {expectedResidue}");
        }

        // ──────────────────────────────────────────────────────────────────────
        // Positional isomer deduplication (fix_009)
        // ──────────────────────────────────────────────────────────────────────

        [Test]
        public static void Digest_RepeatedMotifCircularProtein_RetainsAllPositionalIsomers()
        {
            // "ARARK" canonical → trypsin cuts after R and K.
            // Two distinct "AR" fragments at different ring positions must both survive.
            var protein = new CircularProtein("ARARK", "REPEAT_TEST");

            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1);

            var products = protein
                .Digest(digestionParams, new List<Modification>(), new List<Modification>())
                .ToList();

            var arProducts = products
                .Where(p => p is PeptideWithSetModifications pwsm && pwsm.BaseSequence == "AR")
                .Cast<PeptideWithSetModifications>()
                .ToList();

            Assert.That(arProducts.Count, Is.GreaterThanOrEqualTo(2),
                "Both 'AR' positional isomers at different ring positions must be retained");

            var startPositions = arProducts
                .Select(p => p.OneBasedStartResidueInProtein)
                .Distinct()
                .ToList();

            Assert.That(startPositions.Count, Is.GreaterThanOrEqualTo(2),
                "The 'AR' products must originate at distinct start positions on the ring");
        }

        [Test]
        public static void Digest_RepeatedMotifCircularProtein_TrueDuplicatesStillCollapsed()
        {
            var protein = new CircularProtein("ARARK", "REPEAT_TEST");

            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 2,
                minPeptideLength: 1);

            var products = protein
                .Digest(digestionParams, new List<Modification>(), new List<Modification>())
                .Cast<PeptideWithSetModifications>()
                .Where(p => p is not CircularPeptideWithSetModifications)
                .ToList();

            // No two products should share both the same start position and full sequence
            var keys = products.Select(p =>
                $"{p.OneBasedStartResidueInProtein}:{p.FullSequence}").ToList();

            Assert.That(keys.Distinct().Count(), Is.EqualTo(keys.Count),
                "True duplicates (same position + same sequence) must still be collapsed");
        }

        [Test]
        public static void Digest_SymmetricRepeatProtein_AllPositionsRetained()
        {
            // "ARARAR" — 3 copies of "AR", trypsin cuts after each R
            var symmetricProtein = new CircularProtein("ARARAR", "SYMMETRIC_TEST");

            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1);

            var products = symmetricProtein
                .Digest(digestionParams, new List<Modification>(), new List<Modification>())
                .OfType<PeptideWithSetModifications>()
                .Where(p => p is not CircularPeptideWithSetModifications)
                .Where(p => p.BaseSequence == "AR")
                .ToList();

            Assert.That(products.Count, Is.EqualTo(3),
                "All three 'AR' positional isomers must be retained in a fully symmetric ring");

            var distinctStarts = products
                .Select(p => p.OneBasedStartResidueInProtein)
                .Distinct()
                .ToList();

            Assert.That(distinctStarts.Count, Is.EqualTo(3),
                "Each 'AR' product must have a unique start position");
        }

        [Test]
        public static void Digest_RepeatedMotif_WithMissedCleavages_PositionalIsomersPreserved()
        {
            var protein = new CircularProtein("ARARK", "REPEAT_TEST");

            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 1,
                minPeptideLength: 1);

            var products = protein
                .Digest(digestionParams, new List<Modification>(), new List<Modification>())
                .OfType<PeptideWithSetModifications>()
                .Where(p => p is not CircularPeptideWithSetModifications)
                .ToList();

            // For any sequence that appears at multiple positions, all positions must be present
            var grouped = products
                .GroupBy(p => p.FullSequence)
                .Where(g => g.Select(p => p.OneBasedStartResidueInProtein).Distinct().Count() > 1)
                .ToList();

            foreach (var group in grouped)
            {
                var distinctPositions = group
                    .Select(p => p.OneBasedStartResidueInProtein)
                    .Distinct()
                    .Count();

                Assert.That(group.Count(), Is.EqualTo(distinctPositions),
                    $"Sequence '{group.Key}' has {distinctPositions} distinct positions " +
                    $"but only {group.Count()} products were emitted");
            }
        }
        // ──────────────────────────────────────────────────────────────────────
        // CircularPeptideWithSetModifications constructor validation (fix_010)
        // ──────────────────────────────────────────────────────────────────────

        [Test]
        public static void CircularPeptide_NullProtein_ThrowsArgumentNullException()
        {
            var digestionParams = new DigestionParams();
            var emptyMods = new Dictionary<int, Modification>();

            var ex = Assert.Throws<ArgumentNullException>(() =>
                new CircularPeptideWithSetModifications(
                    protein: null,
                    digestionParams: digestionParams,
                    oneBasedStartResidueInProtein: 1,
                    oneBasedEndResidueInProtein: 5,
                    cleavageSpecificity: CleavageSpecificity.Full,
                    peptideDescription: "test",
                    missedCleavages: 0,
                    allModsOneIsNterminus: emptyMods,
                    numFixedMods: 0));

            Assert.That(ex.ParamName, Is.EqualTo("protein"));
            Assert.That(ex.Message, Does.Contain("non-null"));
        }

        [Test]
        public static void CircularPeptide_NonCircularProtein_ThrowsArgumentException()
        {
            var linearProtein = new Protein("ACDEFG", "LINEAR1");
            var digestionParams = new DigestionParams();
            var emptyMods = new Dictionary<int, Modification>();

            var ex = Assert.Throws<ArgumentException>(() =>
                new CircularPeptideWithSetModifications(
                    protein: linearProtein,
                    digestionParams: digestionParams,
                    oneBasedStartResidueInProtein: 1,
                    oneBasedEndResidueInProtein: 6,
                    cleavageSpecificity: CleavageSpecificity.Full,
                    peptideDescription: "test",
                    missedCleavages: 0,
                    allModsOneIsNterminus: emptyMods,
                    numFixedMods: 0));

            Assert.That(ex.ParamName, Is.EqualTo("protein"));
            Assert.That(ex.Message, Does.Contain("CircularProtein"));
        }

        [Test]
        public static void CircularPeptide_ValidCircularProtein_SetsCircularParent()
        {
            var circularProtein = new CircularProtein("ACDEFG", "CIRC1");
            var digestionParams = new DigestionParams();
            var emptyMods = new Dictionary<int, Modification>();

            var peptide = new CircularPeptideWithSetModifications(
                protein: circularProtein,
                digestionParams: digestionParams,
                oneBasedStartResidueInProtein: 1,
                oneBasedEndResidueInProtein: 6,
                cleavageSpecificity: CleavageSpecificity.Full,
                peptideDescription: "valid",
                missedCleavages: 0,
                allModsOneIsNterminus: emptyMods,
                numFixedMods: 0);

            Assert.That(peptide.CircularParent, Is.Not.Null);
            Assert.That(peptide.CircularParent, Is.SameAs(circularProtein));
        }
        // ──────────────────────────────────────────────────────────────────────
        // Wrap-around fragmentation (fix_013)
        // ──────────────────────────────────────────────────────────────────────

        private static CircularPeptideWithSetModifications BuildFullRingPeptide(string sequence)
        {
            var protein = new CircularProtein(sequence, "TEST");
            int n = protein.BaseSequence.Length;
            return new CircularPeptideWithSetModifications(
                protein,
                new DigestionParams(),
                oneBasedStartResidueInProtein: 1,
                oneBasedEndResidueInProtein: n,
                cleavageSpecificity: CleavageSpecificity.Full,
                peptideDescription: "full-ring",
                missedCleavages: 0,
                allModsOneIsNterminus: new Dictionary<int, Modification>(),
                numFixedMods: 0,
                baseSequence: protein.BaseSequence);
        }

        [Test]
        public static void FragmentInternally_FullRingLength5_CountEqualsNTimesNMinusMinLength()
        {
            // N=5, minLength=2: expected = 5 * (5-2) = 15
            var peptide = BuildFullRingPeptide("ACDEF");
            var products = new List<Product>();
            peptide.FragmentInternally(DissociationType.HCD, 2, products);

            Assert.That(products.Count, Is.EqualTo(15),
                "Full-ring peptide of length 5 should produce N*(N-2) = 15 internal fragments");
        }

        [Test]
        public static void FragmentInternally_FullRing_ContainsWrapAroundFragment()
        {
            var peptide = BuildFullRingPeptide("ACDEF");
            var products = new List<Product>();
            peptide.FragmentInternally(DissociationType.HCD, 2, products);

            // Fragment starting at last residue wrapping to first
            bool hasWrapFragment = products.Any(p =>
                p.FragmentNumber == 5 && p.ResiduePosition == 2);

            Assert.That(hasWrapFragment, Is.True,
                "Should contain a wrap-around fragment starting at residue 5 with length 2");
        }

        [Test]
        public static void FragmentInternally_FullRingLength4_TotalCountIsCorrect()
        {
            // N=4, minLength=2: count = 4 * (4-2) = 8
            var peptide = BuildFullRingPeptide("ACDE");
            var products = new List<Product>();
            peptide.FragmentInternally(DissociationType.HCD, 2, products);

            Assert.That(products.Count, Is.EqualTo(8));
        }

        [Test]
        public static void FragmentInternally_FullRing_AllMassesPositive()
        {
            var peptide = BuildFullRingPeptide("ACDEF");
            var products = new List<Product>();
            peptide.FragmentInternally(DissociationType.HCD, 2, products);

            Assert.That(products.All(p => p.NeutralMass > 0), Is.True,
                "Every fragment mass should be positive, including wrap-around fragments");
        }

        [Test]
        public static void FragmentInternally_FullRing_NoFullLengthFragments()
        {
            var peptide = BuildFullRingPeptide("ACDEF");
            var products = new List<Product>();
            peptide.FragmentInternally(DissociationType.HCD, 2, products);

            Assert.That(products.Any(p => p.ResiduePosition == 5), Is.False,
                "No fragment should span the entire ring (length == N)");
        }

        [Test]
        public static void FragmentInternally_SubPeptide_DoesNotWrapAround()
        {
            var protein = new CircularProtein("ACDEFG", "TEST");
            var subPeptide = new CircularPeptideWithSetModifications(
                protein,
                new DigestionParams(),
                oneBasedStartResidueInProtein: 1,
                oneBasedEndResidueInProtein: 4,
                cleavageSpecificity: CleavageSpecificity.Full,
                peptideDescription: "sub-ring",
                missedCleavages: 0,
                allModsOneIsNterminus: new Dictionary<int, Modification>(),
                numFixedMods: 0,
                baseSequence: "ACDE");

            var products = new List<Product>();
            subPeptide.FragmentInternally(DissociationType.HCD, 2, products);

            // Sub-peptide of length 4, minLength 2: (4-2)(4-2+1)/2 = 3
            Assert.That(products.Count, Is.EqualTo(3),
                "Sub-peptide should not include wrap-around fragments");
        }
    }
}
