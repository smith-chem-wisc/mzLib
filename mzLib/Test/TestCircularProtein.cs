using Chemistry;
using NUnit.Framework;
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
            // "ACDEKFGHIK" has two K's at positions 5 and 10.
            // maxMissedCleavages: 2 allows the linear engine to produce peptides
            // with up to 2 internal missed cleavage sites.
            //
            // Products:
            //   0 cuts (two actual cuts in proxy → sub-peptides, length < N):
            //     "ACDEK" [1-5]   → CircularPeptideWithSetModifications
            //     "FGHIK" [6-10]  → CircularPeptideWithSetModifications
            //
            //   1 cut (one actual cut → full-length linear, length == N):
            //     "ACDEKFGHIK" [1-10]   → PeptideWithSetModifications  (cut at K10, MC=1 internal K5)
            //     "FGHIKACDEK" [6-15]   → PeptideWithSetModifications  (cut at K5,  MC=1 internal K10)
            //
            // There is no "0 cut" CircularPeptideWithSetModifications of length N here
            // because the proxy digestion always finds at least one cleavage site.
            var circular = new CircularProtein("ACDEKFGHIK", "acc_type2");
            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 2,
                minPeptideLength: 1);

            var peptides = circular.Digest(digestionParams, [], [])
                .Cast<PeptideWithSetModifications>()
                .ToList();

            var linearProducts = peptides.Where(p => p is not CircularPeptideWithSetModifications).ToList();
            var circularProducts = peptides.OfType<CircularPeptideWithSetModifications>().ToList();

            Assert.Multiple(() =>
            {
                // Sub-peptides: two cuts, length < N
                Assert.That(circularProducts.Select(p => p.BaseSequence),
                    Is.EquivalentTo(new[] { "ACDEK", "FGHIK" }));
                Assert.That(circularProducts.All(p => p.Length < circular.Length), Is.True);

                // Linear ring-opening products: one cut, length == N
                Assert.That(linearProducts.Select(p => p.BaseSequence),
                    Is.EquivalentTo(new[] { "ACDEKFGHIK", "FGHIKACDEK" }));
                Assert.That(linearProducts.All(p => p.Length == circular.Length), Is.True);
                Assert.That(linearProducts.All(p => p is not CircularPeptideWithSetModifications), Is.True);

                // All bound to the CircularProtein parent
                Assert.That(peptides.All(p => ReferenceEquals(p.Protein, circular)), Is.True);
            });
        }

        [Test]
        public static void Digest_AllPeptidesAreBoundToCircularProteinParent()
        {
            // Every peptide yielded by CircularProtein.Digest() — regardless of type —
            // must be bound to the originating CircularProtein as its parent.
            // Sub-peptides (length < N) are CircularPeptideWithSetModifications.
            // Full-length products (length == N) are linear PeptideWithSetModifications.
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

                // Sub-peptides (two cuts) must be CircularPeptideWithSetModifications
                Assert.That(subPeptides, Is.Not.Empty);
                Assert.That(subPeptides.All(p => p is CircularPeptideWithSetModifications), Is.True,
                    "Sub-peptides (length < N, two cuts) must be CircularPeptideWithSetModifications");

                // Full-length products (one cut) must be plain PeptideWithSetModifications
                Assert.That(fullLengthPeptides, Is.Not.Empty);
                Assert.That(fullLengthPeptides.All(p => p is not CircularPeptideWithSetModifications), Is.True,
                    "Full-length products (length == N, one cut) must be linear PeptideWithSetModifications");

                // All products are bound to the CircularProtein parent
                Assert.That(peptides.All(p => ReferenceEquals(p.Protein, circular)), Is.True,
                    "Every product from CircularProtein.Digest() must reference the CircularProtein parent");
            });
        }

        [Test]
        public static void Digest_CircularPeptide_CircularParentIsTheSameInstance()
        {
            // The CircularParent property of every CircularPeptideWithSetModifications
            // must reference the exact same CircularProtein instance that was digested.
            // Linear ring-opening products (length == N) are PeptideWithSetModifications
            // and do not have CircularParent, but their Protein property must also
            // reference the same CircularProtein instance.
            var circular = new CircularProtein("ACDEKFGHIK", "acc_type4");
            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 1,
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
                Assert.That(circularPeptides, Is.Not.Empty);
                Assert.That(circularPeptides.All(p => ReferenceEquals(p.CircularParent, circular)), Is.True,
                    "CircularParent must be the same instance as the digested CircularProtein");

                // Linear ring-opening products: Protein must be this instance
                Assert.That(linearPeptides, Is.Not.Empty);
                Assert.That(linearPeptides.All(p => ReferenceEquals(p.Protein, circular)), Is.True,
                    "Linear ring-opening products must also reference the CircularProtein as their Protein");
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
        public static void Digest_TwoLysines_ZeroMissedCleavages_ProducesOnlyCircularPeptides()
        {
            // Two cuts produce two sub-peptides — neither spans the full ring.
            // Both must be CircularPeptideWithSetModifications (they retain circular identity).
            var circular = new CircularProtein("ACDEKFGHIK", "acc_linear4");
            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1);

            var peptides = circular.Digest(digestionParams, [], [])
                .Cast<PeptideWithSetModifications>().ToList();

            Assert.Multiple(() =>
            {
                Assert.That(peptides, Has.Count.EqualTo(2));
                Assert.That(peptides.All(p => p is CircularPeptideWithSetModifications), Is.True,
                    "Sub-peptides from a circular digest must be CircularPeptideWithSetModifications");
                Assert.That(peptides.All(p => p.Length < circular.Length), Is.True);
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
        // Missed cleavages within the linear product do not change the cut count.
        // ──────────────────────────────────────────────────────────────────────

        [Test]
        public static void Digest_OneCut_YieldsLinearPeptideOfLengthN_KAtEnd()
        {
            // K at position 9 (last). One cut → linear "ACDEFGHIK", start=1, end=9.
            // maxMissedCleavages=1 makes no difference — still one cut was made.
            var circular = new CircularProtein("ACDEFGHIK", "acc_cuts1a");
            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 1,
                minPeptideLength: 1);

            var peptides = circular.Digest(digestionParams, [], [])
                .Cast<PeptideWithSetModifications>()
                .ToList();

            Assert.Multiple(() =>
            {
                Assert.That(peptides, Has.Count.EqualTo(1),
                    "One cleavage site produces exactly one linear product");
                Assert.That(peptides[0], Is.Not.InstanceOf<CircularPeptideWithSetModifications>(),
                    "One cut → linear PeptideWithSetModifications");
                Assert.That(peptides[0].Length, Is.EqualTo(circular.Length));
                Assert.That(peptides[0].BaseSequence, Is.EqualTo("ACDEFGHIK"));
                Assert.That(peptides[0].OneBasedStartResidueInProtein, Is.EqualTo(1));
                Assert.That(peptides[0].OneBasedEndResidueInProtein, Is.EqualTo(9));
                Assert.That(ReferenceEquals(peptides[0].Protein, circular), Is.True);
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
        // Two or more cuts: all products are PeptideWithSetModifications, length < N
        // ──────────────────────────────────────────────────────────────────────

        [Test]
        public static void Digest_TwoCuts_YieldsOnlySubPeptides_AllCircularType()
        {
            // Two K's → two cuts → two sub-peptides, both length < N.
            // Both are CircularPeptideWithSetModifications.
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
                    "Two cuts → all products shorter than N");
                Assert.That(peptides.All(p => p is CircularPeptideWithSetModifications), Is.True,
                    "Two cuts → all products are CircularPeptideWithSetModifications");
                Assert.That(peptides.Select(p => p.BaseSequence),
                    Is.EquivalentTo(new[] { "ACDEFK", "GHIK" }));
            });
        }

        [Test]
        public static void Digest_TwoCuts_WithMissedCleavages_LinearProductsAreStillLengthN()
        {
            // Two K's with maxMissedCleavages=1 yields both sub-peptides (two cuts)
            // AND full-length linear products (one cut each, internal K is the missed cleavage).
            // Full-length products are PeptideWithSetModifications;
            // sub-peptides are CircularPeptideWithSetModifications.
            var circular = new CircularProtein("ACDEFKGHIK", "acc_cuts2mc");
            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 1,
                minPeptideLength: 1);

            var peptides = circular.Digest(digestionParams, [], [])
                .Cast<PeptideWithSetModifications>()
                .ToList();

            var linear = peptides.Where(p => p is not CircularPeptideWithSetModifications).ToList();
            var circular_ = peptides.OfType<CircularPeptideWithSetModifications>().ToList();

            Assert.Multiple(() =>
            {
                // Sub-peptides from two cuts
                Assert.That(circular_, Has.Count.EqualTo(2));
                Assert.That(circular_.All(p => p.Length < circular.Length), Is.True);
                Assert.That(circular_.Select(p => p.BaseSequence),
                    Is.EquivalentTo(new[] { "ACDEFK", "GHIK" }));

                // Full-length linear products from one cut each
                Assert.That(linear, Has.Count.EqualTo(2));
                Assert.That(linear.All(p => p.Length == circular.Length), Is.True);
                Assert.That(linear.Select(p => p.BaseSequence),
                    Is.EquivalentTo(new[] { "ACDEFKGHIK", "GHIKACDEFK" }));

                // Canonical numbering for linear products
                var fromK10 = linear.Single(p => p.BaseSequence == "ACDEFKGHIK");
                Assert.That(fromK10.OneBasedStartResidueInProtein, Is.EqualTo(1));
                Assert.That(fromK10.OneBasedEndResidueInProtein, Is.EqualTo(10));

                var fromK6 = linear.Single(p => p.BaseSequence == "GHIKACDEFK");
                Assert.That(fromK6.OneBasedStartResidueInProtein, Is.EqualTo(7));
                Assert.That(fromK6.OneBasedEndResidueInProtein, Is.EqualTo(16));

                // All bound to the CircularProtein parent
                Assert.That(peptides.All(p => ReferenceEquals(p.Protein, circular)), Is.True);
            });
        }
    }
}
