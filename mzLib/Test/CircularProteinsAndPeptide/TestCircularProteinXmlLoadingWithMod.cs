using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test.CircularSearch
{
    /// <summary>
    /// Tests that a UniProt-format protein XML database containing a
    /// <c>&lt;feature type="modified residue" description="Phosphoserine"&gt;</c>
    /// annotation at position 4 can be read via
    /// <see cref="ProteinDbLoader.LoadProteinXML"/> and correctly converted to a
    /// <see cref="CircularProtein"/> via <see cref="CircularProtein.FromProtein"/>,
    /// with the modification position remapped to the canonical rotation.
    ///
    /// DATABASE: simpleCircularProteinWithMod.xml
    /// ------------------------------------------
    /// Accession : P05067
    /// Sequence  : LPGSALLLLA  (length 10, S at one-based position 4)
    /// Feature   : Phosphoserine at position 4
    ///
    /// CANONICALIZATION
    /// ----------------
    /// The lexicographically smallest rotation of "LPGSALLLLA" is "ALLLLALPGS"
    /// (offset 4).
    ///
    ///   Original:  L P G S A L L L L A   (positions 1-10)
    ///   Canonical: A L L L L A L P G S   (positions 1-10)
    ///
    /// RemapModifications formula:
    ///   newPos = (oldPos - 1 - offset + N) % N + 1
    ///   newPos = (4 - 1 - 4 + 10) % 10 + 1 = 9 % 10 + 1 = 10
    ///
    /// The phosphoserine must appear at canonical position 10, which is the S
    /// residue at the end of "ALLLLALPGS".
    ///
    /// MOD MATCHING
    /// ------------
    /// <see cref="ProteinDbLoader.LoadProteinXML"/> matches the XML feature
    /// description ("Phosphoserine") against the keys of <c>IdToPossibleMods</c>,
    /// which are built by stripping " on X" from each known mod's OriginalId.
    /// Therefore the mod passed as <c>allKnownModifications</c> must use
    /// <c>_originalId: "Phosphoserine"</c> (not "Phosphorylation on S") so that
    /// the stripped key "Phosphoserine" matches the XML description exactly.
    /// </summary>
    [TestFixture]
    public static class TestCircularProteinXmlLoadingWithMod
    {
        // ── Constants ─────────────────────────────────────────────────────────

        private const string OriginalSequence = "LPGSALLLLA";
        private const string CanonicalSequence = "ALLLLALPGS";
        private const int CanonicalOffset = 4;

        // One-based position of S in the original sequence.
        private const int PhosphoPositionOriginal = 4;

        // One-based position of S in the canonical sequence after remapping:
        //   (4 - 1 - 4 + 10) % 10 + 1 = 10
        private const int PhosphoPositionCanonical = 10;

        private const string Accession = "P05067";
        private const double PhosphoMass = 79.966331;
        private const double MassTolerance = 1e-5;

        [OneTimeSetUp]
        public static void OneTimeSetup()
        {
            Loaders.LoadElements();
        }

        // ── Helper ────────────────────────────────────────────────────────────

        /// <summary>
        /// Builds the phosphoserine <see cref="Modification"/> whose OriginalId
        /// ("Phosphoserine") matches the XML feature description exactly, enabling
        /// <see cref="ProteinDbLoader.LoadProteinXML"/> to resolve it via
        /// <c>IdToPossibleMods</c>.
        /// </summary>
        private static Modification BuildPhosphoSerineMod()
        {
            ModificationMotif.TryGetMotif("S", out ModificationMotif motifS);
            return new Modification(
                _originalId: "Phosphoserine",
                _modificationType: "UniProt",
                _target: motifS,
                _locationRestriction: "Anywhere.",
                _monoisotopicMass: PhosphoMass);
        }

        private static string XmlPath => Path.Combine(
            TestContext.CurrentContext.TestDirectory,
            "DataFiles",
            "simpleCircularProteinWithMod.xml");

        // ── Tests ─────────────────────────────────────────────────────────────

        /// <summary>
        /// Verifies that the XML database loads without error, that the single
        /// entry converts to a <see cref="CircularProtein"/>, and that the
        /// phosphoserine modification is remapped from original position 4 to
        /// canonical position 10.
        ///
        /// Canonical sequence "ALLLLALPGS":
        ///   Position:           1234567890
        ///
        /// Position 10 = 'S' — the only serine in the ring.
        /// Position  4 = 'L' — must carry no modification after remapping.
        /// </summary>
        [Test]
        public static void LoadProteinXmlWithMod_PhosphoSerineRemappedToCanonicalPosition()
        {
            List<Protein> proteins = ProteinDbLoader.LoadProteinXML(
                XmlPath,
                generateTargets: true,
                decoyType: DecoyType.None,
                allKnownModifications: new List<Modification> { BuildPhosphoSerineMod() },
                isContaminant: false,
                modTypesToExclude: new List<string>(),
                unknownModifications: out _);

            Assert.That(proteins, Is.Not.Null.And.Not.Empty,
                "Expected at least one protein to be loaded from the XML database.");

            Protein source = proteins.Single(p => p.Accession == Accession);

            // ── Convert to CircularProtein ────────────────────────────────────
            CircularProtein circular = CircularProtein.FromProtein(source);

            Assert.Multiple(() =>
            {
                // ── Sequence ──────────────────────────────────────────────────
                Assert.That(circular.BaseSequence, Is.EqualTo(CanonicalSequence),
                    "BaseSequence must be the canonical rotation of the original XML sequence.");

                Assert.That(circular.Length, Is.EqualTo(OriginalSequence.Length),
                    "Circular protein length must equal the original sequence length.");

                // ── Modification is present after canonicalization ────────────
                Assert.That(circular.OneBasedPossibleLocalizedModifications, Is.Not.Empty,
                    "The phosphoserine modification must survive the XML read and canonicalization.");

                // ── Remapped to canonical position 10 ─────────────────────────
                Assert.That(
                    circular.OneBasedPossibleLocalizedModifications.ContainsKey(PhosphoPositionCanonical),
                    Is.True,
                    $"Phosphoserine must be at canonical position {PhosphoPositionCanonical} " +
                    $"(the S in \"{CanonicalSequence}\"), not at the original position {PhosphoPositionOriginal}.");

                // ── The residue at the remapped position is indeed S ──────────
                Assert.That(circular.BaseSequence[PhosphoPositionCanonical - 1], Is.EqualTo('S'),
                    $"Canonical position {PhosphoPositionCanonical} must be residue 'S'.");

                // ── Modification mass ─────────────────────────────────────────
                var modsAtCanonical = circular.OneBasedPossibleLocalizedModifications[PhosphoPositionCanonical];
                Assert.That(modsAtCanonical, Is.Not.Empty,
                    $"Modification list at canonical position {PhosphoPositionCanonical} must not be empty.");

                double? mass = modsAtCanonical.First().MonoisotopicMass;
                Assert.That(mass, Is.Not.Null,
                    "Phosphoserine must have a non-null monoisotopic mass.");
                Assert.That(mass.Value, Is.EqualTo(PhosphoMass).Within(MassTolerance),
                    $"Phosphoserine mass must be ~{PhosphoMass} Da.");

                // ── Original position 4 carries no modification ───────────────
                // After remapping, canonical position 4 = 'L', which is not a
                // serine and must not bear the phospho mod.
                Assert.That(
                    circular.OneBasedPossibleLocalizedModifications.ContainsKey(PhosphoPositionOriginal),
                    Is.False,
                    $"No modification should remain at canonical position {PhosphoPositionOriginal} " +
                    $"(residue '{circular.BaseSequence[PhosphoPositionOriginal - 1]}') after remapping.");

                // ── Metadata ──────────────────────────────────────────────────
                Assert.That(circular.Accession, Is.EqualTo(Accession));
                Assert.That(circular.FullName, Is.EqualTo("Amyloid-beta precursor protein"));
                Assert.That(circular.Organism, Is.EqualTo("Homo sapiens"));
                Assert.That(circular.GeneNames.Any(g => g.Item2 == "APP"), Is.True);
                Assert.That(circular.IsDecoy, Is.False);
                Assert.That(circular.IsContaminant, Is.False);
            });
        }

        /// <summary>
        /// Guards against silent mod-resolution failure. If the OriginalId of the
        /// known mod does not match the XML feature description, the parser silently
        /// places the feature in <c>unknownModifications</c> and the protein carries
        /// no modification at all — a failure mode that is invisible without this test.
        /// </summary>
        [Test]
        public static void LoadProteinXmlWithMod_PhosphoSerineNotInUnknownModifications()
        {
            ProteinDbLoader.LoadProteinXML(
                XmlPath,
                generateTargets: true,
                decoyType: DecoyType.None,
                allKnownModifications: new List<Modification> { BuildPhosphoSerineMod() },
                isContaminant: false,
                modTypesToExclude: new List<string>(),
                unknownModifications: out var unknownMods);

            Assert.That(
                unknownMods.ContainsKey("Phosphoserine"),
                Is.False,
                "Phosphoserine must be resolved against allKnownModifications and must NOT " +
                "appear in unknownModifications. If it does, the OriginalId does not match " +
                "the XML feature description.");
        }

        /// <summary>
        /// Verifies that when the circular protein is digested with phosphoserine
        /// as a variable modification, at least one peptide form carries a mass
        /// shift of ~79.97 Da, and that the modification always lands on an S residue.
        ///
        /// "ALLLLALPGS" has no K or R, so trypsin makes zero cuts and the ring is
        /// returned as one CircularPeptideWithSetModifications. With phosphoserine
        /// as a variable mod we expect exactly two forms: unmodified and phosphorylated.
        /// </summary>
        [Test]
        public static void LoadProteinXmlWithMod_VariableDigestionProducesPhosphorylatedForm()
        {
            var phosphoMod = BuildPhosphoSerineMod();

            List<Protein> proteins = ProteinDbLoader.LoadProteinXML(
                XmlPath,
                generateTargets: true,
                decoyType: DecoyType.None,
                allKnownModifications: new List<Modification> { phosphoMod },
                isContaminant: false,
                modTypesToExclude: new List<string>(),
                unknownModifications: out _);

            CircularProtein circular = CircularProtein.FromProtein(
                proteins.Single(p => p.Accession == Accession));

            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 2,
                minPeptideLength: 1);

            var allPeptides = circular
                .Digest(digestionParams,
                        allKnownFixedModifications: new List<Modification>(),
                        variableModifications: new List<Modification> { phosphoMod })
                .Cast<PeptideWithSetModifications>()
                .ToList();

            Assert.That(allPeptides, Is.Not.Empty,
                "Digestion must yield at least one peptide.");

            // All products must be parented to this CircularProtein.
            Assert.That(allPeptides.All(p => ReferenceEquals(p.Protein, circular)), Is.True,
                "Every digestion product must reference the CircularProtein as its parent.");

            // At least one form must carry the phospho mass shift.
            var phosphoPeptides = allPeptides
                .Where(p => p.AllModsOneIsNterminus.Values.Any(m =>
                    m.MonoisotopicMass.HasValue &&
                    Math.Abs(m.MonoisotopicMass.Value - PhosphoMass) < MassTolerance))
                .ToList();

            Assert.That(phosphoPeptides, Is.Not.Empty,
                "At least one peptide must carry the phosphoserine modification.");

            // In every phosphorylated form the mod must sit on 'S'.
            foreach (var pep in phosphoPeptides)
            {
                foreach (var kvp in pep.AllModsOneIsNterminus)
                {
                    if (!kvp.Value.MonoisotopicMass.HasValue) continue;
                    if (Math.Abs(kvp.Value.MonoisotopicMass.Value - PhosphoMass) >= MassTolerance) continue;

                    // AllModsOneIsNterminus key convention:
                    //   key 1       → N-terminus
                    //   key 2..N+1  → residues (key - 2 = 0-based index into BaseSequence)
                    //   key N+2     → C-terminus
                    int zeroBasedIndex = kvp.Key - 2;
                    char residue = pep.BaseSequence[zeroBasedIndex];

                    Assert.That(residue, Is.EqualTo('S'),
                        $"Phosphoserine in peptide \"{pep.BaseSequence}\" at AllMods key {kvp.Key} " +
                        $"must be on residue 'S', not '{residue}'.");
                }
            }
        }
    }
}