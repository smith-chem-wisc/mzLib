using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test.CircularSearch
{
    [TestFixture]
    public static class TestCircularProteinXmlLoading
    {
        [OneTimeSetUp]
        public static void OneTimeSetup()
        {
            Loaders.LoadElements();
        }

        /// <summary>
        /// Verifies that a UniProt-format protein XML database can be read via
        /// <see cref="ProteinDbLoader.LoadProteinXML"/> and that each resulting
        /// <see cref="Protein"/> entry can be successfully converted into a
        /// <see cref="CircularProtein"/> via <see cref="CircularProtein.FromProtein"/>.
        ///
        /// The test database (simpleCircularProtein.xml) contains one entry:
        ///   Accession : P05067
        ///   Sequence  : LPGLALLLLA  (length 10)
        ///   Full name : Amyloid-beta precursor protein
        ///   Organism  : Homo sapiens
        ///   Gene      : APP
        ///
        /// The expected canonical rotation of "LPGLALLLLA" is "ALLLLALPGL"
        /// (lexicographically smallest rotation, offset 4).
        /// </summary>
        [Test]
        public static void LoadProteinXml_SingleEntry_ConvertsToCircularProtein()
        {
            string xmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "simpleCircularProtein.xml");

            // ── Load proteins from XML ────────────────────────────────────────
            List<Protein> proteins = ProteinDbLoader.LoadProteinXML(
                xmlPath,
                generateTargets: true,
                decoyType: DecoyType.None,
                allKnownModifications: new List<Modification>(),
                isContaminant: false,
                modTypesToExclude: new List<string>(),
                unknownModifications: out _);

            Assert.That(proteins, Is.Not.Null);
            Assert.That(proteins, Is.Not.Empty,
                "Expected at least one protein to be loaded from the XML database.");

            // ── Convert each Protein to CircularProtein ───────────────────────
            List<CircularProtein> circularProteins = proteins
                .Select(p => CircularProtein.FromProtein(p))
                .ToList();

            Assert.That(circularProteins, Has.Count.EqualTo(proteins.Count),
                "Every loaded Protein should convert to a CircularProtein without loss.");

            Assert.That(circularProteins.All(cp => cp is CircularProtein), Is.True,
                "FromProtein must return CircularProtein instances.");

            // ── Validate the single entry ─────────────────────────────────────
            CircularProtein entry = circularProteins.Single(cp => cp.Accession == "P05067");

            Assert.Multiple(() =>
            {
                // Sequence must be the canonical (lexicographically smallest) rotation
                // of the original "LPGLALLLLA".
                // Rotations and their canonical: offset 4 → "ALLLLALPGL"
                Assert.That(entry.BaseSequence, Is.EqualTo("ALLLLALPGL"),
                    "BaseSequence must be the canonical rotation of the original XML sequence.");

                Assert.That(entry.CircularSequence, Is.EqualTo(entry.BaseSequence),
                    "CircularSequence must equal BaseSequence for a CircularProtein.");

                Assert.That(entry.Length, Is.EqualTo(10),
                    "Circular protein length must equal the original sequence length.");

                Assert.That(entry.FullName, Is.EqualTo("Amyloid-beta precursor protein"),
                    "FullName must be read from the XML <fullName> element.");

                Assert.That(entry.Organism, Is.EqualTo("Homo sapiens"),
                    "Organism must be read from the XML <organism> element.");

                Assert.That(entry.GeneNames.Any(g => g.Item2 == "APP"), Is.True,
                    "Gene name 'APP' must be read from the XML <gene> element.");

                Assert.That(entry.DatabaseFilePath, Does.EndWith("simpleCircularProtein.xml"),
                    "DatabaseFilePath must reflect the source XML file.");

                Assert.That(entry.IsDecoy, Is.False,
                    "A standard target entry must not be flagged as a decoy.");

                Assert.That(entry.IsContaminant, Is.False,
                    "The protein was loaded as a non-contaminant.");

                // CyclicMonoisotopicMass must be positive and finite
                Assert.That(entry.CyclicMonoisotopicMass, Is.GreaterThan(0),
                    "CyclicMonoisotopicMass must be a positive value.");

                Assert.That(double.IsNaN(entry.CyclicMonoisotopicMass), Is.False,
                    "CyclicMonoisotopicMass must not be NaN.");

                Assert.That(double.IsInfinity(entry.CyclicMonoisotopicMass), Is.False,
                    "CyclicMonoisotopicMass must not be infinite.");
            });
        }

        /// <summary>
        /// Verifies that the circular protein loaded from XML can be digested
        /// without throwing, and that all resulting peptides reference the
        /// correct <see cref="CircularProtein"/> parent.
        /// </summary>
        [Test]
        public static void LoadProteinXml_CircularProtein_CanBeDigested()
        {
            string xmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "simpleCircularProtein.xml");

            List<Protein> proteins = ProteinDbLoader.LoadProteinXML(
                xmlPath,
                generateTargets: true,
                decoyType: DecoyType.None,
                allKnownModifications: new List<Modification>(),
                isContaminant: false,
                modTypesToExclude: new List<string>(),
                unknownModifications: out _);

            CircularProtein circular = CircularProtein.FromProtein(
                proteins.Single(p => p.Accession == "P05067"));

            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 2,
                minPeptideLength: 1);

            var peptides = circular
                .Digest(digestionParams, new List<Modification>(), new List<Modification>())
                .Cast<PeptideWithSetModifications>()
                .ToList();

            Assert.That(peptides, Is.Not.Empty,
                "Digestion of the XML-loaded CircularProtein must yield at least one peptide.");

            Assert.That(peptides.All(p => ReferenceEquals(p.Protein, circular)), Is.True,
                "Every digestion product must reference the CircularProtein as its parent.");
        }
    }
}