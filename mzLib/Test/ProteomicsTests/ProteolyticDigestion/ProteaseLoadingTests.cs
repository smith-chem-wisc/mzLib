using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using MzLibUtil;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class ProteaseLoadingTests
    {
        private static Stopwatch Stopwatch { get; set; }

        [SetUp]
        public static void Setup()
        {
            Stopwatch = new Stopwatch();
            Stopwatch.Start();
        }

        [TearDown]
        public static void TearDown()
        {
            Console.WriteLine($"Analysis time: {Stopwatch.Elapsed.Hours}h {Stopwatch.Elapsed.Minutes}m {Stopwatch.Elapsed.Seconds}s");
        }

        [Test]
        public static void ProteaseLoader()
        {
            string path1 = Path.Combine(TestContext.CurrentContext.TestDirectory, "ProteomicsTests", "ProteaseFilesForLoadingTests", "TestProteases_badMod.tsv");
            string path2 = Path.Combine(TestContext.CurrentContext.TestDirectory, "ProteomicsTests", "ProteaseFilesForLoadingTests", "TestProteases_badMod_dupName.tsv");
            string path3 = Path.Combine(TestContext.CurrentContext.TestDirectory, "ProteomicsTests", "ProteaseFilesForLoadingTests", "TestProteases_dupName.tsv");
            string path4 = Path.Combine(TestContext.CurrentContext.TestDirectory, "ProteomicsTests", "ProteaseFilesForLoadingTests", "TestProteases_Mod_dupName.tsv");
            var proteaseMods = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", "ProteaseMods.txt"), out var errors).ToList();

            NUnit.Framework.Assert.Throws<MzLibUtil.MzLibException>(() => ProteaseDictionary.LoadProteaseDictionary(path1, proteaseMods));
            NUnit.Framework.Assert.Throws<MzLibUtil.MzLibException>(() => ProteaseDictionary.LoadProteaseDictionary(path2, proteaseMods));
            NUnit.Framework.Assert.Throws<MzLibUtil.MzLibException>(() => ProteaseDictionary.LoadProteaseDictionary(path3, proteaseMods));
            NUnit.Framework.Assert.Throws<MzLibUtil.MzLibException>(() => ProteaseDictionary.LoadProteaseDictionary(path4, proteaseMods));
        }

        /// <summary>
        /// Tests that the embedded proteases.tsv resource exists, has the correct naming pattern,
        /// and can be loaded successfully with expected protease definitions.
        /// </summary>
        [Test]
        public static void LoadProteaseDictionary_EmbeddedResource_ExistsAndLoads()
        {
            // Verify the embedded resource exists with the expected naming pattern
            var assembly = Assembly.GetAssembly(typeof(ProteaseDictionary));
            var resourceNames = assembly.GetManifestResourceNames();
            Assert.That(resourceNames, Contains.Item("Proteomics.ProteolyticDigestion.proteases.tsv"),
                $"Expected embedded resource not found. Available resources: {string.Join(", ", resourceNames)}");

            // Verify it loads successfully and contains expected proteases via ResetToDefaults
            ProteaseDictionary.ResetToDefaults();
            var dictionary = ProteaseDictionary.Dictionary;
            Assert.That(dictionary, Is.Not.Null);
            Assert.That(dictionary.Count, Is.GreaterThan(0));

            // Verify well-known proteases exist with expected properties
            Assert.That(dictionary.ContainsKey("trypsin|P"), Is.True);
            Assert.That(dictionary["trypsin|P"].CleavageSpecificity, Is.EqualTo(CleavageSpecificity.Full));
            Assert.That(dictionary["trypsin|P"].DigestionMotifs.Count, Is.EqualTo(2)); // K[P]| and R[P]|
        }

        /// <summary>
        /// Tests backward compatibility for old-style protease names.
        /// Names like "chymotrypsin (don't cleave before proline)" should automatically
        /// resolve to "chymotrypsin|P".
        /// </summary>
        [Test]
        public static void GetProtease_OldStyleName_ResolvesToNewFormat()
        {
            // Reset to defaults to ensure clean state
            ProteaseDictionary.ResetToDefaults();

            // Test various old-style naming patterns
            var testCases = new[]
            {
                ("chymotrypsin (don't cleave before proline)", "chymotrypsin|P"),
                ("trypsin (don't cleave before proline)", "trypsin|P"),
                ("Lys-C (don't cleave before proline)", "Lys-C|P"),
            };

            foreach (var (oldName, expectedNewName) in testCases)
            {
                // Verify normalization
                string normalizedName = ProteaseDictionary.NormalizeProteaseName(oldName);
                Assert.That(normalizedName, Is.EqualTo(expectedNewName),
                    $"Failed to normalize '{oldName}' to '{expectedNewName}'");

                // Verify GetProtease works with old name
                var protease = ProteaseDictionary.GetProtease(oldName);
                Assert.That(protease, Is.Not.Null, $"GetProtease failed for '{oldName}'");
                Assert.That(protease.Name, Is.EqualTo(expectedNewName),
                    $"GetProtease returned wrong protease for '{oldName}'");

                // Verify TryGetProtease works with old name
                bool found = ProteaseDictionary.TryGetProtease(oldName, out var protease2);
                Assert.That(found, Is.True, $"TryGetProtease failed for '{oldName}'");
                Assert.That(protease2.Name, Is.EqualTo(expectedNewName));
            }
        }

        /// <summary>
        /// Tests that GetProtease still works with exact new-style names.
        /// </summary>
        [Test]
        public static void GetProtease_NewStyleName_WorksDirectly()
        {
            ProteaseDictionary.ResetToDefaults();

            var protease = ProteaseDictionary.GetProtease("trypsin|P");
            Assert.That(protease, Is.Not.Null);
            Assert.That(protease.Name, Is.EqualTo("trypsin|P"));

            bool found = ProteaseDictionary.TryGetProtease("chymotrypsin|P", out var protease2);
            Assert.That(found, Is.True);
            Assert.That(protease2.Name, Is.EqualTo("chymotrypsin|P"));
        }

        /// <summary>
        /// Tests that GetProtease throws appropriate exception for unknown protease.
        /// </summary>
        [Test]
        public static void GetProtease_UnknownProtease_ThrowsKeyNotFoundException()
        {
            ProteaseDictionary.ResetToDefaults();

            Assert.Throws<KeyNotFoundException>(() => ProteaseDictionary.GetProtease("nonexistent protease"));

            bool found = ProteaseDictionary.TryGetProtease("nonexistent protease", out var protease);
            Assert.That(found, Is.False);
            Assert.That(protease, Is.Null);
        }

        /// <summary>
        /// Tests that NormalizeProteaseName returns the original name when no pattern matches.
        /// </summary>
        [Test]
        public static void NormalizeProteaseName_NoMatch_ReturnsOriginal()
        {
            Assert.That(ProteaseDictionary.NormalizeProteaseName("trypsin"), Is.EqualTo("trypsin"));
            Assert.That(ProteaseDictionary.NormalizeProteaseName("trypsin|P"), Is.EqualTo("trypsin|P"));
            Assert.That(ProteaseDictionary.NormalizeProteaseName("custom protease"), Is.EqualTo("custom protease"));
            Assert.That(ProteaseDictionary.NormalizeProteaseName(""), Is.EqualTo(""));
            Assert.That(ProteaseDictionary.NormalizeProteaseName(null), Is.Null);
        }

        /// <summary>
        /// Tests that loading a protease definition with insufficient fields throws an appropriate exception.
        /// The parser requires at least 3 fields: Name, Motif, and Specificity.
        /// </summary>
        [Test]
        public static void LoadProteaseDictionary_InsufficientFields_ThrowsWithHelpfulMessage()
        {
            string testFile = Path.Combine(TestContext.CurrentContext.TestDirectory, "test_insufficient_fields.tsv");
            string[] lines =
            {
                "Name\tMotif\tSpecificity\tPSI-MS Accession\tPSI-MS Name\tCleavage Modification",
                "ValidProtease\tK|\tfull\t\t\t",
                "InvalidProtease\tK|"  // Missing Specificity field - only 2 columns
            };
            File.WriteAllLines(testFile, lines);

            var exception = Assert.Throws<MzLibException>(() => ProteaseDictionary.LoadProteaseDictionary(testFile));

            Assert.That(exception.Message, Does.Contain("has only 2 field(s)"));
            Assert.That(exception.Message, Does.Contain("extend to column 3"));
            Assert.That(exception.Message, Does.Contain("InvalidProtease"));

            File.Delete(testFile);
        }
        /// <summary>
        /// Tests that protease definitions with only required fields (Name, Motif, Specificity) 
        /// are parsed correctly, with optional fields defaulting to empty strings.
        /// </summary>
        [Test]
        public static void LoadProteaseDictionary_MinimalFields_DefaultsOptionalFieldsToEmpty()
        {
            string testFile = Path.Combine(TestContext.CurrentContext.TestDirectory, "test_minimal_fields.tsv");
            string[] lines =
            {
                "Name\tMotif\tSpecificity\tPSI-MS Accession\tPSI-MS Name\tCleavage Modification",
                "MinimalProtease\tK|\tfull"  // Only 3 required fields, no optional fields
            };
            File.WriteAllLines(testFile, lines);

            var dictionary = ProteaseDictionary.LoadProteaseDictionary(testFile);

            Assert.That(dictionary.ContainsKey("MinimalProtease"), Is.True);
            var protease = dictionary["MinimalProtease"];
            Assert.That(protease.Name, Is.EqualTo("MinimalProtease"));
            Assert.That(protease.CleavageSpecificity, Is.EqualTo(CleavageSpecificity.Full));
            Assert.That(protease.PsiMsAccessionNumber, Is.EqualTo(string.Empty));
            Assert.That(protease.PsiMsName, Is.EqualTo(string.Empty));
            Assert.That(protease.CleavageMod, Is.Null);

            File.Delete(testFile);
        }

        /// <summary>
        /// Tests the custom protease dictionary functionality including:
        /// - Loading custom proteases from a file and merging into the main dictionary
        /// - Overwriting existing proteases with custom definitions
        /// - Adding new proteases not in the default set
        /// - Using custom proteases for protein digestion
        /// - Resetting to default proteases
        /// 
        /// Custom protease files use TSV format with columns:
        /// Name, Motif, Specificity, PSI-MS Accession, PSI-MS Name, Cleavage Modification
        /// 
        /// Merge rules:
        /// - If protease name matches existing entry: OVERWRITES the built-in definition
        /// - If protease name is new: ADDS to the dictionary
        /// </summary>
        [Test]
        public static void LoadAndMergeCustomProteases_OverwritesAndAddsProteases()
        {
            // Arrange - capture initial state
            int initialProteaseCount = ProteaseDictionary.Dictionary.Count;
            var originalTrypsin = ProteaseDictionary.Dictionary["trypsin|P"];
            Assert.That(originalTrypsin.DigestionMotifs.Count, Is.EqualTo(2)); // K[P]| and R[P]|
            // Verify original trypsin|P cleaves after K and R (not before P)
            Assert.That(originalTrypsin.DigestionMotifs.Any(m => m.InducingCleavage == "K"), Is.True);
            Assert.That(originalTrypsin.DigestionMotifs.Any(m => m.InducingCleavage == "R"), Is.True);

            // Create a custom protease file that:
            // 1. Overrides trypsin|P with a completely different (nonsense) cleavage rule: cleave after L unless followed by P
            // 2. Adds a completely new custom protease
            string customProteaseFile = Path.Combine(TestContext.CurrentContext.TestDirectory, "test_custom_proteases.tsv");
            string[] lines =
            {
                "Name\tMotif\tSpecificity\tPSI-MS Accession\tPSI-MS Name\tCleavage Modification",
                "trypsin|P\tL[P]|\tfull\tMS:1001313\tTrypsin\t",  // Override: change from K[P]|,R[P]| to L[P]| (nonsense rule for testing)
                "MyLabProtease\tE|\tfull\t\tCustom Glu-C variant\t"  // New: cleaves after glutamate
            };
            File.WriteAllLines(customProteaseFile, lines);

            try
            {
                // Act - merge custom proteases
                var addedOrUpdated = ProteaseDictionary.LoadAndMergeCustomProteases(customProteaseFile);

                // Assert - verify merge results
                Assert.That(addedOrUpdated.Count, Is.EqualTo(2));
                Assert.That(addedOrUpdated, Contains.Item("trypsin|P"));
                Assert.That(addedOrUpdated, Contains.Item("MyLabProtease"));
                Assert.That(ProteaseDictionary.Dictionary.Count, Is.EqualTo(initialProteaseCount + 1)); // Only one new protease added (MyLabProtease); trypsin|P was overwritten

                // Verify trypsin|P was overwritten with the new L[P]| motif
                var customTrypsin = ProteaseDictionary.Dictionary["trypsin|P"];
                Assert.That(customTrypsin.DigestionMotifs.Count, Is.EqualTo(1)); // Now only L[P]|
                Assert.That(customTrypsin.DigestionMotifs[0].InducingCleavage, Is.EqualTo("L"));
                Assert.That(customTrypsin.DigestionMotifs[0].PreventingCleavage, Is.EqualTo("P"));

                // Verify new protease was added
                Assert.That(ProteaseDictionary.Dictionary.ContainsKey("MyLabProtease"), Is.True);
                var myLabProtease = ProteaseDictionary.Dictionary["MyLabProtease"];
                Assert.That(myLabProtease.DigestionMotifs.Count, Is.EqualTo(1));
                Assert.That(myLabProtease.DigestionMotifs[0].InducingCleavage, Is.EqualTo("E"));

                // Verify custom proteases work for digestion
                // Protein with L's for testing the overwritten trypsin|P
                // Note: L at position 8 is followed by 'E' (not P), so cleavage will occur there
                var protein = new Protein("PEPTIDELEPEPTIDER", null);

                // Custom trypsin|P should now cleave after L (unless followed by P)
                var customTrypsinParams = new DigestionParams(
                    protease: "trypsin|P",
                    maxMissedCleavages: 0,
                    minPeptideLength: 1);
                var customDigest = protein.Digest(customTrypsinParams, new List<Modification>(), new List<Modification>()).ToList();
                // Should cleave after L at position 8 (L is followed by E, not P), producing: PEPTIDEL, EPEPTIDER
                Assert.That(customDigest.Count, Is.EqualTo(2));
                Assert.That(customDigest[0].BaseSequence, Is.EqualTo("PEPTIDEL"));
                Assert.That(customDigest[1].BaseSequence, Is.EqualTo("EPEPTIDER"));

                // New protease should work
                var myLabParams = new DigestionParams(
                    protease: "MyLabProtease",
                    maxMissedCleavages: 0,
                    minPeptideLength: 1);
                var myLabDigest = protein.Digest(myLabParams, new List<Modification>(), new List<Modification>()).ToList();
                Assert.That(myLabDigest.Count, Is.EqualTo(6)); // PE, PTIDE, LE, PE, PTIDE, R (cleaves after each E)

                // Act - reset to defaults
                ProteaseDictionary.ResetToDefaults();

                // Assert - verify reset worked
                Assert.That(ProteaseDictionary.Dictionary.Count, Is.EqualTo(initialProteaseCount));
                Assert.That(ProteaseDictionary.Dictionary.ContainsKey("MyLabProtease"), Is.False);

                // Verify trypsin|P is back to original behavior (K[P]| and R[P]|)
                var restoredTrypsin = ProteaseDictionary.Dictionary["trypsin|P"];
                Assert.That(restoredTrypsin.DigestionMotifs.Count, Is.EqualTo(2));
                Assert.That(restoredTrypsin.DigestionMotifs.Any(m => m.InducingCleavage == "K"), Is.True);
                Assert.That(restoredTrypsin.DigestionMotifs.Any(m => m.InducingCleavage == "R"), Is.True);
            }
            finally
            {
                // Cleanup - ensure dictionary is reset even if test fails
                ProteaseDictionary.ResetToDefaults();
                if (File.Exists(customProteaseFile))
                {
                    File.Delete(customProteaseFile);
                }
            }
        }
    }
}
