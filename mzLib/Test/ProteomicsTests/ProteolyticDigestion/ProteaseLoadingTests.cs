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

            Assert.Throws<MzLibUtil.MzLibException>(() => ProteaseDictionary.LoadProteaseDictionary(path1, proteaseMods));
            Assert.Throws<MzLibUtil.MzLibException>(() => ProteaseDictionary.LoadProteaseDictionary(path2, proteaseMods));
            Assert.Throws<MzLibUtil.MzLibException>(() => ProteaseDictionary.LoadProteaseDictionary(path3, proteaseMods));
            Assert.Throws<MzLibUtil.MzLibException>(() => ProteaseDictionary.LoadProteaseDictionary(path4, proteaseMods));
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
            Assert.That(dictionary.ContainsKey("trypsin"), Is.True);
            Assert.That(dictionary["trypsin"].CleavageSpecificity, Is.EqualTo(CleavageSpecificity.Full));
            Assert.That(dictionary["trypsin"].DigestionMotifs.Count, Is.EqualTo(2)); // K[P]| and R[P]|
        }

        /// <summary>
        /// Tests backward compatibility for old-style protease names.
        /// Names like "chymotrypsin (don't cleave before proline)" should automatically
        /// resolve to "chymotrypsin".
        /// </summary>
        [Test]
        public static void GetProtease_OldStyleName_ResolvesToNewFormat()
        {
            // Reset to defaults to ensure clean state
            ProteaseDictionary.ResetToDefaults();

            // Test various old-style naming patterns
            var testCases = new[]
            {
                ("chymotrypsin (don't cleave before proline)", "chymotrypsin"),
                ("trypsin (don't cleave before proline)", "trypsin"),
                ("Lys-C (don't cleave before proline)", "Lys-C"),
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

            var protease = ProteaseDictionary.GetProtease("trypsin");
            Assert.That(protease, Is.Not.Null);
            Assert.That(protease.Name, Is.EqualTo("trypsin"));

            bool found = ProteaseDictionary.TryGetProtease("chymotrypsin", out var protease2);
            Assert.That(found, Is.True);
            Assert.That(protease2.Name, Is.EqualTo("chymotrypsin"));
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
            Assert.That(ProteaseDictionary.NormalizeProteaseName("trypsin"), Is.EqualTo("trypsin"));
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
        /// - Skipping (not overwriting) existing built-in proteases
        /// - Adding new proteases not in the default set
        /// - Using custom proteases for protein digestion
        /// - Resetting to default proteases
        /// 
        /// Custom protease files use TSV format with columns:
        /// Name, Motif, Specificity, PSI-MS Accession, PSI-MS Name, Cleavage Modification
        /// 
        /// Merge rules:
        /// - If protease name matches existing entry: SKIPS the custom definition (built-in is retained)
        /// - If protease name is new: ADDS to the dictionary
        /// </summary>
        [Test]
        public static void LoadAndMergeCustomProteases_SkipsExistingAndAddsNew()
        {
            // Arrange - capture initial state
            ProteaseDictionary.ResetToDefaults();
            int initialProteaseCount = ProteaseDictionary.Dictionary.Count;
            var originalTrypsin = ProteaseDictionary.Dictionary["trypsin"];
            Assert.That(originalTrypsin.DigestionMotifs.Count, Is.EqualTo(2)); // K[P]| and R[P]|
            Assert.That(originalTrypsin.DigestionMotifs.Any(m => m.InducingCleavage == "K"), Is.True);
            Assert.That(originalTrypsin.DigestionMotifs.Any(m => m.InducingCleavage == "R"), Is.True);

            // Create a custom protease file that:
            // 1. Attempts to override trypsin — should be SKIPPED (built-in wins)
            // 2. Adds a completely new custom protease
            string customProteaseFile = Path.Combine(TestContext.CurrentContext.TestDirectory, "test_custom_proteases.tsv");
            string[] lines =
            {
                "Name\tMotif\tSpecificity\tPSI-MS Accession\tPSI-MS Name\tCleavage Modification",
                "trypsin\tL[P]|\tfull\tMS:1001313\tTrypsin\t",  // Collision: same name as built-in, will be skipped
                "MyLabProtease\tE|\tfull\t\tCustom Glu-C variant\t"  // New: cleaves after glutamate
            };
            File.WriteAllLines(customProteaseFile, lines);

            try
            {
                // Act - merge custom proteases
                var result = ProteaseDictionary.LoadAndMergeCustomProteases(customProteaseFile);

                // Assert - verify merge results
                // trypsin collides with a built-in → Skipped; MyLabProtease is new → Added
                Assert.That(result.Added.Count, Is.EqualTo(1));
                Assert.That(result.Added, Contains.Item("MyLabProtease"));
                Assert.That(result.Skipped.Count, Is.EqualTo(1));
                Assert.That(result.Skipped, Contains.Item("trypsin"));
                Assert.That(ProteaseDictionary.Dictionary.Count, Is.EqualTo(initialProteaseCount + 1));

                // Verify trypsin was NOT overwritten — original K[P]|, R[P]| motifs intact
                var trypsinAfterMerge = ProteaseDictionary.Dictionary["trypsin"];
                Assert.That(trypsinAfterMerge.DigestionMotifs.Count, Is.EqualTo(2));
                Assert.That(trypsinAfterMerge.DigestionMotifs.Any(m => m.InducingCleavage == "K"), Is.True);
                Assert.That(trypsinAfterMerge.DigestionMotifs.Any(m => m.InducingCleavage == "R"), Is.True);

                // Verify new protease was added
                Assert.That(ProteaseDictionary.Dictionary.ContainsKey("MyLabProtease"), Is.True);
                var myLabProtease = ProteaseDictionary.Dictionary["MyLabProtease"];
                Assert.That(myLabProtease.DigestionMotifs.Count, Is.EqualTo(1));
                Assert.That(myLabProtease.DigestionMotifs[0].InducingCleavage, Is.EqualTo("E"));

                // Verify the new protease works for digestion (cleaves after each E)
                var protein = new Protein("PEPTIDELEPEPTIDER", null);
                var myLabParams = new DigestionParams(
                    protease: "MyLabProtease",
                    maxMissedCleavages: 0,
                    minPeptideLength: 1);
                var myLabDigest = protein.Digest(myLabParams, new List<Modification>(), new List<Modification>()).ToList();
                Assert.That(myLabDigest.Count, Is.EqualTo(6)); // PE, PTIDE, LE, PE, PTIDE, R

                // Act - reset to defaults
                ProteaseDictionary.ResetToDefaults();

                // Assert - verify reset worked
                Assert.That(ProteaseDictionary.Dictionary.Count, Is.EqualTo(initialProteaseCount));
                Assert.That(ProteaseDictionary.Dictionary.ContainsKey("MyLabProtease"), Is.False);

                // Verify trypsin is still the built-in definition
                var restoredTrypsin = ProteaseDictionary.Dictionary["trypsin"];
                Assert.That(restoredTrypsin.DigestionMotifs.Count, Is.EqualTo(2));
                Assert.That(restoredTrypsin.DigestionMotifs.Any(m => m.InducingCleavage == "K"), Is.True);
                Assert.That(restoredTrypsin.DigestionMotifs.Any(m => m.InducingCleavage == "R"), Is.True);
            }
            finally
            {
                ProteaseDictionary.ResetToDefaults();
                if (File.Exists(customProteaseFile))
                    File.Delete(customProteaseFile);
            }
        }

        #region Proline naming convention

        /// <summary>
        /// The bare name is the proline-restricted enzyme and the "/P" suffix is the variant that
        /// ignores the restriction, matching MaxQuant, Mascot, Comet and the PSI-MS cleavage agent
        /// terms of the same spelling. Between January 2026 and the restoration of this convention
        /// the file said the reverse, so this pins the direction rather than merely the motifs.
        /// </summary>
        [Test]
        public static void ProlineSuffixFollowsTheCommunityConvention()
        {
            ProteaseDictionary.ResetToDefaults();

            var restricted = ProteaseDictionary.GetProtease("trypsin");
            var unrestricted = ProteaseDictionary.GetProtease("trypsin/P");

            Assert.That(restricted.DigestionMotifs.All(m => m.PreventingCleavage == "P"), Is.True,
                "bare 'trypsin' must carry the proline restriction (K[P]|,R[P]|)");
            Assert.That(unrestricted.DigestionMotifs.All(m => string.IsNullOrEmpty(m.PreventingCleavage)), Is.True,
                "'trypsin/P' must not carry the proline restriction (K|,R|)");
        }

        /// <summary>
        /// The convention is only real if it changes digestion. KPEPTIDER has a proline
        /// immediately after K1, so the restricted enzyme must read through it.
        /// </summary>
        [Test]
        public static void BareNameDoesNotCleaveBeforeProlineButSlashPDoes()
        {
            ProteaseDictionary.ResetToDefaults();
            var protein = new Protein("KPEPTIDER", "PROLINE_TEST");
            var noMods = new List<Modification>();

            var restricted = protein
                .Digest(new DigestionParams("trypsin", maxMissedCleavages: 0, minPeptideLength: 1), noMods, noMods)
                .Select(p => p.BaseSequence).ToList();
            var unrestricted = protein
                .Digest(new DigestionParams("trypsin/P", maxMissedCleavages: 0, minPeptideLength: 1), noMods, noMods)
                .Select(p => p.BaseSequence).ToList();

            Assert.That(restricted, Is.EquivalentTo(new[] { "KPEPTIDER" }),
                "trypsin must not cut K1|P2, leaving the sequence intact");
            Assert.That(unrestricted, Is.EquivalentTo(new[] { "K", "PEPTIDER" }),
                "trypsin/P must cut K1|P2");
        }

        /// <summary>
        /// Every proline-restricted enzyme drops the suffix. Nothing in the shipped file may
        /// reintroduce a "|P" name, because "|P" now means the opposite of what it used to.
        /// </summary>
        [Test]
        public static void NoShippedProteaseUsesThePipePSuffix()
        {
            ProteaseDictionary.ResetToDefaults();

            var offenders = ProteaseDictionary.Dictionary.Keys
                .Where(k => k.EndsWith("|P", StringComparison.OrdinalIgnoreCase)).ToList();

            Assert.That(offenders, Is.Empty,
                $"'|P' named the restricted enzyme and now means the opposite; found: {string.Join(", ", offenders)}");
        }

        /// <summary>
        /// Settings saved under either historical spelling must resolve to the same enzyme they
        /// named when they were written, so no saved toml silently changes meaning.
        /// </summary>
        [Test]
        [TestCase("trypsin|P", "trypsin")]
        [TestCase("chymotrypsin|P", "chymotrypsin")]
        [TestCase("Lys-C|P", "Lys-C")]
        [TestCase("elastase|P", "elastase")]
        [TestCase("subtilisin|P", "subtilisin")]
        [TestCase("trypsin (don't cleave before proline)", "trypsin")]
        [TestCase("chymotrypsin (don't cleave before proline)", "chymotrypsin")]
        [TestCase("Lys-C (don't cleave before proline)", "Lys-C")]
        [TestCase("trypsin (cleave before proline)", "trypsin/P")]
        [TestCase("chymotrypsin (cleave before proline)", "chymotrypsin/P")]
        [TestCase("Lys-C (cleave before proline)", "Lys-C/P")]
        public static void HistoricalNamesResolveToTheSameEnzyme(string legacyName, string currentName)
        {
            ProteaseDictionary.ResetToDefaults();

            Assert.That(ProteaseDictionary.NormalizeProteaseName(legacyName), Is.EqualTo(currentName));

            Assert.That(ProteaseDictionary.TryGetProtease(legacyName, out var viaLegacy), Is.True,
                $"'{legacyName}' no longer resolves; settings saved under it would fail to load");
            Assert.That(viaLegacy.Name, Is.EqualTo(currentName));

            // DigestionParams is the path a toml actually travels; it must accept the legacy name too.
            Assert.That(new DigestionParams(legacyName).Protease.Name, Is.EqualTo(currentName));
        }

        /// <summary>
        /// Bare "trypsin" is the one name whose meaning changed, and it cannot be detected at lookup
        /// time because both spellings are live keys. Pinned so the trade-off stays deliberate.
        /// </summary>
        [Test]
        public static void BareTrypsinIsNotRemapped()
        {
            Assert.That(ProteaseDictionary.NormalizeProteaseName("trypsin"), Is.EqualTo("trypsin"));
            Assert.That(ProteaseDictionary.NormalizeProteaseName("trypsin/P"), Is.EqualTo("trypsin/P"));
        }


        /// <summary>
        /// mzLib #1005 dropped the unrestricted chymotrypsin and Lys-C variants when it moved to the
        /// "|P" naming, leaving "chymotrypsin (cleave before proline)" and "Lys-C (cleave before
        /// proline)" -- names that really did ship -- with nothing to resolve to. They are restored
        /// here under the "/P" spelling, which is what makes the parenthetical mapping above real.
        /// </summary>
        [Test]
        [TestCase("chymotrypsin/P", "F", "W", "Y", "L")]
        [TestCase("Lys-C/P", "K")]
        public static void RestoredUnrestrictedVariantsIgnoreProline(string name, params string[] residues)
        {
            ProteaseDictionary.ResetToDefaults();

            var protease = ProteaseDictionary.GetProtease(name);
            Assert.That(protease.DigestionMotifs.Select(m => m.InducingCleavage), Is.EquivalentTo(residues));
            Assert.That(protease.DigestionMotifs.All(m => string.IsNullOrEmpty(m.PreventingCleavage)), Is.True,
                $"'{name}' must not carry the proline restriction");

            // And the restricted sibling still does, so the pair is a real pair.
            var restricted = ProteaseDictionary.GetProtease(name.Substring(0, name.Length - 2));
            Assert.That(restricted.DigestionMotifs.All(m => m.PreventingCleavage == "P"), Is.True);
        }

        #endregion
    }
}
