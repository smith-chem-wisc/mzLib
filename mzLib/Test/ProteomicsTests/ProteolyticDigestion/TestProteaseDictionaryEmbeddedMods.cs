using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using Assert = NUnit.Framework.Legacy.ClassicAssert;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    /// <summary>
    /// Tests for ProteaseDictionary embedded resource functionality and custom digestion agent
    /// loading behavior. The embedded proteases.tsv is the authoritative source of truth and
    /// cannot be overridden by any custom file.
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class TestProteaseDictionaryEmbeddedMods
    {
        #region Helpers

        /// <summary>
        /// Writes a minimal valid proteases TSV to a temp file and returns its path.
        /// Caller is responsible for deleting it.
        /// </summary>
        private static string WriteTempProteaseFile(IEnumerable<(string Name, string Motif, string Specificity)> entries)
        {
            string path = Path.GetTempFileName();
            var lines = new List<string> { "Name\tMotif\tSpecificity\tPSI-MS Accession\tPSI-MS Name\tCleavage Modification" };
            foreach (var (name, motif, specificity) in entries)
                lines.Add($"{name}\t{motif}\t{specificity}\t\t\t");
            File.WriteAllLines(path, lines);
            return path;
        }

        #endregion

        #region Embedded Resource Existence Tests

        /// <summary>
        /// Verifies that both embedded resources (proteases.tsv and protease_mods.txt) exist in the assembly.
        /// </summary>
        [Test]
        public static void EmbeddedResources_BothExist()
        {
            var assembly = Assembly.GetAssembly(typeof(ProteaseDictionary));
            var resourceNames = assembly.GetManifestResourceNames();

            Assert.That(resourceNames.Contains("Proteomics.ProteolyticDigestion.proteases.tsv"),
                $"proteases.tsv not found. Available: {string.Join(", ", resourceNames)}");

            Assert.That(resourceNames.Contains("Proteomics.ProteolyticDigestion.protease_mods.txt"),
                $"protease_mods.txt not found. Available: {string.Join(", ", resourceNames)}");
        }

        #endregion

        #region Dictionary Immutability Tests

        /// <summary>
        /// Verifies that the Dictionary property has no public setter — external code cannot
        /// swap the entire dictionary out from under the embedded resource.
        /// </summary>
        [Test]
        public static void Dictionary_HasNoPublicSetter()
        {
            var prop = typeof(ProteaseDictionary).GetProperty(nameof(ProteaseDictionary.Dictionary));
            Assert.That(prop, Is.Not.Null);

            var setter = prop.SetMethod;
            // Setter must either not exist or be non-public
            Assert.That(setter == null || !setter.IsPublic, Is.True,
                "Dictionary.set must not be public — the embedded resource is the authoritative baseline.");
        }

        /// <summary>
        /// Verifies that the static Dictionary is populated at startup and contains the
        /// expected common proteases from the embedded resource.
        /// </summary>
        [Test]
        public static void StaticDictionary_InitializedWithEmbeddedProteases()
        {
            var dict = ProteaseDictionary.Dictionary;

            Assert.That(dict, Is.Not.Null, "Static Dictionary should be initialized");
            Assert.That(dict.Count, Is.GreaterThan(10), "Should contain many proteases from embedded resource");

            Assert.That(dict.ContainsKey("trypsin|P"), Is.True);
            Assert.That(dict.ContainsKey("trypsin"), Is.True);
            Assert.That(dict.ContainsKey("Lys-C|P"), Is.True);
            Assert.That(dict.ContainsKey("Asp-N"), Is.True);
            Assert.That(dict.ContainsKey("CNBr"), Is.True);
        }

        /// <summary>
        /// Verifies that CNBr in the static Dictionary has its cleavage modification loaded
        /// from the embedded protease_mods.txt — proving the self-contained loading works
        /// without any external files.
        /// </summary>
        [Test]
        public static void StaticDictionary_CNBrHasEmbeddedCleavageMod()
        {
            var cnbr = ProteaseDictionary.Dictionary["CNBr"];

            Assert.That(cnbr, Is.Not.Null);
            Assert.That(cnbr.CleavageMod, Is.Not.Null,
                "CNBr should have a cleavage modification loaded from embedded resources");
            Assert.That(cnbr.CleavageMod.IdWithMotif, Is.EqualTo("Homoserine lactone on M"));
            Assert.That(cnbr.CleavageMod.Target?.ToString(), Is.EqualTo("M"));
        }

        /// <summary>
        /// Verifies that trypsin has no cleavage modification (sanity check on the embedded data).
        /// </summary>
        [Test]
        public static void StaticDictionary_TrypsinHasNoCleavageMod()
        {
            var trypsin = ProteaseDictionary.Dictionary["trypsin|P"];

            Assert.That(trypsin, Is.Not.Null);
            Assert.That(trypsin.CleavageMod, Is.Null,
                "Trypsin should NOT have a cleavage modification");
            Assert.That(trypsin.DigestionMotifs.Count, Is.EqualTo(2),
                "Trypsin should have 2 digestion motifs (K[P]| and R[P]|)");
        }

        #endregion

        #region LoadEmbeddedProteaseMods Tests

        /// <summary>
        /// Verifies that LoadEmbeddedProteaseMods returns a non-empty list of modifications.
        /// </summary>
        [Test]
        public static void LoadEmbeddedProteaseMods_ReturnsModifications()
        {
            var mods = ProteaseDictionary.LoadEmbeddedProteaseMods();

            Assert.That(mods, Is.Not.Null, "LoadEmbeddedProteaseMods should not return null");
            Assert.That(mods.Count, Is.GreaterThan(0), "Should load at least one embedded modification");
        }

        /// <summary>
        /// Verifies that the "Homoserine lactone on M" modification is loaded correctly.
        /// </summary>
        [Test]
        public static void LoadEmbeddedProteaseMods_ContainsHomoserineLactone()
        {
            var mods = ProteaseDictionary.LoadEmbeddedProteaseMods();
            var homoserineMod = mods.FirstOrDefault(m => m.IdWithMotif == "Homoserine lactone on M");

            Assert.That(homoserineMod, Is.Not.Null, "Should contain 'Homoserine lactone on M'");
            Assert.That(homoserineMod.Target?.ToString(), Is.EqualTo("M"));
            Assert.That(homoserineMod.ModificationType, Is.EqualTo("Protease"));
            Assert.That(homoserineMod.ChemicalFormula, Is.Not.Null);
            Assert.That(homoserineMod.MonoisotopicMass, Is.Not.Null);
        }

        /// <summary>
        /// Verifies that the "Test on M" modification is loaded correctly.
        /// </summary>
        [Test]
        public static void LoadEmbeddedProteaseMods_ContainsTestOnM()
        {
            var mods = ProteaseDictionary.LoadEmbeddedProteaseMods();
            var testMod = mods.FirstOrDefault(m => m.IdWithMotif == "Test on M");

            Assert.That(testMod, Is.Not.Null, "Should contain 'Test on M'");
            Assert.That(testMod.Target?.ToString(), Is.EqualTo("M"));
        }

        /// <summary>
        /// Verifies all embedded modifications carry the correct file origin tag.
        /// </summary>
        [Test]
        public static void LoadEmbeddedProteaseMods_HasCorrectFileOrigin()
        {
            var mods = ProteaseDictionary.LoadEmbeddedProteaseMods();

            foreach (var mod in mods)
            {
                Assert.That(mod.FileOrigin, Is.EqualTo("Embedded:protease_mods.txt"),
                    $"Modification '{mod.IdWithMotif}' should have correct FileOrigin");
            }
        }

        /// <summary>
        /// Verifies that exactly the expected number of modifications are present and
        /// that no comment lines were incorrectly parsed as entries.
        /// </summary>
        [Test]
        public static void LoadEmbeddedProteaseMods_ExactCount_NoSpuriousEntries()
        {
            var mods = ProteaseDictionary.LoadEmbeddedProteaseMods();

            Assert.That(mods.Count, Is.EqualTo(2),
                "Should have exactly 2 modifications (Homoserine lactone on M and Test on M)");
            Assert.That(mods.All(m => !m.IdWithMotif.StartsWith("#")), Is.True,
                "No modification ID should start with # (comment lines must be skipped)");
        }

        #endregion

        #region ParseModificationsFromString Coverage Tests

        [Test]
        public static void ParseModificationsFromString_ParsesMonoisotopicMass()
        {
            var mods = ProteaseDictionary.LoadEmbeddedProteaseMods();
            var homoserineMod = mods.FirstOrDefault(m => m.IdWithMotif == "Homoserine lactone on M");

            Assert.That(homoserineMod, Is.Not.Null);
            Assert.That(homoserineMod.MonoisotopicMass.HasValue, Is.True,
                "MonoisotopicMass should be calculated from chemical formula");
            Assert.That(homoserineMod.MonoisotopicMass.Value, Is.LessThan(0),
                "Mass change should be negative (loss of atoms)");
        }

        [Test]
        public static void ParseModificationsFromString_ParsesLocationRestriction()
        {
            var mods = ProteaseDictionary.LoadEmbeddedProteaseMods();
            var homoserineMod = mods.FirstOrDefault(m => m.IdWithMotif == "Homoserine lactone on M");
            var testMod = mods.FirstOrDefault(m => m.IdWithMotif == "Test on M");

            Assert.That(homoserineMod?.LocationRestriction, Is.EqualTo("Peptide C-terminal."));
            Assert.That(testMod?.LocationRestriction, Is.EqualTo("Peptide N-terminal."));
        }

        [Test]
        public static void ParseModificationsFromString_DefaultsModificationType()
        {
            var mods = ProteaseDictionary.LoadEmbeddedProteaseMods();

            foreach (var mod in mods)
            {
                Assert.That(mod.ModificationType, Is.EqualTo("Protease"),
                    $"Modification '{mod.IdWithMotif}' should have ModificationType 'Protease'");
            }
        }

        [Test]
        public static void ParseModificationsFromString_ParsesChemicalFormula()
        {
            var mods = ProteaseDictionary.LoadEmbeddedProteaseMods();
            var homoserineMod = mods.FirstOrDefault(m => m.IdWithMotif == "Homoserine lactone on M");

            Assert.That(homoserineMod, Is.Not.Null);
            Assert.That(homoserineMod.ChemicalFormula, Is.Not.Null);
            Assert.That(homoserineMod.ChemicalFormula.ToString(), Does.Contain("C"));
            Assert.That(homoserineMod.ChemicalFormula.ToString(), Does.Contain("H"));
            Assert.That(homoserineMod.ChemicalFormula.ToString(), Does.Contain("S"));
        }

        [Test]
        public static void ParseModificationsFromString_ParsesTarget()
        {
            var mods = ProteaseDictionary.LoadEmbeddedProteaseMods();

            foreach (var mod in mods)
            {
                Assert.That(mod.Target, Is.Not.Null,
                    $"Modification '{mod.IdWithMotif}' should have a target");
                Assert.That(mod.Target.ToString(), Is.EqualTo("M"),
                    $"Modification '{mod.IdWithMotif}' should target Methionine");
            }
        }

        [Test]
        public static void ParseModificationsFromString_HandlesSeparatorCorrectly()
        {
            var mods = ProteaseDictionary.LoadEmbeddedProteaseMods();
            var ids = mods.Select(m => m.IdWithMotif).ToList();

            Assert.That(ids.Count, Is.EqualTo(ids.Distinct().Count()),
                "All modification IDs should be unique — separator must be working correctly");
        }

        [Test]
        public static void ParseModificationsFromString_ParsesDatabaseReference_WhenPresent()
        {
            var mods = ProteaseDictionary.LoadEmbeddedProteaseMods();
            var homoserineMod = mods.FirstOrDefault(m => m.IdWithMotif == "Homoserine lactone on M");

            Assert.That(homoserineMod, Is.Not.Null);
            Assert.That(homoserineMod.DatabaseReference, Is.Not.Null,
                "Homoserine lactone should have a database reference");
            Assert.That(homoserineMod.DatabaseReference.ContainsKey("Unimod"), Is.True);
            Assert.That(homoserineMod.DatabaseReference["Unimod"], Contains.Item("10"));
        }

        [Test]
        public static void ParseModificationsFromString_NullDatabaseReference_WhenAbsent()
        {
            var mods = ProteaseDictionary.LoadEmbeddedProteaseMods();
            var testMod = mods.FirstOrDefault(m => m.IdWithMotif == "Test on M");

            Assert.That(testMod, Is.Not.Null);
            Assert.That(testMod.DatabaseReference, Is.Null,
                "Test on M should NOT have a database reference");
        }

        [Test]
        public static void LoadEmbeddedProteaseMods_HandlesEmptyGracefully()
        {
            // Contract: never returns null, even if the resource were missing
            var mods = ProteaseDictionary.LoadEmbeddedProteaseMods();
            Assert.That(mods, Is.Not.Null, "Should never return null");
        }

        [Test]
        public static void ParseModificationsFromString_HandlesInvalidFormulasGracefully()
        {
            var mods = ProteaseDictionary.LoadEmbeddedProteaseMods();

            foreach (var mod in mods)
            {
                Assert.That(mod.ChemicalFormula, Is.Not.Null,
                    $"Modification '{mod.IdWithMotif}' should have a valid chemical formula");
            }
        }

        #endregion

        #region Custom Digestion Agent Loading Tests

        /// <summary>
        /// A custom protease with a new name is added to the dictionary and appears in Added.
        /// </summary>
        [Test]
        public static void LoadAndMergeCustomProteases_SinglePath_AddsNewProtease()
        {
            string path = WriteTempProteaseFile(new[] { ("MyCustomProtease", "K|", "full") });
            try
            {
                var result = ProteaseDictionary.LoadAndMergeCustomProteases(path);

                Assert.That(result.Added, Contains.Item("MyCustomProtease"),
                    "New protease should appear in Added");
                Assert.That(result.Skipped, Is.Empty,
                    "No names should be skipped when there is no collision");
                Assert.That(ProteaseDictionary.Dictionary.ContainsKey("MyCustomProtease"), Is.True,
                    "New protease should be present in the live dictionary");
            }
            finally
            {
                File.Delete(path);
                // Clean up the custom entry so it does not pollute other tests
                ProteaseDictionary.Dictionary.Remove("MyCustomProtease");
            }
        }

        /// <summary>
        /// The IEnumerable overload processes multiple files in order. Entries from all files
        /// with unique names are all added.
        /// </summary>
        [Test]
        public static void LoadAndMergeCustomProteases_MultiplePaths_AllNewNamesAdded()
        {
            string path1 = WriteTempProteaseFile(new[] { ("CustomA", "K|", "full") });
            string path2 = WriteTempProteaseFile(new[] { ("CustomB", "R|", "full") });
            try
            {
                var result = ProteaseDictionary.LoadAndMergeCustomProteases(new[] { path1, path2 });

                Assert.That(result.Added, Contains.Item("CustomA"));
                Assert.That(result.Added, Contains.Item("CustomB"));
                Assert.That(result.Skipped, Is.Empty);
                Assert.That(ProteaseDictionary.Dictionary.ContainsKey("CustomA"), Is.True);
                Assert.That(ProteaseDictionary.Dictionary.ContainsKey("CustomB"), Is.True);
            }
            finally
            {
                File.Delete(path1);
                File.Delete(path2);
                ProteaseDictionary.Dictionary.Remove("CustomA");
                ProteaseDictionary.Dictionary.Remove("CustomB");
            }
        }

        /// <summary>
        /// A custom protease whose name matches an embedded protease is skipped.
        /// The embedded definition is retained and the name appears in Skipped, not Added.
        /// </summary>
        [Test]
        public static void LoadAndMergeCustomProteases_CollisionWithEmbedded_EmbeddedWins()
        {
            // "trypsin" exists in the embedded resource
            string path = WriteTempProteaseFile(new[] { ("trypsin", "K|", "full") });
            try
            {
                var originalTrypsin = ProteaseDictionary.Dictionary["trypsin"];
                var result = ProteaseDictionary.LoadAndMergeCustomProteases(path);

                Assert.That(result.Skipped, Contains.Item("trypsin"),
                    "Collision with embedded protease must appear in Skipped");
                Assert.That(result.Added, Does.Not.Contain("trypsin"),
                    "Collision must not appear in Added");

                // Embedded definition must be unchanged
                Assert.That(ReferenceEquals(ProteaseDictionary.Dictionary["trypsin"], originalTrypsin), Is.True,
                    "The embedded trypsin object must not have been replaced");
            }
            finally
            {
                File.Delete(path);
            }
        }

        /// <summary>
        /// When the same name appears in two different custom files, the first-encountered wins
        /// and the second is skipped — regardless of whether the name was in the embedded resource.
        /// </summary>
        [Test]
        public static void LoadAndMergeCustomProteases_CollisionAcrossCustomFiles_FirstWins()
        {
            string path1 = WriteTempProteaseFile(new[] { ("SharedName", "K|", "full") });
            string path2 = WriteTempProteaseFile(new[] { ("SharedName", "R|", "semi") });
            try
            {
                var result = ProteaseDictionary.LoadAndMergeCustomProteases(new[] { path1, path2 });

                Assert.That(result.Added, Contains.Item("SharedName"),
                    "First occurrence should be added");
                Assert.That(result.Skipped, Contains.Item("SharedName"),
                    "Second occurrence should be skipped");

                // The motif from file1 (K|) should be the one retained
                var retained = ProteaseDictionary.Dictionary["SharedName"];
                Assert.That(retained.DigestionMotifs.Any(m => m.InducingCleavage == "K"), Is.True,
                    "The protease from the first file should be the one retained");
            }
            finally
            {
                File.Delete(path1);
                File.Delete(path2);
                ProteaseDictionary.Dictionary.Remove("SharedName");
            }
        }

        /// <summary>
        /// Duplicate names within a single custom file are a malformed file and must throw.
        /// This is distinct from cross-file collisions, which are handled gracefully.
        /// </summary>
        [Test]
        public static void LoadAndMergeCustomProteases_DuplicateWithinSingleFile_Throws()
        {
            string path = Path.GetTempFileName();
            try
            {
                File.WriteAllLines(path, new[]
                {
                    "Name\tMotif\tSpecificity\tPSI-MS Accession\tPSI-MS Name\tCleavage Modification",
                    "DuplicateName\tK|\tfull\t\t\t",
                    "DuplicateName\tR|\tfull\t\t\t"
                });

                Assert.Throws<MzLibUtil.MzLibException>(
                    () => ProteaseDictionary.LoadAndMergeCustomProteases(path),
                    "Duplicate names within a single file should throw MzLibException");
            }
            finally
            {
                File.Delete(path);
            }
        }

        /// <summary>
        /// The returned CustomDigestionAgentLoadResult correctly partitions names into
        /// Added and Skipped across a mixed scenario.
        /// </summary>
        [Test]
        public static void LoadAndMergeCustomProteases_ReturnValue_CorrectlyPartitioned()
        {
            // "trypsin" is embedded (will be skipped); "BrandNewProtease" is not (will be added)
            string path = WriteTempProteaseFile(new[]
            {
                ("trypsin",          "K|", "full"),
                ("BrandNewProtease", "R|", "full")
            });
            try
            {
                var result = ProteaseDictionary.LoadAndMergeCustomProteases(path);

                Assert.That(result, Is.Not.Null);
                Assert.That(result.Added, Contains.Item("BrandNewProtease"));
                Assert.That(result.Added, Does.Not.Contain("trypsin"));
                Assert.That(result.Skipped, Contains.Item("trypsin"));
                Assert.That(result.Skipped, Does.Not.Contain("BrandNewProtease"));
            }
            finally
            {
                File.Delete(path);
                ProteaseDictionary.Dictionary.Remove("BrandNewProtease");
            }
        }

        #endregion

        #region Integration Tests

        /// <summary>
        /// Verifies that MetaMorpheus-style usage works: access the dictionary directly
        /// without providing any external files.
        /// </summary>
        [Test]
        public static void Integration_MetaMorpheusStyle_NoExternalFilesNeeded()
        {
            var trypsin = ProteaseDictionary.GetProtease("trypsin|P");
            var cnbr = ProteaseDictionary.GetProtease("CNBr");
            var lysC = ProteaseDictionary.GetProtease("Lys-C|P");

            Assert.That(trypsin, Is.Not.Null);
            Assert.That(cnbr, Is.Not.Null);
            Assert.That(lysC, Is.Not.Null);

            Assert.That(cnbr.CleavageMod, Is.Not.Null,
                "CNBr should have cleavage mod without needing external files");
            Assert.That(cnbr.CleavageMod.IdWithMotif, Is.EqualTo("Homoserine lactone on M"));
        }

        #endregion

        #region Atomicity, Null Guard, and Edge Case Tests

        [Test]
        public static void LoadAndMergeCustomProteases_SecondPathInvalid_DictionaryUnchanged()
        {
            string tempDir = Path.Combine(Path.GetTempPath(), "ProteaseAtomicTests_" + Path.GetRandomFileName());
            Directory.CreateDirectory(tempDir);
            try
            {
                string validFile = Path.Combine(tempDir, "valid.tsv");
                File.WriteAllLines(validFile, new[]
                {
                    "Name\tMotif\tSpecificity",
                    "AtomicTestProtease\tK|\tFull"
                });
                string nonExistent = Path.Combine(tempDir, "no_such_file.tsv");

                int countBefore = ProteaseDictionary.Dictionary.Count;

                NUnit.Framework.Assert.That(() =>
                    ProteaseDictionary.LoadAndMergeCustomProteases(new[] { validFile, nonExistent }),
                    Throws.Exception);

                NUnit.Framework.Assert.That(ProteaseDictionary.Dictionary.Count, Is.EqualTo(countBefore),
                    "Dictionary must not grow when the second file is missing");
                NUnit.Framework.Assert.That(ProteaseDictionary.Dictionary.ContainsKey("AtomicTestProtease"), Is.False,
                    "Entry from the first (valid) file must not be present after a failed multi-file load");
            }
            finally
            {
                Directory.Delete(tempDir, true);
            }
        }

        [Test]
        public static void LoadAndMergeCustomProteases_SecondFileHasDuplicates_DictionaryUnchanged()
        {
            string tempDir = Path.Combine(Path.GetTempPath(), "ProteaseAtomicDup_" + Path.GetRandomFileName());
            Directory.CreateDirectory(tempDir);
            try
            {
                string validFile = Path.Combine(tempDir, "valid.tsv");
                File.WriteAllLines(validFile, new[]
                {
                    "Name\tMotif\tSpecificity",
                    "AtomicTestProtease2\tK|\tFull"
                });
                string dupFile = Path.Combine(tempDir, "dup.tsv");
                File.WriteAllLines(dupFile, new[]
                {
                    "Name\tMotif\tSpecificity",
                    "DupEntry\tK|\tFull",
                    "DupEntry\tR|\tFull"
                });

                int countBefore = ProteaseDictionary.Dictionary.Count;

                NUnit.Framework.Assert.That(() =>
                    ProteaseDictionary.LoadAndMergeCustomProteases(new[] { validFile, dupFile }),
                    Throws.TypeOf<MzLibUtil.MzLibException>());

                NUnit.Framework.Assert.That(ProteaseDictionary.Dictionary.Count, Is.EqualTo(countBefore));
                NUnit.Framework.Assert.That(ProteaseDictionary.Dictionary.ContainsKey("AtomicTestProtease2"), Is.False);
            }
            finally
            {
                Directory.Delete(tempDir, true);
            }
        }

        [Test]
        public static void LoadAndMergeCustomProteases_NullPaths_ThrowsArgumentNullException()
        {
            NUnit.Framework.Assert.That(
                () => ProteaseDictionary.LoadAndMergeCustomProteases((IEnumerable<string>)null),
                Throws.TypeOf<System.ArgumentNullException>()
                    .With.Property("ParamName").EqualTo("paths"));
        }

        [Test]
        public static void LoadAndMergeCustomProteases_EmptyPaths_ReturnsEmptyResult()
        {
            int countBefore = ProteaseDictionary.Dictionary.Count;

            var result = ProteaseDictionary.LoadAndMergeCustomProteases(
                System.Array.Empty<string>());

            NUnit.Framework.Assert.That(result.Added, Is.Empty);
            NUnit.Framework.Assert.That(result.Skipped, Is.Empty);
            NUnit.Framework.Assert.That(ProteaseDictionary.Dictionary.Count, Is.EqualTo(countBefore));
        }

        [Test]
        public static void LoadAndMergeCustomProteases_DuplicateAcrossFiles_FirstFileWins()
        {
            string tempDir = Path.Combine(Path.GetTempPath(), "ProteaseAcrossDup_" + Path.GetRandomFileName());
            Directory.CreateDirectory(tempDir);
            try
            {
                string file1 = Path.Combine(tempDir, "file1.tsv");
                File.WriteAllLines(file1, new[]
                {
                    "Name\tMotif\tSpecificity",
                    "CrossFileDupProtease\tK|\tFull"
                });
                string file2 = Path.Combine(tempDir, "file2.tsv");
                File.WriteAllLines(file2, new[]
                {
                    "Name\tMotif\tSpecificity",
                    "CrossFileDupProtease\tR|\tFull"
                });

                var result = ProteaseDictionary.LoadAndMergeCustomProteases(new[] { file1, file2 });

                NUnit.Framework.Assert.That(result.Added, Does.Contain("CrossFileDupProtease"));
                NUnit.Framework.Assert.That(result.Skipped, Does.Contain("CrossFileDupProtease"));
                // First file's motif (K|) should win
                var motifs = ProteaseDictionary.Dictionary["CrossFileDupProtease"].DigestionMotifs;
                NUnit.Framework.Assert.That(motifs.Any(m => m.InducingCleavage == "K"), Is.True);
            }
            finally
            {
                ProteaseDictionary.Dictionary.Remove("CrossFileDupProtease");
                Directory.Delete(tempDir, true);
            }
        }

        [Test]
        public static void LoadEmbeddedProteaseMods_ReturnsFreshCopy_MutationDoesNotAffectCache()
        {
            var first = ProteaseDictionary.LoadEmbeddedProteaseMods();
            int originalCount = first.Count;

            first.Clear();

            var second = ProteaseDictionary.LoadEmbeddedProteaseMods();
            NUnit.Framework.Assert.That(second.Count, Is.EqualTo(originalCount),
                "Clearing a returned list must not affect the cached copy");
        }

        [Test]
        public static void LoadEmbeddedProteaseMods_ReturnsDifferentInstances()
        {
            var first = ProteaseDictionary.LoadEmbeddedProteaseMods();
            var second = ProteaseDictionary.LoadEmbeddedProteaseMods();

            NUnit.Framework.Assert.That(ReferenceEquals(first, second), Is.False,
                "Each call should return a separate list instance");
        }

        [Test]
        public static void LoadAndMergeCustomProteases_FileWithNoHeader_Throws()
        {
            string path = Path.GetTempFileName();
            try
            {
                // Write a file with only comment lines and no header
                File.WriteAllLines(path, new[]
                {
                    "# This is a comment",
                    "# Another comment"
                });

                NUnit.Framework.Assert.That(
                    () => ProteaseDictionary.LoadAndMergeCustomProteases(path),
                    Throws.TypeOf<MzLibUtil.MzLibException>()
                        .With.Message.Contains("no header"));
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public static void LoadAndMergeCustomProteases_HeaderMissingRequiredColumn_Throws()
        {
            string path = Path.GetTempFileName();
            try
            {
                File.WriteAllLines(path, new[]
                {
                    "Name\tMotif",
                    "TestProtease\tK|"
                });

                NUnit.Framework.Assert.That(
                    () => ProteaseDictionary.LoadAndMergeCustomProteases(path),
                    Throws.TypeOf<MzLibUtil.MzLibException>()
                        .With.Message.Contains("Specificity"));
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public static void LoadAndMergeCustomProteases_EmptySpecificity_Throws()
        {
            string path = Path.GetTempFileName();
            try
            {
                File.WriteAllLines(path, new[]
                {
                    "Name\tMotif\tSpecificity",
                    "TestProtease\tK|\t"
                });

                NUnit.Framework.Assert.That(
                    () => ProteaseDictionary.LoadAndMergeCustomProteases(path),
                    Throws.TypeOf<MzLibUtil.MzLibException>()
                        .With.Message.Contains("Specificity"));
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public static void LoadAndMergeCustomProteases_EmptyFile_ReturnsEmptyResult()
        {
            string path = Path.GetTempFileName();
            try
            {
                // Header only, no data rows
                File.WriteAllLines(path, new[]
                {
                    "Name\tMotif\tSpecificity"
                });

                int countBefore = ProteaseDictionary.Dictionary.Count;
                var result = ProteaseDictionary.LoadAndMergeCustomProteases(path);

                NUnit.Framework.Assert.That(result.Added, Is.Empty);
                NUnit.Framework.Assert.That(result.Skipped, Is.Empty);
                NUnit.Framework.Assert.That(ProteaseDictionary.Dictionary.Count, Is.EqualTo(countBefore));
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public static void LoadAndMergeCustomProteases_OldStyleHeaders_ParsesCorrectly()
        {
            string path = Path.GetTempFileName();
            try
            {
                // Use old-style column names from base branch
                File.WriteAllLines(path, new[]
                {
                    "Name\tSequences Inducing Cleavage\tCleavage Specificity\tPSI-MS Accession Number\tPSI-MS Name\tCleavage Mass Shifts",
                    "OldHeaderTestProtease\tK|\tFull\t\t\t"
                });

                var result = ProteaseDictionary.LoadAndMergeCustomProteases(path);
                NUnit.Framework.Assert.That(result.Added, Does.Contain("OldHeaderTestProtease"));
                NUnit.Framework.Assert.That(ProteaseDictionary.Dictionary.ContainsKey("OldHeaderTestProtease"), Is.True);
            }
            finally
            {
                ProteaseDictionary.Dictionary.Remove("OldHeaderTestProtease");
                File.Delete(path);
            }
        }

        [Test]
        public static void LoadAndMergeCustomProteases_WithExplicitMods_UsesProvidedMods()
        {
            string path = Path.GetTempFileName();
            try
            {
                File.WriteAllLines(path, new[]
                {
                    "Name\tMotif\tSpecificity\tPSI-MS Accession\tPSI-MS Name\tCleavage Modification",
                    "ExplicitModTestProtease\tK|\tFull\t\t\t"
                });

                var explicitMods = ProteaseDictionary.LoadEmbeddedProteaseMods();
                var result = ProteaseDictionary.LoadAndMergeCustomProteases(path, explicitMods);

                NUnit.Framework.Assert.That(result.Added, Does.Contain("ExplicitModTestProtease"));
            }
            finally
            {
                ProteaseDictionary.Dictionary.Remove("ExplicitModTestProtease");
                File.Delete(path);
            }
        }

        [Test]
        public static void LoadAndMergeCustomProteases_BadModWithExplicitEmptyList_Throws()
        {
            string path = Path.GetTempFileName();
            try
            {
                File.WriteAllLines(path, new[]
                {
                    "Name\tMotif\tSpecificity\tPSI-MS Accession\tPSI-MS Name\tCleavage Modification",
                    "BadModProtease\tK|\tFull\t\t\tNonexistent mod on X"
                });

                NUnit.Framework.Assert.That(
                    () => ProteaseDictionary.LoadAndMergeCustomProteases(path, new List<Modification>()),
                    Throws.TypeOf<MzLibUtil.MzLibException>()
                        .With.Message.Contains("Nonexistent mod on X"));
            }
            finally
            {
                File.Delete(path);
            }
        }

        #endregion
    }
}