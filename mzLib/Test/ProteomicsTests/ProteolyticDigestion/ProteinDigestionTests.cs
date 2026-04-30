using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics;
using Omics.BioPolymer;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Omics.Modifications.IO;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases;
using static Chemistry.PeriodicTable;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class ProteinDigestionTests
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

        /// <summary>
        /// Tests that proteases with invalid or duplicate-within-file definitions throw
        /// when loaded via LoadAndMergeCustomProteases. The bad-mod case (specifying a
        /// cleavage modification name that does not exist in the provided mods list)
        /// should throw. The duplicate-name-within-file cases should throw.
        /// </summary>
        [Test]
        public static void ProteaseLoader_InvalidFiles_Throw()
        {
            string path_badMod = Path.Combine(TestContext.CurrentContext.TestDirectory, "ProteomicsTests", "ProteaseFilesForLoadingTests", "TestProteases_badMod.tsv");
            string path_badModDupName = Path.Combine(TestContext.CurrentContext.TestDirectory, "ProteomicsTests", "ProteaseFilesForLoadingTests", "TestProteases_badMod_dupName.tsv");
            string path_dupName = Path.Combine(TestContext.CurrentContext.TestDirectory, "ProteomicsTests", "ProteaseFilesForLoadingTests", "TestProteases_dupName.tsv");
            string path_modDupName = Path.Combine(TestContext.CurrentContext.TestDirectory, "ProteomicsTests", "ProteaseFilesForLoadingTests", "TestProteases_Mod_dupName.tsv");

            // A file referencing a modification name that doesn't exist should throw
            NUnit.Framework.Assert.Throws<MzLibException>(
                () => ProteaseDictionary.LoadAndMergeCustomProteases(path_badMod));

            // A file with both a bad mod reference and a duplicate name should throw
            NUnit.Framework.Assert.Throws<MzLibException>(
                () => ProteaseDictionary.LoadAndMergeCustomProteases(path_badModDupName));

            // A file with duplicate names within a single file should throw
            NUnit.Framework.Assert.Throws<MzLibException>(
                () => ProteaseDictionary.LoadAndMergeCustomProteases(path_dupName));

            NUnit.Framework.Assert.Throws<MzLibException>(
                () => ProteaseDictionary.LoadAndMergeCustomProteases(path_modDupName));
        }

        /// <summary>
        /// Tests CNBr digestion behavior using proteases loaded from a custom file.
        /// CNBr (with cleavage mod), CNBr_old (without), and CNBr_N (N-terminal cleavage
        /// with mod) are loaded as custom entries — they use names that do NOT collide
        /// with the embedded resource so they are added cleanly.
        /// </summary>
        [Test]
        public static void CNBrProteinDigestion()
        {
            var proteaseMods = ModificationLoader.ReadModsFromFile(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", "ProteaseMods.txt"),
                out var errors).ToList();

            var prot = new Protein("PEPTIDEMPEPTIDEM", null);
            var prot2 = new Protein("MPEPTIDEMPEPTIDE", null);

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "ProteomicsTests", "ProteaseFilesForLoadingTests", "DoubleProtease.tsv");
            NUnit.Framework.Assert.That(File.Exists(path));

            // Load the custom file — these proteases must have names not in the embedded resource
            // (e.g. "CNBr_custom", "CNBr_old_custom", "CNBr_N_custom") to be accepted.
            // If the file uses names that collide with embedded entries they will be skipped.
            var result = ProteaseDictionary.LoadAndMergeCustomProteases(path, proteaseMods);

            // Verify the custom proteases were actually added (not silently skipped due to a name
            // collision with the embedded resource). This guards against regressions in the merge
            // semantics — without this, a fallback to the embedded definition would let later
            // assertions pass for the wrong reason.
            Assert.That(result.Added, Contains.Item("CNBr_custom"));
            Assert.That(result.Added, Contains.Item("CNBr_old_custom"));
            Assert.That(result.Added, Contains.Item("CNBr_N_custom"));

            var protease1 = ProteaseDictionary.Dictionary["CNBr_custom"];
            var protease2 = ProteaseDictionary.Dictionary["CNBr_old_custom"];
            var protease3 = ProteaseDictionary.Dictionary["CNBr_N_custom"];

            DigestionParams digestionParams1 = new DigestionParams(
                protease: protease1.Name, maxMissedCleavages: 0, minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            DigestionParams digestionParams2 = new DigestionParams(
                protease: protease2.Name, maxMissedCleavages: 0, minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            DigestionParams digestionParams3 = new DigestionParams(
                protease: protease3.Name, maxMissedCleavages: 0, minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            var peps1 = prot.Digest(digestionParams1, new List<Modification>(), new List<Modification>()).ToList();
            var peps2 = prot.Digest(digestionParams2, new List<Modification>(), new List<Modification>()).ToList();
            var peps3 = prot2.Digest(digestionParams3, new List<Modification>(), new List<Modification>()).ToList();

            Assert.AreNotEqual(null, protease3.CleavageMod);
            Assert.AreEqual("M", protease3.CleavageMod.Target.ToString());
            Assert.AreNotEqual(peps3[0].MonoisotopicMass, peps3[1].MonoisotopicMass);

            Assert.AreEqual(882.39707781799996, peps3[1].MonoisotopicMass);
            Assert.AreEqual(930.400449121, peps3[0].MonoisotopicMass);

            Assert.AreEqual(null, protease2.CleavageMod);
            Assert.AreNotEqual(null, protease1.CleavageMod);
            Assert.AreEqual("M", protease1.CleavageMod.Target.ToString());

            Assert.AreEqual(peps1[1].MonoisotopicMass, peps2[1].MonoisotopicMass);
            Assert.AreEqual(peps1[1].MonoisotopicMass, peps2[0].MonoisotopicMass);
            Assert.AreEqual(peps2[0].MonoisotopicMass, peps2[1].MonoisotopicMass);
            Assert.AreNotEqual(peps1[0].MonoisotopicMass, peps1[1].MonoisotopicMass);
            Assert.AreNotEqual(peps1[0].MonoisotopicMass, peps2[0].MonoisotopicMass);
            Assert.AreNotEqual(peps1[0].MonoisotopicMass, peps2[1].MonoisotopicMass);

            Assert.AreEqual(882.39707781799996, peps1[0].MonoisotopicMass);
            Assert.AreEqual(930.400449121, peps1[1].MonoisotopicMass);

            // Clean up any custom entries that were added
            foreach (var name in result.Added)
                ProteaseDictionary.Dictionary.Remove(name);
        }

        /// <summary>
        /// Verifies that the embedded proteases.tsv resource exists with the correct
        /// assembly resource name, and that the dictionary is populated at startup.
        /// </summary>
        [Test]
        public static void LoadProteaseDictionary_EmbeddedResource_ExistsAndLoads()
        {
            var assembly = Assembly.GetAssembly(typeof(ProteaseDictionary));
            var resourceNames = assembly.GetManifestResourceNames();
            Assert.That(resourceNames, Contains.Item("Proteomics.ProteolyticDigestion.proteases.tsv"),
                $"Expected embedded resource not found. Available: {string.Join(", ", resourceNames)}");

            var dictionary = ProteaseDictionary.Dictionary;
            Assert.That(dictionary, Is.Not.Null);
            Assert.That(dictionary.Count, Is.GreaterThan(0));

            Assert.That(dictionary.ContainsKey("trypsin|P"), Is.True);
            Assert.That(dictionary["trypsin|P"].CleavageSpecificity, Is.EqualTo(CleavageSpecificity.Full));
            Assert.That(dictionary["trypsin|P"].DigestionMotifs.Count, Is.EqualTo(2)); // K[P]| and R[P]|
        }

        /// <summary>
        /// Tests backward compatibility for old-style protease names.
        /// Names like "chymotrypsin (don't cleave before proline)" should automatically
        /// resolve to "chymotrypsin|P" via NormalizeProteaseName.
        /// </summary>
        [Test]
        public static void GetProtease_OldStyleName_ResolvesToNewFormat()
        {
            var testCases = new[]
            {
                ("chymotrypsin (don't cleave before proline)", "chymotrypsin|P"),
                ("trypsin (don't cleave before proline)",      "trypsin|P"),
                ("Lys-C (don't cleave before proline)",        "Lys-C|P"),
            };

            foreach (var (oldName, expectedNewName) in testCases)
            {
                string normalizedName = ProteaseDictionary.NormalizeProteaseName(oldName);
                Assert.That(normalizedName, Is.EqualTo(expectedNewName),
                    $"Failed to normalize '{oldName}' to '{expectedNewName}'");

                var protease = ProteaseDictionary.GetProtease(oldName);
                Assert.That(protease, Is.Not.Null, $"GetProtease failed for '{oldName}'");
                Assert.That(protease.Name, Is.EqualTo(expectedNewName),
                    $"GetProtease returned wrong protease for '{oldName}'");

                bool found = ProteaseDictionary.TryGetProtease(oldName, out var protease2);
                Assert.That(found, Is.True, $"TryGetProtease failed for '{oldName}'");
                Assert.That(protease2.Name, Is.EqualTo(expectedNewName));
            }
        }

        /// <summary>
        /// Tests that GetProtease works with exact new-style names.
        /// </summary>
        [Test]
        public static void GetProtease_NewStyleName_WorksDirectly()
        {
            var protease = ProteaseDictionary.GetProtease("trypsin|P");
            Assert.That(protease, Is.Not.Null);
            Assert.That(protease.Name, Is.EqualTo("trypsin|P"));

            bool found = ProteaseDictionary.TryGetProtease("chymotrypsin|P", out var protease2);
            Assert.That(found, Is.True);
            Assert.That(protease2.Name, Is.EqualTo("chymotrypsin|P"));
        }

        /// <summary>
        /// Tests that GetProtease throws for unknown names and TryGetProtease returns false.
        /// </summary>
        [Test]
        public static void GetProtease_UnknownProtease_ThrowsKeyNotFoundException()
        {
            Assert.Throws<KeyNotFoundException>(
                () => ProteaseDictionary.GetProtease("nonexistent protease"));

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
        /// Tests that a custom file with insufficient fields throws with a helpful message.
        /// This exercises the TSV parser's field-count validation in LoadAndMergeCustomProteases.
        /// </summary>
        [Test]
        public static void LoadAndMergeCustomProteases_InsufficientFields_ThrowsWithHelpfulMessage()
        {
            string testFile = Path.Combine(TestContext.CurrentContext.TestDirectory, "test_insufficient_fields.tsv");
            string[] lines =
            {
                "Name\tMotif\tSpecificity\tPSI-MS Accession\tPSI-MS Name\tCleavage Modification",
                "ValidCustomProtease\tK|\tfull\t\t\t",
                "InvalidProtease\tK|"  // Missing Specificity — only 2 columns
            };
            File.WriteAllLines(testFile, lines);

            try
            {
                var exception = Assert.Throws<MzLibException>(
                    () => ProteaseDictionary.LoadAndMergeCustomProteases(testFile));

                Assert.That(exception.Message, Does.Contain("has only 2 field(s)"));
                Assert.That(exception.Message, Does.Contain("extend to column 3"));
                Assert.That(exception.Message, Does.Contain("InvalidProtease"));
            }
            finally
            {
                File.Delete(testFile);
                // ValidCustomProtease will not have been added because the parse threw mid-file
            }
        }

        /// <summary>
        /// Tests that custom proteases with only the three required fields are accepted,
        /// with optional fields defaulting to empty.
        /// </summary>
        [Test]
        public static void LoadAndMergeCustomProteases_MinimalFields_DefaultsOptionalFieldsToEmpty()
        {
            string testFile = Path.Combine(TestContext.CurrentContext.TestDirectory, "test_minimal_fields.tsv");
            string[] lines =
            {
                "Name\tMotif\tSpecificity\tPSI-MS Accession\tPSI-MS Name\tCleavage Modification",
                "MinimalCustomProtease\tK|\tfull"  // Only 3 required fields
            };
            File.WriteAllLines(testFile, lines);

            try
            {
                var result = ProteaseDictionary.LoadAndMergeCustomProteases(testFile);

                Assert.That(result.Added, Contains.Item("MinimalCustomProtease"));
                Assert.That(result.Skipped, Is.Empty);

                var protease = ProteaseDictionary.Dictionary["MinimalCustomProtease"];
                Assert.That(protease.Name, Is.EqualTo("MinimalCustomProtease"));
                Assert.That(protease.CleavageSpecificity, Is.EqualTo(CleavageSpecificity.Full));
                Assert.That(protease.PsiMsAccessionNumber, Is.EqualTo(string.Empty));
                Assert.That(protease.PsiMsName, Is.EqualTo(string.Empty));
                Assert.That(protease.CleavageMod, Is.Null);
            }
            finally
            {
                File.Delete(testFile);
                ProteaseDictionary.Dictionary.Remove("MinimalCustomProtease");
            }
        }

        [Test]
        public static void TestGoodPeptide()
        {
            var prot = new Protein("MNNNKQQQQ", null);
            var motifList = DigestionMotif.ParseDigestionMotifsFromString("K|");
            var protease = new Protease("CustomizedProtease", CleavageSpecificity.Full, null, null, motifList);
            ProteaseDictionary.Dictionary[protease.Name] = protease;
            try
            {
                DigestionParams digestionParams = new DigestionParams(
                    protease: protease.Name, maxMissedCleavages: 0, minPeptideLength: 1,
                    initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
                var ye = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();

                Assert.AreEqual(2, ye.Count);

                var pep1 = ye[0];
                Assert.IsTrue(pep1.MonoisotopicMass > 0);

                var test = new List<Product>();
                pep1.Fragment(DissociationType.HCD, FragmentationTerminus.Both, test);
                foreach (var huh in test)
                    Assert.IsTrue(huh.NeutralMass > 0);

                var pep2 = ye[1];
                pep1.Fragment(DissociationType.HCD, FragmentationTerminus.Both, test);
                Assert.IsTrue(pep2.MonoisotopicMass > 0);
                foreach (var huh in test)
                    Assert.IsTrue(huh.NeutralMass > 0);
            }
            finally
            {
                ProteaseDictionary.Dictionary.Remove("CustomizedProtease");
            }
        }

        [Test]
        public static void TestNoCleavage()
        {
            List<Modification> fixedModifications = new List<Modification>();
            var prot = new Protein("MNNNKQQQQ", null, null, null,
                new Dictionary<int, List<Modification>>(),
                new List<TruncationProduct> { new TruncationProduct(5, 6, "lala") });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 5);
            var ye = prot.Digest(digestionParams, fixedModifications, new List<Modification>()).ToList();
            Assert.AreEqual(3, ye.Count);
        }

        [Test]
        public static void TestBadPeptide()
        {
            var prot = new Protein("MNNNKQQXQ", null);
            var motifList = DigestionMotif.ParseDigestionMotifsFromString("K|");
            var protease = new Protease("Custom Protease7", CleavageSpecificity.Full, null, null, motifList);
            ProteaseDictionary.Dictionary[protease.Name] = protease;
            try
            {
                DigestionParams digestionParams = new DigestionParams(
                    protease: protease.Name, maxMissedCleavages: 0, minPeptideLength: 1,
                    initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
                var ye = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();

                Assert.AreEqual(2, ye.Count);
                var pep1 = ye[0];
                Assert.IsTrue(pep1.MonoisotopicMass > 0);

                var fragments = new List<Product>();
                pep1.Fragment(DissociationType.HCD, FragmentationTerminus.Both, fragments);
                foreach (var huh in fragments)
                    Assert.IsTrue(huh.NeutralMass > 0);

                var pep2 = ye[1];
                Assert.IsNaN(pep2.MonoisotopicMass);
                var cool = new List<Product>();
                pep2.Fragment(DissociationType.HCD, FragmentationTerminus.Both, cool);
                Assert.IsTrue(cool[0].NeutralMass > 0);
                Assert.IsTrue(cool[1].NeutralMass > 0);
                Assert.IsTrue(cool[2].NeutralMass > 0);
                Assert.IsTrue(double.IsNaN(cool[3].NeutralMass));
                Assert.IsTrue(double.IsNaN(cool[4].NeutralMass));
                Assert.IsTrue(double.IsNaN(cool[5].NeutralMass));
                Assert.IsTrue(cool.Count == 6);
            }
            finally
            {
                ProteaseDictionary.Dictionary.Remove("Custom Protease7");
            }
        }

        [Test]
        public static void TestPeptideWithSetModifications()
        {
            var prot = new Protein("M", null);
            DigestionParams digestionParams = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1, maxModsForPeptides: 3);
            List<Modification> variableModifications = new List<Modification>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif);

            variableModifications.Add(new Modification(_originalId: "ProtNmod", _target: motif, _locationRestriction: "N-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            variableModifications.Add(new Modification(_originalId: "pepNmod", _target: motif, _locationRestriction: "Peptide N-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            variableModifications.Add(new Modification(_originalId: "resMod", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            variableModifications.Add(new Modification(_originalId: "PepCmod", _target: motif, _locationRestriction: "Peptide C-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            variableModifications.Add(new Modification(_originalId: "ProtCmod", _target: motif, _locationRestriction: "C-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));

            var ye = prot.Digest(digestionParams, new List<Modification>(), variableModifications).ToList();
            Assert.AreEqual(3 * 2 * 3, ye.Count);
            Assert.AreEqual("[H]M[H]-[H]", ye.Last().SequenceWithChemicalFormulas);

            double m1 = 5 * GetElement("H").PrincipalIsotope.AtomicMass + Residue.ResidueMonoisotopicMass['M'] + GetElement("O").PrincipalIsotope.AtomicMass;
            m1 = Math.Round(m1, 9, MidpointRounding.AwayFromZero);
            double m2 = ye.Last().MonoisotopicMass;
            double m3 = m1 - m2;
            Assert.IsTrue(m3 < 1e-9);
        }

        [Test]
        public static void TestPeptideDigestion_FixedModifications_ProtModsOverwritePepMods()
        {
            string baseSequence = "M";
            var prot = new Protein(baseSequence, null);
            DigestionParams digestionParams = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1, maxModsForPeptides: 3);
            List<Modification> fixedMods = new List<Modification>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif);
            fixedMods.Add(new Modification(_originalId: "ProtNmod", _target: motif, _locationRestriction: "N-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "pepNmod", _target: motif, _locationRestriction: "Peptide N-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "resMod", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "ProtCmod", _target: motif, _locationRestriction: "C-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "PepCmod", _target: motif, _locationRestriction: "Peptide C-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            var ok = prot.Digest(digestionParams, fixedMods, new List<Modification>()).ToList();

            Assert.AreEqual(1, ok.Count);

            string expectedFullSequence = "[:ProtNmod on M]M[:resMod on M]-[:ProtCmod on M]";
            Assert.AreEqual(expectedFullSequence, ok.First().FullSequence);
            var mods = ok.First().AllModsOneIsNterminus;

            NUnit.Framework.Assert.That(IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(expectedFullSequence), Is.EqualTo(baseSequence));
            NUnit.Framework.Assert.That(IBioPolymerWithSetMods.DetermineFullSequence(baseSequence, mods), Is.EqualTo(expectedFullSequence));
            NUnit.Framework.Assert.That(ok.First().DetermineFullSequence(), Is.EqualTo(expectedFullSequence));

            Assert.AreEqual("[H]M[H]-[H]", ok.First().SequenceWithChemicalFormulas);
            Assert.AreEqual(5 * GetElement("H").PrincipalIsotope.AtomicMass + Residue.ResidueMonoisotopicMass['M'] + GetElement("O").PrincipalIsotope.AtomicMass, ok.Last().MonoisotopicMass, 1e-9);
        }

        [Test]
        public static void TestPeptideDigestion_FixedModifications_ProtModsOverwritePepMods_RandomizedModOrder()
        {
            var rand = new Random(42);
            string baseSequence = "M";
            var prot = new Protein(baseSequence, null);
            DigestionParams digestionParams = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1, maxModsForPeptides: 3);
            List<Modification> fixedMods = new List<Modification>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif);
            fixedMods.Add(new Modification(_originalId: "ProtNmod", _target: motif, _locationRestriction: "N-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "pepNmod", _target: motif, _locationRestriction: "Peptide N-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "resMod", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "ProtCmod", _target: motif, _locationRestriction: "C-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "PepCmod", _target: motif, _locationRestriction: "Peptide C-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));

            int expectedDigestionProducts = 1;
            string expectedFullSequence = "[:ProtNmod on M]M[:resMod on M]-[:ProtCmod on M]";
            string expectedSequenceWithChemicalFormulas = "[H]M[H]-[H]";
            double expectedMonoisotopicMass = 5 * GetElement("H").PrincipalIsotope.AtomicMass + Residue.ResidueMonoisotopicMass['M'] + GetElement("O").PrincipalIsotope.AtomicMass;



            // randomly scramble all mods, digest, and ensure the answer is correct. 
            for (int i = 0; i < 10; i++)
            {
                var shuffledFixedMods = fixedMods.OrderBy(a => rand.Next()).ToList();
                var ok = prot.Digest(digestionParams, shuffledFixedMods, new List<Modification>()).ToList();
                var mods = ok.First().AllModsOneIsNterminus;

                NUnit.Framework.Assert.That(IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(expectedFullSequence), Is.EqualTo(baseSequence));
                NUnit.Framework.Assert.That(IBioPolymerWithSetMods.DetermineFullSequence(baseSequence, mods), Is.EqualTo(expectedFullSequence));
                NUnit.Framework.Assert.That(ok.First().DetermineFullSequence(), Is.EqualTo(expectedFullSequence));

                Assert.AreEqual(expectedDigestionProducts, ok.Count);
                Assert.AreEqual(expectedFullSequence, ok.First().FullSequence);
                Assert.AreEqual(expectedSequenceWithChemicalFormulas, ok.First().SequenceWithChemicalFormulas);
                Assert.AreEqual(expectedMonoisotopicMass, ok.Last().MonoisotopicMass, 1e-9);
            }
        }

        [Test]
        public static void TestPeptideDigestion_FixedModifications_ProtModsOverwritePepMods_TwoProducts()
        {
            var prot = new Protein("MKM", null);
            DigestionParams digestionParams = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1, maxModsForPeptides: 3, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            List<Modification> fixedMods = new List<Modification>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif mMotif);
            ModificationMotif.TryGetMotif("K", out ModificationMotif kMotif);

            fixedMods.Add(new Modification(_originalId: "ProtNmod", _target: mMotif, _locationRestriction: "N-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "ProtNmod", _target: kMotif, _locationRestriction: "N-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "pepNmod", _target: mMotif, _locationRestriction: "Peptide N-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "pepNmod", _target: kMotif, _locationRestriction: "Peptide N-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "resMod", _target: mMotif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "ProtCmod", _target: mMotif, _locationRestriction: "C-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "ProtCmod", _target: kMotif, _locationRestriction: "C-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "PepCmod", _target: mMotif, _locationRestriction: "Peptide C-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new Modification(_originalId: "PepCmod", _target: kMotif, _locationRestriction: "Peptide C-terminal.", _chemicalFormula: ChemicalFormula.ParseFormula("H"), _monoisotopicMass: GetElement(1).PrincipalIsotope.AtomicMass));

            var ok = prot.Digest(digestionParams, fixedMods, new List<Modification>()).ToList();

            Assert.AreEqual(2, ok.Count);
            Assert.AreEqual("[:ProtNmod on M]M[:resMod on M]K-[:PepCmod on K]", ok.First().FullSequence);
            Assert.AreEqual("[:pepNmod on M]M[:resMod on M]-[:ProtCmod on M]", ok.Skip(1).First().FullSequence);
            Assert.AreEqual("[H]M[H]K-[H]", ok.First().SequenceWithChemicalFormulas);
            Assert.AreEqual("[H]M[H]-[H]", ok.Skip(1).First().SequenceWithChemicalFormulas);
            Assert.AreEqual(5 * GetElement("H").PrincipalIsotope.AtomicMass + Residue.ResidueMonoisotopicMass['M'] + GetElement("O").PrincipalIsotope.AtomicMass, ok.Last().MonoisotopicMass, 1e-9);
        }

        [Test]
        public static void TestDigestIndices()
        {
            ModificationMotif.TryGetMotif("N", out ModificationMotif motifN);
            ModificationMotif.TryGetMotif("R", out ModificationMotif motifR);
            Modification modN = new Modification("myMod", null, "myModType", null, motifN, "Anywhere.", null, 10, null, null, null, null, null, null);
            Modification modR = new Modification("myMod", null, "myModType", null, motifR, "Anywhere.", null, 10, null, null, null, null, null, null);
            IDictionary<int, List<Modification>> modDict = new Dictionary<int, List<Modification>>
            {
                {2, new List<Modification> { modN } },
                {8, new List<Modification> { modR } }
            };
            var prot = new Protein("MNNNNKRRRRR", null, null, null, modDict, isDecoy: true);

            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 5, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            var digestedList = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();
            var ok1 = digestedList[1];
            var ok2 = digestedList[3];

            Assert.AreEqual(1, ok1.NumMods);
            Assert.IsTrue(ok1.AllModsOneIsNterminus.ContainsKey(3));
            Assert.AreEqual(1, ok2.NumMods);
            Assert.IsTrue(ok2.AllModsOneIsNterminus.ContainsKey(3));
        }

        [Test]
        public static void TestDigestDecoy()
        {
            ModificationMotif.TryGetMotif("N", out ModificationMotif motifN);
            ModificationMotif.TryGetMotif("R", out ModificationMotif motifR);
            Modification modN = new Modification("myMod", null, "myModType", null, motifN, "Anywhere.", null, 10, null, null, null, null, null, null);
            Modification modR = new Modification("myMod", null, "myModType", null, motifR, "Anywhere.", null, 10, null, null, null, null, null, null);
            IDictionary<int, List<Modification>> modDict = new Dictionary<int, List<Modification>>
            {
                {2, new List<Modification> { modN } },
                {8, new List<Modification> { modR } }
            };
            var prot = new Protein("MNNNNKRRRRR", null, null, null, modDict, isDecoy: true);

            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 5, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            var digestedList = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();
            var ok1 = digestedList[1];
            var ok2 = digestedList[3];

            Assert.AreEqual(1, ok1.NumMods);
            Assert.IsTrue(ok1.AllModsOneIsNterminus.ContainsKey(3));
            Assert.AreEqual(1, ok2.NumMods);
            Assert.IsTrue(ok2.AllModsOneIsNterminus.ContainsKey(3));

            prot = new Protein("MNNNNKRRRRR", null, null, null, modDict);
            ok1 = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).First();
            ok2 = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).Last();

            Assert.AreEqual(0, ok1.NumMods);
            Assert.IsFalse(ok1.AllModsOneIsNterminus.Any());
            Assert.AreEqual(2, ok2.NumMods);
            Assert.IsTrue(ok2.AllModsOneIsNterminus.Any());
        }

        [Test]
        public static void TestGoodPeptideWithLength()
        {
            var prot = new Protein("MNNNKQQQQMNNNKQQQQ", null);

            DigestionParams digestionParams = new DigestionParams("trypsin|P", maxMissedCleavages: 0, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var ye = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();
            digestionParams = new DigestionParams("trypsin|P", maxMissedCleavages: 0, minPeptideLength: 5, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var ye1 = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();
            digestionParams = new DigestionParams("trypsin|P", maxMissedCleavages: 0, minPeptideLength: 1, maxPeptideLength: 5, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var ye2 = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();
            digestionParams = new DigestionParams("trypsin|P", maxMissedCleavages: 0, minPeptideLength: 5, maxPeptideLength: 8, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var ye3 = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();
            Assert.AreEqual(3, ye.Count);
            Assert.AreEqual(2, ye1.Count);
            Assert.AreEqual(2, ye2.Count);
            Assert.AreEqual(1, ye3.Count);
        }

        [Test]
        public static void Test_ProteinDigest()
        {
            DigestionParams d = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 5, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            ModificationMotif.TryGetMotif("D", out ModificationMotif motif);
            Modification mod = new Modification(_originalId: "mod1", _modificationType: "mt", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 10);

            IDictionary<int, List<Modification>> oneBasedModification = new Dictionary<int, List<Modification>>
            {
                { 3, new List<Modification>{ mod } }
            };

            Protein prot1 = new Protein("MEDEEK", "prot1", oneBasedModifications: oneBasedModification);

            var pep1 = prot1.Digest(d, new List<Modification>(), new List<Modification>()).First();
            var pep2 = prot1.Digest(d, new List<Modification>(), new List<Modification>()).Last();

            Assert.AreEqual("MEDEEK", pep1.FullSequence);
            Assert.AreEqual("MED[mt:mod1 on D]EEK", pep2.FullSequence);
        }

        [Test]
        public static void TestDigestionOfSameProteinFromDifferentXmls()
        {
            DigestionParams digestionParams = new DigestionParams("trypsin|P", maxMissedCleavages: 2, minPeptideLength: 7, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            ModificationMotif.TryGetMotif("C", out ModificationMotif motif);
            Modification carbamidomethylOnC = new Modification(_originalId: "Carbamidomethyl on C", _modificationType: "Common Fixed", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("C2H3NO"));
            var fixedModifications = new List<Modification> { carbamidomethylOnC };
            ModificationMotif.TryGetMotif("M", out ModificationMotif motifM);
            Modification oxidationOnM = new Modification(_originalId: "Oxidation on M", _modificationType: "Common Variable", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("O"));
            var variableModifications = new List<Modification> { oxidationOnM };

            var dbFive = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "SingleEntry_ModOrder1.xml");
            var dbSix = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "SingleEntry_ModOrder2.xml");

            var proteins5 = ProteinDbLoader.LoadProteinXML(dbFive, true, DecoyType.None, null, false, null, out var unknownModificationsFive);
            var proteins6 = ProteinDbLoader.LoadProteinXML(dbSix, true, DecoyType.None, null, false, null, out var unknownModificationsSix);

            var fiveMods = ProteinDbLoader.GetPtmListFromProteinXml(dbFive);
            var sixMods = ProteinDbLoader.GetPtmListFromProteinXml(dbSix);

            Assert.AreEqual(fiveMods.Count, sixMods.Count);
            CollectionAssert.AreEquivalent(fiveMods, sixMods);

            var peptides5 = proteins5.First().Digest(digestionParams, fixedModifications, variableModifications).ToList();
            var peptides6 = proteins6.First().Digest(digestionParams, fixedModifications, variableModifications).ToList();
            Assert.AreEqual(peptides5.Count, peptides6.Count);
            CollectionAssert.AreEqual(peptides5, peptides6);
        }

        [Test]
        [TestCase("cRAP_databaseGPTMD.xml")]
        [TestCase("uniprot_aifm1.fasta")]
        public static void TestDecoyScramblingIsReproducible(string fileName)
        {
            var dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", fileName);
            DecoyType decoyType = DecoyType.Reverse;
            List<Protein> proteins1 = null;
            List<Protein> proteins2 = null;
            if (fileName.Contains(".xml"))
            {
                proteins1 = ProteinDbLoader.LoadProteinXML(dbPath, true, decoyType, null, false, null, out var unknownModifications);
                proteins2 = ProteinDbLoader.LoadProteinXML(dbPath, true, decoyType, null, false, null, out unknownModifications);
            }
            else if (fileName.Contains(".fasta"))
            {
                proteins1 = ProteinDbLoader.LoadProteinFasta(dbPath, true, decoyType, false, out var unknownModifications);
                proteins2 = ProteinDbLoader.LoadProteinFasta(dbPath, true, decoyType, false, out unknownModifications);
            }
            else
            {
                NUnit.Framework.Assert.Fail("Unknown file type");
            }

            DigestionParams d = new DigestionParams(maxMissedCleavages: 1, minPeptideLength: 5, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            var pepsToReplace = proteins1.Where(p => !p.IsDecoy)
                .SelectMany(p => p.Digest(d, new List<Modification>(), new List<Modification>()).ToList())
                .Select(pep => pep.BaseSequence)
                .ToHashSet();

            var singleDecoyPeptides = proteins1
                .Where(p => p.IsDecoy)
                .Select(p => p.Digest(d, new List<Modification>(), new List<Modification>()).Skip(2).Take(1))
                .Select(pwsm => pwsm.First().BaseSequence)
                .ToHashSet();

            pepsToReplace.UnionWith(singleDecoyPeptides);

            List<Protein> decoys1 = new();
            foreach (var protein in proteins1.Where(p => p.IsDecoy))
                decoys1.Add(DecoySequenceValidator.ScrambleDecoyBioPolymer(protein, d, pepsToReplace));

            List<Protein> decoys2 = new();
            foreach (var protein in proteins2.Where(p => p.IsDecoy))
                decoys2.Add(DecoySequenceValidator.ScrambleDecoyBioPolymer(protein, d, pepsToReplace));

            Assert.AreEqual(decoys1.Count, decoys2.Count);
            foreach (var decoyPair in decoys1.Concat(decoys2).GroupBy(p => p.Accession))
            {
                Assert.AreEqual(2, decoyPair.Count());
                Assert.AreEqual(decoyPair.First().BaseSequence, decoyPair.Last().BaseSequence);
            }
        }

        [Test]
        public static void TestDecoyScramblerReplacesPeptides()
        {
            DigestionParams d = new DigestionParams(maxMissedCleavages: 1, minPeptideLength: 5, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            Protein target = new Protein("MEDEEKFVGYKYGVFK", "target");
            Protein decoy = new Protein("EEDEMKYGVFKFVGYK", "decoy");

            var targetPep = target.Digest(d, new List<Modification>(), new List<Modification>());
            var decoyPep = decoy.Digest(d, new List<Modification>(), new List<Modification>());

            HashSet<string> targetPepSeqs = targetPep.Select(p => p.FullSequence).ToHashSet();
            var offendingDecoys = decoyPep.Where(p => targetPepSeqs.Contains(p.FullSequence)).Select(d => d.FullSequence).ToList();

            Assert.AreEqual(2, offendingDecoys.Count);

            Protein scrambledDecoy = DecoySequenceValidator.ScrambleDecoyBioPolymer(decoy, d, targetPepSeqs, offendingDecoys);
            var scrambledPep = scrambledDecoy.Digest(d, new List<Modification>(), new List<Modification>());

            Assert.AreEqual(decoyPep.Count(), scrambledPep.Count());
            Assert.IsFalse(scrambledPep.Any(p => offendingDecoys.Contains(p.FullSequence)));

            scrambledDecoy = DecoySequenceValidator.ScrambleDecoyBioPolymer(decoy, d, targetPepSeqs);
            scrambledPep = scrambledDecoy.Digest(d, new List<Modification>(), new List<Modification>());

            Assert.AreEqual(decoyPep.Count(), scrambledPep.Count());
            Assert.IsFalse(scrambledPep.Any(p => offendingDecoys.Contains(p.FullSequence)));
        }

        [Test]
        public static void TestReadPeptideFromString()
        {
            ModificationMotif.TryGetMotif("T", out ModificationMotif target);
            Modification carbamidomethylOnC = new Modification(_originalId: "Carbamidomethyl on C", _modificationType: "Common Fixed", _target: target, _chemicalFormula: ChemicalFormula.ParseFormula("C2H3NO"));
            string sequence = "HQVC[Common Fixed:Carbamidomethyl on C]TPGGTTIAGLC[Common Fixed:Carbamidomethyl on C]VMEEK";

            PeptideWithSetModifications peptide = new PeptideWithSetModifications(sequence, new Dictionary<string, Modification> { { carbamidomethylOnC.IdWithMotif, carbamidomethylOnC } });

            NUnit.Framework.Assert.That(peptide.BaseSequence == "HQVCTPGGTTIAGLCVMEEK");
            NUnit.Framework.Assert.That(peptide.FullSequence == sequence);
            NUnit.Framework.Assert.That(Math.Round(peptide.MonoisotopicMass, 5) == 2187.01225);
            NUnit.Framework.Assert.That(peptide.AllModsOneIsNterminus.First().Value.IdWithMotif == "Carbamidomethyl on C");
            NUnit.Framework.Assert.That(peptide.AllModsOneIsNterminus.Count == 2);
            NUnit.Framework.Assert.That(new HashSet<int>(peptide.AllModsOneIsNterminus.Keys).SetEquals(new HashSet<int>() { 5, 16 }));

            List<Product> theoreticalFragments = new List<Product>();
            peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, theoreticalFragments);
            NUnit.Framework.Assert.That(theoreticalFragments.Count > 0);
        }

        [Test]
        public static void TestGlycoPeptide()
        {
            var prot = new Protein("MNNNYTKQQQQKS", null);
            var motifList = DigestionMotif.ParseDigestionMotifsFromString("K|");
            var protease = new Protease("CustomizedProtease_diffname", CleavageSpecificity.Full, null, null, motifList);
            ProteaseDictionary.Dictionary[protease.Name] = protease;
            try
            {
                DigestionParams digestionParams = new DigestionParams(
                    protease: protease.Name, maxMissedCleavages: 0, minPeptideLength: 1,
                    initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain,
                    keepNGlycopeptide: true, keepOGlycopeptide: true);
                var ye = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();
                Assert.AreEqual(2, ye.Count);
            }
            finally
            {
                ProteaseDictionary.Dictionary.Remove("CustomizedProtease_diffname");
            }
        }

        [Test]
        public static void TestDigestionParamsClone()
        {
            DigestionParams digestionParams = new DigestionParams(
                protease: "top-down", maxMissedCleavages: 0, minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain,
                keepNGlycopeptide: true, keepOGlycopeptide: true);

            DigestionParams digestionParamsClone = (DigestionParams)digestionParams.Clone();
            Assert.AreEqual(digestionParams, digestionParamsClone);
            Assert.AreEqual(digestionParams.InitiatorMethionineBehavior, digestionParamsClone.InitiatorMethionineBehavior);
            Assert.AreEqual(digestionParams.MaxMissedCleavages, digestionParamsClone.MaxMissedCleavages);
            Assert.AreEqual(digestionParams.MaxModificationIsoforms, digestionParamsClone.MaxModificationIsoforms);
            Assert.AreEqual(digestionParams.MinLength, digestionParamsClone.MinLength);
            Assert.AreEqual(digestionParams.MaxLength, digestionParamsClone.MaxLength);
            Assert.AreEqual(digestionParams.MaxMods, digestionParamsClone.MaxMods);
            Assert.AreEqual(digestionParams.Protease, digestionParamsClone.Protease);
            Assert.AreEqual(digestionParams.SearchModeType, digestionParamsClone.SearchModeType);
            Assert.AreEqual(digestionParams.FragmentationTerminus, digestionParamsClone.FragmentationTerminus);
            Assert.AreEqual(digestionParams.GeneratehUnlabeledProteinsForSilac, digestionParamsClone.GeneratehUnlabeledProteinsForSilac);
            Assert.AreEqual(digestionParams.KeepNGlycopeptide, digestionParamsClone.KeepNGlycopeptide);
            Assert.AreEqual(digestionParams.KeepOGlycopeptide, digestionParamsClone.KeepOGlycopeptide);
            Assert.AreEqual(digestionParams.SpecificProtease, digestionParamsClone.SpecificProtease);
            NUnit.Framework.Assert.That(!ReferenceEquals(digestionParams, digestionParamsClone));

            digestionParams = new DigestionParams(
                protease: "top-down", maxMissedCleavages: 0, minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain,
                keepNGlycopeptide: true, keepOGlycopeptide: true,
                maxModificationIsoforms: 5, maxModsForPeptides: 6, maxPeptideLength: 7,
                searchModeType: CleavageSpecificity.None,
                fragmentationTerminus: FragmentationTerminus.C,
                generateUnlabeledProteinsForSilac: false);

            digestionParamsClone = (DigestionParams)digestionParams.Clone();
            Assert.AreEqual(digestionParams, digestionParamsClone);
            Assert.AreEqual(digestionParams.InitiatorMethionineBehavior, digestionParamsClone.InitiatorMethionineBehavior);
            Assert.AreEqual(digestionParams.MaxMissedCleavages, digestionParamsClone.MaxMissedCleavages);
            Assert.AreEqual(digestionParams.MaxModificationIsoforms, digestionParamsClone.MaxModificationIsoforms);
            Assert.AreEqual(digestionParams.MinLength, digestionParamsClone.MinLength);
            Assert.AreEqual(digestionParams.MaxLength, digestionParamsClone.MaxLength);
            Assert.AreEqual(digestionParams.MaxMods, digestionParamsClone.MaxMods);
            Assert.AreEqual(digestionParams.Protease, digestionParamsClone.Protease);
            Assert.AreEqual(digestionParams.SearchModeType, digestionParamsClone.SearchModeType);
            Assert.AreEqual(digestionParams.FragmentationTerminus, digestionParamsClone.FragmentationTerminus);
            Assert.AreEqual(digestionParams.GeneratehUnlabeledProteinsForSilac, digestionParamsClone.GeneratehUnlabeledProteinsForSilac);
            Assert.AreEqual(digestionParams.KeepNGlycopeptide, digestionParamsClone.KeepNGlycopeptide);
            Assert.AreEqual(digestionParams.KeepOGlycopeptide, digestionParamsClone.KeepOGlycopeptide);
            Assert.AreEqual(digestionParams.SpecificProtease, digestionParamsClone.SpecificProtease);
            NUnit.Framework.Assert.That(!ReferenceEquals(digestionParams, digestionParamsClone));
        }

        [Test]
        public static void TestDigestionParamsCloneWithNewTerminus()
        {
            DigestionParams digestionParams = new DigestionParams(
                protease: "top-down", maxMissedCleavages: 0, minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain,
                keepNGlycopeptide: true, keepOGlycopeptide: true);

            DigestionParams digestionParamsClone = (DigestionParams)digestionParams.Clone(FragmentationTerminus.N);
            Assert.AreNotEqual(digestionParams, digestionParamsClone);
            Assert.AreEqual(digestionParams.InitiatorMethionineBehavior, digestionParamsClone.InitiatorMethionineBehavior);
            Assert.AreEqual(digestionParams.MaxMissedCleavages, digestionParamsClone.MaxMissedCleavages);
            Assert.AreEqual(digestionParams.MaxModificationIsoforms, digestionParamsClone.MaxModificationIsoforms);
            Assert.AreEqual(digestionParams.MinLength, digestionParamsClone.MinLength);
            Assert.AreEqual(digestionParams.MaxLength, digestionParamsClone.MaxLength);
            Assert.AreEqual(digestionParams.MaxMods, digestionParamsClone.MaxMods);
            Assert.AreEqual(digestionParams.Protease, digestionParamsClone.Protease);
            Assert.AreEqual(digestionParams.SearchModeType, digestionParamsClone.SearchModeType);
            Assert.AreEqual(FragmentationTerminus.N, digestionParamsClone.FragmentationTerminus);
            Assert.AreEqual(digestionParams.GeneratehUnlabeledProteinsForSilac, digestionParamsClone.GeneratehUnlabeledProteinsForSilac);
            Assert.AreEqual(digestionParams.KeepNGlycopeptide, digestionParamsClone.KeepNGlycopeptide);
            Assert.AreEqual(digestionParams.KeepOGlycopeptide, digestionParamsClone.KeepOGlycopeptide);
            Assert.AreEqual(digestionParams.SpecificProtease, digestionParamsClone.SpecificProtease);
            NUnit.Framework.Assert.That(!ReferenceEquals(digestionParams, digestionParamsClone));

            digestionParams = new DigestionParams(
                protease: "top-down", maxMissedCleavages: 0, minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain,
                keepNGlycopeptide: true, keepOGlycopeptide: true,
                maxModificationIsoforms: 5, maxModsForPeptides: 6, maxPeptideLength: 7,
                searchModeType: CleavageSpecificity.None,
                fragmentationTerminus: FragmentationTerminus.None,
                generateUnlabeledProteinsForSilac: false);

            digestionParamsClone = (DigestionParams)digestionParams.Clone(FragmentationTerminus.N);
            Assert.AreNotEqual(digestionParams, digestionParamsClone);
            Assert.AreEqual(digestionParams.InitiatorMethionineBehavior, digestionParamsClone.InitiatorMethionineBehavior);
            Assert.AreEqual(digestionParams.MaxMissedCleavages, digestionParamsClone.MaxMissedCleavages);
            Assert.AreEqual(digestionParams.MaxModificationIsoforms, digestionParamsClone.MaxModificationIsoforms);
            Assert.AreEqual(digestionParams.MinLength, digestionParamsClone.MinLength);
            Assert.AreEqual(digestionParams.MaxLength, digestionParamsClone.MaxLength);
            Assert.AreEqual(digestionParams.MaxMods, digestionParamsClone.MaxMods);
            Assert.AreEqual(ProteaseDictionary.Dictionary["singleN"], digestionParamsClone.Protease);
            Assert.AreEqual(digestionParams.SearchModeType, digestionParamsClone.SearchModeType);
            Assert.AreEqual(FragmentationTerminus.N, digestionParamsClone.FragmentationTerminus);
            Assert.AreEqual(digestionParams.GeneratehUnlabeledProteinsForSilac, digestionParamsClone.GeneratehUnlabeledProteinsForSilac);
            Assert.AreEqual(digestionParams.KeepNGlycopeptide, digestionParamsClone.KeepNGlycopeptide);
            Assert.AreEqual(digestionParams.KeepOGlycopeptide, digestionParamsClone.KeepOGlycopeptide);
            Assert.AreEqual(digestionParams.SpecificProtease, digestionParamsClone.SpecificProtease);
            NUnit.Framework.Assert.That(!ReferenceEquals(digestionParams, digestionParamsClone));

            digestionParamsClone = (DigestionParams)digestionParams.Clone(FragmentationTerminus.C);
            Assert.AreNotEqual(digestionParams, digestionParamsClone);
            Assert.AreEqual(ProteaseDictionary.Dictionary["singleC"], digestionParamsClone.Protease);
            Assert.AreEqual(FragmentationTerminus.C, digestionParamsClone.FragmentationTerminus);
            NUnit.Framework.Assert.That(!ReferenceEquals(digestionParams, digestionParamsClone));
        }

        [Test]
        public static void TestWhenFixedModIsSamePositionAsUniProtModWithDigestion()
        {
            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml"));
            Dictionary<string, int> formalChargesDictionary = Loaders.GetFormalChargesDictionary(psiModDeserialized);
            List<Modification> UniProtPtms = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt"), formalChargesDictionary).ToList();

            DigestionParams digestionParams = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1, maxModsForPeptides: 3);
            List<Modification> fixedMods = new List<Modification>();
            ModificationMotif.TryGetMotif("S", out ModificationMotif serineMotif);
            ChemicalFormula ohFormula = ChemicalFormula.ParseFormula("OH");
            double ohMass = GetElement("O").PrincipalIsotope.AtomicMass + GetElement("H").PrincipalIsotope.AtomicMass;
            fixedMods.Add(new Modification(_originalId: "serineOhMod", _target: serineMotif, _locationRestriction: "Anywhere.", _chemicalFormula: ohFormula, _monoisotopicMass: ohMass));

            List<Protein> dbProteins = ProteinDbLoader.LoadProteinXML(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"xml.xml"),
                true, DecoyType.Reverse, UniProtPtms.Concat(fixedMods), false,
                new List<string>(), out Dictionary<string, Modification> un);

            Protein prot = dbProteins.First();
            var digestionProducts = prot.Digest(digestionParams, fixedMods, new List<Modification>()).ToList();
            var firstPeptideModifiedForms = digestionProducts.Where(p => p.BaseSequence == "MSGR").ToList();
            List<string> fullSequences = firstPeptideModifiedForms.Select(p => p.FullSequence).ToList();
            List<string> expectedFullSequences = new List<string>
            {
                "MS[:serineOhMod on S]GR",
                "MS[:serineOhMod on S]GR[UniProt:Asymmetric dimethylarginine on R]",
                "MS[:serineOhMod on S]GR[UniProt:Citrulline on R]",
                "MS[:serineOhMod on S]GR[UniProt:Omega-N-methylarginine on R]",
                "MS[:serineOhMod on S]GR[UniProt:Symmetric dimethylarginine on R]",
                "MS[UniProt:Phosphoserine on S]GR",
                "MS[UniProt:Phosphoserine on S]GR[UniProt:Asymmetric dimethylarginine on R]",
                "MS[UniProt:Phosphoserine on S]GR[UniProt:Citrulline on R]",
                "MS[UniProt:Phosphoserine on S]GR[UniProt:Omega-N-methylarginine on R]",
                "MS[UniProt:Phosphoserine on S]GR[UniProt:Symmetric dimethylarginine on R]"
            };
            CollectionAssert.AreEquivalent(expectedFullSequences, fullSequences);
        }

        /// <summary>
        /// Tests that custom proteases with new names are added, that embedded names are
        /// skipped with a warning (not a crash), and that the returned
        /// CustomDigestionAgentLoadResult correctly partitions Added vs Skipped.
        /// </summary>
        [Test]
        public static void LoadAndMergeCustomProteases_AddsNewAndSkipsEmbeddedCollisions()
        {
            int initialCount = ProteaseDictionary.Dictionary.Count;

            // "trypsin|P" exists in embedded resource → must be skipped
            // "MyLabProtease" is genuinely new → must be added
            string customProteaseFile = Path.Combine(TestContext.CurrentContext.TestDirectory, "test_custom_proteases.tsv");
            string[] lines =
            {
                "Name\tMotif\tSpecificity\tPSI-MS Accession\tPSI-MS Name\tCleavage Modification",
                "trypsin|P\tL[P]|\tfull\tMS:1001313\tTrypsin\t",
                "MyLabProtease\tE|\tfull\t\tCustom Glu-C variant\t"
            };
            File.WriteAllLines(customProteaseFile, lines);

            try
            {
                var originalTrypsin = ProteaseDictionary.Dictionary["trypsin|P"];

                var result = ProteaseDictionary.LoadAndMergeCustomProteases(customProteaseFile);

                // trypsin|P must be skipped, not added
                Assert.That(result.Skipped, Contains.Item("trypsin|P"),
                    "Embedded protease name must appear in Skipped");
                Assert.That(result.Added, Does.Not.Contain("trypsin|P"),
                    "Embedded protease name must not appear in Added");

                // MyLabProtease must be added
                Assert.That(result.Added, Contains.Item("MyLabProtease"),
                    "New protease name must appear in Added");
                Assert.That(result.Skipped, Does.Not.Contain("MyLabProtease"));

                // Dictionary count increases by exactly 1 (only MyLabProtease added)
                Assert.That(ProteaseDictionary.Dictionary.Count, Is.EqualTo(initialCount + 1));

                // Embedded trypsin|P definition must be unchanged
                Assert.That(ReferenceEquals(ProteaseDictionary.Dictionary["trypsin|P"], originalTrypsin), Is.True,
                    "The embedded trypsin|P object must not have been replaced");
                Assert.That(ProteaseDictionary.Dictionary["trypsin|P"].DigestionMotifs.Count, Is.EqualTo(2),
                    "Embedded trypsin|P should still have its original 2 motifs");

                // MyLabProtease should work for digestion
                var protein = new Protein("PEPTIDEEPEPTIDER", null);
                var myLabParams = new DigestionParams(
                    protease: "MyLabProtease", maxMissedCleavages: 0, minPeptideLength: 1);
                var myLabDigest = protein.Digest(myLabParams, new List<Modification>(), new List<Modification>()).ToList();
                Assert.That(myLabDigest.Count, Is.GreaterThan(0));
            }
            finally
            {
                File.Delete(customProteaseFile);
                ProteaseDictionary.Dictionary.Remove("MyLabProtease");
            }
        }
    }
}
