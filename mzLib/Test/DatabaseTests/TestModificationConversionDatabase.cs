using Chemistry;
using NUnit.Framework;
using Omics.BioPolymer;
using Omics.Modifications;
using Proteomics;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class TestModificationConversionDatabase
    {
        [Test]
        public static void TestConvertModsOnProtein_OneBasedPossibleLocalizedModifications()
        {
            // Create modifications in MetaMorpheus convention
            ModificationMotif.TryGetMotif("S", out var motifS);
            ModificationMotif.TryGetMotif("T", out var motifT);
            ModificationMotif.TryGetMotif("K", out var motifK);
            
            var phosphoMM = new Modification(
                _originalId: "Phosphorylation",
                _modificationType: "Common Biological",
                _target: motifS,
                _locationRestriction: "Anywhere.",
                _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));

            var phosphoTMM = new Modification(
                _originalId: "Phosphorylation",
                _modificationType: "Common Biological",
                _target: motifT,
                _locationRestriction: "Anywhere.",
                _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));
            
            var acetylMM = new Modification(
                _originalId: "Acetylation",
                _modificationType: "Common Biological",
                _target: motifK,
                _locationRestriction: "Anywhere.",
                _chemicalFormula: ChemicalFormula.ParseFormula("C2H2O1"));

            // Create protein with modifications
            var oneBasedMods = new Dictionary<int, List<Modification>>
            {
                { 3, new List<Modification> { phosphoMM } },      // S at position 3
                { 5, new List<Modification> { phosphoTMM } },     // T at position 5
                { 7, new List<Modification> { acetylMM } }        // K at position 7
            };

            var protein = new Protein(
                "MASATDKE",
                "TestProtein",
                oneBasedModifications: oneBasedMods);

            // Verify original mods are MetaMorpheus style
            Assert.That(protein.OneBasedPossibleLocalizedModifications[3][0].ModificationType, 
                Is.EqualTo("Common Biological"));
            Assert.That(protein.OneBasedPossibleLocalizedModifications[5][0].ModificationType, 
                Is.EqualTo("Common Biological"));
            Assert.That(protein.OneBasedPossibleLocalizedModifications[7][0].ModificationType, 
                Is.EqualTo("Common Biological"));

            // Convert to UniProt convention
            protein.ConvertMods(ModificationNamingConvention.UniProt);

            // Verify conversions
            Assert.That(protein.OneBasedPossibleLocalizedModifications[3][0].ModificationType, 
                Is.EqualTo("UniProt"));
            Assert.That(protein.OneBasedPossibleLocalizedModifications[5][0].ModificationType, 
                Is.EqualTo("UniProt"));
            Assert.That(protein.OneBasedPossibleLocalizedModifications[7][0].ModificationType, 
                Is.EqualTo("UniProt"));

            // Verify chemical formulas are preserved
            Assert.That(protein.OneBasedPossibleLocalizedModifications[3][0].ChemicalFormula.Equals(
                ChemicalFormula.ParseFormula("H1O3P1")), Is.True);
            Assert.That(protein.OneBasedPossibleLocalizedModifications[5][0].ChemicalFormula.Equals(
                ChemicalFormula.ParseFormula("H1O3P1")), Is.True);
            Assert.That(protein.OneBasedPossibleLocalizedModifications[7][0].ChemicalFormula.Equals(
                ChemicalFormula.ParseFormula("C2H2O1")), Is.True);
        }

        [Test]
        public static void TestConvertModsOnProtein_SequenceVariations()
        {
            // Create modifications
            ModificationMotif.TryGetMotif("S", out var motifS);
            
            var phosphoMM = new Modification(
                _originalId: "Phosphorylation",
                _modificationType: "Common Biological",
                _target: motifS,
                _locationRestriction: "Anywhere.",
                _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));

            // Create sequence variation with modification
            var variantMods = new Dictionary<int, List<Modification>>
            {
                { 1, new List<Modification> { phosphoMM } }  // Mod on the variant sequence
            };

            var sequenceVariation = new SequenceVariation(
                oneBasedBeginPosition: 3,
                oneBasedEndPosition: 3,
                originalSequence: "A",
                variantSequence: "S",
                description: "A3S",
                oneBasedModifications: variantMods);

            var protein = new Protein(
                "MAAADE",
                "TestProtein",
                sequenceVariations: new List<SequenceVariation> { sequenceVariation });

            // Apply the variation
            protein = protein.GetVariantBioPolymers().Skip(1).First();

            // Verify original mod is MetaMorpheus style
            Assert.That(protein.SequenceVariations.First().OneBasedModifications[1][0].ModificationType,
                Is.EqualTo("Common Biological"));

            // Convert to UniProt convention
            protein.ConvertMods(ModificationNamingConvention.UniProt);

            // Verify conversion in sequence variation
            Assert.That(protein.SequenceVariations.First().OneBasedModifications[1][0].ModificationType,
                Is.EqualTo("UniProt"));

            // Also verify it converted in AppliedSequenceVariations
            Assert.That(protein.AppliedSequenceVariations.First().OneBasedModifications[1][0].ModificationType,
                Is.EqualTo("UniProt"));
        }

        [Test]
        public static void TestConvertModsOnProtein_OriginalNonVariantModifications()
        {
            // This tests the OriginalNonVariantModifications dictionary conversion
            ModificationMotif.TryGetMotif("K", out var motifK);
            
            var acetylMM = new Modification(
                _originalId: "Acetylation",
                _modificationType: "Common Biological",
                _target: motifK,
                _locationRestriction: "Anywhere.",
                _chemicalFormula: ChemicalFormula.ParseFormula("C2H2O1"));

            // Create modification dictionaries
            var oneBasedMods = new Dictionary<int, List<Modification>>
            {
                { 4, new List<Modification> { acetylMM } }
            };

            // Apply a variant that doesn't affect the modification position
            var variant = new SequenceVariation(2, "A", "V", "A2V");
            var proteinWithVariant = new Protein(
                "MVP KDE",
                "TestProteinWithVariant",
                sequenceVariations: new List<SequenceVariation> { variant },
                oneBasedModifications: oneBasedMods);

            proteinWithVariant = proteinWithVariant.GetVariantBioPolymers().Skip(1).First();

            // The original mods should be stored in OriginalNonVariantModifications
            // Verify it's MetaMorpheus style initially
            if (proteinWithVariant.OriginalNonVariantModifications.Any())
            {
                Assert.That(proteinWithVariant.OriginalNonVariantModifications.First().Value[0].ModificationType,
                    Is.EqualTo("Common Biological"));
            }

            // Convert to UniProt
            proteinWithVariant.ConvertMods(ModificationNamingConvention.UniProt);

            // Verify conversion in OriginalNonVariantModifications if present
            if (proteinWithVariant.OriginalNonVariantModifications.Any())
            {
                Assert.That(proteinWithVariant.OriginalNonVariantModifications.First().Value[0].ModificationType,
                    Is.EqualTo("UniProt"));
            }
        }

        [Test]
        public static void TestConvertModsOnProteinWithMultipleModsPerSite()
        {
            // Test conversion when multiple modifications are possible at the same site
            ModificationMotif.TryGetMotif("S", out var motifS);
            
            var phosphoMM = new Modification(
                _originalId: "Phosphorylation",
                _modificationType: "Common Biological",
                _target: motifS,
                _locationRestriction: "Anywhere.",
                _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));

            var sulfoMM = new Modification(
                _originalId: "Sulfation",
                _modificationType: "Common Biological",
                _target: motifS,
                _locationRestriction: "Anywhere.",
                _chemicalFormula: ChemicalFormula.ParseFormula("O3S1"));

            // Multiple mods at same position
            var oneBasedMods = new Dictionary<int, List<Modification>>
            {
                { 3, new List<Modification> { phosphoMM, sulfoMM } }
            };

            var protein = new Protein(
                "MASIDE",
                "TestProtein",
                oneBasedModifications: oneBasedMods);

            // Verify both mods are MetaMorpheus style
            Assert.That(protein.OneBasedPossibleLocalizedModifications[3][0].ModificationType,
                Is.EqualTo("Common Biological"));
            Assert.That(protein.OneBasedPossibleLocalizedModifications[3][1].ModificationType,
                Is.EqualTo("Common Biological"));

            // Convert to UniProt
            protein.ConvertMods(ModificationNamingConvention.UniProt);

            // Verify both mods were converted
            Assert.That(protein.OneBasedPossibleLocalizedModifications[3][0].ModificationType,
                Is.EqualTo("UniProt"));
            Assert.That(protein.OneBasedPossibleLocalizedModifications[3][1].ModificationType,
                Is.EqualTo("UniProt"));
        }

        [Test]
        public static void TestProteinConversionRoundTrip()
        {
            // Test that converting back and forth preserves chemistry
            ModificationMotif.TryGetMotif("S", out var motifS);
            
            var phosphoMM = new Modification(
                _originalId: "Phosphorylation",
                _modificationType: "Common Biological",
                _target: motifS,
                _locationRestriction: "Anywhere.",
                _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));

            var oneBasedMods = new Dictionary<int, List<Modification>>
            {
                { 3, new List<Modification> { phosphoMM } }
            };

            var protein = new Protein(
                "MASIDE",
                "TestProtein",
                oneBasedModifications: oneBasedMods);

            var originalFormula = protein.OneBasedPossibleLocalizedModifications[3][0].ChemicalFormula;
            var originalTarget = protein.OneBasedPossibleLocalizedModifications[3][0].Target.ToString();

            // Convert to UniProt and back
            protein.ConvertMods(ModificationNamingConvention.UniProt);
            var uniprotFormula = protein.OneBasedPossibleLocalizedModifications[3][0].ChemicalFormula;
            
            protein.ConvertMods(ModificationNamingConvention.MetaMorpheus);
            var finalFormula = protein.OneBasedPossibleLocalizedModifications[3][0].ChemicalFormula;
            var finalTarget = protein.OneBasedPossibleLocalizedModifications[3][0].Target.ToString();

            // Verify chemistry is preserved
            Assert.That(originalFormula.Equals(uniprotFormula), Is.True);
            Assert.That(originalFormula.Equals(finalFormula), Is.True);
            Assert.That(originalTarget, Is.EqualTo(finalTarget));
        }

        [Test]
        public static void TestDatabaseWriteReadWithModificationConversion()
        {
            // Create dummy proteins with MetaMorpheus-style modifications
            ModificationMotif.TryGetMotif("S", out var motifS);
            ModificationMotif.TryGetMotif("T", out var motifT);
            ModificationMotif.TryGetMotif("K", out var motifK);
            
            var phosphoMM = new Modification(
                _originalId: "Phosphorylation",
                _modificationType: "Common Biological",
                _target: motifS,
                _locationRestriction: "Anywhere.",
                _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));

            var acetylMM = new Modification(
                _originalId: "Acetylation",
                _modificationType: "Common Biological",
                _target: motifK,
                _locationRestriction: "Anywhere.",
                _chemicalFormula: ChemicalFormula.ParseFormula("C2H2O1"));

            // Create proteins with mods
            var protein1Mods = new Dictionary<int, List<Modification>>
            {
                { 2, new List<Modification> { acetylMM } },
                { 5, new List<Modification> { phosphoMM } }
            };

            var protein2Mods = new Dictionary<int, List<Modification>>
            {
                { 3, new List<Modification> { phosphoMM } },
                { 7, new List<Modification> { acetylMM } }
            };

            var protein1 = new Protein(
                "MKPSIDE",
                "Protein1",
                oneBasedModifications: protein1Mods);

            var protein2 = new Protein(
                "MASKTDE",
                "Protein2",
                oneBasedModifications: protein2Mods);

            var proteins = new List<Protein> { protein1, protein2 };

            // Write MetaMorpheus proteins to XML
            string mmXmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "test_metamorpheus_mods.xml");
            ProteinDbWriter.WriteXmlDatabase(
                new Dictionary<string, HashSet<System.Tuple<int, Modification>>>(),
                proteins,
                mmXmlPath);

            // Read back MetaMorpheus proteins
            var mmProteins = ProteinDbLoader.LoadProteinXML(
                mmXmlPath,
                true,
                DecoyType.None,
                new List<Modification>(),
                false,
                new List<string>(),
                out var unknownMods1);

            // Collect all modification strings from MetaMorpheus file
            var mmModStrings = new List<string>();
            var mmFileLines = File.ReadAllLines(mmXmlPath);
            foreach (var line in mmFileLines)
            {
                if (line.Contains("modified residue"))
                {
                    mmModStrings.Add(line.Trim());
                }
            }

            // Convert proteins to UniProt naming convention
            foreach (var protein in proteins)
            {
                protein.ConvertMods(ModificationNamingConvention.UniProt);
            }

            // Write converted proteins to second XML file
            string uniprotXmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "test_uniprot_mods.xml");
            ProteinDbWriter.WriteXmlDatabase(
                new Dictionary<string, HashSet<System.Tuple<int, Modification>>>(),
                proteins,
                uniprotXmlPath);

            // Read back UniProt proteins
            var uniprotProteins = ProteinDbLoader.LoadProteinXML(
                uniprotXmlPath,
                true,
                DecoyType.None,
                new List<Modification>(),
                false,
                new List<string>(),
                out var unknownMods2);

            // Collect all modification strings from UniProt file
            var uniprotModStrings = new List<string>();
            var uniprotFileLines = File.ReadAllLines(uniprotXmlPath);
            foreach (var line in uniprotFileLines)
            {
                if (line.Contains("modified residue"))
                {
                    uniprotModStrings.Add(line.Trim());
                }
            }

            // Verify that modification naming conventions are different
            var mmHasPhosphorylation = mmModStrings.Any(s => s.Contains("Phosphorylation")) && !mmModStrings.Any(s => s.Contains("Phosphoserine"));
            var uniprotHasPhosphoserine = uniprotModStrings.Any(s => s.Contains("Phosphoserine")) && !uniprotModStrings.Any(s => s.Contains("Phosphorylation"));

            Assert.That(mmHasPhosphorylation, Is.True, "MetaMorpheus file should have 'Phosphorylation' mod type");
            Assert.That(uniprotHasPhosphoserine, Is.True, "UniProt file should have 'UniProt' mod type");

            // Verify both sets of proteins have the same number
            Assert.That(mmProteins.Count, Is.EqualTo(uniprotProteins.Count));

            // Verify chemical equivalence for each protein
            for (int i = 0; i < mmProteins.Count; i++)
            {
                // Same sequence
                Assert.That(mmProteins[i].BaseSequence, Is.EqualTo(uniprotProteins[i].BaseSequence));

                // Same modification positions
                var mmModPositions = mmProteins[i].OneBasedPossibleLocalizedModifications.Keys.OrderBy(k => k).ToList();
                var upModPositions = uniprotProteins[i].OneBasedPossibleLocalizedModifications.Keys.OrderBy(k => k).ToList();
                Assert.That(mmModPositions, Is.EqualTo(upModPositions));

                // Chemical formulas match at each position
                foreach (var position in mmModPositions)
                {
                    var mmMod = mmProteins[i].OneBasedPossibleLocalizedModifications[position][0];
                    var upMod = uniprotProteins[i].OneBasedPossibleLocalizedModifications[position][0];

                    if (mmMod.ChemicalFormula != null && upMod.ChemicalFormula != null)
                    {
                        Assert.That(mmMod.ChemicalFormula.Equals(upMod.ChemicalFormula), Is.True,
                            $"Chemical formulas should match for protein {i} at position {position}");
                    }

                    // Targets should match
                    Assert.That(mmMod.Target.ToString(), Is.EqualTo(upMod.Target.ToString()),
                        $"Targets should match for protein {i} at position {position}");
                }
            }

            // Verify the modification type strings are actually different
            var mmFirstModType = mmProteins[0].OneBasedPossibleLocalizedModifications.First().Value[0].ModificationType;
            var upFirstModType = uniprotProteins[0].OneBasedPossibleLocalizedModifications.First().Value[0].ModificationType;
            Assert.That(mmFirstModType, Is.Not.EqualTo(upFirstModType),
                "Modification types should be different between conventions");

            // Clean up
            File.Delete(mmXmlPath);
            File.Delete(uniprotXmlPath);
        }

        [Test]
        public static void TestDatabaseWriteReadWithSequenceVariationMods()
        {
            // Test that sequence variation mods survive write/read cycle
            ModificationMotif.TryGetMotif("S", out var motifS);
            
            var phosphoMM = new Modification(
                _originalId: "Phosphorylation",
                _modificationType: "Common Biological",
                _target: motifS,
                _locationRestriction: "Anywhere.",
                _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));

            var variantMods = new Dictionary<int, List<Modification>>
            {
                { 1, new List<Modification> { phosphoMM } }
            };

            var sequenceVariation = new SequenceVariation(
                oneBasedBeginPosition: 3,
                oneBasedEndPosition: 3,
                originalSequence: "A",
                variantSequence: "S",
                description: "A3S",
                oneBasedModifications: variantMods);

            var protein = new Protein(
                "MAAADE",
                "TestProteinWithVariant",
                sequenceVariations: new List<SequenceVariation> { sequenceVariation });

            var proteins = new List<Protein> { protein };

            // Write to XML
            string xmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "test_variant_mods.xml");
            ProteinDbWriter.WriteXmlDatabase(
                new Dictionary<string, HashSet<System.Tuple<int, Modification>>>(),
                proteins,
                xmlPath);

            // Read back
            var readProteins = ProteinDbLoader.LoadProteinXML(
                xmlPath,
                true,
                DecoyType.None,
                new List<Modification>(),
                false,
                new List<string>(),
                out var unknownMods);

            // Verify sequence variation and its mods were preserved
            Assert.That(readProteins[0].SequenceVariations.Count, Is.EqualTo(1));
            Assert.That(readProteins[0].SequenceVariations.First().OneBasedModifications.Count, Is.EqualTo(1));
            
            var readMod = readProteins[0].SequenceVariations.First().OneBasedModifications[1][0];
            Assert.That(readMod.ChemicalFormula.Equals(ChemicalFormula.ParseFormula("H1O3P1")), Is.True);

            // Now test conversion
            readProteins[0].ConvertMods(ModificationNamingConvention.UniProt);
            Assert.That(readProteins[0].SequenceVariations.First().OneBasedModifications[1][0].ModificationType,
                Is.EqualTo("UniProt"));

            // Clean up
            File.Delete(xmlPath);
        }

        [Test]
        public static void TestDatabaseWriteReadPreservesModificationDetails()
        {
            // Test that all mod details (neutral losses, diagnostic ions) are preserved
            ModificationMotif.TryGetMotif("S", out var motifS);
            
            var phosphoWithNL = new Modification(
                _originalId: "Phosphorylation",
                _modificationType: "Common Biological",
                _target: motifS,
                _locationRestriction: "Anywhere.",
                _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"),
                _neutralLosses: new Dictionary<MassSpectrometry.DissociationType, List<double>>
                {
                    { MassSpectrometry.DissociationType.HCD, new List<double> { 97.976896 } }
                });

            var oneBasedMods = new Dictionary<int, List<Modification>>
            {
                { 3, new List<Modification> { phosphoWithNL } }
            };

            var protein = new Protein(
                "MASIDE",
                "TestProteinWithNL",
                oneBasedModifications: oneBasedMods);

            var proteins = new List<Protein> { protein };

            // Write to XML
            string xmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "test_nl_mods.xml");
            ProteinDbWriter.WriteXmlDatabase(
                new Dictionary<string, HashSet<System.Tuple<int, Modification>>>(),
                proteins,
                xmlPath);

            // Read back
            var readProteins = ProteinDbLoader.LoadProteinXML(
                xmlPath,
                true,
                DecoyType.None,
                new List<Modification>(),
                false,
                new List<string>(),
                out var unknownMods);

            // Verify neutral losses were preserved
            var readMod = readProteins[0].OneBasedPossibleLocalizedModifications[3][0];
            Assert.That(readMod.NeutralLosses, Is.Not.Null);
            Assert.That(readMod.NeutralLosses.ContainsKey(MassSpectrometry.DissociationType.HCD), Is.True);

            readProteins.ForEach( p => p.ConvertMods(ModificationNamingConvention.UniProt) );
            Assert.That(readProteins[0].OneBasedPossibleLocalizedModifications[3][0].ModificationType,
                Is.EqualTo("UniProt"));

            // Clean up
            File.Delete(xmlPath);
        }

        [Test]
        public static void TestConversionWithEmptyProtein()
        {
            // Test that conversion works on protein with no modifications
            var protein = new Protein("MAPSIDE", "TestProteinNoMods");

            // Should not throw
            Assert.DoesNotThrow(() => protein.ConvertMods(ModificationNamingConvention.UniProt));
            Assert.DoesNotThrow(() => protein.ConvertMods(ModificationNamingConvention.MetaMorpheus));
        }

        [Test]
        public static void Dummy()
        {
            string inPath = @"B:\Users\Nic\SharedWithMe\ClaireMulti_PrunedDB\1-5-25-DBs";
            string db = @"B:\Users\Nic\SharedWithMe\ClaireMulti_PrunedDB\1-5-25-DBs\uniprotkb_taxonomy_id_9606_AND_reviewed_2025_11_22-HIV_WTJB474_updated_JustGag-JustPolGPTMDproteinPruned-noCitCarbox.xml";

            //foreach (var db in Directory.GetFiles(inPath, "*.xml"))
            //{
                //if (db.EndsWith("_uniprot.xml")) 
                //    continue;

                var inProteins = ProteinDbLoader.LoadProteinXML(
                    db,
                    true,
                    DecoyType.None,
                    new List<Modification>(),
                    false,
                    new List<string>(),
                    out var unknownMods);
                
                foreach (var prot in inProteins)
                {
                    prot.ConvertMods(ModificationNamingConvention.UniProt);
                }

                var oneBased = inProteins.SelectMany(p => p.OneBasedPossibleLocalizedModifications.Values).SelectMany(m => m).Distinct().ToList();
                var original = inProteins.SelectMany(p => p.OriginalNonVariantModifications.Values).SelectMany(m => m).Distinct().ToList();
                var seqVars = inProteins.SelectMany(p => p.SequenceVariations.SelectMany(p => p.OneBasedModifications.Values).SelectMany(m => m)).Distinct().ToList();
                var appliedVars = inProteins.SelectMany(p => p.AppliedSequenceVariations.SelectMany(p => p.OneBasedModifications.Values).SelectMany(m => m)).Distinct().ToList();

                string outPath = Path.Combine(Path.GetDirectoryName(db)!, Path.GetFileNameWithoutExtension(db) + "_uniprot.xml");
                ProteinDbWriter.WriteXmlDatabase(
                    new Dictionary<string, HashSet<System.Tuple<int, Modification>>>(),
                    inProteins,
                    outPath, namingConvention: ModificationNamingConvention.UniProt);
            //}

            var cache = ModificationConverter.ExtractCache(ModificationNamingConvention.UniProt);
        }
    }
}
