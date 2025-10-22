using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Xml;
using Omics;
using UsefulProteomicsDatabases;
using Omics.BioPolymer;
using Transcriptomics;

namespace Test.DatabaseTests.VariantTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    internal class TestProteinXmlWriteVariants
    {
        private static Modification NewMod(string originalId)
        {
            ModificationMotif.TryGetMotif("X", out var motifAny);
            // IdWithMotif will be computed internally as "<OriginalId> on X" when appropriate for this codebase
            var m = new Modification(
                _originalId: originalId,
                _accession: null,
                _modificationType: "mt",
                _featureType: null,
                _target: motifAny,
                _locationRestriction: "Anywhere.",
                _chemicalFormula: null,
                _monoisotopicMass: 1,
                _databaseReference: null,
                _taxonomicRange: null,
                _keywords: null,
                _neutralLosses: null,
                _diagnosticIons: null,
                _fileOrigin: null);
            return m;
        }

        private static Protein BuildConsensusProtein(out SequenceVariation sv, out Modification baseA, out Modification baseZ, out Modification svMod)
        {
            // Base sequence ACDE; variant D3->E (point substitution)
            baseA = NewMod("A1 on X");
            baseZ = NewMod("Z9 on X");
            svMod = NewMod("VarMod on X");

            var baseMods = new Dictionary<int, List<Modification>>
            {
                { 2, new List<Modification> { baseZ, baseA } } // Intentional order to verify sorting
            };

            sv = new SequenceVariation(
                oneBasedBeginPosition: 3,
                oneBasedEndPosition: 3,
                originalSequence: "D",
                variantSequence: "E",
                description: null,
                variantCallFormatDataString: null,
                oneBasedModifications: new Dictionary<int, List<Modification>>
                {
                    { 3, new List<Modification> { svMod } }
                });

            // Unsorted DatabaseReference properties; writer should sort by type then value
            var dbRef = new DatabaseReference(
                "Xref",
                "ID",
                new[]
                {
                    Tuple.Create("z", "2"),
                    Tuple.Create("a", "2"),
                    Tuple.Create("a", "1")
                });

            var prot = new Protein(
                sequence: "ACDE",
                accession: "PBASE",
                organism: "Org",
                geneNames: new List<Tuple<string, string>> { Tuple.Create("primary", "GENE") },
                oneBasedModifications: baseMods,
                proteolysisProducts: null,
                name: "Name",
                fullName: "Full",
                isDecoy: false,
                isContaminant: false,
                databaseReferences: new List<DatabaseReference> { dbRef },
                sequenceVariations: new List<SequenceVariation> { sv },
                appliedSequenceVariations: null,
                sampleNameForVariants: null,
                disulfideBonds: null,
                spliceSites: null,
                databaseFilePath: null);

            return prot;
        }

        private static Protein BuildAppliedVariantProtein(Protein consensus, SequenceVariation sv)
        {
            // Apply the variant (D3->E): ACDE -> ACEE
            var applied = new Protein(
                variantBaseSequence: "ACEE",
                protein: consensus,
                appliedSequenceVariations: new[] { sv },
                applicableProteolysisProducts: null,
                oneBasedModifications: new Dictionary<int, List<Modification>>(), // no extra base mods
                sampleNameForVariants: "sampleX");

            return applied;
        }

        private static XmlDocument LoadXml(string path)
        {
            var doc = new XmlDocument();
            doc.Load(path);
            return doc;
        }

        private static XmlElement FindEntryByAccession(XmlDocument doc, string accession)
        {
            foreach (XmlElement entry in doc.GetElementsByTagName("entry"))
            {
                var acc = entry.GetElementsByTagName("accession").OfType<XmlElement>().FirstOrDefault();
                if (acc != null && string.Equals(acc.InnerText, accession, StringComparison.Ordinal))
                {
                    return entry;
                }
            }
            return null;
        }

        [Test]
        public void ProteinXml_AppliedVariantEntries_And_ModCatalog_And_Sorting()
        {
            // Arrange consensus + applied
            var consensus = BuildConsensusProtein(out var sv, out var baseA, out var baseZ, out var svMod);
            var applied = BuildAppliedVariantProtein(consensus, sv);

            // Additional mods: 2 new at positions 1 and 4 (counted twice), and 1 duplicate of base at pos 2 (not counted)
            var extraNew = NewMod("ExtraMod on X");
            var extraDup = NewMod("A1 on X"); // duplicate id; should not increment NewModResEntries

            // Variant-specific additional mod keyed to the applied accession
            var varExtra = NewMod("VarExtra on X");

            string outPath = Path.Combine(TestContext.CurrentContext.WorkDirectory, "prot_variant_write.xml");

            // includeAppliedVariantEntries = true ? both entries written
            var additional = new Dictionary<string, HashSet<Tuple<int, Modification>>>(StringComparer.Ordinal)
            {
                {
                    consensus.Accession,
                    new HashSet<Tuple<int, Modification>>
                    {
                        Tuple.Create(1, extraNew),
                        Tuple.Create(4, extraNew),
                        Tuple.Create(2, extraDup)
                    }
                },
                {
                    applied.Accession,
                    new HashSet<Tuple<int, Modification>>
                    {
                        Tuple.Create(3, varExtra)
                    }
                }
            };

            try
            {
                // Act
                var newCounts = ProteinDbWriter.WriteXmlDatabase(
                    additionalModsToAddToProteins: additional,
                    proteinList: new List<Protein> { consensus, applied },
                    outputFileName: outPath,
                    updateTimeStamp: true,
                    includeAppliedVariantEntries: true,
                    includeAppliedVariantFeatures: true);

                // Assert: file created
                Assert.That(File.Exists(outPath), Is.True);

                // Assert: NewModResEntries counts (2 new positions on base + 1 variant-extra on applied)
                Assert.That(newCounts, Contains.Key("ExtraMod on X"));
                Assert.That(newCounts["ExtraMod on X"], Is.EqualTo(2), "ExtraMod should be counted for two new positions on base accession.");
                Assert.That(newCounts, Contains.Key("VarExtra on X"));
                Assert.That(newCounts["VarExtra on X"], Is.EqualTo(1), "VarExtra should be counted once on the applied accession.");
                Assert.That(newCounts.ContainsKey("A1 on X"), Is.False);

                // Parse
                var doc = LoadXml(outPath);

                // Two entries expected: base + applied
                var baseEntry = FindEntryByAccession(doc, consensus.Accession);
                var varEntry = FindEntryByAccession(doc, applied.Accession);
                Assert.That(baseEntry, Is.Not.Null, "Base entry not found.");
                Assert.That(varEntry, Is.Not.Null, "Applied variant entry not found.");

                // Applied entry should be annotated and have updated modified date
                Assert.That(varEntry.HasAttribute("variant"), Is.True, "Applied entry missing 'variant' attribute.");
                Assert.That(varEntry.GetAttribute("variant"), Is.EqualTo("true"));
                var modifiedAttr = varEntry.GetAttribute("modified");
                Assert.That(modifiedAttr, Does.Match(@"^\d{4}-\d{2}-\d{2}$"), "Modified date missing/invalid.");

                // Base entry: candidate "sequence variant" features present
                var baseSeqVarFeatures = baseEntry.GetElementsByTagName("feature").OfType<XmlElement>()
                    .Where(e => e.HasAttribute("type") && e.GetAttribute("type") == "sequence variant")
                    .ToList();
                Assert.That(baseSeqVarFeatures.Count, Is.GreaterThanOrEqualTo(1), "Expected candidate sequence variant feature(s) on base entry.");

                // Applied entry: no sequence variant features should be written
                var appliedSeqVarFeatures = varEntry.GetElementsByTagName("feature").OfType<XmlElement>()
                    .Where(e => e.HasAttribute("type") && e.GetAttribute("type") == "sequence variant")
                    .ToList();
                Assert.That(appliedSeqVarFeatures.Count, Is.EqualTo(0), "Applied entries must not contain sequence variant features.");

                // Variant-specific subfeatures exist under base entry's variant features (svMod at pos3)
                var baseAnySubfeatureMod = baseSeqVarFeatures
                    .SelectMany(f => f.GetElementsByTagName("subfeature").OfType<XmlElement>())
                    .Any(sf => sf.HasAttribute("type") && sf.GetAttribute("type") == "modified residue");
                Assert.That(baseAnySubfeatureMod, Is.True, "Expected variant-specific modified residue subfeature(s) on base entry.");

                // Base entry: base mods + additional mods features exist; mod IDs at same position sorted lexicographically
                var baseFeatures = baseEntry.GetElementsByTagName("feature").OfType<XmlElement>()
                    .Where(e => e.HasAttribute("type") && e.GetAttribute("type") == "modified residue")
                    .ToList();
                Assert.That(baseFeatures.Count, Is.GreaterThanOrEqualTo(3), "Expected at least 3 modified residue features (2 base at pos2 + extras).");

                // Extract modified residue descriptions for position 2 and validate order
                var pos2ModDescs = baseEntry
                    .GetElementsByTagName("feature").OfType<XmlElement>()
                    .Where(e => e.HasAttribute("type") && e.GetAttribute("type") == "modified residue")
                    .Where(e =>
                    {
                        var loc = e.GetElementsByTagName("location").OfType<XmlElement>().FirstOrDefault();
                        var pos = loc?.GetElementsByTagName("position").OfType<XmlElement>().FirstOrDefault();
                        return pos?.GetAttribute("position") == "2";
                    })
                    .Select(e => e.GetAttribute("description"))
                    .ToList();

                var expectedOrder = new[] { "A1 on X", "Z9 on X" };
                Assert.That(pos2ModDescs, Is.EquivalentTo(expectedOrder));
                Assert.That(pos2ModDescs, Is.EqualTo(expectedOrder), "Mod IDs at same position should be ordered lexicographically.");

                // DatabaseReference property sorting: expect ("a","1"), ("a","2"), ("z","2")
                var dbRef = baseEntry.GetElementsByTagName("dbReference").OfType<XmlElement>().FirstOrDefault(e => e.HasAttribute("type") && e.GetAttribute("type") == "Xref");
                Assert.That(dbRef, Is.Not.Null, "dbReference 'Xref' not found.");
                var props = dbRef!.GetElementsByTagName("property").OfType<XmlElement>()
                    .Select(p => (type: p.GetAttribute("type"), value: p.GetAttribute("value")))
                    .ToList();
                Assert.That(props.Count, Is.EqualTo(3));
                Assert.That(props[0], Is.EqualTo(("a", "1")));
                Assert.That(props[1], Is.EqualTo(("a", "2")));
                Assert.That(props[2], Is.EqualTo(("z", "2")));

                // Modification catalog: baseA, baseZ, svMod, extraNew, varExtra
                var modCatalog = doc.GetElementsByTagName("modification").OfType<XmlElement>().ToList();
                var expectedUnique = new HashSet<string>(StringComparer.Ordinal)
                {
                    baseA.IdWithMotif, baseZ.IdWithMotif, svMod.IdWithMotif, extraNew.IdWithMotif, varExtra.IdWithMotif
                };
                Assert.That(modCatalog.Count, Is.EqualTo(expectedUnique.Count), "Modification catalog unique count mismatch.");

                // Global entry ordering by accession (ascending)
                var entryAccOrder = doc.GetElementsByTagName("entry").OfType<XmlElement>()
                    .Select(e => e.GetElementsByTagName("accession").OfType<XmlElement>().First().InnerText)
                    .ToList();
                var sorted = entryAccOrder.OrderBy(a => a, StringComparer.Ordinal).ToList();
                Assert.That(entryAccOrder, Is.EqualTo(sorted), "Entries should be ordered by accession.");
            }
            finally
            {
                if (File.Exists(outPath))
                    File.Delete(outPath);
            }
        }

        [Test]
        public void ProteinXml_AppliedVariantFeatures_Toggle()
        {
            // Arrange
            var consensus = BuildConsensusProtein(out var sv, out _, out _, out _);
            var applied = BuildAppliedVariantProtein(consensus, sv);

            string outPathTrue = Path.Combine(TestContext.CurrentContext.WorkDirectory, "prot_variant_features_true.xml");
            string outPathFalse = Path.Combine(TestContext.CurrentContext.WorkDirectory, "prot_variant_features_false.xml");

            try
            {
                // includeAppliedVariantFeatures = true
                ProteinDbWriter.WriteXmlDatabase(
                    additionalModsToAddToProteins: new Dictionary<string, HashSet<Tuple<int, Modification>>>(),
                    proteinList: new List<Protein> { consensus, applied },
                    outputFileName: outPathTrue,
                    updateTimeStamp: false,
                    includeAppliedVariantEntries: true,
                    includeAppliedVariantFeatures: true);

                // includeAppliedVariantFeatures = false
                ProteinDbWriter.WriteXmlDatabase(
                    additionalModsToAddToProteins: new Dictionary<string, HashSet<Tuple<int, Modification>>>(),
                    proteinList: new List<Protein> { consensus, applied },
                    outputFileName: outPathFalse,
                    updateTimeStamp: false,
                    includeAppliedVariantEntries: true,
                    includeAppliedVariantFeatures: false);

                var docTrue = LoadXml(outPathTrue);
                var docFalse = LoadXml(outPathFalse);

                var baseEntryTrue = FindEntryByAccession(docTrue, consensus.Accession);
                var varEntryTrue = FindEntryByAccession(docTrue, applied.Accession);
                var baseEntryFalse = FindEntryByAccession(docFalse, consensus.Accession);
                var varEntryFalse = FindEntryByAccession(docFalse, applied.Accession);

                Assert.That(baseEntryTrue, Is.Not.Null);
                Assert.That(varEntryTrue, Is.Not.Null);
                Assert.That(baseEntryFalse, Is.Not.Null);
                Assert.That(varEntryFalse, Is.Not.Null);

                // True ? base has sequence variant features; applied has none
                var baseFeaturesTrue = baseEntryTrue!.GetElementsByTagName("feature").OfType<XmlElement>()
                    .Where(e => e.HasAttribute("type") && e.GetAttribute("type") == "sequence variant")
                    .ToList();
                Assert.That(baseFeaturesTrue.Count, Is.GreaterThanOrEqualTo(1), "Expected sequence variant features on consensus when enabled.");

                var appliedFeaturesTrue = varEntryTrue!.GetElementsByTagName("feature").OfType<XmlElement>()
                    .Where(e => e.HasAttribute("type") && e.GetAttribute("type") == "sequence variant")
                    .ToList();
                Assert.That(appliedFeaturesTrue.Count, Is.EqualTo(0), "Applied entries must not contain sequence variant features.");

                // False ? no sequence variant features anywhere
                var baseFeaturesFalse = baseEntryFalse!.GetElementsByTagName("feature").OfType<XmlElement>()
                    .Where(e => e.HasAttribute("type") && e.GetAttribute("type") == "sequence variant")
                    .ToList();
                Assert.That(baseFeaturesFalse.Count, Is.EqualTo(0), "Consensus entry should not have sequence variant features when disabled.");

                var appliedFeaturesFalse = varEntryFalse!.GetElementsByTagName("feature").OfType<XmlElement>()
                    .Where(e => e.HasAttribute("type") && e.GetAttribute("type") == "sequence variant")
                    .ToList();
                Assert.That(appliedFeaturesFalse.Count, Is.EqualTo(0), "Applied entries must not contain sequence variant features.");
            }
            finally
            {
                if (File.Exists(outPathTrue)) File.Delete(outPathTrue);
                if (File.Exists(outPathFalse)) File.Delete(outPathFalse);
            }
        }

        [Test]
        public void ProteinXml_AdditionalMods_NewCounts_And_Catalog_Filter_When_No_Applied_Entries()
        {
            // Arrange
            var consensus = BuildConsensusProtein(out var sv, out _, out _, out var svMod);
            var applied = BuildAppliedVariantProtein(consensus, sv);

            var extraNew = NewMod("ExtraMod on X");
            var varExtra = NewMod("VarExtra on X");

            string outPath = Path.Combine(TestContext.CurrentContext.WorkDirectory, "prot_variant_no_applied.xml");

            // includeAppliedVariantEntries = false ? applied entry not written
            var additional = new Dictionary<string, HashSet<Tuple<int, Modification>>>(StringComparer.Ordinal)
            {
                { consensus.Accession, new HashSet<Tuple<int, Modification>> { Tuple.Create(1, extraNew) } },
                { applied.Accession,   new HashSet<Tuple<int, Modification>> { Tuple.Create(3, varExtra) } } // should be ignored entirely
            };

            try
            {
                var counts = ProteinDbWriter.WriteXmlDatabase(
                    additionalModsToAddToProteins: additional,
                    proteinList: new List<Protein> { consensus, applied },
                    outputFileName: outPath,
                    updateTimeStamp: false,
                    includeAppliedVariantEntries: false,
                    includeAppliedVariantFeatures: true);

                // Assert: counts only reflect the base accession addition; variant-keyed additional mod not counted
                Assert.That(counts, Contains.Key("ExtraMod on X"));
                Assert.That(counts["ExtraMod on X"], Is.EqualTo(1));
                Assert.That(counts.ContainsKey("VarExtra on X"), Is.False, "Variant-keyed additional mod should not be counted when applied entries are not written.");

                var doc = LoadXml(outPath);

                // Only base entry present
                Assert.That(FindEntryByAccession(doc, consensus.Accession), Is.Not.Null);
                Assert.That(FindEntryByAccession(doc, applied.Accession), Is.Null);

                // Modification catalog should include: base mods + candidate variant mod + base additional; not variant-keyed additional
                var modCatalog = doc.GetElementsByTagName("modification").OfType<XmlElement>().ToList();
                Assert.That(modCatalog.Count, Is.EqualTo(4), "Catalog should exclude variant-keyed additional mod when applied entries are not written.");
            }
            finally
            {
                if (File.Exists(outPath)) File.Delete(outPath);
            }
        }

        [Test]
        public void WriteXmlDatabase_Dispatch_By_IBioPolymer_For_Protein_And_RNA()
        {
            // Protein path (use concrete overload)
            var consensus = BuildConsensusProtein(out _, out _, out _, out _);
            string outProt = Path.Combine(TestContext.CurrentContext.WorkDirectory, "dispatch_protein.xml");
            string outRna = Path.Combine(TestContext.CurrentContext.WorkDirectory, "dispatch_rna.xml");
            try
            {
                var retProt = ProteinDbWriter.WriteXmlDatabase(
                    additionalModsToAddToProteins: new Dictionary<string, HashSet<Tuple<int, Modification>>>(),
                    proteinList: new List<Protein> { consensus },
                    outputFileName: outProt);
                Assert.That(File.Exists(outProt), Is.True);
                Assert.That(retProt, Is.Not.Null);
            }
            finally
            {
                if (File.Exists(outProt)) File.Delete(outProt);
            }

            // RNA path (use concrete overload)
            var rna = new RNA(
                sequence: "AUGC",
                accession: "RNA001",
                oneBasedPossibleModifications: null,
                fivePrimeTerminus: null,
                threePrimeTerminus: null,
                name: "rna1",
                organism: "org",
                databaseFilePath: null,
                isContaminant: false,
                isDecoy: false,
                geneNames: new List<Tuple<string, string>> { Tuple.Create("primary", "GENE1") },
                databaseAdditionalFields: null,
                truncationProducts: null,
                sequenceVariations: new List<SequenceVariation>
                {
                    new SequenceVariation(oneBasedBeginPosition: 2, oneBasedEndPosition: 2, originalSequence: "U", variantSequence: "C", description: null, variantCallFormatDataString: null, oneBasedModifications: null)
                },
                appliedSequenceVariations: null,
                sampleNameForVariants: null,
                fullName: "full");

            try
            {
                var retRna = ProteinDbWriter.WriteXmlDatabase(
                    additionalModsToAddToNucleicAcids: new Dictionary<string, HashSet<Tuple<int, Modification>>>(),
                    nucleicAcidList: new List<RNA> { rna },
                    outputFileName: outRna);
                Assert.That(File.Exists(outRna), Is.True);
                Assert.That(retRna, Is.Not.Null);
            }
            finally
            {
                if (File.Exists(outRna)) File.Delete(outRna);
            }
        }
    }
}