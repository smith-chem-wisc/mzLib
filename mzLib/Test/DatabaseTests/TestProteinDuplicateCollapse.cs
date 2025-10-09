using NUnit.Framework;
using Proteomics;
using UsefulProteomicsDatabases;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Xml;
using Omics.Modifications;
using Omics.BioPolymer;

namespace Test.DatabaseTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    internal class TestProteinDuplicateCollapse
    {
        private static Modification NewMod(string originalId)
        {
            ModificationMotif.TryGetMotif("X", out var motifAny);
            return new Modification(
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
        }

        private static Protein BuildConsensusProtein(out SequenceVariation sv, out Modification baseMod)
        {
            // Base: ACDE; variant D3->E
            baseMod = NewMod("BaseMod on X");
            var baseMods = new Dictionary<int, List<Modification>>
            {
                { 1, new List<Modification> { baseMod } }
            };

            sv = new SequenceVariation(
                oneBasedBeginPosition: 3,
                oneBasedEndPosition: 3,
                originalSequence: "D",
                variantSequence: "E",
                description: "D3E",
                variantCallFormatDataString: null,
                oneBasedModifications: new Dictionary<int, List<Modification>>
                {
                    { 3, new List<Modification> { NewMod("VarMod on X") } }
                });

            return new Protein(
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
                databaseReferences: null,
                sequenceVariations: new List<SequenceVariation> { sv },
                appliedSequenceVariations: null,
                sampleNameForVariants: null,
                disulfideBonds: null,
                spliceSites: null,
                databaseFilePath: null);
        }

        private static Protein BuildAppliedVariantProtein(Protein consensus, SequenceVariation sv, out Modification appliedOnlyMod)
        {
            // Apply variant D3->E => ACEE; add an applied-only mod at pos 4 to be merged
            appliedOnlyMod = NewMod("AppliedOnly on X");
            var appliedMods = new Dictionary<int, List<Modification>>
            {
                { 4, new List<Modification> { appliedOnlyMod } }
            };

            return new Protein(
                variantBaseSequence: "ACEE",
                protein: consensus,
                appliedSequenceVariations: new[] { sv },
                applicableProteolysisProducts: null,
                oneBasedModifications: appliedMods,
                sampleNameForVariants: "sampleX");
        }

        private static T InvokeInternalStatic<T>(Type type, string method, params object[] args)
        {
            var mi = type.GetMethod(method, BindingFlags.NonPublic | BindingFlags.Static);
            Assert.That(mi, Is.Not.Null, $"Internal method {type.Name}.{method} not found.");
            return (T)mi.Invoke(null, args);
        }
        [Test]
        public void Loader_Collapses_Duplicate_AppliedVariant_FromConsensusExpansion()
        {
            // Arrange: build consensus + a pre-existing applied entry (same accession/sequence as expansion will generate)
            var consensus = BuildConsensusProtein(out var sv, out var baseMod);
            var applied = BuildAppliedVariantProtein(consensus, sv, out var appliedOnlyMod);

            string outPath = Path.Combine(TestContext.CurrentContext.WorkDirectory, "dup_collapse.xml");

            try
            {
                // Write both entries
                ProteinDbWriter.WriteXmlDatabase(
                    additionalModsToAddToProteins: new Dictionary<string, HashSet<Tuple<int, Modification>>>(),
                    proteinList: new List<Protein> { consensus, applied },
                    outputFileName: outPath,
                    updateTimeStamp: false,
                    includeAppliedVariantEntries: true,
                    includeAppliedVariantFeatures: true);

                // Act: read and expand variants (LoadProteinXML auto-collapses duplicates)
                var options = new ProteinDbLoader.ProteinXmlLoadOptions
                {
                    GenerateTargets = true,
                    DecoyType = DecoyType.None,
                    AllKnownModifications = Array.Empty<Modification>(),
                    IsContaminant = false,
                    ModTypesToExclude = Array.Empty<string>(),
                    MaxThreads = -1,
                    MaxSequenceVariantsPerIsoform = 4,
                    MinAlleleDepth = 1,
                    MaxSequenceVariantIsoforms = 1,
                    AddTruncations = false,
                    DecoyIdentifier = "DECOY"
                };
                _ = ProteinDbLoader.LoadProteinXML(outPath, options, out var unknownMods);

                // Assert: unknowns empty
                Assert.That(unknownMods, Is.Not.Null);
                Assert.That(unknownMods.Count, Is.EqualTo(0));

                // Re-read with GenerateTargets true to get the expanded list (again)
                var proteins = ProteinDbLoader.LoadProteinXML(outPath, options, out _);

                // There should be exactly one applied entry with the variant accession (duplicate collapsed)
                var appliedAccession = applied.Accession;
                var applieds = proteins.Where(p => p.Accession == appliedAccession && p.BaseSequence == applied.BaseSequence).ToList();
                Assert.That(applieds.Count, Is.EqualTo(1), "Duplicate applied entries should be collapsed.");

                var mergedApplied = applieds[0];

                // Applied entry should NOT inherit consensus base mods; it should include applied-only mods
                var pos1HasBaseOnApplied = mergedApplied.OneBasedPossibleLocalizedModifications.TryGetValue(1, out var a1)
                                           && a1.Any(m => string.Equals(m.IdWithMotif, baseMod.IdWithMotif, StringComparison.Ordinal));
                var pos4HasAppliedOnly = mergedApplied.OneBasedPossibleLocalizedModifications.TryGetValue(4, out var a4)
                                           && a4.Any(m => string.Equals(m.IdWithMotif, appliedOnlyMod.IdWithMotif, StringComparison.Ordinal));
                Assert.That(pos1HasBaseOnApplied, Is.False, "Applied entry should not inherit base mods from consensus.");
                Assert.That(pos4HasAppliedOnly, Is.True, "Merged applied entry should include applied-only mod from prewritten applied.");

                // Consensus entry should still contain its base mod
                var consensusEntry = proteins.FirstOrDefault(p => p.Accession == consensus.Accession && p.BaseSequence == consensus.BaseSequence);
                Assert.That(consensusEntry, Is.Not.Null, "Consensus entry was not found after load/expand.");
                var pos1HasBaseOnConsensus = consensusEntry!.OneBasedPossibleLocalizedModifications.TryGetValue(1, out var c1)
                                             && c1.Any(m => string.Equals(m.IdWithMotif, baseMod.IdWithMotif, StringComparison.Ordinal));
                Assert.That(pos1HasBaseOnConsensus, Is.True, "Consensus entry should retain its base modification(s).");

                // Applied proteoform identity and variant application should be reflected in accession and base sequence
                Assert.That(mergedApplied.Accession, Is.EqualTo(applied.Accession), "Applied accession should be preserved.");
                Assert.That(mergedApplied.BaseSequence, Is.EqualTo("ACEE"), "Applied base sequence should reflect the applied variant.");
            }
            finally
            {
                if (File.Exists(outPath)) File.Delete(outPath);
            }
        }
        [Test]
        public void Internal_FindDuplicateGroups_Discovers_Duplicates_By_Accession_And_BaseSequence()
        {
            var consensus = BuildConsensusProtein(out var sv, out _);
            var appliedA = BuildAppliedVariantProtein(consensus, sv, out _);
            // Create a synthetic duplicate applied with same accession/base sequence (no mods)
            var appliedB = new Protein(
                sequence: appliedA.BaseSequence,
                accession: appliedA.Accession,
                organism: appliedA.Organism,
                geneNames: new List<Tuple<string, string>>(appliedA.GeneNames),
                oneBasedModifications: new Dictionary<int, List<Modification>>(),
                proteolysisProducts: null,
                name: appliedA.Name,
                fullName: appliedA.FullName,
                isDecoy: appliedA.IsDecoy,
                isContaminant: appliedA.IsContaminant,
                databaseReferences: new List<DatabaseReference>(appliedA.DatabaseReferences),
                sequenceVariations: new List<SequenceVariation>(),
                appliedSequenceVariations: new List<SequenceVariation>(appliedA.AppliedSequenceVariations),
                sampleNameForVariants: appliedA.SampleNameForVariants,
                disulfideBonds: new List<DisulfideBond>(appliedA.DisulfideBonds),
                spliceSites: new List<SpliceSite>(appliedA.SpliceSites),
                databaseFilePath: appliedA.DatabaseFilePath);

            var proteins = new List<Protein> { consensus, appliedA, appliedB };

            var groups = InvokeInternalStatic<IEnumerable<IGrouping<(string accession, string baseSequence), Protein>>>(
                typeof(ProteinDbLoader),
                "FindDuplicateGroupsByAccessionAndBaseSequence",
                proteins);

            var dupGroup = groups.FirstOrDefault(g => g.Key.accession == appliedA.Accession && g.Key.baseSequence == appliedA.BaseSequence);
            Assert.That(dupGroup, Is.Not.Null);
            Assert.That(dupGroup.Count(), Is.EqualTo(2));
        }
        [Test]
        public void Internal_Collapse_Merges_Unique_Mods_And_DeDuplicates()
        {
            var consensus = BuildConsensusProtein(out var sv, out var baseMod);
            var appliedA = BuildAppliedVariantProtein(consensus, sv, out var appliedOnlyA);
            var appliedB = BuildAppliedVariantProtein(consensus, sv, out var appliedOnlyB);

            // Put different unique mods in A and B at different positions; also duplicate one id in both
            var common = NewMod("Common on X");
            appliedA.OneBasedPossibleLocalizedModifications[2] = new List<Modification> { common };
            appliedB.OneBasedPossibleLocalizedModifications[2] = new List<Modification> { common };

            // Use a valid position within ACEE (length 4); previously used 5 which is invalid and gets filtered out
            appliedB.OneBasedPossibleLocalizedModifications[1] = new List<Modification> { NewMod("BOnly on X") };

            var collapsed = InvokeInternalStatic<List<Protein>>(
                typeof(ProteinDbLoader),
                "CollapseDuplicateProteinsByAccessionAndBaseSequence",
                new List<Protein> { consensus, appliedA, appliedB });

            // Exactly one applied in collapsed set
            var merged = collapsed.Where(p => p.Accession == appliedA.Accession && p.BaseSequence == appliedA.BaseSequence).Single();

            // Check union across applied duplicates:
            // - pos4 from appliedA
            // - pos1 from appliedB (valid position within ACEE)
            // - pos2 common de-duplicated
            Assert.That(merged.OneBasedPossibleLocalizedModifications.ContainsKey(4), "Missing applied-only A mod position 4.");
            Assert.That(merged.OneBasedPossibleLocalizedModifications.ContainsKey(1), "Missing applied-only B mod position 1.");
            Assert.That(merged.OneBasedPossibleLocalizedModifications.ContainsKey(2), "Missing common mod position 2.");
            var commons = merged.OneBasedPossibleLocalizedModifications[2].Where(m => m.IdWithMotif == common.IdWithMotif).ToList();
            Assert.That(commons.Count, Is.EqualTo(1), "Common mod should be de-duplicated.");
        }
        [Test]
        public void Internal_Collapse_Does_Not_Collapse_When_BaseSequence_Diff()
        {
            var p1 = new Protein(sequence: "AAAA", accession: "SAME", organism: "o",
                geneNames: new List<Tuple<string, string>>(), oneBasedModifications: null, proteolysisProducts: null);
            var p2 = new Protein(sequence: "AAAB", accession: "SAME", organism: "o",
                geneNames: new List<Tuple<string, string>>(), oneBasedModifications: null, proteolysisProducts: null);

            var collapsed = InvokeInternalStatic<List<Protein>>(
                typeof(ProteinDbLoader),
                "CollapseDuplicateProteinsByAccessionAndBaseSequence",
                new List<Protein> { p1, p2 });

            // Both remain because BaseSequence differs
            Assert.That(collapsed.Count(p => p.Accession == "SAME"), Is.EqualTo(2));
        }
    }
}