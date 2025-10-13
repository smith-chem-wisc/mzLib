using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Proteomics;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class TestProteinXmlDiagnostics
    {
        [Test]
        public static void Diagnose_SequenceVariation_CoordinateErrors_DecoyStepwise()
        {
            string xmlPath = @"E:\Projects\Mann_11cell_lines\A549\A549_1\uniprotkb_taxonomy_id_9606_AND_reviewed_2024_10_07.xml";
            string logPath = Path.Combine(Path.GetDirectoryName(xmlPath), "diagnostic_sequence_variation_errors_decoy_stepwise.txt");
            var log = new List<string>();

            List<Protein> proteins = null;
            try
            {
                proteins = ProteinDbLoader.LoadProteinXML(
                    xmlPath,
                    generateTargets: true,
                    decoyType: DecoyType.None, // Only targets, so we can step through decoy generation
                    allKnownModifications: null,
                    isContaminant: false,
                    modTypesToExclude: null,
                    unknownModifications: out var _,
                    maxThreads: -1,
                    maxSequenceVariantsPerIsoform: 0,
                    minAlleleDepth: 0,
                    maxSequenceVariantIsoforms: 1,
                    addTruncations: false);
            }
            catch (Exception ex)
            {
                log.Add($"[LoadProteinXML (targets only) EXCEPTION] {ex}");
                File.WriteAllLines(logPath, log);
                Assert.Fail($"LoadProteinXML (targets only) failed: {ex.Message}");
                return;
            }

            if (proteins == null)
            {
                log.Add("[LoadProteinXML] Protein list is null");
                File.WriteAllLines(logPath, log);
                Assert.Fail("Protein list is null");
                return;
            }

            // Try to generate decoys one by one, logging any errors
            foreach (var protein in proteins)
            {
                try
                {
                    // This is a simplified version of what DecoyProteinGenerator.GenerateDecoys does for Reverse decoys
                    var decoys = UsefulProteomicsDatabases.DecoyProteinGenerator.GenerateDecoys(
                        new List<Protein> { protein },
                        DecoyType.Reverse,
                        maxThreads: 1, // or -1 for default
                        decoyIdentifier: "DECOY"
                    );
                    // Optionally, check sequence variations in the decoy
                    foreach (var decoy in decoys)
                    {
                        if (decoy.SequenceVariations != null)
                        {
                            foreach (var sv in decoy.SequenceVariations)
                            {
                                var isValid = sv.GetType().GetMethod("AreValid")?.Invoke(sv, null);
                                if (isValid is bool valid && !valid)
                                {
                                    log.Add($"[DECOY INVALID] Accession: {decoy.Accession}, SequenceVariation: {DescribeVariation(sv)}");
                                }
                            }
                        }
                    }
                }
                catch (Exception ex)
                {
                    log.Add($"[DECOY EXCEPTION] Target Accession: {protein.Accession}, Error: {ex}");
                    // Optionally, log the target's sequence variations for context
                    if (protein.SequenceVariations != null)
                    {
                        foreach (var sv in protein.SequenceVariations)
                        {
                            log.Add($"[TARGET VARIATION] Accession: {protein.Accession}, SequenceVariation: {DescribeVariation(sv)}");
                        }
                    }
                }
            }

            File.WriteAllLines(logPath, log);
            Assert.Pass($"Diagnostics complete. See {logPath} for details.");
        }

        private static string DescribeVariation(object sv)
        {
            try
            {
                var begin = sv.GetType().GetProperty("OneBasedBeginPosition")?.GetValue(sv);
                var end = sv.GetType().GetProperty("OneBasedEndPosition")?.GetValue(sv);
                var orig = sv.GetType().GetProperty("OriginalSequence")?.GetValue(sv);
                var varSeq = sv.GetType().GetProperty("VariantSequence")?.GetValue(sv);
                var desc = sv.GetType().GetProperty("Description")?.GetValue(sv);
                return $"[{begin}-{end}] {orig}?{varSeq} ({desc})";
            }
            catch
            {
                return sv?.ToString() ?? "(null)";
            }
        }
    }
}