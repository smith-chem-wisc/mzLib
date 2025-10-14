using NUnit.Framework;
using Omics.BioPolymer;
using Omics.Modifications;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public class TestProteinXmlDiagnostics
    {
        [Test]
        public void DiagnosticTest_LoadProteinXML_WithDecoyValidation()
        {
            // Path to the large XML file on your hard drive
            string proteinDbLocation = @"E:\Projects\Mann_11cell_lines\A549\A549_1\uniprotkb_taxonomy_id_9606_AND_reviewed_2024_10_07.xml";

            // Path to the log file where results will be written
            string logFilePath = @"E:\Projects\Mann_11cell_lines\A549\A549_1\DiagnosticResults.log";

            // Path to the second log file for decoy generation results
            string decoyLogFilePath = @"E:\Projects\Mann_11cell_lines\A549\A549_1\DecoyGenerationResults.log";

            // Ensure the log files are empty before starting
            if (File.Exists(logFilePath))
            {
                File.Delete(logFilePath);
            }
            if (File.Exists(decoyLogFilePath))
            {
                File.Delete(decoyLogFilePath);
            }

            // Options for loading the protein database
            bool generateTargets = true;
            DecoyType decoyType = DecoyType.Reverse; // Do NOT generate decoys initially
            IEnumerable<Modification> allKnownModifications = new List<Modification>(); // No predefined modifications
            bool isContaminant = false;
            IEnumerable<string> modTypesToExclude = new List<string>();
            int maxThreads = Environment.ProcessorCount; // Use all available threads
            int maxSequenceVariantsPerIsoform = 1;
            int minAlleleDepth = 0;
            int maxSequenceVariantIsoforms = 10;
            bool addTruncations = false;
            string decoyIdentifier = "DECOY";

            try
            {
                // Load the protein database (targets only, no decoys)
                Dictionary<string, Modification> unknownModifications;
                List<Protein> proteins = ProteinDbLoader.LoadProteinXML(
                    proteinDbLocation,
                    generateTargets,
                    decoyType,
                    allKnownModifications,
                    isContaminant,
                    modTypesToExclude,
                    out unknownModifications,
                    maxThreads,
                    maxSequenceVariantsPerIsoform,
                    minAlleleDepth,
                    maxSequenceVariantIsoforms,
                    addTruncations,
                    decoyIdentifier
                );

                // Count target proteins with AppliedVariations
                int targetWithAppliedVariations = proteins.Count(p => p.AppliedSequenceVariations != null && p.AppliedSequenceVariations.Count > 0);

                // Log the results of protein loading
                using (StreamWriter logWriter = new StreamWriter(logFilePath, append: true))
                {
                    logWriter.WriteLine($"Diagnostic Test Results for {proteinDbLocation}");
                    logWriter.WriteLine($"Total Proteins Loaded: {proteins.Count}");
                    logWriter.WriteLine($"Unknown Modifications: {unknownModifications.Count}");
                    logWriter.WriteLine($"Target Proteins with Applied Variations: {targetWithAppliedVariations}");

                    foreach (var kvp in unknownModifications)
                    {
                        logWriter.WriteLine($"Unknown Modification: {kvp.Key} - {kvp.Value}");
                    }

                    foreach (Protein protein in proteins)
                    {
                        logWriter.WriteLine($"Successfully loaded Protein Accession: {protein.Accession}");
                    }
                }

                // Attempt to create decoys for each protein and log the results
                using (StreamWriter decoyLogWriter = new StreamWriter(decoyLogFilePath, append: true))
                {
                    decoyLogWriter.WriteLine($"Decoy Generation Results for {proteinDbLocation}");

                    int decoyWithAppliedVariations = 0;

                    foreach (Protein protein in proteins)
                    {
                        try
                        {
                            // Attempt to create a decoy for the current protein
                            List<Protein> decoys = DecoyProteinGenerator.GenerateDecoys(
                                new List<Protein> { protein },
                                DecoyType.Reverse // Generate reverse decoys
                            );

                            // Count decoys with AppliedVariations
                            decoyWithAppliedVariations += decoys.Count(d => d.AppliedSequenceVariations != null && d.AppliedSequenceVariations.Count > 0);

                            // Log success
                            decoyLogWriter.WriteLine($"Successfully created decoy for Protein Accession: {protein.Accession}");
                        }
                        catch (Exception ex)
                        {
                            // Log failure
                            decoyLogWriter.WriteLine($"Failed to create decoy for Protein Accession: {protein.Accession}");
                            decoyLogWriter.WriteLine($"Error Message: {ex.Message}");
                        }
                    }

                    // Log the count of decoys with AppliedVariations
                    decoyLogWriter.WriteLine($"Decoy Proteins with Applied Variations: {decoyWithAppliedVariations}");
                }
            }
            catch (Exception ex)
            {
                // Log any critical errors that prevent the test from running
                File.AppendAllText(logFilePath, $"Critical Error: {ex.Message}{Environment.NewLine}");
            }
        }

        [Test]
        public void DiagnosticTest_LoadProteinXML_WithoutDecoys()
        {
            // Path to the large XML file on your hard drive
            string proteinDbLocation = @"E:\Projects\Mann_11cell_lines\A549\A549_1\uniprotkb_taxonomy_id_9606_AND_reviewed_2024_10_07.xml";

            // Path to the log file where results will be written
            string logFilePath = @"E:\Projects\Mann_11cell_lines\A549\A549_1\DiagnosticResults_NoDecoys.log";

            // Ensure the log file is empty before starting
            if (File.Exists(logFilePath))
            {
                File.Delete(logFilePath);
            }

            // Options for loading the protein database
            bool generateTargets = true;
            DecoyType decoyType = DecoyType.None; // Do NOT generate decoys
            IEnumerable<Modification> allKnownModifications = new List<Modification>(); // No predefined modifications
            bool isContaminant = false;
            IEnumerable<string> modTypesToExclude = new List<string>();
            int maxThreads = Environment.ProcessorCount; // Use all available threads
            int maxSequenceVariantsPerIsoform = 1;
            int minAlleleDepth = 0;
            int maxSequenceVariantIsoforms = 2;
            bool addTruncations = false;
            string decoyIdentifier = "DECOY";

            try
            {
                // Load the protein database (targets only, no decoys)
                Dictionary<string, Modification> unknownModifications;
                List<Protein> proteins = ProteinDbLoader.LoadProteinXML(
                    proteinDbLocation,
                    generateTargets,
                    decoyType,
                    allKnownModifications,
                    isContaminant,
                    modTypesToExclude,
                    out unknownModifications,
                    maxThreads,
                    maxSequenceVariantsPerIsoform,
                    minAlleleDepth,
                    maxSequenceVariantIsoforms,
                    addTruncations,
                    decoyIdentifier
                );

                // Log the results
                using (StreamWriter logWriter = new StreamWriter(logFilePath, append: true))
                {
                    logWriter.WriteLine($"Diagnostic Test Results for {proteinDbLocation}");
                    logWriter.WriteLine($"Total Proteins Loaded: {proteins.Count}");
                    logWriter.WriteLine($"Unknown Modifications: {unknownModifications.Count}");

                    foreach (var kvp in unknownModifications)
                    {
                        logWriter.WriteLine($"Unknown Modification: {kvp.Key} - {kvp.Value}");
                    }

                    foreach (Protein protein in proteins)
                    {
                        try
                        {
                            // Validate the protein
                            if (protein.BaseSequence == null || protein.BaseSequence.Length == 0)
                            {
                                throw new Exception("Protein has an empty or null base sequence.");
                            }

                            // Log successful protein creation
                            logWriter.WriteLine($"Successfully loaded Protein Accession: {protein.Accession}");
                        }
                        catch (Exception ex)
                        {
                            // Log the error for this protein
                            logWriter.WriteLine($"Error with Protein Accession: {protein.Accession}");
                            logWriter.WriteLine($"Error Message: {ex.Message}");
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                // Log any critical errors that prevent the test from running
                File.AppendAllText(logFilePath, $"Critical Error: {ex.Message}{Environment.NewLine}");
            }
        }

        [Test]
        public void DiagnosticTest_SingleProteinDecoyValidation()
        {
            // Path to the single protein XML file
            string proteinDbLocation = @"E:\Projects\Mann_11cell_lines\A549\A549_1\O76039.xml";

            // Path to the log file where results will be written
            string logFilePath = @"E:\Projects\Mann_11cell_lines\A549\A549_1\SingleProteinDiagnosticResults.log";

            // Path to the second log file for decoy generation results
            string decoyLogFilePath = @"E:\Projects\Mann_11cell_lines\A549\A549_1\SingleProteinDecoyGenerationResults.log";

            // Ensure the log files are empty before starting
            if (File.Exists(logFilePath))
            {
                File.Delete(logFilePath);
            }
            if (File.Exists(decoyLogFilePath))
            {
                File.Delete(decoyLogFilePath);
            }

            // Options for loading the protein database
            bool generateTargets = true;
            DecoyType decoyType = DecoyType.None; // Do NOT generate decoys initially
            IEnumerable<Modification> allKnownModifications = new List<Modification>(); // No predefined modifications
            bool isContaminant = false;
            IEnumerable<string> modTypesToExclude = new List<string>();
            int maxThreads = 1; // Single-threaded for simplicity
            int maxSequenceVariantsPerIsoform = 0;
            int minAlleleDepth = 0;
            int maxSequenceVariantIsoforms = 1;
            bool addTruncations = false;
            string decoyIdentifier = "DECOY";

            try
            {
                // Load the single protein (targets only, no decoys)
                Dictionary<string, Modification> unknownModifications;
                List<Protein> proteins = ProteinDbLoader.LoadProteinXML(
                    proteinDbLocation,
                    generateTargets,
                    decoyType,
                    allKnownModifications,
                    isContaminant,
                    modTypesToExclude,
                    out unknownModifications,
                    maxThreads,
                    maxSequenceVariantsPerIsoform,
                    minAlleleDepth,
                    maxSequenceVariantIsoforms,
                    addTruncations,
                    decoyIdentifier
                );

                // Log the results of protein loading
                using (StreamWriter logWriter = new StreamWriter(logFilePath, append: true))
                {
                    logWriter.WriteLine($"Diagnostic Test Results for {proteinDbLocation}");
                    logWriter.WriteLine($"Total Proteins Loaded: {proteins.Count}");
                    logWriter.WriteLine($"Unknown Modifications: {unknownModifications.Count}");

                    foreach (var kvp in unknownModifications)
                    {
                        logWriter.WriteLine($"Unknown Modification: {kvp.Key} - {kvp.Value}");
                    }

                    foreach (Protein protein in proteins)
                    {
                        logWriter.WriteLine($"Successfully loaded Protein Accession: {protein.Accession}");
                        logWriter.WriteLine($"Protein Base Sequence: {protein.BaseSequence}");
                        logWriter.WriteLine($"Protein Length: {protein.Length}");
                        logWriter.WriteLine($"Protein Modifications: {protein.OneBasedPossibleLocalizedModifications.Count}");
                        foreach (var mod in protein.OneBasedPossibleLocalizedModifications)
                        {
                            logWriter.WriteLine($"  Modification at Position {mod.Key}: {string.Join(", ", mod.Value)}");
                        }
                    }
                }

                // Attempt to create a decoy for the single protein and log the results
                using (StreamWriter decoyLogWriter = new StreamWriter(decoyLogFilePath, append: true))
                {
                    decoyLogWriter.WriteLine($"Decoy Generation Results for {proteinDbLocation}");

                    foreach (Protein protein in proteins)
                    {
                        try
                        {
                            // Attempt to create a decoy for the current protein
                            List<Protein> decoys = DecoyProteinGenerator.GenerateDecoys(
                                new List<Protein> { protein },
                                DecoyType.Reverse // Generate reverse decoys
                            );

                            // Log success
                            foreach (var decoy in decoys)
                            {
                                decoyLogWriter.WriteLine($"Successfully created decoy for Protein Accession: {protein.Accession}");
                                decoyLogWriter.WriteLine($"Decoy Base Sequence: {decoy.BaseSequence}");
                                decoyLogWriter.WriteLine($"Decoy Length: {decoy.Length}");
                                decoyLogWriter.WriteLine($"Decoy Modifications: {decoy.OneBasedPossibleLocalizedModifications.Count}");
                                foreach (var mod in decoy.OneBasedPossibleLocalizedModifications)
                                {
                                    decoyLogWriter.WriteLine($"  Modification at Position {mod.Key}: {string.Join(", ", mod.Value)}");
                                }
                            }
                        }
                        catch (Exception ex)
                        {
                            // Log failure
                            decoyLogWriter.WriteLine($"Failed to create decoy for Protein Accession: {protein.Accession}");
                            decoyLogWriter.WriteLine($"Error Message: {ex.Message}");
                            decoyLogWriter.WriteLine($"Stack Trace: {ex.StackTrace}");
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                // Log any critical errors that prevent the test from running
                File.AppendAllText(logFilePath, $"Critical Error: {ex.Message}{Environment.NewLine}");
                File.AppendAllText(logFilePath, $"Stack Trace: {ex.StackTrace}{Environment.NewLine}");
            }
        }
        [Test]
        public void DiagnosticTest_SingleProteinSequenceVariantDecoyValidation()
        {
            // Path to the single protein XML file
            string proteinDbLocation = @"E:\Projects\Mann_11cell_lines\A549\A549_1\O76039.xml";

            // Path to the log file where results will be written
            string logFilePath = @"E:\Projects\Mann_11cell_lines\A549\A549_1\SingleProteinVariantDiagnosticResults.log";

            // Path to the second log file for decoy generation results
            string decoyLogFilePath = @"E:\Projects\Mann_11cell_lines\A549\A549_1\SingleProteinVariantDecoyGenerationResults.log";

            // Ensure the log files are empty before starting
            if (File.Exists(logFilePath))
            {
                File.Delete(logFilePath);
            }
            if (File.Exists(decoyLogFilePath))
            {
                File.Delete(decoyLogFilePath);
            }

            // Options for loading the protein database
            bool generateTargets = true;
            DecoyType decoyType = DecoyType.None; // Do NOT generate decoys initially
            IEnumerable<Modification> allKnownModifications = new List<Modification>(); // No predefined modifications
            bool isContaminant = false;
            IEnumerable<string> modTypesToExclude = new List<string>();
            int maxThreads = 1; // Single-threaded for simplicity
            int maxSequenceVariantsPerIsoform = 0;
            int minAlleleDepth = 0;
            int maxSequenceVariantIsoforms = 1;
            bool addTruncations = false;
            string decoyIdentifier = "DECOY";

            try
            {
                // Load the single protein (targets only, no decoys)
                Dictionary<string, Modification> unknownModifications;
                List<Protein> proteins = ProteinDbLoader.LoadProteinXML(
                    proteinDbLocation,
                    generateTargets,
                    decoyType,
                    allKnownModifications,
                    isContaminant,
                    modTypesToExclude,
                    out unknownModifications,
                    maxThreads,
                    maxSequenceVariantsPerIsoform,
                    minAlleleDepth,
                    maxSequenceVariantIsoforms,
                    addTruncations,
                    decoyIdentifier
                );

                // Log the results of protein loading
                using (StreamWriter logWriter = new StreamWriter(logFilePath, append: true))
                {
                    logWriter.WriteLine($"Diagnostic Test Results for {proteinDbLocation}");
                    logWriter.WriteLine($"Total Proteins Loaded: {proteins.Count}");
                    logWriter.WriteLine($"Unknown Modifications: {unknownModifications.Count}");

                    foreach (var kvp in unknownModifications)
                    {
                        logWriter.WriteLine($"Unknown Modification: {kvp.Key} - {kvp.Value}");
                    }

                    foreach (Protein protein in proteins)
                    {
                        logWriter.WriteLine($"Successfully loaded Protein Accession: {protein.Accession}");
                        logWriter.WriteLine($"Protein Base Sequence: {protein.BaseSequence}");
                        logWriter.WriteLine($"Protein Length: {protein.Length}");
                        logWriter.WriteLine($"Protein Modifications: {protein.OneBasedPossibleLocalizedModifications.Count}");
                        foreach (var mod in protein.OneBasedPossibleLocalizedModifications)
                        {
                            logWriter.WriteLine($"  Modification at Position {mod.Key}: {string.Join(", ", mod.Value)}");
                        }

                        // Process each sequence variation individually
                        var sequenceVariations = protein.SequenceVariations;
                        logWriter.WriteLine($"Total Sequence Variations: {sequenceVariations.Count}");

                        foreach (var variation in sequenceVariations)
                        {
                            logWriter.WriteLine($"Sequence Variation Details:");
                            logWriter.WriteLine($"  Begin Position: {variation.OneBasedBeginPosition}");
                            logWriter.WriteLine($"  End Position: {variation.OneBasedEndPosition}");
                            logWriter.WriteLine($"  Original Sequence: {variation.OriginalSequence}");
                            logWriter.WriteLine($"  Variant Sequence: {variation.VariantSequence}");
                            logWriter.WriteLine($"  Description: {variation.Description}");
                            if (variation.OneBasedModifications != null)
                            {
                                foreach (var mod in variation.OneBasedModifications)
                                {
                                    logWriter.WriteLine($"    Modification at Position {mod.Key}: {string.Join(", ", mod.Value)}");
                                }
                            }
                        }
                    }
                }

                // Attempt to create decoys for each sequence variation and log the results
                using (StreamWriter decoyLogWriter = new StreamWriter(decoyLogFilePath, append: true))
                {
                    decoyLogWriter.WriteLine($"Decoy Generation Results for {proteinDbLocation}");

                    foreach (Protein protein in proteins)
                    {
                        foreach (var variation in protein.SequenceVariations)
                        {
                            try
                            {
                                // Create a new target protein with only this sequence variation
                                var targetProteinWithVariation = new Protein(
                                    variantBaseSequence: variation.VariantSequence,
                                    protein: protein,
                                    appliedSequenceVariations: new List<SequenceVariation> { variation },
                                    applicableProteolysisProducts: protein.TruncationProducts,
                                    oneBasedModifications: protein.OneBasedPossibleLocalizedModifications,
                                    sampleNameForVariants: protein.SampleNameForVariants
                                );

                                // Attempt to create a decoy for the new target protein
                                List<Protein> decoys = DecoyProteinGenerator.GenerateDecoys(
                                    new List<Protein> { targetProteinWithVariation },
                                    DecoyType.Reverse // Generate reverse decoys
                                );

                                // Log success
                                foreach (var decoy in decoys)
                                {
                                    decoyLogWriter.WriteLine($"Successfully created decoy for Protein Accession: {protein.Accession}");
                                    decoyLogWriter.WriteLine($"  Sequence Variation:");
                                    decoyLogWriter.WriteLine($"    Begin Position: {variation.OneBasedBeginPosition}");
                                    decoyLogWriter.WriteLine($"    End Position: {variation.OneBasedEndPosition}");
                                    decoyLogWriter.WriteLine($"    Original Sequence: {variation.OriginalSequence}");
                                    decoyLogWriter.WriteLine($"    Variant Sequence: {variation.VariantSequence}");
                                    decoyLogWriter.WriteLine($"    Description: {variation.Description}");
                                    decoyLogWriter.WriteLine($"  Decoy Base Sequence: {decoy.BaseSequence}");
                                    decoyLogWriter.WriteLine($"  Decoy Length: {decoy.Length}");
                                    decoyLogWriter.WriteLine($"  Decoy Modifications: {decoy.OneBasedPossibleLocalizedModifications.Count}");
                                    foreach (var mod in decoy.OneBasedPossibleLocalizedModifications)
                                    {
                                        decoyLogWriter.WriteLine($"    Modification at Position {mod.Key}: {string.Join(", ", mod.Value)}");
                                    }
                                }
                            }
                            catch (Exception ex)
                            {
                                // Log failure
                                decoyLogWriter.WriteLine($"Failed to create decoy for Protein Accession: {protein.Accession}");
                                decoyLogWriter.WriteLine($"  Sequence Variation:");
                                decoyLogWriter.WriteLine($"    Begin Position: {variation.OneBasedBeginPosition}");
                                decoyLogWriter.WriteLine($"    End Position: {variation.OneBasedEndPosition}");
                                decoyLogWriter.WriteLine($"    Original Sequence: {variation.OriginalSequence}");
                                decoyLogWriter.WriteLine($"    Variant Sequence: {variation.VariantSequence}");
                                decoyLogWriter.WriteLine($"    Description: {variation.Description}");
                                decoyLogWriter.WriteLine($"  Error Message: {ex.Message}");
                                decoyLogWriter.WriteLine($"  Stack Trace: {ex.StackTrace}");
                            }
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                // Log any critical errors that prevent the test from running
                File.AppendAllText(logFilePath, $"Critical Error: {ex.Message}{Environment.NewLine}");
                File.AppendAllText(logFilePath, $"Stack Trace: {ex.StackTrace}{Environment.NewLine}");
            }
        }
        [Test]
        public void DiagnosticTest_LoadProteinXML_SequenceVariantValidation()
        {
            // Path to the large XML file on your hard drive
            string proteinDbLocation = @"E:\Projects\Mann_11cell_lines\A549\A549_1\uniprotkb_taxonomy_id_9606_AND_reviewed_2024_10_07.xml";

            // Path to the log file where results will be written
            string logFilePath = @"E:\Projects\Mann_11cell_lines\A549\A549_1\SequenceVariantDiagnosticResults.log";

            // Path to the second log file for decoy generation results
            string decoyLogFilePath = @"E:\Projects\Mann_11cell_lines\A549\A549_1\SequenceVariantDecoyResults.log";

            // Ensure the log files are empty before starting
            if (File.Exists(logFilePath))
            {
                File.Delete(logFilePath);
            }
            if (File.Exists(decoyLogFilePath))
            {
                File.Delete(decoyLogFilePath);
            }

            // Options for loading the protein database
            bool generateTargets = true;
            DecoyType decoyType = DecoyType.Reverse; // Generate reverse decoys
            IEnumerable<Modification> allKnownModifications = new List<Modification>(); // No predefined modifications
            bool isContaminant = false;
            IEnumerable<string> modTypesToExclude = new List<string>();
            int maxThreads = Environment.ProcessorCount; // Use all available threads
            int maxSequenceVariantsPerIsoform = 1;
            int minAlleleDepth = 0;
            int maxSequenceVariantIsoforms = 10;
            bool addTruncations = false;
            string decoyIdentifier = "DECOY";

            try
            {
                // Load the protein database (targets only, no decoys)
                Dictionary<string, Modification> unknownModifications;
                List<Protein> proteins = ProteinDbLoader.LoadProteinXML(
                    proteinDbLocation,
                    generateTargets,
                    decoyType,
                    allKnownModifications,
                    isContaminant,
                    modTypesToExclude,
                    out unknownModifications,
                    maxThreads,
                    maxSequenceVariantsPerIsoform,
                    minAlleleDepth,
                    maxSequenceVariantIsoforms,
                    addTruncations,
                    decoyIdentifier
                );

                // Log the results of protein loading
                using (StreamWriter logWriter = new StreamWriter(logFilePath, append: true))
                {
                    logWriter.WriteLine($"Sequence Variant Diagnostic Test Results for {proteinDbLocation}");
                    logWriter.WriteLine($"Total Proteins Loaded: {proteins.Count}");
                    logWriter.WriteLine($"Unknown Modifications: {unknownModifications.Count}");

                    foreach (var kvp in unknownModifications)
                    {
                        logWriter.WriteLine($"Unknown Modification: {kvp.Key} - {kvp.Value}");
                    }

                    foreach (Protein protein in proteins)
                    {
                        int validSequenceVariants = protein.SequenceVariations.Count(v => v.AreValid());
                        logWriter.WriteLine($"Protein Accession: {protein.Accession}");
                        logWriter.WriteLine($"  Total Sequence Variants: {protein.SequenceVariations.Count}");
                        logWriter.WriteLine($"  Valid Sequence Variants: {validSequenceVariants}");
                    }
                }

                // Attempt to create decoys for each protein and log the results
                using (StreamWriter decoyLogWriter = new StreamWriter(decoyLogFilePath, append: true))
                {
                    decoyLogWriter.WriteLine($"Sequence Variant Decoy Results for {proteinDbLocation}");

                    foreach (Protein protein in proteins)
                    {
                        try
                        {
                            // Attempt to create a decoy for the current protein
                            List<Protein> decoys = DecoyProteinGenerator.GenerateDecoys(
                                new List<Protein> { protein },
                                DecoyType.Reverse // Generate reverse decoys
                            );

                            foreach (Protein decoy in decoys)
                            {
                                int validSequenceVariants = decoy.SequenceVariations.Count(v => v.AreValid());
                                decoyLogWriter.WriteLine($"Decoy Protein Accession: {decoy.Accession}");
                                decoyLogWriter.WriteLine($"  Total Sequence Variants: {decoy.SequenceVariations.Count}");
                                decoyLogWriter.WriteLine($"  Valid Sequence Variants: {validSequenceVariants}");
                            }
                        }
                        catch (Exception ex)
                        {
                            // Log failure
                            decoyLogWriter.WriteLine($"Failed to create decoy for Protein Accession: {protein.Accession}");
                            decoyLogWriter.WriteLine($"Error Message: {ex.Message}");
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                // Log any critical errors that prevent the test from running
                File.AppendAllText(logFilePath, $"Critical Error: {ex.Message}{Environment.NewLine}");
            }
        }
        [Test]
        public void DiagnosticTest_SingleProteinXML_SequenceVariantValidation()
        {
            // Path to the single protein XML file
            string proteinDbLocation = @"E:\Projects\Mann_11cell_lines\A549\A549_1\O76039.xml";

            // Path to the log file where results will be written
            string logFilePath = @"E:\Projects\Mann_11cell_lines\A549\A549_1\SingleProteinSequenceVariantDiagnosticResults.log";

            // Path to the second log file for decoy generation results
            string decoyLogFilePath = @"E:\Projects\Mann_11cell_lines\A549\A549_1\SingleProteinSequenceVariantDecoyResults.log";

            // Ensure the log files are empty before starting
            if (File.Exists(logFilePath))
            {
                File.Delete(logFilePath);
            }
            if (File.Exists(decoyLogFilePath))
            {
                File.Delete(decoyLogFilePath);
            }

            // Options for loading the protein database
            bool generateTargets = true;
            DecoyType decoyType = DecoyType.Reverse; // Generate reverse decoys
            IEnumerable<Modification> allKnownModifications = new List<Modification>(); // No predefined modifications
            bool isContaminant = false;
            IEnumerable<string> modTypesToExclude = new List<string>();
            int maxThreads = 1; // Single-threaded for simplicity
            int maxSequenceVariantsPerIsoform = 1;
            int minAlleleDepth = 0;
            int maxSequenceVariantIsoforms = 1;
            bool addTruncations = false;
            string decoyIdentifier = "DECOY";

            try
            {
                // Load the single protein (targets only, no decoys)
                Dictionary<string, Modification> unknownModifications;
                List<Protein> proteins = ProteinDbLoader.LoadProteinXML(
                    proteinDbLocation,
                    generateTargets,
                    decoyType,
                    allKnownModifications,
                    isContaminant,
                    modTypesToExclude,
                    out unknownModifications,
                    maxThreads,
                    maxSequenceVariantsPerIsoform,
                    minAlleleDepth,
                    maxSequenceVariantIsoforms,
                    addTruncations,
                    decoyIdentifier
                );

                // Log the results of protein loading
                using (StreamWriter logWriter = new StreamWriter(logFilePath, append: true))
                {
                    logWriter.WriteLine($"Single Protein Sequence Variant Diagnostic Test Results for {proteinDbLocation}");
                    logWriter.WriteLine($"Total Proteins Loaded: {proteins.Count}");
                    logWriter.WriteLine($"Unknown Modifications: {unknownModifications.Count}");

                    foreach (var kvp in unknownModifications)
                    {
                        logWriter.WriteLine($"Unknown Modification: {kvp.Key} - {kvp.Value}");
                    }

                    foreach (Protein protein in proteins)
                    {
                        int validSequenceVariants = protein.SequenceVariations.Count(v => v.AreValid());
                        logWriter.WriteLine($"Protein Accession: {protein.Accession}");
                        logWriter.WriteLine($"  Total Sequence Variants: {protein.SequenceVariations.Count}");
                        logWriter.WriteLine($"  Valid Sequence Variants: {validSequenceVariants}");
                    }
                }

                // Attempt to create a decoy for the single protein and log the results
                using (StreamWriter decoyLogWriter = new StreamWriter(decoyLogFilePath, append: true))
                {
                    decoyLogWriter.WriteLine($"Single Protein Sequence Variant Decoy Results for {proteinDbLocation}");

                    foreach (Protein protein in proteins)
                    {
                        try
                        {
                            // Attempt to create a decoy for the current protein
                            List<Protein> decoys = DecoyProteinGenerator.GenerateDecoys(
                                new List<Protein> { protein },
                                DecoyType.Reverse // Generate reverse decoys
                            );

                            foreach (Protein decoy in decoys)
                            {
                                int validSequenceVariants = decoy.SequenceVariations.Count(v => v.AreValid());
                                decoyLogWriter.WriteLine($"Decoy Protein Accession: {decoy.Accession}");
                                decoyLogWriter.WriteLine($"  Total Sequence Variants: {decoy.SequenceVariations.Count}");
                                decoyLogWriter.WriteLine($"  Valid Sequence Variants: {validSequenceVariants}");
                            }
                        }
                        catch (Exception ex)
                        {
                            // Log failure
                            decoyLogWriter.WriteLine($"Failed to create decoy for Protein Accession: {protein.Accession}");
                            decoyLogWriter.WriteLine($"Error Message: {ex.Message}");
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                // Log any critical errors that prevent the test from running
                File.AppendAllText(logFilePath, $"Critical Error: {ex.Message}{Environment.NewLine}");
            }
        }
    }
}