using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using UsefulProteomicsDatabases;
using Assert = NUnit.Framework.Legacy.ClassicAssert;

namespace Test.DatabaseTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class TestProteoformMapping
    {
        private static StreamWriter _logWriter;
        private static string _logPath;

        private static void Log(string message)
        {
            TestContext.WriteLine(message);
            _logWriter?.WriteLine(message);
        }

        [Test]
        public static void MapObservedProteoformsToTheoretical()
        {
            // === INPUT FILE PATHS ===
            string observedFile1 = @"E:\Projects\ProteoformLists\MM-AllProteoforms.psmtsv";
            string observedFile2 = @"E:\Projects\ProteoformLists\ProteomeDiscovererProteoforms.csv";
            string observedFile3 = @"E:\Projects\ProteoformLists\toppic.tsv";

            string xmlDatabase1 = @"E:\Projects\ProteoformLists\TrEMBLonlyPrune.xml";
            string xmlDatabase2 = @"E:\Projects\ProteoformLists\uniprotkb_9606_reviewed-HIV_WTJB474_updated_JustGag-JustPolGPTMDproteinPruned-noCitCarbox-notrunc.xml";

            string outputPath = @"E:\Projects\ProteoformLists\output.txt";
            _logPath = @"E:\Projects\ProteoformLists\mapping_log.txt";

            try
            {
                _logWriter = new StreamWriter(_logPath) { AutoFlush = true };
            }
            catch (Exception ex)
            {
                TestContext.WriteLine($"WARNING: Could not create log file: {ex.Message}");
                _logWriter = null;
            }

            try
            {
                Log("=== STARTING PROTEOFORM MAPPING TEST ===");
                Log($"Timestamp: {DateTime.Now:yyyy-MM-dd HH:mm:ss}");
                Log($"Log file: {_logPath}");

                Log("\n=== TERMINOLOGY ===");
                Log("Full-Length: Observed sequence EXACTLY matches a theoretical proteoform (chain, processed form, etc.)");
                Log("Full-Length (minus Met): Same as Full-Length, but parent protein starts with M and proteoform starts at position 2");
                Log("Partial Match: Observed sequence is a SUBSTRING of a theoretical proteoform (missing N-terminal, C-terminal, or both residues)");
                Log("No Match: Observed sequence not found in any theoretical proteoform or protein");

                // === CHECK FILE EXISTENCE ===
                Log("\n=== CHECKING INPUT FILES ===");
                Log($"observedFile1: {observedFile1} - Exists: {File.Exists(observedFile1)}");
                Log($"observedFile2: {observedFile2} - Exists: {File.Exists(observedFile2)}");
                Log($"observedFile3: {observedFile3} - Exists: {File.Exists(observedFile3)}");
                Log($"xmlDatabase1: {xmlDatabase1} - Exists: {File.Exists(xmlDatabase1)}");
                Log($"xmlDatabase2: {xmlDatabase2} - Exists: {File.Exists(xmlDatabase2)}");

                // === LOAD PROTEIN DATABASES ===
                Log("\n=== LOADING PROTEIN DATABASES ===");
                var allProteins = new List<Protein>();

                if (File.Exists(xmlDatabase1))
                {
                    try
                    {
                        var proteins1 = ProteinDbLoader.LoadProteinXML(xmlDatabase1, true, DecoyType.None,
                            null, false, null, out var unknownMods1);
                        allProteins.AddRange(proteins1);
                        Log($"Loaded {proteins1.Count} proteins from {Path.GetFileName(xmlDatabase1)}");
                    }
                    catch (Exception ex)
                    {
                        Log($"ERROR loading {xmlDatabase1}: {ex.Message}");
                    }
                }

                if (File.Exists(xmlDatabase2))
                {
                    try
                    {
                        var proteins2 = ProteinDbLoader.LoadProteinXML(xmlDatabase2, true, DecoyType.None,
                            null, false, null, out var unknownMods2);
                        allProteins.AddRange(proteins2);
                        Log($"Loaded {proteins2.Count} proteins from {Path.GetFileName(xmlDatabase2)}");
                    }
                    catch (Exception ex)
                    {
                        Log($"ERROR loading {xmlDatabase2}: {ex.Message}");
                    }
                }

                Log($"Total proteins loaded: {allProteins.Count}");

                if (allProteins.Count == 0)
                {
                    Log("ERROR: No proteins loaded from XML databases");
                    Assert.Fail("No proteins loaded from XML databases");
                    return;
                }

                // === GENERATE THEORETICAL PROTEOFORMS ===
                Log("\n=== GENERATING THEORETICAL PROTEOFORMS ===");
                Log("Using 'top-down' protease which generates proteoforms from chains, signal peptides, etc.");

                var digestionParams = new DigestionParams(
                    protease: "top-down",
                    maxMissedCleavages: 0,
                    minPeptideLength: 2,
                    maxPeptideLength: int.MaxValue);

                var emptyMods = new List<Modification>();
                var theoreticalProteoforms = new List<(string Accession, string Sequence, int StartResidue, int EndResidue, Protein ParentProtein)>();

                int proteinCount = 0;
                foreach (var protein in allProteins)
                {
                    proteinCount++;
                    try
                    {
                        var digested = protein.Digest(digestionParams, emptyMods, emptyMods);
                        foreach (var proteoform in digested)
                        {
                            theoreticalProteoforms.Add((
                                protein.Accession,
                                proteoform.BaseSequence,
                                proteoform.OneBasedStartResidueInProtein,
                                proteoform.OneBasedEndResidueInProtein,
                                protein));
                        }
                    }
                    catch (Exception ex)
                    {
                        Log($"ERROR digesting protein {protein.Accession}: {ex.Message}");
                    }

                    if (proteinCount % 100 == 0)
                        Log($"  Processed {proteinCount} proteins, generated {theoreticalProteoforms.Count} proteoforms...");
                }

                Log($"Generated {theoreticalProteoforms.Count} theoretical proteoforms from {allProteins.Count} proteins");

                // Build lookup dictionary
                var sequenceToProteoforms = new Dictionary<string, List<(string Accession, string Sequence, int StartResidue, int EndResidue, Protein ParentProtein)>>();
                foreach (var proteoform in theoreticalProteoforms)
                {
                    if (!sequenceToProteoforms.ContainsKey(proteoform.Sequence))
                        sequenceToProteoforms[proteoform.Sequence] = new();
                    sequenceToProteoforms[proteoform.Sequence].Add(proteoform);
                }

                Log($"Unique theoretical proteoform sequences: {sequenceToProteoforms.Count}");

                // === READ OBSERVED PROTEOFORMS ===
                Log("\n=== READING OBSERVED PROTEOFORMS ===");
                var observedProteoforms = new List<(string FileName, string Accession, string Sequence, string OriginalSequence, int LineNumber)>();
                var observedFiles = new[] { observedFile1, observedFile2, observedFile3 }
                    .Where(f => !string.IsNullOrEmpty(f) && File.Exists(f))
                    .ToList();

                Log($"Found {observedFiles.Count} existing observed files");

                foreach (var filePath in observedFiles)
                {
                    try
                    {
                        var (accessionCol, sequenceCol, sequenceType, delimiter) = ParseHeader(filePath);

                        Log($"\nFile: {Path.GetFileName(filePath)}");
                        Log($"  Detected delimiter: '{(delimiter == '\t' ? "TAB" : delimiter == ',' ? "COMMA" : delimiter.ToString())}'");
                        Log($"  Accession col: {accessionCol}, Sequence col: {sequenceCol}, Type: {sequenceType}");

                        if (sequenceCol < 0)
                        {
                            string header = File.ReadLines(filePath).FirstOrDefault() ?? "";
                            Log($"  WARNING: Could not find sequence column");
                            Log($"  Full Header: {header}");
                            continue;
                        }

                        int lineNum = 0;
                        int proteoformCount = 0;
                        foreach (var line in File.ReadLines(filePath))
                        {
                            lineNum++;
                            if (lineNum == 1) continue;

                            var parts = line.Split(delimiter);
                            if (parts.Length <= sequenceCol) continue;

                            string accession = accessionCol >= 0 && parts.Length > accessionCol
                                ? parts[accessionCol].Trim()
                                : "";
                            string originalSequence = parts[sequenceCol].Trim();
                            string sequence = ProcessSequence(originalSequence, sequenceType);

                            if (string.IsNullOrWhiteSpace(sequence)) continue;

                            observedProteoforms.Add((Path.GetFileName(filePath), accession, sequence, originalSequence, lineNum));
                            proteoformCount++;

                            if (proteoformCount <= 3)
                                Log($"  Sample {proteoformCount}: Accession='{accession}', Seq='{sequence.Substring(0, Math.Min(30, sequence.Length))}...'");
                        }
                        Log($"  Read {proteoformCount} proteoforms");
                    }
                    catch (Exception ex)
                    {
                        Log($"ERROR reading {filePath}: {ex.Message}");
                    }
                }

                Log($"\nTotal observed proteoforms: {observedProteoforms.Count}");

                if (observedProteoforms.Count == 0)
                {
                    Log("ERROR: No observed proteoforms loaded");
                    Assert.Fail("No observed proteoforms loaded");
                    return;
                }

                // === MAP OBSERVED TO THEORETICAL PROTEOFORMS ===
                Log("\n=== MAPPING OBSERVED TO THEORETICAL ===");
                var results = new List<ProteoformMappingResult>();

                int processed = 0;
                foreach (var (fileName, accession, observedSequence, originalSequence, lineNum) in observedProteoforms)
                {
                    processed++;
                    var result = new ProteoformMappingResult
                    {
                        SourceFile = fileName,
                        LineNumber = lineNum,
                        ObservedAccession = accession,
                        ObservedSequence = observedSequence,
                        OriginalSequence = originalSequence
                    };

                    // Check for EXACT sequence match in theoretical proteoforms
                    if (sequenceToProteoforms.TryGetValue(observedSequence, out var exactMatches))
                    {
                        var shortestMatch = exactMatches
                            .OrderBy(m => m.Sequence.Length)
                            .ThenBy(m => m.ParentProtein.BaseSequence.Length)
                            .First();

                        bool isMinusMet = shortestMatch.StartResidue == 2 &&
                                          shortestMatch.ParentProtein.BaseSequence.StartsWith("M");

                        result.MatchType = isMinusMet ? "Full-Length (minus Met)" : "Full-Length";
                        result.MatchedAccessions.Add(shortestMatch.Accession);
                        result.TheoreticalLength = shortestMatch.ParentProtein.BaseSequence.Length;
                        result.ProteoformStart = shortestMatch.StartResidue;
                        result.ProteoformEnd = shortestMatch.EndResidue;
                        result.ProteoformLength = shortestMatch.Sequence.Length;
                        result.TotalMatchCount = exactMatches.Count;
                        result.TruncationInfo = $"{shortestMatch.StartResidue}-{shortestMatch.EndResidue} of {shortestMatch.ParentProtein.BaseSequence.Length}";
                    }
                    else
                    {
                        var containingProteoforms = theoreticalProteoforms
                            .Where(p => p.Sequence.Contains(observedSequence))
                            .OrderBy(p => p.Sequence.Length)
                            .ThenBy(p => p.ParentProtein.BaseSequence.Length)
                            .ToList();

                        if (containingProteoforms.Any())
                        {
                            var shortestMatch = containingProteoforms.First();

                            int startInProteoform = shortestMatch.Sequence.IndexOf(observedSequence);
                            int endInProteoform = startInProteoform + observedSequence.Length - 1;
                            int missingFromNterm = startInProteoform;
                            int missingFromCterm = shortestMatch.Sequence.Length - 1 - endInProteoform;

                            string missingDescription = "";
                            if (missingFromNterm > 0 && missingFromCterm > 0)
                                missingDescription = $"missing {missingFromNterm} N-term and {missingFromCterm} C-term residues";
                            else if (missingFromNterm > 0)
                                missingDescription = $"missing {missingFromNterm} N-term residues";
                            else if (missingFromCterm > 0)
                                missingDescription = $"missing {missingFromCterm} C-term residues";

                            result.MatchType = "Partial Match";
                            result.MatchedAccessions.Add(shortestMatch.Accession);
                            result.TheoreticalLength = shortestMatch.ParentProtein.BaseSequence.Length;
                            result.ProteoformLength = shortestMatch.Sequence.Length;
                            result.TotalMatchCount = containingProteoforms.Count;

                            int startInProtein = shortestMatch.StartResidue + startInProteoform;
                            int endInProtein = startInProtein + observedSequence.Length - 1;
                            result.TruncationInfo = $"Observed {startInProtein}-{endInProtein}, Proteoform {shortestMatch.StartResidue}-{shortestMatch.EndResidue} ({missingDescription})";
                            result.ProteoformStart = startInProtein;
                            result.ProteoformEnd = endInProtein;
                        }
                        else
                        {
                            var containingProteins = allProteins
                                .Where(p => p.BaseSequence.Contains(observedSequence))
                                .OrderBy(p => p.BaseSequence.Length)
                                .ToList();

                            if (containingProteins.Any())
                            {
                                var shortestProtein = containingProteins.First();
                                result.MatchType = "Partial Match (no matching proteoform)";
                                result.MatchedAccessions.Add(shortestProtein.Accession);
                                result.TheoreticalLength = shortestProtein.BaseSequence.Length;
                                result.ProteoformLength = shortestProtein.BaseSequence.Length;
                                result.TotalMatchCount = containingProteins.Count;

                                int startIndex = shortestProtein.BaseSequence.IndexOf(observedSequence) + 1;
                                int endIndex = startIndex + observedSequence.Length - 1;
                                result.TruncationInfo = $"Found in protein at {startIndex}-{endIndex} of {shortestProtein.BaseSequence.Length}, but no matching proteoform";
                                result.ProteoformStart = startIndex;
                                result.ProteoformEnd = endIndex;
                            }
                            else
                            {
                                result.MatchType = "No Match";
                            }
                        }
                    }

                    results.Add(result);

                    if (processed % 100 == 0)
                        Log($"  Mapped {processed} of {observedProteoforms.Count} proteoforms...");
                }

                // === OUTPUT RESULTS ===
                Log("\n=== WRITING OUTPUT ===");
                try
                {
                    using (var writer = new StreamWriter(outputPath))
                    {
                        writer.WriteLine("SourceFile\tLineNumber\tObservedAccession\tOriginalSequence\tProcessedSequence\tObservedLength\tMatchType\tMatchedAccession\tProteoformStart\tProteoformEnd\tProteoformLength\tTheoreticalProteinLength\tDetails\tTotalMatchCount");

                        foreach (var result in results)
                        {
                            writer.WriteLine($"{result.SourceFile}\t" +
                                             $"{result.LineNumber}\t" +
                                             $"{result.ObservedAccession}\t" +
                                             $"{result.OriginalSequence}\t" +
                                             $"{result.ObservedSequence}\t" +
                                             $"{result.ObservedSequence.Length}\t" +
                                             $"{result.MatchType}\t" +
                                             $"{string.Join(";", result.MatchedAccessions)}\t" +
                                             $"{result.ProteoformStart}\t" +
                                             $"{result.ProteoformEnd}\t" +
                                             $"{result.ProteoformLength}\t" +
                                             $"{result.TheoreticalLength}\t" +
                                             $"{result.TruncationInfo}\t" +
                                             $"{result.TotalMatchCount}");
                        }
                    }
                    Log($"Results written to: {outputPath}");
                }
                catch (Exception ex)
                {
                    Log($"ERROR writing output: {ex.Message}");
                }

                // === SUMMARY ===
                Log("\n=== MAPPING SUMMARY ===");
                Log($"Total observed proteoforms: {results.Count}");
                Log($"Full-Length matches: {results.Count(r => r.MatchType == "Full-Length")}");
                Log($"Full-Length (minus Met) matches: {results.Count(r => r.MatchType == "Full-Length (minus Met)")}");
                Log($"Partial matches: {results.Count(r => r.MatchType.StartsWith("Partial"))}");
                Log($"No matches: {results.Count(r => r.MatchType == "No Match")}");

                Log("\n=== SUMMARY BY SOURCE FILE ===");
                foreach (var group in results.GroupBy(r => r.SourceFile))
                {
                    Log($"\n{group.Key}:");
                    Log($"  Total: {group.Count()}");
                    Log($"  Full-Length: {group.Count(r => r.MatchType == "Full-Length")}");
                    Log($"  Full-Length (minus Met): {group.Count(r => r.MatchType == "Full-Length (minus Met)")}");
                    Log($"  Partial: {group.Count(r => r.MatchType.StartsWith("Partial"))}");
                    Log($"  No Match: {group.Count(r => r.MatchType == "No Match")}");
                }

                Log($"\n=== TEST COMPLETE ===");
                Log($"Timestamp: {DateTime.Now:yyyy-MM-dd HH:mm:ss}");
                Log($"Log file saved to: {_logPath}");

                Assert.IsTrue(results.Count > 0, "No proteoforms were processed");
                Assert.IsTrue(File.Exists(outputPath), "Output file was not created");
            }
            finally
            {
                _logWriter?.Close();
                _logWriter?.Dispose();
                _logWriter = null;
            }
        }

        /// <summary>
        /// Parses the header to find accession and sequence columns.
        /// </summary>
        private static (int accessionCol, int sequenceCol, string sequenceType, char delimiter) ParseHeader(string filePath)
        {
            string header = File.ReadLines(filePath).FirstOrDefault() ?? "";

            // Determine delimiter based on file extension
            char delimiter;
            if (filePath.EndsWith(".csv", StringComparison.OrdinalIgnoreCase))
                delimiter = ',';
            else if (filePath.EndsWith(".tsv", StringComparison.OrdinalIgnoreCase) ||
                     filePath.EndsWith(".psmtsv", StringComparison.OrdinalIgnoreCase))
                delimiter = '\t';
            else
            {
                int tabCount = header.Count(c => c == '\t');
                int commaCount = header.Count(c => c == ',');
                delimiter = tabCount >= commaCount ? '\t' : ',';
            }

            var columns = header.Split(delimiter);
            int accessionCol = -1;
            int sequenceCol = -1;
            string sequenceType = "Unknown";

            // Track candidates for sequence column by priority
            int baseSequenceCol = -1;
            int exactSequenceCol = -1;
            int proteoformCol = -1;

            // Log all columns for debugging
            Log($"  Parsing header with {columns.Length} columns:");
            for (int i = 0; i < columns.Length; i++)
            {
                string col = columns[i].Trim();
                string colLower = col.ToLower();
                Log($"    Col {i}: '{col}'");

                // Find accession column
                if (accessionCol < 0 &&
                    (colLower.Contains("accession") || (colLower.Contains("protein") && !colLower.Contains("proteoform"))))
                {
                    accessionCol = i;
                    Log($"      -> Selected as ACCESSION column");
                }

                // Collect sequence column candidates (don't select yet - collect all first)
                // Priority 1: "base sequence" (exact or contains)
                if (colLower == "base sequence" || colLower.Contains("base sequence"))
                {
                    baseSequenceCol = i;
                    Log($"      -> Candidate for SEQUENCE column (BaseSequence)");
                }
                // Priority 2: exact match "sequence" only
                else if (colLower == "sequence")
                {
                    exactSequenceCol = i;
                    Log($"      -> Candidate for SEQUENCE column (Sequence)");
                }
                // Priority 3: "proteoform" but NOT if it contains "confidence", "score", "characterization", etc.
                else if (colLower.Contains("proteoform") &&
                         !colLower.Contains("confidence") &&
                         !colLower.Contains("score") &&
                         !colLower.Contains("characterization") &&
                         !colLower.Contains("count") &&
                         !colLower.Contains("number"))
                {
                    proteoformCol = i;
                    Log($"      -> Candidate for SEQUENCE column (Proteoform)");
                }
            }

            // Select sequence column by priority: BaseSequence > Sequence > Proteoform
            if (baseSequenceCol >= 0)
            {
                sequenceCol = baseSequenceCol;
                sequenceType = "BaseSequence";
                Log($"  Selected column {sequenceCol} as SEQUENCE (BaseSequence)");
            }
            else if (exactSequenceCol >= 0)
            {
                sequenceCol = exactSequenceCol;
                sequenceType = "Sequence";
                Log($"  Selected column {sequenceCol} as SEQUENCE (Sequence)");
            }
            else if (proteoformCol >= 0)
            {
                sequenceCol = proteoformCol;
                sequenceType = "Proteoform";
                Log($"  Selected column {sequenceCol} as SEQUENCE (Proteoform)");
            }

            return (accessionCol, sequenceCol, sequenceType, delimiter);
        }

        private static string ProcessSequence(string sequence, string sequenceType)
        {
            if (string.IsNullOrWhiteSpace(sequence))
                return "";

            if (sequenceType == "BaseSequence" || sequenceType == "Sequence")
                return sequence.Trim();

            if (sequenceType == "Proteoform")
            {
                string processed = sequence.Trim();

                // Remove brackets and any content inside (e.g., "PEP[123.45]TIDE" -> "PEPTIDE")
                processed = Regex.Replace(processed, @"\[[^\]]*\]", "");

                // Remove prefix: 0-1 letter followed by '.' (e.g., "M.PEPTIDE" or ".PEPTIDE" -> "PEPTIDE")
                processed = Regex.Replace(processed, @"^[A-Za-z]?\.", "");

                // Remove suffix: '.' followed by 0-1 letter (e.g., "PEPTIDE.M" or "PEPTIDE." -> "PEPTIDE")
                processed = Regex.Replace(processed, @"\.[A-Za-z]?$", "");

                // Remove parentheses
                processed = processed.Replace("(", "").Replace(")", "");

                return processed.Trim();
            }

            return sequence.Trim();
        }

        private class ProteoformMappingResult
        {
            public string SourceFile { get; set; } = "";
            public int LineNumber { get; set; }
            public string ObservedAccession { get; set; } = "";
            public string ObservedSequence { get; set; } = "";
            public string OriginalSequence { get; set; } = "";
            public string MatchType { get; set; } = "";
            public List<string> MatchedAccessions { get; set; } = new();
            public int TheoreticalLength { get; set; }
            public int ProteoformStart { get; set; }
            public int ProteoformEnd { get; set; }
            public int ProteoformLength { get; set; }
            public string TruncationInfo { get; set; } = "";
            public int TotalMatchCount { get; set; }
        }
    }
}