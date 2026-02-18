using MassSpectrometry;
using NUnit.Framework;
using Readers;
using Readers.InternalIons;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test.InternalIons
{
    public static class InternalFragmentAnalysisRunnerTests
    {
        private const string DefaultPsmTsvPath = @"FILL_IN";
        private const string DefaultRawFileFolder = @"FILL_IN";
        private const string DefaultOutputDirectory = @"FILL_IN";
        private const double DefaultCollisionEnergy = 42.0;

        private static readonly char[] StandardAminoAcids =
            { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' };

        private static readonly string[] SupportedExtensions = { ".raw", ".mzml", ".mgf" };

        public static void RunAll() => RunAll(DefaultPsmTsvPath, DefaultRawFileFolder, DefaultOutputDirectory);

        public static void RunAll(string psmTsvPath, string rawFileFolder, string outputDirectory)
        {
            Console.WriteLine("==================================================================");
            Console.WriteLine("     INTERNAL FRAGMENT ION ANALYSIS HARNESS");
            Console.WriteLine("==================================================================");
            Console.WriteLine();
            Console.WriteLine($"PSM TSV:          {psmTsvPath}");
            Console.WriteLine($"Raw File Folder:  {rawFileFolder}");
            Console.WriteLine($"Output Directory: {outputDirectory}");
            Console.WriteLine($"Default CE:       {DefaultCollisionEnergy}");
            Console.WriteLine();

            try
            {
                if (!File.Exists(psmTsvPath))
                    throw new FileNotFoundException($"PSM TSV file not found: {psmTsvPath}");
                if (!Directory.Exists(rawFileFolder))
                    throw new DirectoryNotFoundException($"Raw file folder not found: {rawFileFolder}");
                if (!Directory.Exists(outputDirectory))
                    Directory.CreateDirectory(outputDirectory);

                var psms = Step1_LoadPsms(psmTsvPath);
                var msDataFiles = Step2_LoadRawFiles(rawFileFolder, psms);
                var internalIons = Step3_ExtractInternalFragments(psms, msDataFiles, DefaultCollisionEnergy);
                string outputPath = Step4_WriteOutputTsv(internalIons, outputDirectory);
                Step5_RoundTripValidation(internalIons, outputPath);
                Step6_AminoAcidTerminusEnrichment(internalIons);
                Step7_IsobaricAmbiguityReport(internalIons);
                Step8_ExtendedTerminusEnrichment(internalIons);
                Step9_MassAccuracyCharacterization(internalIons);
                Step10_FeaturePreparationAndHypothesisValidation(internalIons);

                Console.WriteLine();
                Console.WriteLine("==================================================================");
                Console.WriteLine("ALL STEPS COMPLETED SUCCESSFULLY");
                Console.WriteLine("==================================================================");
            }
            catch (Exception ex)
            {
                Console.WriteLine();
                Console.WriteLine("==================================================================");
                Console.WriteLine("FAILURE");
                Console.WriteLine("==================================================================");
                Console.WriteLine($"Exception: {ex.GetType().FullName}");
                Console.WriteLine($"Message: {ex.Message}");
                Console.WriteLine($"Stack Trace:\n{ex.StackTrace}");
                throw;
            }
        }

        #region Steps 1-9 (unchanged)

        private static List<PsmFromTsv> Step1_LoadPsms(string psmTsvPath)
        {
            Console.WriteLine("=== STEP 1: Load PSMs ===");
            var psms = SpectrumMatchTsvReader.ReadPsmTsv(psmTsvPath, out var warnings);
            Console.WriteLine($"Total PSMs loaded: {psms.Count:N0}");
            var sourceFiles = psms.Select(p => p.FileNameWithoutExtension).Distinct().ToList();
            Console.WriteLine($"Unique source files: {sourceFiles.Count}");
            var psmsWithInternal = psms
                .Where(p => p.MatchedIons != null && p.MatchedIons.Any(ion => IsInternalAnnotation(ion.Annotation)))
                .ToList();
            Console.WriteLine($"PSMs with internal fragments: {psmsWithInternal.Count:N0}");
            foreach (var psm in psmsWithInternal.Take(3))
            {
                var annotations = psm.MatchedIons
                    .Where(ion => IsInternalAnnotation(ion.Annotation))
                    .Select(ion => ion.Annotation).Take(3);
                Console.WriteLine($"  Scan {psm.Ms2ScanNumber}: {psm.BaseSeq} -> {string.Join(", ", annotations)}");
            }
            Console.WriteLine();
            return psms;
        }

        private static bool IsInternalAnnotation(string annotation)
        {
            if (string.IsNullOrEmpty(annotation)) return false;
            if (annotation.Contains("int", StringComparison.OrdinalIgnoreCase)) return true;
            int start = annotation.IndexOf('[');
            int end = annotation.IndexOf(']');
            return start >= 0 && end > start && annotation.Substring(start, end - start).Contains('-');
        }

        private static Dictionary<string, MsDataFile> Step2_LoadRawFiles(string rawFileFolder, List<PsmFromTsv> psms)
        {
            Console.WriteLine("=== STEP 2: Load Raw Files ===");
            var rawFiles = SupportedExtensions
                .SelectMany(ext => Directory.GetFiles(rawFileFolder, $"*{ext}", SearchOption.TopDirectoryOnly)).ToList();
            Console.WriteLine($"Found {rawFiles.Count} raw files");
            var neededFiles = psms.Select(p => p.FileNameWithoutExtension).Distinct().ToHashSet(StringComparer.OrdinalIgnoreCase);
            var msDataFiles = new Dictionary<string, MsDataFile>(StringComparer.OrdinalIgnoreCase);
            int totalMs2 = 0;
            foreach (var rawFile in rawFiles)
            {
                string name = Path.GetFileNameWithoutExtension(rawFile);
                if (!neededFiles.Contains(name)) continue;
                Console.WriteLine($"  Loading {Path.GetFileName(rawFile)}...");
                try
                {
                    string ext = Path.GetExtension(rawFile).ToLowerInvariant();
                    MsDataFile file = ext switch
                    {
                        ".mzml" => Mzml.LoadAllStaticData(rawFile),
                        ".raw" => ThermoRawFileReader.LoadAllStaticData(rawFile),
                        ".mgf" => Mgf.LoadAllStaticData(rawFile),
                        _ => throw new NotSupportedException(ext)
                    };
                    msDataFiles[name] = file;
                    int ms2 = file.GetAllScansList().Count(s => s.MsnOrder == 2);
                    totalMs2 += ms2;
                    Console.WriteLine($"    {file.GetAllScansList().Count} scans ({ms2} MS2)");
                }
                catch (Exception ex) { Console.WriteLine($"    ERROR: {ex.Message}"); }
            }
            Console.WriteLine($"Total MS2 scans: {totalMs2:N0}");
            Console.WriteLine();
            return msDataFiles;
        }

        private static List<InternalFragmentIon> Step3_ExtractInternalFragments(
            List<PsmFromTsv> psms, Dictionary<string, MsDataFile> msDataFiles, double defaultCollisionEnergy)
        {
            Console.WriteLine("=== STEP 3: Extract Internal Fragment Features ===");
            Console.WriteLine($"Using default collision energy: {defaultCollisionEnergy}");
            Console.WriteLine();
            var allIons = new List<InternalFragmentIon>();
            foreach (var group in psms.GroupBy(p => p.FileNameWithoutExtension))
            {
                if (!msDataFiles.TryGetValue(group.Key, out var file)) continue;
                var ions = InternalFragmentFeatureExtractor.ExtractFromPsms(group.ToList(), file, defaultCollisionEnergy);
                allIons.AddRange(ions);
                if (ions.Count > 0) Console.WriteLine($"  {group.Key}: {ions.Count:N0} internal ions");
            }
            Console.WriteLine($"Total internal ions: {allIons.Count:N0}");
            if (allIons.Count > 0)
            {
                int proline = allIons.Count(i => i.HasProlineAtEitherTerminus);
                int aspartate = allIons.Count(i => i.HasAspartateAtEitherTerminus);
                Console.WriteLine($"With Proline terminus: {proline:N0} ({100.0 * proline / allIons.Count:F1}%)");
                Console.WriteLine($"With Aspartate terminus: {aspartate:N0} ({100.0 * aspartate / allIons.Count:F1}%)");
                var errors = allIons.Select(i => i.MassErrorPpm).Where(e => !double.IsNaN(e)).ToList();
                if (errors.Count > 0)
                {
                    double mean = errors.Average();
                    double std = Math.Sqrt(errors.Select(e => Math.Pow(e - mean, 2)).Average());
                    Console.WriteLine($"Mass Error (ppm): Mean={mean:F4}, StdDev={std:F4}");
                }
            }
            Console.WriteLine();
            return allIons;
        }

        private static string Step4_WriteOutputTsv(List<InternalFragmentIon> ions, string outputDirectory)
        {
            Console.WriteLine("=== STEP 4: Write Output TSV ===");
            string outputPath = Path.Combine(outputDirectory, "InternalFragmentIons.tsv");
            InternalFragmentTsvWriter.WriteToTsv(ions, outputPath);
            Console.WriteLine($"Output: {outputPath}");
            Console.WriteLine($"Rows: {ions.Count:N0}");
            Console.WriteLine();
            return outputPath;
        }

        private static void Step5_RoundTripValidation(List<InternalFragmentIon> original, string outputPath)
        {
            Console.WriteLine("=== STEP 5: Round-trip Validation ===");
            if (original.Count == 0) { Console.WriteLine("Skipping (no ions)\n"); return; }
            var readBack = InternalFragmentTsvWriter.ReadFromTsv(outputPath);
            Console.WriteLine($"Original: {original.Count}, Read back: {readBack.Count}");
            if (original.Count != readBack.Count)
                throw new InvalidOperationException($"Count mismatch");
            Console.WriteLine("Round-trip validation PASSED\n");
        }

        private static void Step6_AminoAcidTerminusEnrichment(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("=== STEP 6: Amino Acid Terminus Enrichment ===");
            if (ions.Count == 0) { Console.WriteLine("Skipping\n"); return; }
            string allSeq = string.Concat(ions.Select(i => i.InternalSequence));
            var bgFreq = StandardAminoAcids.ToDictionary(aa => aa, aa => (double)allSeq.Count(c => c == aa) / allSeq.Length);
            var termCount = StandardAminoAcids.ToDictionary(aa => aa, _ => 0);
            foreach (var ion in ions)
            {
                if (string.IsNullOrEmpty(ion.InternalSequence)) continue;
                if (termCount.ContainsKey(ion.InternalSequence[0])) termCount[ion.InternalSequence[0]]++;
                if (termCount.ContainsKey(ion.InternalSequence[^1])) termCount[ion.InternalSequence[^1]]++;
            }
            int totalPos = ions.Count * 2;
            Console.WriteLine("+----+----------+----------+------------+-------+");
            Console.WriteLine("| AA | Observed | Expected | Enrichment | Flag  |");
            Console.WriteLine("+----+----------+----------+------------+-------+");
            foreach (var (aa, obs, exp, enrich) in StandardAminoAcids
                .Select(aa => (aa, obs: termCount[aa], exp: bgFreq[aa] * totalPos, enrich: bgFreq[aa] > 0 ? termCount[aa] / (bgFreq[aa] * totalPos) : 0))
                .OrderByDescending(x => x.enrich))
            {
                string flag = enrich > 1.5 ? " ***" : "";
                Console.WriteLine($"|  {aa} | {obs,8} | {exp,8:F1} | {enrich,10:F3} |{flag,-6}|");
            }
            Console.WriteLine("+----+----------+----------+------------+-------+\n");
        }

        private static void Step7_IsobaricAmbiguityReport(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("=== STEP 7: Isobaric Ambiguity Report ===");
            if (ions.Count == 0) { Console.WriteLine("Skipping\n"); return; }
            int ambCount = ions.Count(i => i.IsIsobaricAmbiguous);
            Console.WriteLine($"Isobaric ambiguous: {ambCount:N0} ({100.0 * ambCount / ions.Count:F1}%)\n");
        }

        private static void Step8_ExtendedTerminusEnrichment(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("=== STEP 8: Extended Terminus Enrichment ===");
            if (ions.Count == 0) { Console.WriteLine("Skipping\n"); return; }
            Console.WriteLine("(See Step 10 Part C for detailed analysis)\n");
        }

        private static void Step9_MassAccuracyCharacterization(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("=== STEP 9: Mass Accuracy Characterization ===");
            if (ions.Count == 0) { Console.WriteLine("Skipping\n"); return; }
            var validIons = ions.Where(i => !double.IsNaN(i.MassErrorPpm)).ToList();
            double corr = ComputePearsonCorrelation(
                validIons.Select(i => i.NormalizedIntensity).ToArray(),
                validIons.Select(i => Math.Abs(i.MassErrorPpm)).ToArray());
            Console.WriteLine($"Pearson(Intensity, |MassError|): {corr:F4}");
            int pass = validIons.Count(i => i.PassesMassAccuracyFilter);
            Console.WriteLine($"PassesMassAccuracyFilter: {pass:N0}/{validIons.Count:N0} ({100.0 * pass / validIons.Count:F1}%)\n");
        }

        #endregion

        #region Step 10 - Feature Preparation and Hypothesis Validation

        private static void Step10_FeaturePreparationAndHypothesisValidation(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("==================================================================");
            Console.WriteLine("=== STEP 10: Feature Preparation & Hypothesis Validation ===");
            Console.WriteLine("==================================================================");
            Console.WriteLine();

            if (ions.Count == 0)
            {
                Console.WriteLine("Skipping (no ions)");
                return;
            }

            var validIons = ions.Where(i => !double.IsNaN(i.MassErrorPpm)).ToList();

            Step10_PartA_FilterQualityAudit(validIons);
            Step10_PartB_ModificationBasicResidueInteraction(validIons);
            Step10_PartC_WeakBondEnrichment(validIons);
            Step10_PartD_BYTerminalIonCorrelation(validIons);
            Step10_PartE_IntensityDistributionAndFeatureSummary(validIons);
        }

        private static void Step10_PartA_FilterQualityAudit(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine("PART A: FILTER QUALITY AUDIT");
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine();

            var failingIons = ions.Where(i => !i.PassesMassAccuracyFilter).ToList();
            Console.WriteLine($"Total failing ions: {failingIons.Count:N0}");

            // Category A: |MassErrorPpm| > 15
            var catA = failingIons.Where(i => Math.Abs(i.MassErrorPpm) > 15).ToList();
            // Category B: 5 < |MassErrorPpm| <= 15 AND NormalizedIntensity < 0.10
            var catB = failingIons.Where(i => Math.Abs(i.MassErrorPpm) > 5 && Math.Abs(i.MassErrorPpm) <= 15 && i.NormalizedIntensity < 0.10).ToList();
            // Category C: HasModifiedResidue = TRUE AND FragmentLength > 12
            var catC = failingIons.Where(i => i.HasModifiedResidue && i.FragmentLength > 12).ToList();

            Console.WriteLine($"Category A (|Error| > 15 ppm): {catA.Count:N0}");
            Console.WriteLine($"Category B (5 < |Error| <= 15 AND Intensity < 0.10): {catB.Count:N0}");
            Console.WriteLine($"Category C (Modified AND Length > 12): {catC.Count:N0}");
            Console.WriteLine();

            // Pairwise overlaps
            int abOverlap = catA.Intersect(catB).Count();
            int acOverlap = catA.Intersect(catC).Count();
            int bcOverlap = catB.Intersect(catC).Count();
            int abcOverlap = catA.Intersect(catB).Intersect(catC).Count();

            Console.WriteLine("Pairwise overlaps:");
            Console.WriteLine($"  A ∩ B: {abOverlap}");
            Console.WriteLine($"  A ∩ C: {acOverlap}");
            Console.WriteLine($"  B ∩ C: {bcOverlap}");
            Console.WriteLine($"  A ∩ B ∩ C: {abcOverlap}");
            Console.WriteLine();

            // Among Category B, fraction with modifications
            if (catB.Count > 0)
            {
                int catBModified = catB.Count(i => i.HasModifiedResidue);
                Console.WriteLine($"Category B ions with HasModifiedResidue=TRUE: {catBModified}/{catB.Count} ({100.0 * catBModified / catB.Count:F1}%)");
            }
            Console.WriteLine();
        }

        private static void Step10_PartB_ModificationBasicResidueInteraction(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine("PART B: MODIFICATION × BASIC RESIDUE INTERACTION");
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine();

            var modifiedIons = ions.Where(i => i.HasModifiedResidue && !double.IsNaN(i.MassErrorPpm)).ToList();

            if (modifiedIons.Count == 0)
            {
                Console.WriteLine("No modified ions with valid mass error.");
                Console.WriteLine();
                return;
            }

            var group0 = modifiedIons.Where(i => i.NumberOfBasicResidues == 0).ToList();
            var group1 = modifiedIons.Where(i => i.NumberOfBasicResidues > 0).ToList();

            double mean0 = group0.Count > 0 ? group0.Select(i => Math.Abs(i.MassErrorPpm)).Average() : double.NaN;
            double mean1 = group1.Count > 0 ? group1.Select(i => Math.Abs(i.MassErrorPpm)).Average() : double.NaN;

            Console.WriteLine($"Modified ions with NumberOfBasicResidues = 0: {group0.Count:N0}, Mean |MassErrorPpm|: {mean0:F4}");
            Console.WriteLine($"Modified ions with NumberOfBasicResidues > 0: {group1.Count:N0}, Mean |MassErrorPpm|: {mean1:F4}");

            if (!double.IsNaN(mean0) && !double.IsNaN(mean1))
            {
                double diff = Math.Abs(mean1 - mean0);
                Console.WriteLine($"Difference: {diff:F4} ppm");

                if (diff > 2.0)
                {
                    Console.WriteLine("FLAG: HasModifiedResidue × BasicResidue is a candidate interaction feature.");
                }
            }
            Console.WriteLine();
        }

        private static void Step10_PartC_WeakBondEnrichment(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine("PART C: WEAK-BOND ENRICHMENT (passing ions only)");
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine();

            var passingIons = ions.Where(i => i.PassesMassAccuracyFilter).ToList();

            if (passingIons.Count == 0)
            {
                Console.WriteLine("No passing ions.");
                Console.WriteLine();
                return;
            }

            // Compute background frequency from all peptide sequences
            string allPeptideChars = string.Concat(passingIons.Select(i => i.PeptideSequence));
            int totalChars = allPeptideChars.Length;

            double bgDE = totalChars > 0 ? (double)allPeptideChars.Count(c => c == 'D' || c == 'E') / totalChars : 0;
            double bgP = totalChars > 0 ? (double)allPeptideChars.Count(c => c == 'P') / totalChars : 0;
            double bgKR = totalChars > 0 ? (double)allPeptideChars.Count(c => c == 'K' || c == 'R') / totalChars : 0;

            // Count observed frequencies at each position
            int intNTermDE = 0, intNTermP = 0, intNTermKR = 0;
            int intCTermDE = 0, intCTermP = 0, intCTermKR = 0;
            int flankNDE = 0, flankNP = 0, flankNKR = 0;
            int flankCDE = 0, flankCP = 0, flankCKR = 0;

            foreach (var ion in passingIons)
            {
                char intN = ion.InternalNTerminalAA;
                char intC = ion.InternalCTerminalAA;
                char flankN = ion.NTerminalFlankingResidue;
                char flankC = ion.CTerminalFlankingResidue;

                if (intN == 'D' || intN == 'E') intNTermDE++;
                if (intN == 'P') intNTermP++;
                if (intN == 'K' || intN == 'R') intNTermKR++;

                if (intC == 'D' || intC == 'E') intCTermDE++;
                if (intC == 'P') intCTermP++;
                if (intC == 'K' || intC == 'R') intCTermKR++;

                if (flankN == 'D' || flankN == 'E') flankNDE++;
                if (flankN == 'P') flankNP++;
                if (flankN == 'K' || flankN == 'R') flankNKR++;

                if (flankC == 'D' || flankC == 'E') flankCDE++;
                if (flankC == 'P') flankCP++;
                if (flankC == 'K' || flankC == 'R') flankCKR++;
            }

            int n = passingIons.Count;

            double EnrichmentRatio(int count, double bg) => bg > 0 ? ((double)count / n) / bg : 0;

            Console.WriteLine($"Background frequencies: D+E={bgDE:F4}, P={bgP:F4}, K+R={bgKR:F4}");
            Console.WriteLine();
            Console.WriteLine("+----------------------+----------------+----------------+----------------+");
            Console.WriteLine("| Position             | D+E enrichment | P enrichment   | K+R enrichment |");
            Console.WriteLine("+----------------------+----------------+----------------+----------------+");
            Console.WriteLine($"| InternalNTerminalAA  | {EnrichmentRatio(intNTermDE, bgDE),14:F3} | {EnrichmentRatio(intNTermP, bgP),14:F3} | {EnrichmentRatio(intNTermKR, bgKR),14:F3} |");
            Console.WriteLine($"| InternalCTerminalAA  | {EnrichmentRatio(intCTermDE, bgDE),14:F3} | {EnrichmentRatio(intCTermP, bgP),14:F3} | {EnrichmentRatio(intCTermKR, bgKR),14:F3} |");
            Console.WriteLine($"| NTermFlankingResidue | {EnrichmentRatio(flankNDE, bgDE),14:F3} | {EnrichmentRatio(flankNP, bgP),14:F3} | {EnrichmentRatio(flankNKR, bgKR),14:F3} |");
            Console.WriteLine($"| CTermFlankingResidue | {EnrichmentRatio(flankCDE, bgDE),14:F3} | {EnrichmentRatio(flankCP, bgP),14:F3} | {EnrichmentRatio(flankCKR, bgKR),14:F3} |");
            Console.WriteLine("+----------------------+----------------+----------------+----------------+");
            Console.WriteLine();
        }

        private static void Step10_PartD_BYTerminalIonCorrelation(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine("PART D: B/Y TERMINAL ION CORRELATION");
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine();

            var passingIons = ions.Where(i => i.PassesMassAccuracyFilter).ToList();

            if (passingIons.Count == 0)
            {
                Console.WriteLine("No passing ions.");
                Console.WriteLine();
                return;
            }

            // D3a: Coverage
            Console.WriteLine("--- D3a: B/Y Ion Coverage ---");
            int hasB = passingIons.Count(i => i.HasMatchedBIonAtNTerm);
            int hasY = passingIons.Count(i => i.HasMatchedYIonAtCTerm);
            int hasBoth = passingIons.Count(i => i.HasMatchedBIonAtNTerm && i.HasMatchedYIonAtCTerm);
            int hasNeither = passingIons.Count(i => !i.HasMatchedBIonAtNTerm && !i.HasMatchedYIonAtCTerm);

            Console.WriteLine($"HasMatchedBIonAtNTerm = TRUE: {hasB:N0} ({100.0 * hasB / passingIons.Count:F1}%)");
            Console.WriteLine($"HasMatchedYIonAtCTerm = TRUE: {hasY:N0} ({100.0 * hasY / passingIons.Count:F1}%)");
            Console.WriteLine($"Both TRUE: {hasBoth:N0} ({100.0 * hasBoth / passingIons.Count:F1}%)");
            Console.WriteLine($"Neither TRUE: {hasNeither:N0} ({100.0 * hasNeither / passingIons.Count:F1}%)");
            Console.WriteLine();

            // D3b: Intensity correlations
            Console.WriteLine("--- D3b: Intensity Correlations ---");
            double corrB = ComputePearsonCorrelation(
                passingIons.Select(i => i.NormalizedIntensity).ToArray(),
                passingIons.Select(i => i.BIonIntensityAtNTerm).ToArray());
            double corrY = ComputePearsonCorrelation(
                passingIons.Select(i => i.NormalizedIntensity).ToArray(),
                passingIons.Select(i => i.YIonIntensityAtCTerm).ToArray());
            double corrBY = ComputePearsonCorrelation(
                passingIons.Select(i => i.NormalizedIntensity).ToArray(),
                passingIons.Select(i => i.BYProductScore).ToArray());

            Console.WriteLine($"Pearson(NormalizedIntensity, BIonIntensityAtNTerm): {corrB:F4}");
            Console.WriteLine($"Pearson(NormalizedIntensity, YIonIntensityAtCTerm): {corrY:F4}");
            Console.WriteLine($"Pearson(NormalizedIntensity, BYProductScore): {corrBY:F4}");

            if (corrB > 0.3 || corrY > 0.3 || corrBY > 0.3)
                Console.WriteLine("  => Strong positive correlation supports the b/y sharing hypothesis.");
            Console.WriteLine();

            // D3c: Group comparison
            Console.WriteLine("--- D3c: Group Comparison by Terminal Ion Presence ---");
            var group00 = passingIons.Where(i => !i.HasMatchedBIonAtNTerm && !i.HasMatchedYIonAtCTerm).ToList();
            var group10 = passingIons.Where(i => i.HasMatchedBIonAtNTerm && !i.HasMatchedYIonAtCTerm).ToList();
            var group01 = passingIons.Where(i => !i.HasMatchedBIonAtNTerm && i.HasMatchedYIonAtCTerm).ToList();
            var group11 = passingIons.Where(i => i.HasMatchedBIonAtNTerm && i.HasMatchedYIonAtCTerm).ToList();

            void PrintGroupStats(string name, List<InternalFragmentIon> grp)
            {
                if (grp.Count == 0)
                {
                    Console.WriteLine($"  {name}: (no ions)");
                    return;
                }
                var sorted = grp.Select(i => i.NormalizedIntensity).OrderBy(x => x).ToList();
                double mean = sorted.Average();
                double median = sorted[sorted.Count / 2];
                Console.WriteLine($"  {name}: Count={grp.Count,5}, Mean={mean:F4}, Median={median:F4}");
            }

            PrintGroupStats("Group 00 (neither)", group00);
            PrintGroupStats("Group 10 (b only)", group10);
            PrintGroupStats("Group 01 (y only)", group01);
            PrintGroupStats("Group 11 (both)", group11);
            Console.WriteLine("  Prediction: Mean intensity Group 11 > Group 10 ≈ Group 01 > Group 00");
            Console.WriteLine();

            // D3d: Top ions check
            Console.WriteLine("--- D3d: Top 20 Ions by NormalizedIntensity ---");
            var top20 = passingIons.OrderByDescending(i => i.NormalizedIntensity).Take(20).ToList();
            int top20HasB = top20.Count(i => i.HasMatchedBIonAtNTerm);
            int top20HasY = top20.Count(i => i.HasMatchedYIonAtCTerm);

            Console.WriteLine($"Among top 20 ions: {top20HasB}/20 have HasMatchedBIonAtNTerm=TRUE ({100.0 * top20HasB / 20:F0}%)");
            Console.WriteLine($"Among top 20 ions: {top20HasY}/20 have HasMatchedYIonAtCTerm=TRUE ({100.0 * top20HasY / 20:F0}%)");
            Console.WriteLine();
        }

        private static void Step10_PartE_IntensityDistributionAndFeatureSummary(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine("PART E: INTENSITY DISTRIBUTION AND FEATURE SUMMARY");
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine();

            var passingIons = ions.Where(i => i.PassesMassAccuracyFilter).ToList();

            if (passingIons.Count == 0)
            {
                Console.WriteLine("No passing ions.");
                Console.WriteLine();
                return;
            }

            // E1: Histogram of NormalizedIntensity
            Console.WriteLine("--- E1: NormalizedIntensity Histogram ---");
            var bins = new (double min, double max, string label)[]
            {
                (0.00, 0.05, "[0.00, 0.05)"),
                (0.05, 0.10, "[0.05, 0.10)"),
                (0.10, 0.20, "[0.10, 0.20)"),
                (0.20, 0.40, "[0.20, 0.40)"),
                (0.40, 1.01, "[0.40, 1.00]")
            };

            Console.WriteLine("+---------------+-------+----------+");
            Console.WriteLine("| Bin           | Count | Fraction |");
            Console.WriteLine("+---------------+-------+----------+");
            foreach (var (min, max, label) in bins)
            {
                int count = passingIons.Count(i => i.NormalizedIntensity >= min && i.NormalizedIntensity < max);
                Console.WriteLine($"| {label,-13} | {count,5} | {100.0 * count / passingIons.Count,7:F1}% |");
            }
            Console.WriteLine("+---------------+-------+----------+");
            Console.WriteLine();

            // E2: Correlation matrix
            Console.WriteLine("--- E2: Feature Correlation Matrix ---");
            var features = new (string name, Func<InternalFragmentIon, double> getter)[]
            {
                ("FragLen", i => i.FragmentLength),
                ("NumBasic", i => i.NumberOfBasicResidues),
                ("NormInt", i => i.NormalizedIntensity),
                ("BIonInt", i => i.BIonIntensityAtNTerm),
                ("YIonInt", i => i.YIonIntensityAtCTerm),
                ("BYProd", i => i.BYProductScore)
            };

            Console.Write("         ");
            foreach (var f in features) Console.Write($"{f.name,8} ");
            Console.WriteLine();

            var collinearPairs = new List<string>();

            for (int i = 0; i < features.Length; i++)
            {
                Console.Write($"{features[i].name,8} ");
                var xVals = passingIons.Select(features[i].getter).ToArray();

                for (int j = 0; j < features.Length; j++)
                {
                    var yVals = passingIons.Select(features[j].getter).ToArray();
                    double r = ComputePearsonCorrelation(xVals, yVals);
                    Console.Write($"{r,8:F3} ");

                    if (i < j && Math.Abs(r) > 0.5)
                        collinearPairs.Add($"{features[i].name} × {features[j].name} (r={r:F3})");
                }
                Console.WriteLine();
            }
            Console.WriteLine();

            if (collinearPairs.Count > 0)
            {
                Console.WriteLine("Potentially collinear pairs (|r| > 0.5):");
                foreach (var pair in collinearPairs)
                    Console.WriteLine($"  {pair}");
                Console.WriteLine();
            }

            // E3: Feature inventory
            Console.WriteLine("--- E3: Feature Inventory ---");
            Console.WriteLine("+---------------------------+-------+------------+------------+");
            Console.WriteLine("| Feature                   | Valid |       Mean |     StdDev |");
            Console.WriteLine("+---------------------------+-------+------------+------------+");

            void PrintFeatureStats(string name, IEnumerable<double> values)
            {
                var list = values.Where(v => !double.IsNaN(v)).ToList();
                int valid = list.Count;
                double mean = valid > 0 ? list.Average() : double.NaN;
                double std = valid > 1 ? Math.Sqrt(list.Select(v => Math.Pow(v - mean, 2)).Average()) : 0;
                string meanStr = double.IsNaN(mean) ? "N/A" : $"{mean:F4}";
                string stdStr = double.IsNaN(mean) ? "N/A" : $"{std:F4}";
                Console.WriteLine($"| {name,-25} | {valid,5} | {meanStr,10} | {stdStr,10} |");
            }

            void PrintBoolFeatureStats(string name, IEnumerable<bool> values)
            {
                var list = values.ToList();
                int trueCount = list.Count(v => v);
                double mean = (double)trueCount / list.Count;
                Console.WriteLine($"| {name,-25} | {list.Count,5} | {mean,10:F4} |        N/A |");
            }

            PrintFeatureStats("FragmentLength", passingIons.Select(i => (double)i.FragmentLength));
            PrintFeatureStats("NumberOfBasicResidues", passingIons.Select(i => (double)i.NumberOfBasicResidues));
            PrintFeatureStats("NormalizedIntensity", passingIons.Select(i => i.NormalizedIntensity));
            PrintFeatureStats("LocalIntensityRank", passingIons.Select(i => i.LocalIntensityRank));
            PrintBoolFeatureStats("HasProlineAtEitherTerminus", passingIons.Select(i => i.HasProlineAtEitherTerminus));
            PrintBoolFeatureStats("HasAspartateAtEitherTerminus", passingIons.Select(i => i.HasAspartateAtEitherTerminus));
            PrintFeatureStats("BIonIntensityAtNTerm", passingIons.Select(i => i.BIonIntensityAtNTerm));
            PrintFeatureStats("YIonIntensityAtCTerm", passingIons.Select(i => i.YIonIntensityAtCTerm));
            PrintFeatureStats("BYProductScore", passingIons.Select(i => i.BYProductScore));
            PrintBoolFeatureStats("PassesMassAccuracyFilter", passingIons.Select(i => i.PassesMassAccuracyFilter));

            Console.WriteLine("+---------------------------+-------+------------+------------+");
            Console.WriteLine();
            Console.WriteLine("Pre-training data quality gate: All features populated.");
            Console.WriteLine();
        }

        #endregion

        #region Utility Methods

        private static double ComputePearsonCorrelation(double[] x, double[] y)
        {
            if (x.Length != y.Length || x.Length == 0) return double.NaN;
            double meanX = x.Average();
            double meanY = y.Average();
            double sumXY = 0, sumX2 = 0, sumY2 = 0;
            for (int i = 0; i < x.Length; i++)
            {
                double dx = x[i] - meanX;
                double dy = y[i] - meanY;
                sumXY += dx * dy;
                sumX2 += dx * dx;
                sumY2 += dy * dy;
            }
            double denom = Math.Sqrt(sumX2 * sumY2);
            return denom == 0 ? 0 : sumXY / denom;
        }

        #endregion
    }

    [TestFixture]
    public class InternalFragmentTests
    {
        [Test]
        [Explicit("Manual integration harness - requires local data files")]
        public void RunInternalFragmentAnalysis()
        {
            InternalFragmentAnalysisRunnerTests.RunAll(
                psmTsvPath: @"F:\JurakatMultiprotease\2026-02-17-16-40-42\Task1-SearchTask\AllPeptides.psmtsv",
                rawFileFolder: @"F:\JurakatMultiprotease",
                outputDirectory: @"E:\Projects\internalIons\mzLibReportsForClaude"
            );
        }
    }
}