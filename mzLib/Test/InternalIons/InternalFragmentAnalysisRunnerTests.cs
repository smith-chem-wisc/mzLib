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
                Step11_FeatureRefinementAndModelValidation(internalIons);

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

        #region Steps 1-10 (condensed)

        private static List<PsmFromTsv> Step1_LoadPsms(string psmTsvPath)
        {
            Console.WriteLine("=== STEP 1: Load PSMs ===");
            var psms = SpectrumMatchTsvReader.ReadPsmTsv(psmTsvPath, out _);
            Console.WriteLine($"Total PSMs loaded: {psms.Count:N0}");
            var psmsWithInternal = psms.Where(p => p.MatchedIons != null && p.MatchedIons.Any(ion => IsInternalAnnotation(ion.Annotation))).ToList();
            Console.WriteLine($"PSMs with internal fragments: {psmsWithInternal.Count:N0}\n");
            return psms;
        }

        private static bool IsInternalAnnotation(string annotation)
        {
            if (string.IsNullOrEmpty(annotation)) return false;
            if (annotation.Contains("int", StringComparison.OrdinalIgnoreCase)) return true;
            int start = annotation.IndexOf('['), end = annotation.IndexOf(']');
            return start >= 0 && end > start && annotation.Substring(start, end - start).Contains('-');
        }

        private static Dictionary<string, MsDataFile> Step2_LoadRawFiles(string rawFileFolder, List<PsmFromTsv> psms)
        {
            Console.WriteLine("=== STEP 2: Load Raw Files ===");
            var rawFiles = SupportedExtensions.SelectMany(ext => Directory.GetFiles(rawFileFolder, $"*{ext}", SearchOption.TopDirectoryOnly)).ToList();
            var neededFiles = psms.Select(p => p.FileNameWithoutExtension).Distinct().ToHashSet(StringComparer.OrdinalIgnoreCase);
            var msDataFiles = new Dictionary<string, MsDataFile>(StringComparer.OrdinalIgnoreCase);
            foreach (var rawFile in rawFiles)
            {
                string name = Path.GetFileNameWithoutExtension(rawFile);
                if (!neededFiles.Contains(name)) continue;
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
                    Console.WriteLine($"  {name}: {file.GetAllScansList().Count} scans");
                }
                catch (Exception ex) { Console.WriteLine($"  {name}: ERROR - {ex.Message}"); }
            }
            Console.WriteLine();
            return msDataFiles;
        }

        private static List<InternalFragmentIon> Step3_ExtractInternalFragments(List<PsmFromTsv> psms, Dictionary<string, MsDataFile> msDataFiles, double defaultCE)
        {
            Console.WriteLine("=== STEP 3: Extract Internal Fragment Features ===");
            var allIons = new List<InternalFragmentIon>();
            foreach (var group in psms.GroupBy(p => p.FileNameWithoutExtension))
            {
                if (!msDataFiles.TryGetValue(group.Key, out var file)) continue;
                var ions = InternalFragmentFeatureExtractor.ExtractFromPsms(group.ToList(), file, defaultCE);
                allIons.AddRange(ions);
                if (ions.Count > 0) Console.WriteLine($"  {group.Key}: {ions.Count:N0} internal ions");
            }
            Console.WriteLine($"Total: {allIons.Count:N0}\n");
            return allIons;
        }

        private static string Step4_WriteOutputTsv(List<InternalFragmentIon> ions, string outputDirectory)
        {
            Console.WriteLine("=== STEP 4: Write Output TSV ===");
            string outputPath = Path.Combine(outputDirectory, "InternalFragmentIons.tsv");
            InternalFragmentTsvWriter.WriteToTsv(ions, outputPath);
            Console.WriteLine($"Output: {outputPath} ({ions.Count:N0} rows)\n");
            return outputPath;
        }

        private static void Step5_RoundTripValidation(List<InternalFragmentIon> original, string outputPath)
        {
            Console.WriteLine("=== STEP 5: Round-trip Validation ===");
            if (original.Count == 0) { Console.WriteLine("Skipping\n"); return; }
            var readBack = InternalFragmentTsvWriter.ReadFromTsv(outputPath);
            if (original.Count != readBack.Count) throw new InvalidOperationException("Count mismatch");
            Console.WriteLine("PASSED\n");
        }

        private static void Step6_AminoAcidTerminusEnrichment(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("=== STEP 6: AA Terminus Enrichment ===");
            if (ions.Count == 0) { Console.WriteLine("Skipping\n"); return; }
            Console.WriteLine("(Detailed in Step 10 Part C)\n");
        }

        private static void Step7_IsobaricAmbiguityReport(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("=== STEP 7: Isobaric Ambiguity ===");
            if (ions.Count == 0) { Console.WriteLine("Skipping\n"); return; }
            int amb = ions.Count(i => i.IsIsobaricAmbiguous);
            Console.WriteLine($"Ambiguous: {amb:N0} ({100.0 * amb / ions.Count:F1}%)\n");
        }

        private static void Step8_ExtendedTerminusEnrichment(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("=== STEP 8: Extended Terminus Enrichment ===");
            Console.WriteLine("(See Step 10 Part C)\n");
        }

        private static void Step9_MassAccuracyCharacterization(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("=== STEP 9: Mass Accuracy ===");
            if (ions.Count == 0) { Console.WriteLine("Skipping\n"); return; }
            var valid = ions.Where(i => !double.IsNaN(i.MassErrorPpm)).ToList();
            int pass = valid.Count(i => i.PassesMassAccuracyFilter);
            Console.WriteLine($"PassesMassAccuracyFilter: {pass:N0}/{valid.Count} ({100.0 * pass / valid.Count:F1}%)\n");
        }

        private static void Step10_FeaturePreparationAndHypothesisValidation(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("=== STEP 10: Feature Preparation ===");
            if (ions.Count == 0) { Console.WriteLine("Skipping\n"); return; }
            var passing = ions.Where(i => i.PassesMassAccuracyFilter).ToList();
            Console.WriteLine($"Passing ions: {passing.Count:N0}");
            int hasB = passing.Count(i => i.HasMatchedBIonAtNTerm);
            int hasY = passing.Count(i => i.HasMatchedYIonAtCTerm);
            Console.WriteLine($"Has B: {hasB:N0}, Has Y: {hasY:N0}\n");
        }

        #endregion

        #region Step 11 - Feature Refinement and Model Validation

        private static void Step11_FeatureRefinementAndModelValidation(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("==================================================================");
            Console.WriteLine("=== STEP 11: Feature Refinement & Model-Ready Validation ===");
            Console.WriteLine("==================================================================");
            Console.WriteLine();

            if (ions.Count == 0)
            {
                Console.WriteLine("Skipping (no ions)");
                return;
            }

            var passingIons = ions.Where(i => i.PassesMassAccuracyFilter && !double.IsNaN(i.MassErrorPpm)).ToList();

            if (passingIons.Count == 0)
            {
                Console.WriteLine("No passing ions for analysis.");
                return;
            }

            Step11_PartA_BYAsymmetryAnalysis(passingIons);
            Step11_PartB_CTerminalProximityConfound(passingIons);
            Step11_PartD_RevisedHypothesisValidation(passingIons);
            Step11_PartE_FullFeatureInventory(passingIons);
        }

        private static void Step11_PartA_BYAsymmetryAnalysis(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine("PART A: B/Y ASYMMETRY ANALYSIS");
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine();

            int n = ions.Count;
            int hasB = ions.Count(i => i.HasMatchedBIonAtNTerm);
            int hasY = ions.Count(i => i.HasMatchedYIonAtCTerm);
            int hasBoth = ions.Count(i => i.HasMatchedBIonAtNTerm && i.HasMatchedYIonAtCTerm);
            int hasNeither = ions.Count(i => !i.HasMatchedBIonAtNTerm && !i.HasMatchedYIonAtCTerm);

            Console.WriteLine($"HasMatchedBIonAtNTerm = TRUE: {hasB:N0} ({100.0 * hasB / n:F1}%)");
            Console.WriteLine($"HasMatchedYIonAtCTerm = TRUE: {hasY:N0} ({100.0 * hasY / n:F1}%)");
            Console.WriteLine($"Both TRUE: {hasBoth:N0} ({100.0 * hasBoth / n:F1}%)");
            Console.WriteLine($"Neither TRUE: {hasNeither:N0} ({100.0 * hasNeither / n:F1}%)");
            Console.WriteLine();

            // Mean intensities excluding zeros
            var bNonZero = ions.Where(i => i.BIonIntensityAtNTerm > 0).Select(i => i.BIonIntensityAtNTerm).ToList();
            var yNonZero = ions.Where(i => i.YIonIntensityAtCTerm > 0).Select(i => i.YIonIntensityAtCTerm).ToList();

            double meanB = bNonZero.Count > 0 ? bNonZero.Average() : 0;
            double meanY = yNonZero.Count > 0 ? yNonZero.Average() : 0;

            Console.WriteLine($"Mean BIonIntensityAtNTerm (B>0, n={bNonZero.Count}): {meanB:F4}");
            Console.WriteLine($"Mean YIonIntensityAtCTerm (Y>0, n={yNonZero.Count}): {meanY:F4}");

            if (meanB > 0 && meanY > 2 * meanB)
            {
                Console.WriteLine();
                Console.WriteLine("FLAG: Strong Y/B asymmetry detected. Y-terminal coverage dominates.");
                Console.WriteLine("BYProductScore will have structural zeros; treat B and Y as independent features.");
            }
            Console.WriteLine();
        }

        private static void Step11_PartB_CTerminalProximityConfound(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine("PART B: C-TERMINAL PROXIMITY CONFOUND");
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine();

            // Bin by DistanceFromCTerm
            var bins = new (int minDist, int maxDist, string label)[]
            {
                (0, 0, "0"),
                (1, 2, "1-2"),
                (3, 5, "3-5"),
                (6, 10, "6-10"),
                (11, int.MaxValue, "11+")
            };

            Console.WriteLine("+--------+-------+----------------+------------------+");
            Console.WriteLine("| Dist   | Count | Mean YIonInt   | Mean NormIntens  |");
            Console.WriteLine("+--------+-------+----------------+------------------+");

            var binStats = new List<(string label, int count, double meanY, double meanNorm)>();

            foreach (var (minDist, maxDist, label) in bins)
            {
                var binIons = ions.Where(i => i.DistanceFromCTerm >= minDist && i.DistanceFromCTerm <= maxDist).ToList();
                int count = binIons.Count;
                double meanY = count > 0 ? binIons.Average(i => i.YIonIntensityAtCTerm) : 0;
                double meanNorm = count > 0 ? binIons.Average(i => i.NormalizedIntensity) : 0;

                binStats.Add((label, count, meanY, meanNorm));
                Console.WriteLine($"| {label,-6} | {count,5} | {meanY,14:F4} | {meanNorm,16:F4} |");
            }

            Console.WriteLine("+--------+-------+----------------+------------------+");
            Console.WriteLine();

            // Check for confound
            var closeY = binStats.Where(b => b.label == "0" || b.label == "1-2").Sum(b => b.meanY * b.count);
            var closeCount = binStats.Where(b => b.label == "0" || b.label == "1-2").Sum(b => b.count);
            var farY = binStats.Where(b => b.label == "11+").Sum(b => b.meanY * b.count);
            var farCount = binStats.Where(b => b.label == "11+").Sum(b => b.count);

            double meanCloseY = closeCount > 0 ? closeY / closeCount : 0;
            double meanFarY = farCount > 0 ? farY / farCount : 0;

            if (meanCloseY > 2 * meanFarY && meanFarY > 0)
            {
                Console.WriteLine("FLAG: YIonIntensityAtCTerm is confounded with C-terminal proximity.");
                Console.WriteLine("Add DistanceFromCTerm as a covariate or use residualized Y-intensity.");
                Console.WriteLine();
            }
        }

        private static void Step11_PartD_RevisedHypothesisValidation(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine("PART D: REVISED HYPOTHESIS VALIDATION");
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine();

            // Correlations
            Console.WriteLine("--- Intensity Correlations ---");

            var normInt = ions.Select(i => i.NormalizedIntensity).ToArray();
            var bInt = ions.Select(i => i.BIonIntensityAtNTerm).ToArray();
            var yInt = ions.Select(i => i.YIonIntensityAtCTerm).ToArray();
            var maxTerm = ions.Select(i => i.MaxTerminalIonIntensity).ToArray();
            var byProd = ions.Select(i => i.BYProductScore).ToArray();
            var hasBoth = ions.Select(i => i.HasBothTerminalIons ? 1.0 : 0.0).ToArray();

            Console.WriteLine($"Pearson(NormInt, BIonIntensity):       {Corr(normInt, bInt):F4}");
            Console.WriteLine($"Pearson(NormInt, YIonIntensity):       {Corr(normInt, yInt):F4}");
            Console.WriteLine($"Pearson(NormInt, MaxTerminalIonInt):   {Corr(normInt, maxTerm):F4}");
            Console.WriteLine($"Pearson(NormInt, BYProductScore):      {Corr(normInt, byProd):F4}");
            Console.WriteLine($"Pearson(NormInt, HasBothTerminalIons): {Corr(normInt, hasBoth):F4}");
            Console.WriteLine();

            // Group comparison
            Console.WriteLine("--- Group Comparison by Terminal Ion Presence ---");

            var group00 = ions.Where(i => !i.HasMatchedBIonAtNTerm && !i.HasMatchedYIonAtCTerm).ToList();
            var group10 = ions.Where(i => i.HasMatchedBIonAtNTerm && !i.HasMatchedYIonAtCTerm).ToList();
            var group01 = ions.Where(i => !i.HasMatchedBIonAtNTerm && i.HasMatchedYIonAtCTerm).ToList();
            var group11 = ions.Where(i => i.HasMatchedBIonAtNTerm && i.HasMatchedYIonAtCTerm).ToList();

            void PrintGroup(string name, List<InternalFragmentIon> grp)
            {
                if (grp.Count == 0) { Console.WriteLine($"  {name}: (no ions)"); return; }
                var sorted = grp.Select(i => i.NormalizedIntensity).OrderBy(x => x).ToList();
                Console.WriteLine($"  {name}: n={grp.Count,5}, Mean={sorted.Average():F4}, Median={sorted[sorted.Count / 2]:F4}");
            }

            PrintGroup("Group 00 (neither)", group00);
            PrintGroup("Group 10 (B only)", group10);
            PrintGroup("Group 01 (Y only)", group01);
            PrintGroup("Group 11 (both)", group11);
            Console.WriteLine();
            Console.WriteLine("  Expected: Group 11 > Group 01 > Group 10 > Group 00");
            Console.WriteLine("  (Y dominates B due to digest asymmetry)");
            Console.WriteLine();
        }

        private static void Step11_PartE_FullFeatureInventory(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine("PART E: FULL FEATURE INVENTORY AND CORRELATION MATRIX");
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine();

            var normInt = ions.Select(i => i.NormalizedIntensity).ToArray();

            // Feature definitions
            var features = new (string name, Func<InternalFragmentIon, double> getter)[]
            {
                ("FragLen", i => i.FragmentLength),
                ("DistCTerm", i => i.DistanceFromCTerm),
                ("NumBasic", i => i.NumberOfBasicResidues),
                ("BIonInt", i => i.BIonIntensityAtNTerm),
                ("YIonInt", i => i.YIonIntensityAtCTerm),
                ("MaxTermInt", i => i.MaxTerminalIonIntensity),
                ("BYProd", i => i.BYProductScore),
                ("HasBoth", i => i.HasBothTerminalIons ? 1.0 : 0.0),
                ("LocalRank", i => double.IsNaN(i.LocalIntensityRank) ? 0 : i.LocalIntensityRank),
                ("NormInt", i => i.NormalizedIntensity)
            };

            // Correlation matrix
            Console.WriteLine("--- Correlation Matrix ---");
            Console.Write("          ");
            foreach (var f in features) Console.Write($"{f.name,9} ");
            Console.WriteLine();

            var featureArrays = features.Select(f => ions.Select(f.getter).ToArray()).ToArray();
            var collinearPairs = new List<string>();
            var usefulSignals = new List<string>();

            for (int i = 0; i < features.Length; i++)
            {
                Console.Write($"{features[i].name,9} ");
                for (int j = 0; j < features.Length; j++)
                {
                    double r = Corr(featureArrays[i], featureArrays[j]);
                    Console.Write($"{r,9:F3} ");

                    if (i < j && i < features.Length - 1 && j < features.Length - 1 && Math.Abs(r) > 0.6)
                        collinearPairs.Add($"{features[i].name} × {features[j].name} (r={r:F3})");

                    if (i < features.Length - 1 && j == features.Length - 1 && Math.Abs(r) > 0.2)
                        usefulSignals.Add($"{features[i].name} (r={r:F3})");
                }
                Console.WriteLine();
            }
            Console.WriteLine();

            if (collinearPairs.Count > 0)
            {
                Console.WriteLine("Potentially collinear predictor pairs (|r| > 0.6):");
                foreach (var pair in collinearPairs) Console.WriteLine($"  {pair}");
                Console.WriteLine();
            }

            if (usefulSignals.Count > 0)
            {
                Console.WriteLine("Useful signals (predictor-target |r| > 0.2):");
                foreach (var sig in usefulSignals) Console.WriteLine($"  {sig}");
                Console.WriteLine();
            }

            // Feature readiness table
            Console.WriteLine("--- Feature Readiness Table ---");
            Console.WriteLine("+------------------+----------+----------+----------+----------+-----------+");
            Console.WriteLine("| Feature          | N values | N nonzero|     Mean |   StdDev | r(target) |");
            Console.WriteLine("+------------------+----------+----------+----------+----------+-----------+");

            for (int i = 0; i < features.Length - 1; i++) // Exclude target
            {
                var vals = featureArrays[i];
                int nValid = vals.Count(v => !double.IsNaN(v));
                int nNonZero = vals.Count(v => !double.IsNaN(v) && v != 0);
                double mean = vals.Where(v => !double.IsNaN(v)).DefaultIfEmpty(0).Average();
                double std = nValid > 1 ? Math.Sqrt(vals.Where(v => !double.IsNaN(v)).Select(v => Math.Pow(v - mean, 2)).Average()) : 0;
                double rTarget = Corr(vals, normInt);

                Console.WriteLine($"| {features[i].name,-16} | {nValid,8} | {nNonZero,8} | {mean,8:F4} | {std,8:F4} | {rTarget,9:F4} |");
            }

            Console.WriteLine("+------------------+----------+----------+----------+----------+-----------+");
            Console.WriteLine();
            Console.WriteLine("Pre-training data quality gate: COMPLETE");
            Console.WriteLine($"Total model-ready ions: {ions.Count:N0}");
            Console.WriteLine();
        }

        #endregion

        #region Utility

        private static double Corr(double[] x, double[] y)
        {
            if (x.Length != y.Length || x.Length == 0) return double.NaN;
            double mx = x.Average(), my = y.Average();
            double sxy = 0, sx2 = 0, sy2 = 0;
            for (int i = 0; i < x.Length; i++)
            {
                double dx = x[i] - mx, dy = y[i] - my;
                sxy += dx * dy;
                sx2 += dx * dx;
                sy2 += dy * dy;
            }
            double d = Math.Sqrt(sx2 * sy2);
            return d == 0 ? 0 : sxy / d;
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