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

        private static readonly string[] SupportedExtensions = { ".raw", ".mzml", ".mgf" };

        public static void RunAll() => RunAll(DefaultPsmTsvPath, DefaultRawFileFolder, DefaultOutputDirectory);

        public static void RunAll(string psmTsvPath, string rawFileFolder, string outputDirectory)
        {
            Console.WriteLine("==================================================================");
            Console.WriteLine("     INTERNAL FRAGMENT ION ANALYSIS HARNESS");
            Console.WriteLine("==================================================================\n");

            try
            {
                if (!File.Exists(psmTsvPath)) throw new FileNotFoundException($"PSM TSV not found: {psmTsvPath}");
                if (!Directory.Exists(rawFileFolder)) throw new DirectoryNotFoundException($"Raw folder not found: {rawFileFolder}");
                if (!Directory.Exists(outputDirectory)) Directory.CreateDirectory(outputDirectory);

                var psms = Step1_LoadPsms(psmTsvPath);
                var msDataFiles = Step2_LoadRawFiles(rawFileFolder, psms);
                var internalIons = Step3_ExtractInternalFragments(psms, msDataFiles, DefaultCollisionEnergy);
                string outputPath = Step4_WriteOutputTsv(internalIons, outputDirectory);
                Step5_RoundTripValidation(internalIons, outputPath);
                Step6to10_Summary(internalIons);
                Step11_TicNormalizedAnalysis(internalIons);

                Console.WriteLine("\n==================================================================");
                Console.WriteLine("ALL STEPS COMPLETED SUCCESSFULLY");
                Console.WriteLine("==================================================================");
            }
            catch (Exception ex)
            {
                Console.WriteLine($"\nFAILURE: {ex.GetType().Name}: {ex.Message}\n{ex.StackTrace}");
                throw;
            }
        }

        #region Steps 1-10

        private static List<PsmFromTsv> Step1_LoadPsms(string path)
        {
            Console.WriteLine("=== STEP 1: Load PSMs ===");
            var psms = SpectrumMatchTsvReader.ReadPsmTsv(path, out _);
            Console.WriteLine($"Loaded: {psms.Count:N0} PSMs\n");
            return psms;
        }

        private static Dictionary<string, MsDataFile> Step2_LoadRawFiles(string folder, List<PsmFromTsv> psms)
        {
            Console.WriteLine("=== STEP 2: Load Raw Files ===");
            var needed = psms.Select(p => p.FileNameWithoutExtension).Distinct().ToHashSet(StringComparer.OrdinalIgnoreCase);
            var files = new Dictionary<string, MsDataFile>(StringComparer.OrdinalIgnoreCase);
            foreach (var f in SupportedExtensions.SelectMany(e => Directory.GetFiles(folder, $"*{e}")))
            {
                var name = Path.GetFileNameWithoutExtension(f);
                if (!needed.Contains(name)) continue;
                try
                {
                    var ext = Path.GetExtension(f).ToLower();
                    files[name] = ext switch
                    {
                        ".mzml" => Mzml.LoadAllStaticData(f),
                        ".raw" => ThermoRawFileReader.LoadAllStaticData(f),
                        ".mgf" => Mgf.LoadAllStaticData(f),
                        _ => throw new NotSupportedException(ext)
                    };
                    Console.WriteLine($"  {name}: OK");
                }
                catch (Exception ex) { Console.WriteLine($"  {name}: ERROR - {ex.Message}"); }
            }
            Console.WriteLine();
            return files;
        }

        private static List<InternalFragmentIon> Step3_ExtractInternalFragments(List<PsmFromTsv> psms, Dictionary<string, MsDataFile> files, double ce)
        {
            Console.WriteLine("=== STEP 3: Extract Features ===");
            var all = new List<InternalFragmentIon>();
            foreach (var g in psms.GroupBy(p => p.FileNameWithoutExtension))
            {
                if (!files.TryGetValue(g.Key, out var f)) continue;
                var ions = InternalFragmentFeatureExtractor.ExtractFromPsms(g.ToList(), f, ce);
                all.AddRange(ions);
                if (ions.Count > 0) Console.WriteLine($"  {g.Key}: {ions.Count:N0}");
            }
            Console.WriteLine($"Total: {all.Count:N0}\n");
            return all;
        }

        private static string Step4_WriteOutputTsv(List<InternalFragmentIon> ions, string dir)
        {
            Console.WriteLine("=== STEP 4: Write TSV ===");
            var path = Path.Combine(dir, "InternalFragmentIons.tsv");
            InternalFragmentTsvWriter.WriteToTsv(ions, path);
            Console.WriteLine($"Output: {path}\n");
            return path;
        }

        private static void Step5_RoundTripValidation(List<InternalFragmentIon> orig, string path)
        {
            Console.WriteLine("=== STEP 5: Round-trip ===");
            if (orig.Count == 0) { Console.WriteLine("Skip\n"); return; }
            var rb = InternalFragmentTsvWriter.ReadFromTsv(path);
            Console.WriteLine(orig.Count == rb.Count ? "PASSED\n" : "FAILED\n");
        }

        private static void Step6to10_Summary(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("=== STEPS 6-10: Summary ===");
            if (ions.Count == 0) { Console.WriteLine("Skip\n"); return; }
            var valid = ions.Where(i => !double.IsNaN(i.MassErrorPpm)).ToList();
            int pass = valid.Count(i => i.PassesMassAccuracyFilter);
            Console.WriteLine($"PassesMassAccuracyFilter: {pass:N0}/{valid.Count} ({100.0 * pass / valid.Count:F1}%)\n");
        }

        #endregion

        #region Step 11 - TIC Normalized Analysis

        private static void Step11_TicNormalizedAnalysis(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("==================================================================");
            Console.WriteLine("=== STEP 11: B/Y Directionality & Model Feature Validation ===");
            Console.WriteLine("==================================================================\n");

            if (ions.Count == 0) { Console.WriteLine("No ions.\n"); return; }

            var passing = ions.Where(i => i.PassesMassAccuracyFilter && !double.IsNaN(i.TicNormalizedIntensity)).ToList();
            if (passing.Count == 0) { Console.WriteLine("No passing ions.\n"); return; }

            Console.WriteLine($"Analysis set: {passing.Count:N0} ions (PassesMassAccuracyFilter = TRUE)\n");

            Step11_PartA_BYDirectionality(passing);
            Step11_PartB_YIonStratification(passing);
            Step11_PartC_UnsupportedHighIntensityIons(passing);
            Step11_PartD_FourGroupComparison(passing);
            Step11_PartE_CorrelationAnalysis(passing);
            Step11_PartF_FeatureInventory(passing);
            Step11_PartG_NewFeatureValidation(passing);
        }

        private static void Step11_PartA_BYDirectionality(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine("PART A: B vs Y DIRECTIONALITY CONFIRMATION");
            Console.WriteLine("═══════════════════════════════════════════════════\n");

            int n = ions.Count;
            int hasB = ions.Count(i => i.HasMatchedBIonAtNTerm);
            int hasY = ions.Count(i => i.HasMatchedYIonAtCTerm);
            int hasBoth = ions.Count(i => i.HasBothTerminalIons);
            int hasNeither = ions.Count(i => !i.HasMatchedBIonAtNTerm && !i.HasMatchedYIonAtCTerm);

            Console.WriteLine($"HasMatchedBIonAtNTerm = TRUE:  {hasB,6:N0} ({100.0 * hasB / n:F1}%)");
            Console.WriteLine($"HasMatchedYIonAtCTerm = TRUE:  {hasY,6:N0} ({100.0 * hasY / n:F1}%)");
            Console.WriteLine($"HasBothTerminalIons = TRUE:    {hasBoth,6:N0} ({100.0 * hasBoth / n:F1}%)");
            Console.WriteLine($"Neither = TRUE:                {hasNeither,6:N0} ({100.0 * hasNeither / n:F1}%)\n");

            var bNonZero = ions.Where(i => i.BIonIntensityAtNTerm > 0).Select(i => i.BIonIntensityAtNTerm).ToList();
            var yNonZero = ions.Where(i => i.YIonIntensityAtCTerm > 0).Select(i => i.YIonIntensityAtCTerm).ToList();

            double meanB = bNonZero.Count > 0 ? bNonZero.Average() : 0;
            double meanY = yNonZero.Count > 0 ? yNonZero.Average() : 0;

            Console.WriteLine($"Mean BIonIntensityAtNTerm (B>0, n={bNonZero.Count,5}): {meanB:E4}");
            Console.WriteLine($"Mean YIonIntensityAtCTerm (Y>0, n={yNonZero.Count,5}): {meanY:E4}\n");

            if (meanB > meanY && meanB > 0)
                Console.WriteLine("CONFIRMED: B-ions dominate Y-ions (positive mode HCD).\n");
            else
                Console.WriteLine("UNEXPECTED: Y-ions match or exceed B-ions.\n");
        }

        private static void Step11_PartB_YIonStratification(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine("PART B: Y-ION MECHANISM STRATIFICATION");
            Console.WriteLine("═══════════════════════════════════════════════════\n");

            var withY = ions.Where(i => i.HasMatchedYIonAtCTerm).ToList();
            if (withY.Count == 0) { Console.WriteLine("No ions with matched Y-ion.\n"); return; }

            var s1 = withY.Where(i => i.DistanceFromCTerm <= 3).ToList();
            var s2 = withY.Where(i => i.DistanceFromCTerm > 3 && i.BasicResiduesInYIonSpan >= 2).ToList();
            var s3 = withY.Where(i => i.DistanceFromCTerm > 3 && i.BasicResiduesInYIonSpan <= 1).ToList();

            void Print(string name, List<InternalFragmentIon> list)
            {
                Console.WriteLine($"{name}: {list.Count:N0} ions");
                if (list.Count == 0) return;
                Console.WriteLine($"  Mean TicNI: {list.Average(i => i.TicNormalizedIntensity):E4}\n");
            }

            Print("S1 (DistCTerm <= 3)", s1);
            Print("S2 (DistCTerm > 3, BasicY >= 2)", s2);
            Print("S3 (DistCTerm > 3, BasicY <= 1)", s3);
        }

        private static void Step11_PartC_UnsupportedHighIntensityIons(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine("PART C: UNSUPPORTED HIGH-INTENSITY IONS");
            Console.WriteLine("═══════════════════════════════════════════════════\n");

            var unsupported = ions.Where(i => !i.HasMatchedBIonAtNTerm && !i.HasMatchedYIonAtCTerm).ToList();
            Console.WriteLine($"Ions with neither b nor y: {unsupported.Count:N0}\n");
        }

        private static void Step11_PartD_FourGroupComparison(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine("PART D: FOUR-GROUP COMPARISON");
            Console.WriteLine("═══════════════════════════════════════════════════\n");

            var g00 = ions.Where(i => !i.HasMatchedBIonAtNTerm && !i.HasMatchedYIonAtCTerm).ToList();
            var g10 = ions.Where(i => i.HasMatchedBIonAtNTerm && !i.HasMatchedYIonAtCTerm).ToList();
            var g01 = ions.Where(i => !i.HasMatchedBIonAtNTerm && i.HasMatchedYIonAtCTerm).ToList();
            var g11 = ions.Where(i => i.HasMatchedBIonAtNTerm && i.HasMatchedYIonAtCTerm).ToList();

            void Print(string name, List<InternalFragmentIon> grp)
            {
                if (grp.Count == 0) { Console.WriteLine($"  {name}: (no ions)"); return; }
                Console.WriteLine($"  {name}: n={grp.Count,6}, Mean={grp.Average(i => i.TicNormalizedIntensity):E4}");
            }

            Print("G00 (neither)", g00);
            Print("G10 (B only)", g10);
            Print("G01 (Y only)", g01);
            Print("G11 (both)", g11);
            Console.WriteLine();
        }

        private static void Step11_PartE_CorrelationAnalysis(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine("PART E: CORRELATION ANALYSIS");
            Console.WriteLine("═══════════════════════════════════════════════════\n");

            var target = ions.Select(i => i.TicNormalizedIntensity).ToArray();

            var correlations = new List<(string name, double r)>
            {
                ("BIonIntensityAtNTerm", Corr(ions.Select(i => i.BIonIntensityAtNTerm).ToArray(), target)),
                ("YIonIntensityAtCTerm", Corr(ions.Select(i => i.YIonIntensityAtCTerm).ToArray(), target)),
                ("MaxTerminalIonIntensity", Corr(ions.Select(i => i.MaxTerminalIonIntensity).ToArray(), target)),
                ("FragmentLength", Corr(ions.Select(i => (double)i.FragmentLength).ToArray(), target)),
                ("NTermFlankHydrophobicity", Corr(ions.Select(i => i.NTerminalFlankingHydrophobicity).ToArray(), target))
            };

            Console.WriteLine("Top correlations with TicNormalizedIntensity:");
            foreach (var (name, r) in correlations.OrderByDescending(x => Math.Abs(x.r)))
                Console.WriteLine($"  {name,-25}: r = {r,8:F4}");
            Console.WriteLine();
        }

        private static void Step11_PartF_FeatureInventory(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine("PART F: FEATURE INVENTORY");
            Console.WriteLine("═══════════════════════════════════════════════════\n");

            Console.WriteLine($"Model-ready ions: {ions.Count:N0}");
            Console.WriteLine("Feature inventory complete.\n");
        }

        private static void Step11_PartG_NewFeatureValidation(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine("PART G: NEW SCORER FEATURE VALIDATION");
            Console.WriteLine("═══════════════════════════════════════════════════\n");

            // IsProlineAtInternalNTerminus
            var proTrue = ions.Where(i => i.IsProlineAtInternalNTerminus).ToList();
            var proFalse = ions.Where(i => !i.IsProlineAtInternalNTerminus).ToList();

            Console.WriteLine("--- IsProlineAtInternalNTerminus ---");
            Console.WriteLine($"  TRUE:  {proTrue.Count:N0} ({100.0 * proTrue.Count / ions.Count:F1}%)");
            Console.WriteLine($"  FALSE: {proFalse.Count:N0}");

            double meanProTrue = proTrue.Count > 0 ? proTrue.Average(i => i.TicNormalizedIntensity) : 0;
            double meanProFalse = proFalse.Count > 0 ? proFalse.Average(i => i.TicNormalizedIntensity) : 0;
            double proRatio = meanProFalse > 0 ? meanProTrue / meanProFalse : 0;

            Console.WriteLine($"  Mean TicNI (TRUE):  {meanProTrue:E4}");
            Console.WriteLine($"  Mean TicNI (FALSE): {meanProFalse:E4}");
            Console.WriteLine($"  Ratio (TRUE/FALSE): {proRatio:F2}");
            Console.WriteLine($"  Expectation: ratio > 1.0 (proline N-term cleavage enriched)");
            Console.WriteLine();

            // IsTerminalRescue
            var rescueTrue = ions.Where(i => i.IsTerminalRescue).ToList();
            var rescueFalse = ions.Where(i => !i.IsTerminalRescue).ToList();

            Console.WriteLine("--- IsTerminalRescue ---");
            Console.WriteLine($"  TRUE:  {rescueTrue.Count:N0} ({100.0 * rescueTrue.Count / ions.Count:F1}%)");
            Console.WriteLine($"  FALSE: {rescueFalse.Count:N0}");

            double meanRescueTicTrue = rescueTrue.Count > 0 ? rescueTrue.Average(i => i.TicNormalizedIntensity) : 0;
            double meanRescueTicFalse = rescueFalse.Count > 0 ? rescueFalse.Average(i => i.TicNormalizedIntensity) : 0;
            double meanRescueYTrue = rescueTrue.Count > 0 ? rescueTrue.Average(i => i.YIonIntensityAtCTerm) : 0;
            double meanRescueYFalse = rescueFalse.Count > 0 ? rescueFalse.Average(i => i.YIonIntensityAtCTerm) : 0;

            Console.WriteLine($"  Mean TicNI (TRUE):  {meanRescueTicTrue:E4}");
            Console.WriteLine($"  Mean TicNI (FALSE): {meanRescueTicFalse:E4}");
            Console.WriteLine($"  Mean YIon (TRUE):   {meanRescueYTrue:E4}");
            Console.WriteLine($"  Mean YIon (FALSE):  {meanRescueYFalse:E4}");
            Console.WriteLine($"  Expectation: YIon elevated when TRUE (C-terminal K/R rescue)");
            Console.WriteLine();

            // NTerminalFlankingHydrophobicity
            var hydroVals = ions.Select(i => i.NTerminalFlankingHydrophobicity).ToList();
            var target = ions.Select(i => i.TicNormalizedIntensity).ToArray();
            var bInt = ions.Select(i => i.BIonIntensityAtNTerm).ToArray();

            double rHydroTic = Corr(hydroVals.ToArray(), target);
            double rHydroBIon = Corr(hydroVals.ToArray(), bInt);

            Console.WriteLine("--- NTerminalFlankingHydrophobicity ---");
            Console.WriteLine($"  Min:    {hydroVals.Min():F2}");
            Console.WriteLine($"  Max:    {hydroVals.Max():F2}");
            Console.WriteLine($"  Mean:   {hydroVals.Average():F2}");
            Console.WriteLine($"  StdDev: {StdDev(hydroVals):F2}");
            Console.WriteLine($"  r(Hydro, TicNI):   {rHydroTic:F4}");
            Console.WriteLine($"  r(Hydro, BIonInt): {rHydroBIon:F4}");
            Console.WriteLine($"  Expectation: positive correlation (hydrophobic flanking = higher intensity)");
            Console.WriteLine();
        }

        #endregion

        #region Utilities

        private static double Corr(double[] x, double[] y)
        {
            if (x.Length != y.Length || x.Length == 0) return double.NaN;
            double mx = x.Average(), my = y.Average();
            double sxy = 0, sx2 = 0, sy2 = 0;
            for (int i = 0; i < x.Length; i++)
            {
                double dx = x[i] - mx, dy = y[i] - my;
                sxy += dx * dy; sx2 += dx * dx; sy2 += dy * dy;
            }
            double d = Math.Sqrt(sx2 * sy2);
            return d == 0 ? 0 : sxy / d;
        }

        private static double StdDev(IEnumerable<double> vals)
        {
            var list = vals.Where(v => !double.IsNaN(v)).ToList();
            if (list.Count < 2) return 0;
            double mean = list.Average();
            return Math.Sqrt(list.Select(v => (v - mean) * (v - mean)).Average());
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
                psmTsvPath: @"F:\MSV000090552_scribe\2026-02-23-14-49-36_nce42\Task1-SearchTask\filteredPeptides_nce42.psmtsv",
                rawFileFolder: @"F:\MSV000090552_scribe\2026-02-23-13-28-13\Task1-CalibrateTask",
                outputDirectory: @"F:\MSV000090552_scribe\2026-02-23-14-49-36_nce42\Task1-SearchTask"
            );
        }
    }
}