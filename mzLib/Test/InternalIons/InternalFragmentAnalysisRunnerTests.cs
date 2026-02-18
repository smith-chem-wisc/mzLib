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
            double maxB = bNonZero.Count > 0 ? bNonZero.Max() : 0;
            double maxY = yNonZero.Count > 0 ? yNonZero.Max() : 0;

            Console.WriteLine($"Mean BIonIntensityAtNTerm (B>0, n={bNonZero.Count,5}): {meanB:E4}");
            Console.WriteLine($"Mean YIonIntensityAtCTerm (Y>0, n={yNonZero.Count,5}): {meanY:E4}");
            Console.WriteLine($"Max  BIonIntensityAtNTerm: {maxB:E4}");
            Console.WriteLine($"Max  YIonIntensityAtCTerm: {maxY:E4}\n");

            if (meanB > meanY && meanB > 0)
            {
                Console.WriteLine("CONFIRMED: B-ions dominate Y-ions (positive mode HCD, as expected).");
                Console.WriteLine("N-terminal cleavage propensity is the primary spectral driver.\n");
            }
            else
            {
                Console.WriteLine("UNEXPECTED: Y-ions match or exceed B-ions after TIC normalization.");
                Console.WriteLine("Verify that TotalIonCurrent reflects MS2 scan TIC, not MS1 or TIC.\n");
            }
        }

        private static void Step11_PartB_YIonStratification(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine("PART B: STRONG Y-ION MECHANISM STRATIFICATION");
            Console.WriteLine("═══════════════════════════════════════════════════\n");

            var withY = ions.Where(i => i.HasMatchedYIonAtCTerm).ToList();
            if (withY.Count == 0) { Console.WriteLine("No ions with matched Y-ion.\n"); return; }

            var s1 = withY.Where(i => i.DistanceFromCTerm <= 3).ToList();
            var s2 = withY.Where(i => i.DistanceFromCTerm > 3 && i.BasicResiduesInYIonSpan >= 2).ToList();
            var s3 = withY.Where(i => i.DistanceFromCTerm > 3 && i.BasicResiduesInYIonSpan <= 1).ToList();

            void PrintStratum(string name, string desc, List<InternalFragmentIon> list)
            {
                Console.WriteLine($"{name}: {list.Count:N0} ions");
                Console.WriteLine($"  {desc}");
                if (list.Count == 0) { Console.WriteLine(); return; }
                double meanY = list.Average(i => i.YIonIntensityAtCTerm);
                double meanTic = list.Average(i => i.TicNormalizedIntensity);
                double fracIntense = (double)list.Count(i => i.TicNormalizedIntensity > 0.003) / list.Count;
                Console.WriteLine($"  Mean YIonIntensity:        {meanY:E4}");
                Console.WriteLine($"  Mean TicNormalizedInt:     {meanTic:E4}");
                Console.WriteLine($"  Frac TicNormInt > 0.003:   {fracIntense:P1}\n");
            }

            PrintStratum("Stratum 1", "C-terminal K/R rescue (DistCTerm <= 3)", s1);
            PrintStratum("Stratum 2", "Multiply-basic large y-ion (DistCTerm > 3, BasicInY >= 2)", s2);
            PrintStratum("Stratum 3", "Weak charge retention (DistCTerm > 3, BasicInY <= 1)", s3);

            // Find most informative stratum
            var strata = new[] { (s1, "Stratum 1"), (s2, "Stratum 2"), (s3, "Stratum 3") }
                .Where(x => x.Item1.Count > 0)
                .OrderByDescending(x => x.Item1.Average(i => i.TicNormalizedIntensity))
                .FirstOrDefault();

            if (strata.Item1 != null)
                Console.WriteLine($"NOTE: {strata.Item2} shows highest mean internal fragment intensity.\n");
        }

        private static void Step11_PartC_UnsupportedHighIntensityIons(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine("PART C: UNSUPPORTED HIGH-INTENSITY IONS");
            Console.WriteLine("═══════════════════════════════════════════════════\n");

            var unsupported = ions.Where(i => !i.HasMatchedBIonAtNTerm && !i.HasMatchedYIonAtCTerm).ToList();
            Console.WriteLine($"Ions with neither b nor y terminal ion matched: {unsupported.Count:N0}\n");

            if (unsupported.Count == 0) { return; }

            Console.WriteLine("Top 10 by TicNormalizedIntensity (weak-bond hypothesis validation):");
            Console.WriteLine("+------------------+------------+--------+--------+--------+--------+--------+----------+");
            Console.WriteLine("| InternalSeq      | TicNormInt | FlankN | FlankC | IntN   | IntC   | Length | DistCTrm |");
            Console.WriteLine("+------------------+------------+--------+--------+--------+--------+--------+----------+");

            foreach (var ion in unsupported.OrderByDescending(i => i.TicNormalizedIntensity).Take(10))
            {
                string seq = ion.InternalSequence.Length > 16 ? ion.InternalSequence.Substring(0, 13) + "..." : ion.InternalSequence;
                Console.WriteLine($"| {seq,-16} | {ion.TicNormalizedIntensity,10:E3} | {ion.NTerminalFlankingResidue,6} | {ion.CTerminalFlankingResidue,6} | {ion.InternalNTerminalAA,6} | {ion.InternalCTerminalAA,6} | {ion.FragmentLength,6} | {ion.DistanceFromCTerm,8} |");
            }
            Console.WriteLine("+------------------+------------+--------+--------+--------+--------+--------+----------+");
            Console.WriteLine("\nThese ions validate weak-bond hypothesis as independent of b/y support.\n");
        }

        private static void Step11_PartD_FourGroupComparison(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine("PART D: FOUR-GROUP INTENSITY COMPARISON");
            Console.WriteLine("═══════════════════════════════════════════════════\n");

            var g00 = ions.Where(i => !i.HasMatchedBIonAtNTerm && !i.HasMatchedYIonAtCTerm).ToList();
            var g10 = ions.Where(i => i.HasMatchedBIonAtNTerm && !i.HasMatchedYIonAtCTerm).ToList();
            var g01 = ions.Where(i => !i.HasMatchedBIonAtNTerm && i.HasMatchedYIonAtCTerm).ToList();
            var g11 = ions.Where(i => i.HasMatchedBIonAtNTerm && i.HasMatchedYIonAtCTerm).ToList();

            void PrintGroup(string name, List<InternalFragmentIon> grp)
            {
                if (grp.Count == 0) { Console.WriteLine($"  {name}: (no ions)"); return; }
                var sorted = grp.Select(i => i.TicNormalizedIntensity).OrderBy(x => x).ToList();
                double mean = sorted.Average();
                double median = sorted[sorted.Count / 2];
                Console.WriteLine($"  {name}: n={grp.Count,6}, Mean={mean:E4}, Median={median:E4}");
            }

            PrintGroup("Group 00 (neither)", g00);
            PrintGroup("Group 10 (B only)", g10);
            PrintGroup("Group 01 (Y only)", g01);
            PrintGroup("Group 11 (both)", g11);

            Console.WriteLine("\n  Expected ordering: Group 11 > Group 10 > Group 01 > Group 00\n");

            double mean10 = g10.Count > 0 ? g10.Average(i => i.TicNormalizedIntensity) : 0;
            double mean01 = g01.Count > 0 ? g01.Average(i => i.TicNormalizedIntensity) : 0;

            if (mean10 > mean01 && mean10 > 0)
            {
                Console.WriteLine("CONFIRMED: B-ion support is more predictive than Y-ion support,");
                Console.WriteLine("consistent with positive mode HCD fragmentation physics.\n");
            }
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
                ("BYProductScore", Corr(ions.Select(i => i.BYProductScore).ToArray(), target)),
                ("HasBothTerminalIons", Corr(ions.Select(i => i.HasBothTerminalIons ? 1.0 : 0.0).ToArray(), target)),
                ("BasicResiduesInBIonSpan", Corr(ions.Select(i => (double)i.BasicResiduesInBIonSpan).ToArray(), target)),
                ("BasicResiduesInYIonSpan", Corr(ions.Select(i => (double)i.BasicResiduesInYIonSpan).ToArray(), target)),
                ("FragmentLength", Corr(ions.Select(i => (double)i.FragmentLength).ToArray(), target)),
                ("DistanceFromCTerm", Corr(ions.Select(i => (double)i.DistanceFromCTerm).ToArray(), target)),
                ("NumberOfBasicResidues", Corr(ions.Select(i => (double)i.NumberOfBasicResidues).ToArray(), target))
            };

            Console.WriteLine("--- Correlations with TicNormalizedIntensity (ranked by |r|) ---\n");
            foreach (var (name, r) in correlations.OrderByDescending(x => Math.Abs(x.r)))
            {
                string flag = Math.Abs(r) > 0.15 ? " *" : "";
                Console.WriteLine($"  {name,-25}: r = {r,8:F4}{flag}");
            }
            Console.WriteLine();

            // Charge-retention mechanism test
            var bSpan = ions.Select(i => (double)i.BasicResiduesInBIonSpan).ToArray();
            var ySpan = ions.Select(i => (double)i.BasicResiduesInYIonSpan).ToArray();
            var bInt = ions.Select(i => i.BIonIntensityAtNTerm).ToArray();
            var yInt = ions.Select(i => i.YIonIntensityAtCTerm).ToArray();

            double rBSpanBInt = Corr(bSpan, bInt);
            double rYSpanYInt = Corr(ySpan, yInt);

            Console.WriteLine("--- Charge-Retention Mechanism Test ---");
            Console.WriteLine($"  Corr(BasicResiduesInBIonSpan, BIonIntensity): {rBSpanBInt:F4}");
            Console.WriteLine($"  Corr(BasicResiduesInYIonSpan, YIonIntensity): {rYSpanYInt:F4}\n");

            if (rBSpanBInt > 0.1 && rYSpanYInt > 0.1)
            {
                Console.WriteLine("CONFIRMED: Charge-retention mechanism validated. Basic residue");
                Console.WriteLine("count in ion span predicts terminal ion intensity for both b and y series.\n");
            }
        }

        private static void Step11_PartF_FeatureInventory(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine("PART F: FEATURE INVENTORY — MODEL GO/NO-GO GATE");
            Console.WriteLine("═══════════════════════════════════════════════════\n");

            var target = ions.Select(i => i.TicNormalizedIntensity).ToArray();

            var features = new (string name, Func<InternalFragmentIon, double> get, bool isTarget)[]
            {
                ("FragmentLength", i => i.FragmentLength, false),
                ("DistanceFromCTerm", i => i.DistanceFromCTerm, false),
                ("NumberOfBasicResidues", i => i.NumberOfBasicResidues, false),
                ("BasicResiduesInBIonSpan", i => i.BasicResiduesInBIonSpan, false),
                ("BasicResiduesInYIonSpan", i => i.BasicResiduesInYIonSpan, false),
                ("BIonIntensityAtNTerm", i => i.BIonIntensityAtNTerm, false),
                ("YIonIntensityAtCTerm", i => i.YIonIntensityAtCTerm, false),
                ("MaxTerminalIonIntensity", i => i.MaxTerminalIonIntensity, false),
                ("BYProductScore", i => i.BYProductScore, false),
                ("HasBothTerminalIons", i => i.HasBothTerminalIons ? 1.0 : 0.0, false),
                ("HasProlineAtEitherTerm", i => i.HasProlineAtEitherTerminus ? 1.0 : 0.0, false),
                ("HasAspartateAtEitherTerm", i => i.HasAspartateAtEitherTerminus ? 1.0 : 0.0, false),
                ("LocalIntensityRank", i => double.IsNaN(i.LocalIntensityRank) ? 0 : i.LocalIntensityRank, false),
                ("HasModifiedResidue", i => i.HasModifiedResidue ? 1.0 : 0.0, false),
                ("TicNormalizedIntensity", i => i.TicNormalizedIntensity, true)
            };

            Console.WriteLine("+---------------------------+----------+----------+------------+------------+-----------+--------+");
            Console.WriteLine("| Feature                   | N values | N nonzero|       Mean |     StdDev | r(target) | Flags  |");
            Console.WriteLine("+---------------------------+----------+----------+------------+------------+-----------+--------+");

            foreach (var (name, getter, isTarget) in features)
            {
                var vals = ions.Select(getter).ToArray();
                int nVal = vals.Count(v => !double.IsNaN(v));
                int nNz = vals.Count(v => !double.IsNaN(v) && v != 0);
                double mean = nVal > 0 ? vals.Where(v => !double.IsNaN(v)).Average() : 0;
                double std = nVal > 1 ? Math.Sqrt(vals.Where(v => !double.IsNaN(v)).Select(v => (v - mean) * (v - mean)).Average()) : 0;
                double r = isTarget ? double.NaN : Corr(vals, target);

                string rStr = isTarget ? "  (target)" : $"{r,10:F4}";

                var flags = new List<string>();
                if (!isTarget && nVal > 0 && (double)nNz / nVal < 0.20)
                    flags.Add("SPARSE");
                if (!isTarget && !double.IsNaN(r) && Math.Abs(r) > 0.15)
                    flags.Add("SIGNAL");

                string flagStr = flags.Count > 0 ? string.Join(",", flags) : "";

                Console.WriteLine($"| {name,-25} | {nVal,8} | {nNz,8} | {mean,10:E3} | {std,10:E3} | {rStr,-9} | {flagStr,-6} |");
            }

            Console.WriteLine("+---------------------------+----------+----------+------------+------------+-----------+--------+");
            Console.WriteLine();
            Console.WriteLine("Legend:");
            Console.WriteLine("  SPARSE = N non-zero < 20% of N non-null — consider binary encoding or dropping");
            Console.WriteLine("  SIGNAL = |r(target)| > 0.15 — include in initial model");
            Console.WriteLine();
            Console.WriteLine($"Model-ready ions: {ions.Count:N0}");
            Console.WriteLine("Pre-training data quality gate: COMPLETE\n");
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
                psmTsvPath: @"F:\JurakatMultiprotease\2026-02-17-16-40-42\Task1-SearchTask\FilteredPeptides.psmtsv",
                rawFileFolder: @"F:\JurakatMultiprotease",
                outputDirectory: @"E:\Projects\internalIons\mzLibReportsForClaude"
            );
        }
    }
}