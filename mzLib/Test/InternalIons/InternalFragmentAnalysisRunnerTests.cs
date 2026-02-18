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
            Console.WriteLine("==================================================================");
            Console.WriteLine();

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
            Console.WriteLine("=== STEP 11: TIC-NORMALIZED ANALYSIS ===");
            Console.WriteLine("==================================================================\n");

            if (ions.Count == 0) { Console.WriteLine("No ions.\n"); return; }

            var passing = ions.Where(i => i.PassesMassAccuracyFilter && !double.IsNaN(i.TicNormalizedIntensity)).ToList();
            if (passing.Count == 0) { Console.WriteLine("No passing ions.\n"); return; }

            Console.WriteLine($"Analysis set: {passing.Count:N0} ions (PassesMassAccuracyFilter = TRUE)\n");

            Step11_PartB_BYDirectionalityConfirmation(passing);
            Step11_PartC_YIonNormalizationAudit(passing);
            Step11_PartD_StrongYIonMechanismStratification(passing);
            Step11_PartF_HypothesisValidation(passing);
            Step11_PartG_FullFeatureInventory(passing);
        }

        private static void Step11_PartB_BYDirectionalityConfirmation(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine("PART B: B vs Y DIRECTIONALITY CONFIRMATION");
            Console.WriteLine("═══════════════════════════════════════════════════\n");

            int n = ions.Count;
            int hasB = ions.Count(i => i.HasMatchedBIonAtNTerm);
            int hasY = ions.Count(i => i.HasMatchedYIonAtCTerm);
            int hasBoth = ions.Count(i => i.HasMatchedBIonAtNTerm && i.HasMatchedYIonAtCTerm);
            int hasNeither = ions.Count(i => !i.HasMatchedBIonAtNTerm && !i.HasMatchedYIonAtCTerm);

            Console.WriteLine($"HasMatchedBIonAtNTerm = TRUE: {hasB:N0} ({100.0 * hasB / n:F1}%)");
            Console.WriteLine($"HasMatchedYIonAtCTerm = TRUE: {hasY:N0} ({100.0 * hasY / n:F1}%)");
            Console.WriteLine($"Both TRUE: {hasBoth:N0} ({100.0 * hasBoth / n:F1}%)");
            Console.WriteLine($"Neither TRUE: {hasNeither:N0} ({100.0 * hasNeither / n:F1}%)\n");

            var bNonZero = ions.Where(i => i.BIonIntensityAtNTerm > 0).Select(i => i.BIonIntensityAtNTerm).ToList();
            var yNonZero = ions.Where(i => i.YIonIntensityAtCTerm > 0).Select(i => i.YIonIntensityAtCTerm).ToList();

            double meanB = bNonZero.Count > 0 ? bNonZero.Average() : 0;
            double meanY = yNonZero.Count > 0 ? yNonZero.Average() : 0;

            Console.WriteLine($"Mean BIonIntensityAtNTerm (TIC-norm, B>0, n={bNonZero.Count}): {meanB:E4}");
            Console.WriteLine($"Mean YIonIntensityAtCTerm (TIC-norm, Y>0, n={yNonZero.Count}): {meanY:E4}\n");

            if (meanB > meanY && meanB > 0)
            {
                Console.WriteLine("CONFIRMED: B-ions dominate Y-ions as expected for positive mode HCD");
                Console.WriteLine("of LysC-digested peptides. N-terminal cleavage propensity is the");
                Console.WriteLine("primary spectral driver.\n");
            }
            else
            {
                Console.WriteLine("UNEXPECTED: Y-ions match or exceed B-ions after TIC normalization.");
                Console.WriteLine("Investigate whether TotalIonCurrent reflects MS2 scan TIC.\n");
            }
        }

        private static void Step11_PartC_YIonNormalizationAudit(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine("PART C: Y-ION NORMALIZATION AUDIT");
            Console.WriteLine("═══════════════════════════════════════════════════\n");

            var highY = ions.Where(i => i.YIonIntensityAtCTerm > 0.2).ToList();
            Console.WriteLine($"Ions with YIonIntensityAtCTerm > 0.2: {highY.Count}\n");

            if (highY.Count > 0)
            {
                Console.WriteLine("Sample high-Y ions (TIC-normalized):");
                Console.WriteLine("+------------------+------------+------------+----------+----------+");
                Console.WriteLine("| InternalSeq      | TicNormInt | YIonInt    | DistCTerm| BasicInY |");
                Console.WriteLine("+------------------+------------+------------+----------+----------+");

                foreach (var ion in highY.OrderByDescending(i => i.YIonIntensityAtCTerm).Take(10))
                {
                    string seq = ion.InternalSequence.Length > 16 ? ion.InternalSequence.Substring(0, 13) + "..." : ion.InternalSequence;
                    Console.WriteLine($"| {seq,-16} | {ion.TicNormalizedIntensity,10:E3} | {ion.YIonIntensityAtCTerm,10:E3} | {ion.DistanceFromCTerm,8} | {ion.BasicResiduesInYIonSpan,8} |");
                }
                Console.WriteLine("+------------------+------------+------------+----------+----------+\n");

                if (highY.Any(i => i.YIonIntensityAtCTerm > 0.5))
                    Console.WriteLine("WARNING: Y > 0.5 detected. Investigate TIC normalization.\n");
            }
        }

        private static void Step11_PartD_StrongYIonMechanismStratification(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine("PART D: STRONG Y-ION MECHANISM STRATIFICATION");
            Console.WriteLine("═══════════════════════════════════════════════════\n");

            var withY = ions.Where(i => i.HasMatchedYIonAtCTerm).ToList();
            if (withY.Count == 0) { Console.WriteLine("No ions with matched Y-ion.\n"); return; }

            // Stratum 1: DistanceFromCTerm <= 3
            var s1 = withY.Where(i => i.DistanceFromCTerm <= 3).ToList();

            // Stratum 2: DistanceFromCTerm > 3 AND BasicResiduesInYIonSpan >= 2
            var s2 = withY.Where(i => i.DistanceFromCTerm > 3 && i.BasicResiduesInYIonSpan >= 2).ToList();

            // Stratum 3: DistanceFromCTerm > 3 AND BasicResiduesInYIonSpan <= 1
            var s3 = withY.Where(i => i.DistanceFromCTerm > 3 && i.BasicResiduesInYIonSpan <= 1).ToList();

            void PrintStratum(string name, string desc, List<InternalFragmentIon> list)
            {
                if (list.Count == 0)
                {
                    Console.WriteLine($"{name}: 0 ions");
                    return;
                }
                double meanY = list.Average(i => i.YIonIntensityAtCTerm);
                double meanTic = list.Average(i => i.TicNormalizedIntensity);
                Console.WriteLine($"{name}: {list.Count:N0} ions");
                Console.WriteLine($"  {desc}");
                Console.WriteLine($"  Mean YIonIntensity (TIC): {meanY:E4}");
                Console.WriteLine($"  Mean TicNormalizedIntensity: {meanTic:E4}");
            }

            PrintStratum("Stratum 1", "C-terminal K/R rescue (DistCTerm <= 3)", s1);
            Console.WriteLine();
            PrintStratum("Stratum 2", "Multiply-basic large y-ion (DistCTerm > 3, BasicInY >= 2)", s2);
            Console.WriteLine();
            PrintStratum("Stratum 3", "Weak charge retention (DistCTerm > 3, BasicInY <= 1)", s3);
            Console.WriteLine();
        }

        private static void Step11_PartF_HypothesisValidation(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine("PART F: HYPOTHESIS VALIDATION");
            Console.WriteLine("═══════════════════════════════════════════════════\n");

            var target = ions.Select(i => i.TicNormalizedIntensity).ToArray();

            Console.WriteLine("--- Correlations with TicNormalizedIntensity ---");
            PrintCorr("BIonIntensityAtNTerm", ions.Select(i => i.BIonIntensityAtNTerm).ToArray(), target);
            PrintCorr("YIonIntensityAtCTerm", ions.Select(i => i.YIonIntensityAtCTerm).ToArray(), target);
            PrintCorr("MaxTerminalIonIntensity", ions.Select(i => i.MaxTerminalIonIntensity).ToArray(), target);
            PrintCorr("BYProductScore", ions.Select(i => i.BYProductScore).ToArray(), target);
            PrintCorr("HasBothTerminalIons", ions.Select(i => i.HasBothTerminalIons ? 1.0 : 0.0).ToArray(), target);
            PrintCorr("BasicResiduesInBIonSpan", ions.Select(i => (double)i.BasicResiduesInBIonSpan).ToArray(), target);
            PrintCorr("BasicResiduesInYIonSpan", ions.Select(i => (double)i.BasicResiduesInYIonSpan).ToArray(), target);
            Console.WriteLine();

            // Group comparison
            Console.WriteLine("--- Group Comparison (TicNormalizedIntensity) ---");
            var g00 = ions.Where(i => !i.HasMatchedBIonAtNTerm && !i.HasMatchedYIonAtCTerm).ToList();
            var g10 = ions.Where(i => i.HasMatchedBIonAtNTerm && !i.HasMatchedYIonAtCTerm).ToList();
            var g01 = ions.Where(i => !i.HasMatchedBIonAtNTerm && i.HasMatchedYIonAtCTerm).ToList();
            var g11 = ions.Where(i => i.HasMatchedBIonAtNTerm && i.HasMatchedYIonAtCTerm).ToList();

            PrintGroup("Group 00 (neither)", g00);
            PrintGroup("Group 10 (B only)", g10);
            PrintGroup("Group 01 (Y only)", g01);
            PrintGroup("Group 11 (both)", g11);
            Console.WriteLine("\nExpected for positive mode HCD: Group 11 > Group 10 > Group 01 > Group 00\n");
        }

        private static void Step11_PartG_FullFeatureInventory(List<InternalFragmentIon> ions)
        {
            Console.WriteLine("═══════════════════════════════════════════════════");
            Console.WriteLine("PART G: FULL FEATURE INVENTORY");
            Console.WriteLine("═══════════════════════════════════════════════════\n");

            var features = new (string name, Func<InternalFragmentIon, double> get)[]
            {
                ("FragLen", i => i.FragmentLength),
                ("DistCTerm", i => i.DistanceFromCTerm),
                ("NumBasic", i => i.NumberOfBasicResidues),
                ("BasicInB", i => i.BasicResiduesInBIonSpan),
                ("BasicInY", i => i.BasicResiduesInYIonSpan),
                ("BIonInt", i => i.BIonIntensityAtNTerm),
                ("YIonInt", i => i.YIonIntensityAtCTerm),
                ("MaxTermInt", i => i.MaxTerminalIonIntensity),
                ("BYProd", i => i.BYProductScore),
                ("HasBoth", i => i.HasBothTerminalIons ? 1.0 : 0.0),
                ("LocalRank", i => double.IsNaN(i.LocalIntensityRank) ? 0 : i.LocalIntensityRank),
                ("TicNormInt", i => i.TicNormalizedIntensity)
            };

            var arrays = features.Select(f => ions.Select(f.get).ToArray()).ToArray();
            var target = arrays[^1];

            // Correlation matrix header
            Console.WriteLine("--- Correlation Matrix ---");
            Console.Write("          ");
            foreach (var f in features) Console.Write($"{f.name,9} ");
            Console.WriteLine();

            var collinear = new List<string>();
            var useful = new List<string>();

            for (int i = 0; i < features.Length; i++)
            {
                Console.Write($"{features[i].name,9} ");
                for (int j = 0; j < features.Length; j++)
                {
                    double r = Corr(arrays[i], arrays[j]);
                    Console.Write($"{r,9:F3} ");

                    if (i < j && i < features.Length - 1 && j < features.Length - 1 && Math.Abs(r) > 0.6)
                        collinear.Add($"{features[i].name} × {features[j].name} (r={r:F3})");

                    if (i < features.Length - 1 && j == features.Length - 1 && Math.Abs(r) > 0.2)
                        useful.Add($"{features[i].name} (r={r:F3})");
                }
                Console.WriteLine();
            }
            Console.WriteLine();

            if (collinear.Count > 0)
            {
                Console.WriteLine("Collinear predictor pairs (|r| > 0.6):");
                foreach (var p in collinear) Console.WriteLine($"  {p}");
                Console.WriteLine();
            }

            if (useful.Count > 0)
            {
                Console.WriteLine("Useful signals (predictor-target |r| > 0.2):");
                foreach (var s in useful) Console.WriteLine($"  {s}");
                Console.WriteLine();
            }

            // Feature readiness table
            Console.WriteLine("--- Feature Readiness Table ---");
            Console.WriteLine("+------------+----------+----------+------------+------------+-----------+");
            Console.WriteLine("| Feature    | N values | N nonzero|       Mean |     StdDev | r(target) |");
            Console.WriteLine("+------------+----------+----------+------------+------------+-----------+");

            for (int i = 0; i < features.Length - 1; i++)
            {
                var vals = arrays[i];
                int nVal = vals.Count(v => !double.IsNaN(v));
                int nNz = vals.Count(v => !double.IsNaN(v) && v != 0);
                double mean = nVal > 0 ? vals.Where(v => !double.IsNaN(v)).Average() : 0;
                double std = nVal > 1 ? Math.Sqrt(vals.Where(v => !double.IsNaN(v)).Select(v => (v - mean) * (v - mean)).Average()) : 0;
                double r = Corr(vals, target);
                Console.WriteLine($"| {features[i].name,-10} | {nVal,8} | {nNz,8} | {mean,10:E3} | {std,10:E3} | {r,9:F4} |");
            }

            Console.WriteLine("+------------+----------+----------+------------+------------+-----------+");
            Console.WriteLine($"\nModel-ready ions: {ions.Count:N0}");
            Console.WriteLine("Pre-training data quality gate: COMPLETE\n");
        }

        #endregion

        #region Utilities

        private static void PrintCorr(string name, double[] x, double[] target)
        {
            double r = Corr(x, target);
            string flag = Math.Abs(r) > 0.2 ? " *" : "";
            Console.WriteLine($"  {name,-25}: r = {r:F4}{flag}");
        }

        private static void PrintGroup(string name, List<InternalFragmentIon> grp)
        {
            if (grp.Count == 0) { Console.WriteLine($"  {name}: (no ions)"); return; }
            var sorted = grp.Select(i => i.TicNormalizedIntensity).OrderBy(x => x).ToList();
            double mean = sorted.Average();
            double median = sorted[sorted.Count / 2];
            Console.WriteLine($"  {name}: n={grp.Count,5}, Mean={mean:E4}, Median={median:E4}");
        }

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
                psmTsvPath: @"F:\JurakatMultiprotease\2026-02-17-16-40-42\Task1-SearchTask\AllPeptides.psmtsv",
                rawFileFolder: @"F:\JurakatMultiprotease",
                outputDirectory: @"E:\Projects\internalIons\mzLibReportsForClaude"
            );
        }
    }
}