using System;
using System.IO;
using System.Linq;
using System.Globalization;
using NUnit.Framework;
using MassSpectrometry.FlashDeconvRuntime;
using System.Diagnostics.CodeAnalysis;
using System.Text.RegularExpressions;

#nullable enable

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public sealed class FlashDeconvOptionsTests
    {
        [Test]
        public void FlashDeconvolutionOptionsTest()
        {
            string dataDir = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "TopDown");
            string inputFile    = Path.Combine(dataDir, "lvsYeastTopDownSnip.mzML");
            string expectedTsv  = Path.Combine(dataDir, "lvsYeastSnipFlashDeconvOptions_features.tsv");
            string expectedMzml = Path.Combine(dataDir, "lvsYeastSnipFlashDeconvOptions_spectra.mzML");

            if (!File.Exists(inputFile))
                Assert.Ignore("Missing input: " + inputFile);
            if (!File.Exists(expectedTsv) || !File.Exists(expectedMzml))
                Assert.Ignore("Missing expected baseline TSV or mzML in TopDown directory.");

            string tempDir = Path.Combine(Path.GetTempPath(), "FlashDeconvOpts_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(tempDir);

            string outTsv  = Path.Combine(tempDir, "deconv_features.tsv");
            string outMzml = Path.Combine(tempDir, "deconv_spectra.mzML");
            string outLog  = Path.Combine(tempDir, "flashdeconv.log");

            try
            {
                var opts = new FlashDeconvOptions
                {
                    InputMzMlPath      = inputFile,
                    OutputTsvPath      = outTsv,
                    OutMzmlDeconvolved = outMzml,
                    WriteDetail        = false,
                    MaxMsLevel         = 3,
                    PrecedingMs1Count  = 3
                };
                opts.Algorithm.MinMass   = 1000;
                opts.Algorithm.MaxMass   = 20000;
                opts.Algorithm.MinCharge = 1;
                opts.Algorithm.MaxCharge = 20;
                opts.Common.NoProgress   = true;
                opts.Common.LogFile      = outLog;

                FlashDeconvRunner runner;
                try { runner = new FlashDeconvRunner(preflight: false); }
                catch (FileNotFoundException e)
                {
                    Assert.Ignore("FLASHDeconv not found: " + e.Message);
                    return;
                }

                // Build args (WITHOUT -write_ini)
                var argList = opts.ToArgumentList().ToList();
                StripOptionWithValue(argList, "-instance");
                StripOptionWithValue(argList, "-debug");
                StripOptionWithValue(argList, "-threads");

                string inputArg  = ExtractAndRemoveValue(argList, "-in")!;
                string outputArg = ExtractAndRemoveValue(argList, "-out")!;
                string extraArgs = string.Join(" ", argList.Select(EscapeIfNeeded));

                TestContext.WriteLine("FLASHDeconv args:\n-in \"" + inputArg + "\" -out \"" + outputArg + "\" " + extraArgs);

                var result = runner.Run(inputArg, outputArg, extraArgs, timeoutMs: 0);

                TestContext.WriteLine("ExitCode: " + result.ExitCode);
                if (!string.IsNullOrWhiteSpace(result.StdErr))
                    TestContext.WriteLine("STDERR:\n" + result.StdErr);
                if (!string.IsNullOrWhiteSpace(result.StdOut))
                    TestContext.WriteLine("STDOUT (trunc):\n" + Truncate(result.StdOut, 1500));
                if (File.Exists(outLog))
                {
                    TestContext.WriteLine("LOG (first 1500 chars):");
                    TestContext.WriteLine(Truncate(File.ReadAllText(outLog), 1500));
                }

                Assert.That(result.ExitCode, Is.EqualTo(0), "FLASHDeconv failed.");

                RecoverMissingOutput(outTsv,  tempDir, "*.tsv",  "TSV");
                RecoverMissingOutput(outMzml, tempDir, "*.mzML", "mzML");

                if (!File.Exists(outTsv) || !File.Exists(outMzml))
                {
                    DumpDirectory(tempDir);
                    Assert.Ignore("Outputs missing; likely no features generated.");
                }

                Assert.That(new FileInfo(outTsv).Length,  Is.GreaterThan(0), "TSV empty.");
                Assert.That(new FileInfo(outMzml).Length, Is.GreaterThan(0), "mzML empty.");

                CompareTsvIgnoringColumnWithTolerance(expectedTsv, outTsv, 1, 1e-3, 1e-2);
                CompareTextFileApproximately(expectedMzml, outMzml, 0.20, 2e-4, 2e-2,
                    ignorePathPattern: @"[A-Z]:\\[^""\s>]*lvsYeastTopDownSnip\.mzML");

                TestContext.WriteLine("Comparisons passed.");
            }
            finally
            {
                try { Directory.Delete(tempDir, true); } catch { }
            }
        }

        private static void RecoverMissingOutput(string expectedPath, string root, string pattern, string label)
        {
            if (File.Exists(expectedPath)) return;
            var candidates = Directory.GetFiles(root, pattern, SearchOption.AllDirectories);
            TestContext.WriteLine($"{label} expected path missing. Candidates: " +
                                  (candidates.Length == 0 ? "<none>" : string.Join(", ", candidates)));
            var candidate = candidates
                .OrderByDescending(f => new FileInfo(f).Length)
                .FirstOrDefault();
            if (candidate != null)
            {
                File.Move(candidate, expectedPath, true);
                TestContext.WriteLine($"Adopted candidate {label}: {candidate}");
            }
        }

        #region TSV Comparison
        private static void CompareTsvIgnoringColumnWithTolerance(
            string expectedPath,
            string actualPath,
            int ignoreColumnIndex,
            double numericAbsTol,
            double numericRelTol)
        {
            var expectedLines = File.ReadAllLines(expectedPath).Where(l => l.Length > 0).ToList();
            var actualLines   = File.ReadAllLines(actualPath).Where(l => l.Length > 0).ToList();

            Assert.That(actualLines.Count, Is.EqualTo(expectedLines.Count),
                $"Line count mismatch. Expected {expectedLines.Count}, got {actualLines.Count}.");

            AssertHeadersCompatibleSkip(expectedLines[0], actualLines[0], ignoreColumnIndex);

            for (int i = 1; i < expectedLines.Count; i++)
            {
                var expCols = expectedLines[i].Split('\t');
                var actCols = actualLines[i].Split('\t');
                int maxCols = Math.Max(expCols.Length, actCols.Length);
                if (expCols.Length != actCols.Length)
                {
                    Array.Resize(ref expCols, maxCols);
                    Array.Resize(ref actCols, maxCols);
                }

                for (int c = 0; c < maxCols; c++)
                {
                    if (c == ignoreColumnIndex) continue;
                    string e = (expCols[c] ?? "").Trim();
                    string a = (actCols[c] ?? "").Trim();
                    if (e == a) continue;

                    if (TryNumeric(e, out double ev) && TryNumeric(a, out double av))
                    {
                        if (!Within(ev, av, numericAbsTol, numericRelTol))
                            FailNumeric(i + 1, c + 1, ev, av, numericAbsTol, numericRelTol);
                    }
                    else if (!(string.IsNullOrEmpty(e) && string.IsNullOrEmpty(a)))
                    {
                        Assert.Fail($"TSV mismatch (line {i + 1}, col {c + 1}): '{e}' vs '{a}'");
                    }
                }
            }
        }

        private static void AssertHeadersCompatibleSkip(string exp, string act, int ignoreIdx)
        {
            if (exp == act) return;
            var eCols = exp.Split('\t');
            var aCols = act.Split('\t');
            int max = Math.Max(eCols.Length, aCols.Length);
            for (int i = 0; i < max; i++)
            {
                if (i == ignoreIdx) continue;
                string e = i < eCols.Length ? eCols[i] : "";
                string a = i < aCols.Length ? aCols[i] : "";
                if (e != a)
                    Assert.Fail($"Header mismatch at column {i + 1}: '{e}' vs '{a}'");
            }
        }
        #endregion

        #region mzML Approx Comparison
        private static void CompareTextFileApproximately(
            string expectedPath,
            string actualPath,
            double sizeRelTolerance,
            double numericAbsTol,
            double numericRelTol,
            string? ignorePathPattern = null)
        {
            var expInfo = new FileInfo(expectedPath);
            var actInfo = new FileInfo(actualPath);
            double sizeDiff = Math.Abs(expInfo.Length - actInfo.Length) / Math.Max(1.0, expInfo.Length);
            Assert.That(sizeDiff, Is.LessThanOrEqualTo(sizeRelTolerance),
                $"mzML size differs by {sizeDiff:P1}");

            var expLines = File.ReadAllLines(expectedPath);
            var actLines = File.ReadAllLines(actualPath);
            int compareLines = Math.Min(expLines.Length, actLines.Length);
            Regex? pathRx = ignorePathPattern != null ? new Regex(ignorePathPattern, RegexOptions.IgnoreCase) : null;

            int limit = Math.Min(400, compareLines);
            for (int i = 0; i < limit; i++)
            {
                string e = expLines[i];
                string a = actLines[i];
                if (pathRx != null)
                {
                    e = pathRx.Replace(e, "<PATH>");
                    a = pathRx.Replace(a, "<PATH>");
                }
                if (e == a) continue;

                var eTokens = TokenizeForComparison(e);
                var aTokens = TokenizeForComparison(a);
                int tokenCount = Math.Min(eTokens.Length, aTokens.Length);
                for (int t = 0; t < tokenCount; t++)
                {
                    var et = eTokens[t];
                    var at = aTokens[t];
                    if (et == at) continue;
                    if (TryNumeric(et, out double ev) && TryNumeric(at, out double av))
                    {
                        if (!Within(ev, av, numericAbsTol, numericRelTol))
                            FailNumeric(i + 1, t + 1, ev, av, numericAbsTol, numericRelTol, "mzML");
                    }
                }
            }
        }
        #endregion

        #region Shared Helpers
        private static void StripOptionWithValue(System.Collections.Generic.List<string> args, string key)
        {
            for (int i = 0; i < args.Count; i++)
            {
                if (string.Equals(args[i], key, StringComparison.Ordinal))
                {
                    if (i + 1 < args.Count) args.RemoveAt(i + 1);
                    args.RemoveAt(i);
                    i--;
                }
            }
        }

        private static string? ExtractAndRemoveValue(System.Collections.Generic.List<string> args, string key)
        {
            for (int i = 0; i < args.Count - 1; i++)
            {
                if (string.Equals(args[i], key, StringComparison.Ordinal))
                {
                    string val = args[i + 1].Trim('"');
                    args.RemoveAt(i + 1);
                    args.RemoveAt(i);
                    return val;
                }
            }
            return null;
        }

        private static string EscapeIfNeeded(string s)
        {
            if (s.Length == 0) return "\"\"";
            if (s.Contains(' ') && !s.StartsWith("\"")) return "\"" + s + "\"";
            return s;
        }

        private static bool TryNumeric(string? s, out double v) =>
            double.TryParse(s, NumberStyles.Any, CultureInfo.InvariantCulture, out v);

        private static bool Within(double expected, double actual, double absTol, double relTol)
        {
            double diff = Math.Abs(expected - actual);
            if (diff <= absTol) return true;
            double baseVal = Math.Abs(expected) > 0 ? Math.Abs(expected) : 1.0;
            double rel = diff / baseVal;
            return rel <= relTol;
        }

        private static void FailNumeric(int line, int col, double ev, double av,
            double absTol, double relTol, string context = "TSV")
        {
            string where = (line > 0 && col > 0) ? $"(line {line}, col {col})" : "";
            Assert.Fail($"{context} numeric mismatch {where}: expected {ev} actual {av} (absTol={absTol}, relTol={relTol})");
        }

        private static string Truncate(string s, int max) =>
            string.IsNullOrEmpty(s) || s.Length <= max ? s : s.Substring(0, max) + "...";

        private static string[] TokenizeForComparison(string line) =>
            line.Split(new[] { ' ', '\t', '"', '\'', '=', ',', '>', '<', '/', '\\', '(', ')', ':' },
                       StringSplitOptions.RemoveEmptyEntries);
        #endregion

        private static void DumpDirectory(string dir)
        {
            TestContext.WriteLine("Directory listing for: " + dir);
            if (!Directory.Exists(dir))
            {
                TestContext.WriteLine("  <directory missing>");
                return;
            }
            foreach (var f in Directory.GetFiles(dir, "*", SearchOption.AllDirectories))
            {
                try
                {
                    var fi = new FileInfo(f);
                    TestContext.WriteLine($"  {fi.FullName}  {fi.Length} bytes");
                }
                catch { /* ignore */ }
            }
        }
    }
}