using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;

namespace Test.Deconvolution
{
    /// <summary>
    /// Tests for RealFLASHDeconvolution (NUnit 4).
    ///
    /// Uses SmallCalibratibleYeast.mzml (already in Test/DataFiles, copied to output dir).
    /// Verified from manual FLASHDeconv run:
    ///   - 3149 MS1 envelopes across 24 MS1 scans
    ///   - ScanNum = 1 for all rows (scan ID extraction failed but data is present)
    ///
    /// Key implementation note: -out_spec requires TWO paths for MS1+MS2 files.
    /// Passing only one path causes FLASHDeconv to write MS2 output there,
    /// leaving MS1 output missing entirely.
    ///
    /// Coverage strategy:
    ///   - ResolveExePath was refactored to accept the well-known-paths list and
    ///     PATH-env value as parameters, so its three branches (well-known found,
    ///     PATH-search hit, not-found-anywhere throw) can be exercised
    ///     deterministically (Section E below).
    ///   - Deconvolute and DeconvoluteFile accept an injected FLASHDeconvRunner
    ///     delegate. Tests pass a stub that writes canned TSV files instead of
    ///     spawning the real exe, so the surrounding orchestration (mzML writing,
    ///     range filtering, TSV parsing, missing-MS1 handling) is unit-testable
    ///     (Section F below).
    ///
    /// Accepted coverage gap (intentionally not unit-tested):
    ///   - RunFLASHDeconvDefault: pure subprocess plumbing -- spawn Process,
    ///     async-drain stdout/stderr, enforce timeout, throw on non-zero exit.
    ///     Marked [ExcludeFromCodeCoverage] so it doesn't drag the line-coverage
    ///     metric. Real-exe coverage is provided by the RealFLASH_RealExe_*
    ///     integration tests below when FLASHDECONV_EXE is set or OpenMS is
    ///     installed at KnownExePath.
    /// </summary>
    [TestFixture]
    public class TestRealFLASHDeconvolutionInfrastructure
    {
        private static readonly string KnownExePath =
            @"C:\Program Files\OpenMS-3.0.0-pre-HEAD-2023-06-17\bin\FLASHDeconv.exe";

        // ── A: Parameters ─────────────────────────────────────────────────────

        [Test]
        public void RealFLASH_Params_DefaultConstruct()
        {
            var p = new RealFLASHDeconvolutionParameters();
            Assert.That(p, Is.Not.Null);
            Assert.That(p.DeconvolutionType, Is.EqualTo(DeconvolutionType.RealFLASHDeconvolution));
        }

        [Test]
        public void RealFLASH_Params_Defaults()
        {
            var p = new RealFLASHDeconvolutionParameters();
            Assert.That(p.MinAssumedChargeState, Is.EqualTo(1));
            Assert.That(p.MaxAssumedChargeState, Is.EqualTo(60));
            Assert.That(p.TolerancePpm, Is.EqualTo(10.0).Within(1e-9));
            Assert.That(p.MinMass, Is.EqualTo(50.0).Within(1e-9));
            Assert.That(p.MaxMass, Is.EqualTo(100_000.0).Within(1e-9));
            Assert.That(p.MinIsotopeCosine, Is.EqualTo(0.85).Within(1e-9));
            Assert.That(p.ProcessTimeoutSeconds, Is.EqualTo(300));
            Assert.That(p.FLASHDeconvExePath, Is.Null);
        }

        [Test]
        public void RealFLASH_Params_ExplicitValues()
        {
            var p = new RealFLASHDeconvolutionParameters(
                minCharge: 2, maxCharge: 30, tolerancePpm: 5.0,
                minMass: 1000, maxMass: 50_000, minIsotopeCosine: 0.75,
                polarity: Polarity.Negative, flashDeconvExePath: KnownExePath);

            Assert.That(p.MinAssumedChargeState, Is.EqualTo(2));
            Assert.That(p.MaxAssumedChargeState, Is.EqualTo(30));
            Assert.That(p.TolerancePpm, Is.EqualTo(5.0).Within(1e-9));
            Assert.That(p.Polarity, Is.EqualTo(Polarity.Negative));
            Assert.That(p.FLASHDeconvExePath, Is.EqualTo(KnownExePath));
        }

        [Test]
        public void RealFLASH_Params_IsDeconvolutionParameters()
        {
            Assert.That(new RealFLASHDeconvolutionParameters(),
                Is.InstanceOf<DeconvolutionParameters>());
        }

        // ── B: Enum / factory ─────────────────────────────────────────────────

        [Test]
        public void RealFLASH_Enum_Exists()
        {
            Assert.That(Enum.IsDefined(typeof(DeconvolutionType),
                DeconvolutionType.RealFLASHDeconvolution), Is.True);
        }

        [Test]
        public void RealFLASH_Factory_EmptySpectrum_ReturnsEmptyWithoutCallingExe()
        {
            var p = new RealFLASHDeconvolutionParameters(
                flashDeconvExePath: @"C:\DoesNotExist\FLASHDeconv.exe");
            var empty = new MzSpectrum(Array.Empty<double>(), Array.Empty<double>(), shouldCopy: false);
            Assert.That(() =>
            {
                var result = Deconvoluter.Deconvolute(empty, p).ToList();
                Assert.That(result.Count, Is.EqualTo(0));
            }, Throws.Nothing);
        }

        [Test]
        public void Deconvolute_EmptySpectrum_RunnerIsNotInvoked()
        {
            // Pins the empty-spectrum short-circuit DIRECTLY via a stub runner that
            // fails the test if invoked, instead of relying on the indirect side
            // effect of an unresolvable exe path. A future refactor that moves or
            // restructures ResolveExePath could quietly cause _runner(...) to be
            // called on empty spectra without breaking the indirect test above.
            string fakeExe = CreateFakeExeFile();
            try
            {
                var p = new RealFLASHDeconvolutionParameters(flashDeconvExePath: fakeExe);
                RealFLASHDeconvolutionAlgorithm.FLASHDeconvRunner failIfCalled =
                    (_, _, _, _, _, _) => Assert.Fail("Runner must not be invoked for an empty spectrum.");

                var algorithm = new RealFLASHDeconvolutionAlgorithm(p, failIfCalled);
                var empty = new MzSpectrum(Array.Empty<double>(), Array.Empty<double>(), shouldCopy: false);

                var result = algorithm.Deconvolute(empty, new MzRange(0, 5000)).ToList();

                Assert.That(result, Is.Empty);
            }
            finally { File.Delete(fakeExe); }
        }

        [Test]
        public void RealFLASH_Factory_BadExePath_ThrowsFileNotFound()
        {
            var p = new RealFLASHDeconvolutionParameters(
                flashDeconvExePath: @"C:\DoesNotExist\FLASHDeconv.exe");
            Assert.That(
                () => Deconvoluter.Deconvolute(BuildSyntheticSpectrum(), p).ToList(),
                Throws.TypeOf<FileNotFoundException>());
        }

        [Test]
        public void DeconvoluteFile_MissingMzml_ThrowsFileNotFound()
        {
            // Pins the input-file guard at the top of DeconvoluteFile. Catching this
            // before launching FLASHDeconv produces a clear error message instead of
            // letting the exe fail with an opaque parse error.
            string nonexistent = Path.Combine(Path.GetTempPath(),
                $"realflash_missing_{Guid.NewGuid():N}.mzML");
            var p = new RealFLASHDeconvolutionParameters(
                flashDeconvExePath: @"C:\unused\FLASHDeconv.exe");

            Assert.That(
                () => RealFLASHDeconvolutionAlgorithm.DeconvoluteFile(nonexistent, p),
                Throws.TypeOf<FileNotFoundException>());
        }

        // ── C: Real-exe tests ─────────────────────────────────────────────────

        private static string GetExeOrSkip()
        {
            string? env = Environment.GetEnvironmentVariable("FLASHDECONV_EXE");
            if (!string.IsNullOrWhiteSpace(env) && File.Exists(env)) return env!;
            if (File.Exists(KnownExePath)) return KnownExePath;
            Assert.Inconclusive("FLASHDeconv not found. Set FLASHDECONV_EXE or install OpenMS at: " + KnownExePath);
            return null!;
        }

        private static string FindMzmlOrSkip(string filename)
        {
            var candidates = new List<string>
            {
                // Copied to test output directory by csproj
                Path.Combine(TestContext.CurrentContext.TestDirectory, filename),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", filename),
            };

            // Walk up to repo root
            string? dir = TestContext.CurrentContext.TestDirectory;
            for (int i = 0; i < 6 && dir != null; i++)
            {
                candidates.Add(Path.Combine(dir, "Test", "DataFiles", filename));
                candidates.Add(Path.Combine(dir, "mzLib", "Test", "DataFiles", filename));
                dir = Path.GetDirectoryName(dir);
            }

            string? found = candidates.FirstOrDefault(File.Exists);
            if (found != null) return found;
            Assert.Inconclusive($"{filename} not found. Searched:\n{string.Join("\n", candidates)}");
            return null!;
        }

        [Test]
        public void RealFLASH_RealExe_WholeFile_SmallYeast_ReturnsEnvelopes()
        {
            string exe = GetExeOrSkip();
            string mzml = FindMzmlOrSkip("SmallCalibratibleYeast.mzml");

            var p = new RealFLASHDeconvolutionParameters(
                flashDeconvExePath: exe,
                minCharge: 1, maxCharge: 10,
                minMass: 50, maxMass: 5000,
                minIsotopeCosine: 0.6);

            var byScan = RealFLASHDeconvolutionAlgorithm.DeconvoluteFile(mzml, p);
            int total = byScan.Values.Sum(l => l.Count);

            Console.WriteLine($"DeconvoluteFile: {byScan.Count} scan keys, {total} total envelopes.");

            Assert.That(total, Is.GreaterThan(0),
                "Expected envelopes from SmallCalibratibleYeast.");
        }

        [Test]
        public void RealFLASH_RealExe_WholeFile_SmallYeast_EnvelopeCountMatchesReference()
        {
            string exe = GetExeOrSkip();
            string mzml = FindMzmlOrSkip("SmallCalibratibleYeast.mzml");

            var p = new RealFLASHDeconvolutionParameters(
                flashDeconvExePath: exe,
                minCharge: 1, maxCharge: 10,
                minMass: 50, maxMass: 5000,
                minIsotopeCosine: 0.6);

            var byScan = RealFLASHDeconvolutionAlgorithm.DeconvoluteFile(mzml, p);
            int total = byScan.Values.Sum(l => l.Count);

            Console.WriteLine($"Total envelopes: {total}");

            // Reference run on OpenMS-3.0.0-pre-HEAD-2023-06-17 produced 3149 envelopes.
            // ~5% tolerance is generous enough to absorb legitimate cross-version drift
            // (minor OpenMS bumps, scoring-threshold tweaks) but tight enough that a
            // regression cutting the output by half would no longer silently pass --
            // matching what the test name "MatchesReference" promises.
            //
            // The 3149 constant is calibrated specifically to OpenMS-3.0.0-pre-HEAD-2023-06-17.
            // Contributors who run this test against a newer OpenMS (3.3.0 / 3.4.0 / 3.5.0,
            // all of which appear in the algorithm's well-known-paths search) may see the
            // count drift outside ±150 due to upstream filter or scoring-threshold tweaks
            // -- this is an environment difference, not an mzLib regression. Either update
            // the constant after verifying the new count is plausible, or run with
            // FLASHDECONV_EXE pointed at the reference version.
            Assert.That(total, Is.EqualTo(3149).Within(150),
                "Envelope count must be within 5% of the reference run (3149) to catch real regressions.");
        }

        [Test]
        public void RealFLASH_RealExe_WholeFile_SmallYeast_EnvelopeFieldsValid()
        {
            string exe = GetExeOrSkip();
            string mzml = FindMzmlOrSkip("SmallCalibratibleYeast.mzml");

            var p = new RealFLASHDeconvolutionParameters(
                flashDeconvExePath: exe,
                minCharge: 1, maxCharge: 10,
                minMass: 50, maxMass: 5000,
                minIsotopeCosine: 0.6);

            var allEnvs = RealFLASHDeconvolutionAlgorithm.DeconvoluteFile(mzml, p)
                .Values.SelectMany(x => x).ToList();

            Assert.That(allEnvs.Count, Is.GreaterThan(0),
                "Regression that empties FLASHDeconv output should be a hard failure, not Inconclusive");

            foreach (var env in allEnvs)
            {
                Assert.That(env.MonoisotopicMass, Is.GreaterThan(0));
                Assert.That(env.Charge, Is.Not.EqualTo(0));
                Assert.That(env.TotalIntensity, Is.GreaterThan(0));
                Assert.That(env.Peaks, Is.Not.Null.And.Not.Empty);
                Assert.That(env.Score, Is.GreaterThanOrEqualTo(0));
            }
        }

        [Test]
        public void RealFLASH_RealExe_WholeFile_SmallYeast_MassesInRange()
        {
            string exe = GetExeOrSkip();
            string mzml = FindMzmlOrSkip("SmallCalibratibleYeast.mzml");

            var p = new RealFLASHDeconvolutionParameters(
                flashDeconvExePath: exe,
                minCharge: 1, maxCharge: 10,
                minMass: 50, maxMass: 5000,
                minIsotopeCosine: 0.6);

            var allEnvs = RealFLASHDeconvolutionAlgorithm.DeconvoluteFile(mzml, p)
                .Values.SelectMany(x => x).ToList();

            Assert.That(allEnvs.Count, Is.GreaterThan(0),
                "Regression that empties FLASHDeconv output should be a hard failure, not Inconclusive");

            // All masses must be within the requested bounds
            Assert.That(allEnvs.All(e => e.MonoisotopicMass >= 50 && e.MonoisotopicMass <= 5000),
                Is.True, "All envelope masses should be within requested min/max mass bounds.");

            var top = allEnvs.MaxBy(e => e.TotalIntensity)!;
            Console.WriteLine($"Top: mass={top.MonoisotopicMass:F2} Da, z={top.Charge}, score={top.Score:F4}");
        }

        // ── D: Parser unit tests (no exe required) ────────────────────────────

        [TestCase(true)]
        [TestCase(false)]
        public void ParseSpecTsvByScan_TrailingTabInHeader_AllRowsParsed(bool trailingTab)
        {
            // Real FLASHDeconv output emits a trailing tab on the header row.
            // The parser strips empty header tokens; if that ever regresses,
            // every data row is rejected because c.Length < header.Length.
            string tsv = Path.Combine(Path.GetTempPath(), $"realflash_trailing_tab_{Guid.NewGuid():N}.tsv");
            try
            {
                string headerSuffix = trailingTab ? "\t" : "";
                File.WriteAllLines(tsv, new[]
                {
                    "ScanNum\tMonoisotopicMass\tSumIntensity\tRepresentativeCharge\t" +
                    "RepresentativeMzStart\tRepresentativeMzEnd\tIsotopeCosine\tMassSNR\tQscore" + headerSuffix,
                    "1\t1000.0\t5.0e5\t5\t201.0\t201.4\t0.95\t12.0\t0.8",
                    "1\t2000.5\t6.0e5\t8\t251.5\t251.9\t0.88\t14.0\t0.7",
                    "2\t1500.0\t4.5e5\t6\t251.0\t251.3\t0.90\t11.0\t0.75",
                });

                var p = new RealFLASHDeconvolutionParameters(
                    flashDeconvExePath: @"C:\unused\FLASHDeconv.exe",
                    minMass: 50, maxMass: 5000);

                var byScan = RealFLASHDeconvolutionAlgorithm.ParseSpecTsvByScan(tsv, p);

                int total = byScan.Values.Sum(l => l.Count);
                Assert.That(total, Is.EqualTo(3),
                    "All 3 data rows should parse regardless of trailing header tab.");
                Assert.That(byScan.Keys, Is.EquivalentTo(new[] { 1, 2 }));
            }
            finally
            {
                if (File.Exists(tsv)) File.Delete(tsv);
            }
        }

        [Test]
        public void WriteSingleScanMzml_DoesNotEmitBomAndSetsInstrumentConfigRef()
        {
            // OpenMS rejects a UTF-8 BOM with a misleading parse error, and refuses
            // <run>/<spectrum> elements that don't reference an instrumentConfiguration.
            // Both invariants are otherwise only catchable via a real-exe run.
            string mzml = Path.Combine(Path.GetTempPath(), $"realflash_writer_{Guid.NewGuid():N}.mzML");
            try
            {
                RealFLASHDeconvolutionAlgorithm.WriteSingleScanMzml(
                    BuildSyntheticSpectrum(), mzml, Polarity.Positive);

                byte[] firstThree = new byte[3];
                using (var fs = File.OpenRead(mzml))
                {
                    int read = fs.Read(firstThree, 0, 3);
                    Assert.That(read, Is.EqualTo(3));
                }
                Assert.That(firstThree, Is.Not.EqualTo(new byte[] { 0xEF, 0xBB, 0xBF }),
                    "Output must not start with a UTF-8 BOM (OpenMS rejects it).");

                string xml = File.ReadAllText(mzml);
                Assert.That(xml, Does.Contain("<run defaultInstrumentConfigurationRef=\"IC1\""),
                    "<run> must declare defaultInstrumentConfigurationRef=\"IC1\".");
                Assert.That(xml, Does.Contain("<spectrum")
                    .And.Contains("defaultInstrumentConfigurationRef=\"IC1\""),
                    "<spectrum> must declare defaultInstrumentConfigurationRef=\"IC1\".");
            }
            finally
            {
                if (File.Exists(mzml)) File.Delete(mzml);
            }
        }

        [Test]
        public void ParseSpecTsvByScan_NegativePolarity_AllChargesNegative()
        {
            string tsv = Path.Combine(Path.GetTempPath(), $"realflash_neg_{Guid.NewGuid():N}.tsv");
            try
            {
                // TSV uses positive z values; parser should sign-flip per Polarity.Negative.
                File.WriteAllLines(tsv, new[]
                {
                    "ScanNum\tMonoisotopicMass\tSumIntensity\tRepresentativeCharge\t" +
                    "RepresentativeMzStart\tRepresentativeMzEnd\tIsotopeCosine\tMassSNR\tQscore",
                    "1\t1000.0\t5.0e5\t5\t201.0\t201.4\t0.95\t12.0\t0.8",
                    "1\t2000.5\t6.0e5\t8\t251.5\t251.9\t0.88\t14.0\t0.7",
                    "2\t1500.0\t4.5e5\t6\t251.0\t251.3\t0.90\t11.0\t0.75",
                });

                var p = new RealFLASHDeconvolutionParameters(
                    polarity: Polarity.Negative,
                    flashDeconvExePath: @"C:\unused\FLASHDeconv.exe",
                    minMass: 50, maxMass: 5000);

                var allEnvs = RealFLASHDeconvolutionAlgorithm.ParseSpecTsvByScan(tsv, p)
                    .Values.SelectMany(x => x).ToList();

                Assume.That(allEnvs, Is.Not.Empty);
                Assert.That(allEnvs.All(e => e.Charge < 0), Is.True,
                    "All envelope charges must be negative when Polarity.Negative is set.");
            }
            finally
            {
                if (File.Exists(tsv)) File.Delete(tsv);
            }
        }

        /// <summary>
        /// Locks in the score-priority contract: Qscore (when parseable) wins, else
        /// MassSNR (when parseable), else IsotopeCosine. The cascade is driven by
        /// successful-parse, NOT by non-zero value — a Qscore that legitimately
        /// reads 0.0 must yield Score=0, not silently fall through to SNR.
        /// Cells written "NA" simulate a missing/unparseable column.
        /// </summary>
        [TestCase("0.7", "5.0", 0.7, TestName = "ScorePriority_QscoreNonZero_UsesQscore")]
        [TestCase("0", "5.0", 0.0, TestName = "ScorePriority_QscoreZero_StillUsesQscore")]
        [TestCase("NA", "5.0", 5.0, TestName = "ScorePriority_QscoreUnparseable_FallsToSnr")]
        [TestCase("NA", "NA", 0.95, TestName = "ScorePriority_QscoreAndSnrUnparseable_FallsToCosine")]
        public void ParseSpecTsvByScan_ScorePriorityContract(string qscoreCell, string snrCell, double expectedScore)
        {
            string tsv = Path.Combine(Path.GetTempPath(), $"realflash_score_{Guid.NewGuid():N}.tsv");
            try
            {
                File.WriteAllLines(tsv, new[]
                {
                    "ScanNum\tMonoisotopicMass\tSumIntensity\tRepresentativeCharge\t" +
                    "RepresentativeMzStart\tRepresentativeMzEnd\tIsotopeCosine\tMassSNR\tQscore",
                    $"1\t1000.0\t5.0e5\t5\t201.0\t201.4\t0.95\t{snrCell}\t{qscoreCell}",
                });

                var p = new RealFLASHDeconvolutionParameters(
                    flashDeconvExePath: @"C:\unused\FLASHDeconv.exe",
                    minMass: 50, maxMass: 5000);

                var allEnvs = RealFLASHDeconvolutionAlgorithm.ParseSpecTsvByScan(tsv, p)
                    .Values.SelectMany(x => x).ToList();

                Assume.That(allEnvs, Has.Count.EqualTo(1));
                Assert.That(allEnvs[0].Score, Is.EqualTo(expectedScore).Within(1e-9));
            }
            finally
            {
                if (File.Exists(tsv)) File.Delete(tsv);
            }
        }

        [Test]
        public void ParseSpecTsvByScan_RequiredColumnMissing_ThrowsMzLibException()
        {
            // Pins the required-column contract in Col(): a header missing any of
            // {ScanNum, MonoisotopicMass, RepresentativeCharge, SumIntensity,
            // RepresentativeMzStart, IsotopeCosine} must fail loudly with the
            // tried-name list -- not silently produce zero envelopes.
            string tsv = Path.Combine(Path.GetTempPath(),
                $"realflash_missingcol_{Guid.NewGuid():N}.tsv");
            try
            {
                // Header omits the required ScanNum column.
                File.WriteAllLines(tsv, new[]
                {
                    "MonoisotopicMass\tSumIntensity\tRepresentativeCharge\t" +
                    "RepresentativeMzStart\tRepresentativeMzEnd\tIsotopeCosine\tMassSNR\tQscore",
                    "1000.0\t5.0e5\t5\t201.0\t201.4\t0.95\t12.0\t0.8",
                });

                var p = new RealFLASHDeconvolutionParameters(
                    flashDeconvExePath: @"C:\unused\FLASHDeconv.exe",
                    minMass: 50, maxMass: 5000);

                Assert.That(
                    () => RealFLASHDeconvolutionAlgorithm.ParseSpecTsvByScan(tsv, p),
                    Throws.TypeOf<MzLibException>()
                          .With.Message.Contains("ScanNum"));
            }
            finally
            {
                if (File.Exists(tsv)) File.Delete(tsv);
            }
        }

        [Test]
        public void ParseSpecTsvByScan_OptionalColumnsAbsent_ParsesAndUsesFallbacks()
        {
            // Pins the optional-column behavior of ColOpt: when the header has none
            // of the alternative names, ColOpt returns -1 and downstream parsing
            // gracefully degrades (mzEnd missing -> repMz = mzStart; both score
            // sources missing -> falls through to IsotopeCosine).
            string tsv = Path.Combine(Path.GetTempPath(),
                $"realflash_nooptional_{Guid.NewGuid():N}.tsv");
            try
            {
                // Header has only the six required columns -- no RepresentativeMzEnd,
                // Qscore, or MassSNR. Exercises ColOpt's "return -1" path for all three.
                File.WriteAllLines(tsv, new[]
                {
                    "ScanNum\tMonoisotopicMass\tSumIntensity\tRepresentativeCharge\t" +
                    "RepresentativeMzStart\tIsotopeCosine",
                    "1\t1000.0\t5.0e5\t5\t201.0\t0.95",
                });

                var p = new RealFLASHDeconvolutionParameters(
                    flashDeconvExePath: @"C:\unused\FLASHDeconv.exe",
                    minMass: 50, maxMass: 5000);

                var allEnvs = RealFLASHDeconvolutionAlgorithm.ParseSpecTsvByScan(tsv, p)
                    .Values.SelectMany(x => x).ToList();

                Assert.That(allEnvs, Has.Count.EqualTo(1));
                // Score falls through to IsotopeCosine when both Qscore and MassSNR
                // columns are absent from the header.
                Assert.That(allEnvs[0].Score, Is.EqualTo(0.95).Within(1e-9));
            }
            finally
            {
                if (File.Exists(tsv)) File.Delete(tsv);
            }
        }

        [Test]
        public void ParseSpecTsvByScan_HeaderOnly_ReturnsEmptyDictionary()
        {
            // Pins the lines.Length < 2 early-return: FLASHDeconv occasionally
            // produces an MS1 TSV with only the header row when filters reject
            // every candidate envelope. The parser must return an empty dict
            // cleanly rather than throw or proceed into the row loop.
            string tsv = Path.Combine(Path.GetTempPath(),
                $"realflash_headeronly_{Guid.NewGuid():N}.tsv");
            try
            {
                File.WriteAllLines(tsv, new[]
                {
                    "ScanNum\tMonoisotopicMass\tSumIntensity\tRepresentativeCharge\t" +
                    "RepresentativeMzStart\tRepresentativeMzEnd\tIsotopeCosine\tMassSNR\tQscore",
                });

                var p = new RealFLASHDeconvolutionParameters(
                    flashDeconvExePath: @"C:\unused\FLASHDeconv.exe",
                    minMass: 50, maxMass: 5000);

                var result = RealFLASHDeconvolutionAlgorithm.ParseSpecTsvByScan(tsv, p);

                Assert.That(result, Is.Empty);
            }
            finally
            {
                if (File.Exists(tsv)) File.Delete(tsv);
            }
        }

        [Test]
        public void ParseSpecTsvByScan_RowsWithUnparseableRequiredCells_AreSilentlySkipped()
        {
            // Pins the per-row silent-skip cascade: real FLASHDeconv output occasionally
            // contains "NA" or other non-numeric cells in required columns (ScanNum,
            // MonoisotopicMass, RepresentativeCharge, ...). The parser's documented contract
            // is "skip the row, return what you can" -- not "throw" and not "include with
            // garbage values." Mix valid + unparseable rows in one file and verify only
            // the valid ones come back.
            string tsv = Path.Combine(Path.GetTempPath(),
                $"realflash_unparseable_{Guid.NewGuid():N}.tsv");
            try
            {
                File.WriteAllLines(tsv, new[]
                {
                    "ScanNum\tMonoisotopicMass\tSumIntensity\tRepresentativeCharge\t" +
                    "RepresentativeMzStart\tRepresentativeMzEnd\tIsotopeCosine\tMassSNR\tQscore",
                    "1\t1000.0\t5.0e5\t5\t201.0\t201.4\t0.95\t12.0\t0.8",       // valid
                    "NA\t2000.0\t6.0e5\t8\t251.5\t251.9\t0.88\t14.0\t0.7",     // bad ScanNum
                    "2\tNA\t6.0e5\t8\t251.5\t251.9\t0.88\t14.0\t0.7",          // bad mono
                    "3\t1500.0\t4.5e5\tNA\t251.0\t251.3\t0.90\t11.0\t0.75",    // bad charge
                    "4\t1800.0\t7.0e5\t6\t301.0\t301.4\t0.92\t13.0\t0.85",     // valid
                });

                var p = new RealFLASHDeconvolutionParameters(
                    flashDeconvExePath: @"C:\unused\FLASHDeconv.exe",
                    minMass: 50, maxMass: 5000);

                var byScan = RealFLASHDeconvolutionAlgorithm.ParseSpecTsvByScan(tsv, p);

                int total = byScan.Values.Sum(l => l.Count);
                Assert.That(total, Is.EqualTo(2),
                    "Only the two well-formed rows should have parsed; the three with unparseable required cells must be silently skipped.");
                Assert.That(byScan.Keys, Is.EquivalentTo(new[] { 1, 4 }));
            }
            finally
            {
                if (File.Exists(tsv)) File.Delete(tsv);
            }
        }

        // ── E0: FlashDeconvExePathRegistry (caching layer) ────────────────────

        [Test]
        public void FlashDeconvExePathRegistry_RegisterValidPath_AddsEntryAndCounts()
        {
            FlashDeconvExePathRegistry.Clear();
            string fake = CreateFakeExeFile();
            try
            {
                Assert.That(FlashDeconvExePathRegistry.Count, Is.EqualTo(0));
                FlashDeconvExePathRegistry.Register(fake);
                Assert.That(FlashDeconvExePathRegistry.Count, Is.EqualTo(1));
            }
            finally
            {
                FlashDeconvExePathRegistry.Clear();
                if (File.Exists(fake)) File.Delete(fake);
            }
        }

        [Test]
        public void FlashDeconvExePathRegistry_RegisterMissingPath_ThrowsFileNotFound()
        {
            FlashDeconvExePathRegistry.Clear();
            string nonexistent = Path.Combine(Path.GetTempPath(),
                $"realflash_missing_{Guid.NewGuid():N}.exe");

            Assert.That(
                () => FlashDeconvExePathRegistry.Register(nonexistent),
                Throws.TypeOf<FileNotFoundException>());
            Assert.That(FlashDeconvExePathRegistry.Count, Is.EqualTo(0),
                "Failed registration must not leak a partial entry into the registry.");
        }

        [Test]
        public void FlashDeconvExePathRegistry_RegisterNullOrEmpty_ThrowsArgumentException()
        {
            Assert.That(() => FlashDeconvExePathRegistry.Register(null!),
                Throws.TypeOf<ArgumentException>());
            Assert.That(() => FlashDeconvExePathRegistry.Register(""),
                Throws.TypeOf<ArgumentException>());
            Assert.That(() => FlashDeconvExePathRegistry.Register("   "),
                Throws.TypeOf<ArgumentException>());
        }

        [Test]
        public void Deconvolute_WithRegisteredPath_SkipsFilesystemValidation()
        {
            // Pin the optimization: once a path is in the registry, the algorithm
            // uses it without re-checking File.Exists. Prove the cache is being
            // consulted by registering the path, deleting the file, and observing
            // that Deconvolute still proceeds to the stub runner -- a non-cached
            // resolve would throw FileNotFoundException at this point.
            FlashDeconvExePathRegistry.Clear();
            string fakeExe = CreateFakeExeFile();
            FlashDeconvExePathRegistry.Register(fakeExe);
            File.Delete(fakeExe);  // cache says it exists; filesystem says otherwise

            try
            {
                var p = new RealFLASHDeconvolutionParameters(
                    flashDeconvExePath: fakeExe,
                    minMass: 50, maxMass: 5000);

                string tsv = string.Join("\n",
                    Ms1TsvHeader,
                    "1\t1000.0\t5.0e5\t5\t201.0\t201.4\t0.95\t12.0\t0.8");

                var algorithm = new RealFLASHDeconvolutionAlgorithm(p, BuildStubRunner(tsv));

                var result = algorithm
                    .Deconvolute(BuildSyntheticSpectrum(), new MzRange(0, 5000))
                    .ToList();

                Assert.That(result, Has.Count.EqualTo(1),
                    "Cached path must be honored without re-validating that the file still exists; the deleted file would otherwise throw before the stub runner runs.");
            }
            finally
            {
                FlashDeconvExePathRegistry.Clear();
            }
        }

        // ── E: ResolveExePath unit tests ──────────────────────────────────────

        [Test]
        public void ResolveExePath_NoExplicit_FoundInWellKnown_ReturnsThatPath()
        {
            // Pins the well-known-search branch: when no explicit path is given
            // and one of the well-known paths exists on disk, that path is returned.
            string fake = CreateFakeExeFile();
            try
            {
                string result = RealFLASHDeconvolutionAlgorithm.ResolveExePath(
                    explicitPath: null,
                    wellKnownPaths: new[] { fake },
                    pathEnv: "");

                Assert.That(result, Is.EqualTo(fake));
            }
            finally { File.Delete(fake); }
        }

        [Test]
        public void ResolveExePath_NoExplicit_NoWellKnown_FoundInPathEnv_ReturnsFromPath()
        {
            // Pins the PATH-search branch: when no explicit path is given, no
            // well-known path matches, but a directory on PATH contains the exe,
            // it is found and returned. Uses an isolated temp directory as the
            // sole PATH entry so the test isn't affected by what the host has.
            string dir = Path.Combine(Path.GetTempPath(), $"realflash_path_{Guid.NewGuid():N}");
            Directory.CreateDirectory(dir);
            string fakeName = OperatingSystem.IsWindows() ? "FLASHDeconv.exe" : "FLASHDeconv";
            string fake = Path.Combine(dir, fakeName);
            File.WriteAllText(fake, "stub");
            try
            {
                string result = RealFLASHDeconvolutionAlgorithm.ResolveExePath(
                    explicitPath: null,
                    wellKnownPaths: Array.Empty<string>(),
                    pathEnv: dir);

                Assert.That(result, Is.EqualTo(fake));
            }
            finally { Directory.Delete(dir, recursive: true); }
        }

        [Test]
        public void ResolveExePath_NoExplicit_NoWellKnown_NoPath_ThrowsFileNotFound()
        {
            // Pins the not-found-anywhere throw with a clear message that points
            // the user at the override knob.
            Assert.That(
                () => RealFLASHDeconvolutionAlgorithm.ResolveExePath(
                    explicitPath: null,
                    wellKnownPaths: Array.Empty<string>(),
                    pathEnv: ""),
                Throws.TypeOf<FileNotFoundException>()
                      .With.Message.Contains("FLASHDeconvExePath"));
        }

        // ── F: Orchestration unit tests with stubbed FLASHDeconvRunner ────────

        [Test]
        public void Deconvolute_StubRunnerWritesValidTsv_ReturnsParsedEnvelopes()
        {
            // Pins the per-scan happy path: WriteSingleScanMzml -> stub runner
            // produces canned MS1 TSV -> ParseSpecTsvByScan -> range filter ->
            // returned envelopes. No real exe, no integration assumption.
            string fakeExe = CreateFakeExeFile();
            try
            {
                var p = new RealFLASHDeconvolutionParameters(
                    flashDeconvExePath: fakeExe,
                    minMass: 50, maxMass: 5000);

                string tsv = string.Join("\n",
                    Ms1TsvHeader,
                    "1\t1000.0\t5.0e5\t5\t201.0\t201.4\t0.95\t12.0\t0.8");

                var algorithm = new RealFLASHDeconvolutionAlgorithm(p, BuildStubRunner(tsv));

                var result = algorithm
                    .Deconvolute(BuildSyntheticSpectrum(), new MzRange(0, 5000))
                    .ToList();

                Assert.That(result, Has.Count.EqualTo(1));
                Assert.That(result[0].MonoisotopicMass, Is.EqualTo(1000.0).Within(1e-6));
            }
            finally { File.Delete(fakeExe); }
        }

        [Test]
        public void Deconvolute_StubRunnerWritesNoMs1_ReturnsEmpty()
        {
            // Pins the missing-MS1-file early-return at the top of the parse step.
            // FLASHDeconv occasionally produces no MS1 TSV when no envelopes pass
            // the configured filters; the algorithm must return empty cleanly
            // rather than throw on a missing file.
            string fakeExe = CreateFakeExeFile();
            try
            {
                var p = new RealFLASHDeconvolutionParameters(flashDeconvExePath: fakeExe);
                var algorithm = new RealFLASHDeconvolutionAlgorithm(
                    p, BuildStubRunner(ms1Content: null));

                var result = algorithm
                    .Deconvolute(BuildSyntheticSpectrum(), new MzRange(0, 5000))
                    .ToList();

                Assert.That(result, Is.Empty);
            }
            finally { File.Delete(fakeExe); }
        }

        [Test]
        public void Deconvolute_StubRunner_FiltersEnvelopesByRange()
        {
            // Pins the range-filter behaviour: even when the runner produces
            // envelopes outside the requested mz range, only those whose
            // representative-mz midpoint is within range are returned.
            string fakeExe = CreateFakeExeFile();
            try
            {
                var p = new RealFLASHDeconvolutionParameters(
                    flashDeconvExePath: fakeExe,
                    minMass: 50, maxMass: 100_000);

                string tsv = string.Join("\n",
                    Ms1TsvHeader,
                    "1\t1000.0\t5.0e5\t5\t201.0\t201.4\t0.95\t12.0\t0.8",
                    "1\t2000.0\t6.0e5\t8\t9999.0\t9999.5\t0.88\t14.0\t0.7");

                var algorithm = new RealFLASHDeconvolutionAlgorithm(p, BuildStubRunner(tsv));

                var result = algorithm
                    .Deconvolute(BuildSyntheticSpectrum(), new MzRange(0, 1000))
                    .ToList();

                Assert.That(result, Has.Count.EqualTo(1));
                Assert.That(result[0].MonoisotopicMass, Is.EqualTo(1000.0).Within(1e-6));
            }
            finally { File.Delete(fakeExe); }
        }

        [Test]
        public void Deconvolute_StubRunner_RangeFilterIsInclusiveOnUpperBound()
        {
            // Pins the boundary semantic: MzRange.Contains is inclusive on both ends
            // (DoubleRange.Contains, "True if the item is within the range (inclusive)"),
            // so an envelope whose representative-mz midpoint sits EXACTLY at the upper
            // bound must be kept, and one that's epsilon past it must be dropped.
            // Without this, a future change to MzRange's endpoint behaviour or to the
            // Where-clause expression could silently shift edge results.
            string fakeExe = CreateFakeExeFile();
            try
            {
                var p = new RealFLASHDeconvolutionParameters(
                    flashDeconvExePath: fakeExe,
                    minMass: 50, maxMass: 100_000);

                // Three rows: midpoint exactly at 1000.0 (must be IN), midpoint just past
                // (1000.05, must be OUT), and clearly in (200.0) as a sanity baseline.
                // Midpoint = (mzStart + mzEnd) / 2.
                string tsv = string.Join("\n",
                    Ms1TsvHeader,
                    "1\t1000.0\t5.0e5\t5\t999.9\t1000.1\t0.95\t12.0\t0.8",  // mid=1000.0, IN (inclusive)
                    "1\t2000.0\t6.0e5\t8\t999.9\t1000.2\t0.88\t14.0\t0.7",  // mid=1000.05, OUT
                    "1\t3000.0\t7.0e5\t6\t199.5\t200.5\t0.92\t13.0\t0.85"); // mid=200.0, IN baseline

                var algorithm = new RealFLASHDeconvolutionAlgorithm(p, BuildStubRunner(tsv));

                var result = algorithm
                    .Deconvolute(BuildSyntheticSpectrum(), new MzRange(0, 1000))
                    .ToList();

                Assert.That(result, Has.Count.EqualTo(2),
                    "Range filter is inclusive on the upper bound: 1000.0 is IN, 1000.05 is OUT, baseline is IN.");
                Assert.That(result.Select(e => e.MonoisotopicMass).OrderBy(m => m),
                    Is.EqualTo(new[] { 1000.0, 3000.0 }).AsCollection,
                    "Expected the at-bound and baseline envelopes; the just-past envelope must be dropped.");
            }
            finally { File.Delete(fakeExe); }
        }

        [Test]
        public void DeconvoluteFile_StubRunnerWritesValidTsv_ReturnsScansKeyedByScanNum()
        {
            // Pins the whole-file orchestration: input mzML existence check ->
            // stub runner writes canned MS1 TSV -> ParseSpecTsvByScan groups
            // by scan number -> returned dictionary. mzML content is irrelevant
            // because the runner is stubbed.
            string fakeExe = CreateFakeExeFile();
            string fakeMzml = Path.Combine(Path.GetTempPath(),
                $"realflash_in_{Guid.NewGuid():N}.mzML");
            File.WriteAllText(fakeMzml, "<mzML/>");
            try
            {
                var p = new RealFLASHDeconvolutionParameters(
                    flashDeconvExePath: fakeExe,
                    minMass: 50, maxMass: 5000);

                string tsv = string.Join("\n",
                    Ms1TsvHeader,
                    "1\t1000.0\t5.0e5\t5\t201.0\t201.4\t0.95\t12.0\t0.8",
                    "2\t2000.0\t6.0e5\t8\t251.5\t251.9\t0.88\t14.0\t0.7");

                var byScan = RealFLASHDeconvolutionAlgorithm.DeconvoluteFile(
                    fakeMzml, p, BuildStubRunner(tsv));

                Assert.That(byScan.Keys, Is.EquivalentTo(new[] { 1, 2 }));
                Assert.That(byScan.Values.Sum(l => l.Count), Is.EqualTo(2));
            }
            finally
            {
                File.Delete(fakeExe);
                File.Delete(fakeMzml);
            }
        }

        [Test]
        public void DeconvoluteFile_StubRunnerWritesNoMs1_ReturnsEmptyDictionary()
        {
            // Pins the missing-MS1-file early-return for the whole-file path.
            string fakeExe = CreateFakeExeFile();
            string fakeMzml = Path.Combine(Path.GetTempPath(),
                $"realflash_in_{Guid.NewGuid():N}.mzML");
            File.WriteAllText(fakeMzml, "<mzML/>");
            try
            {
                var p = new RealFLASHDeconvolutionParameters(flashDeconvExePath: fakeExe);

                var byScan = RealFLASHDeconvolutionAlgorithm.DeconvoluteFile(
                    fakeMzml, p, BuildStubRunner(ms1Content: null));

                Assert.That(byScan, Is.Empty);
            }
            finally
            {
                File.Delete(fakeExe);
                File.Delete(fakeMzml);
            }
        }

        [Test]
        public void Deconvolute_WrongParameterType_ThrowsMzLibException()
        {
            // Pins the up-front parameter-type guard at the top of Deconvolute:
            // constructing the algorithm with a parameter class meant for a
            // different deconvolution algorithm must throw clearly rather than
            // proceed and fail deeper in.
            var wrongParams = new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0);
            var algorithm = new RealFLASHDeconvolutionAlgorithm(wrongParams);

            Assert.That(
                () => algorithm.Deconvolute(BuildSyntheticSpectrum(), new MzRange(0, 5000)).ToList(),
                Throws.TypeOf<MzLibException>().With.Message.Contains("do not match"));
        }

        // ── Stub-runner helpers ───────────────────────────────────────────────

        private const string Ms1TsvHeader =
            "ScanNum\tMonoisotopicMass\tSumIntensity\tRepresentativeCharge\t" +
            "RepresentativeMzStart\tRepresentativeMzEnd\tIsotopeCosine\tMassSNR\tQscore";

        private static string CreateFakeExeFile()
        {
            string path = Path.Combine(Path.GetTempPath(),
                $"realflash_fake_exe_{Guid.NewGuid():N}.exe");
            File.WriteAllText(path, "stub");
            return path;
        }

        private static RealFLASHDeconvolutionAlgorithm.FLASHDeconvRunner BuildStubRunner(
            string? ms1Content)
            => (exe, inMzml, outFeat, outMs1, outMs2, p) =>
            {
                if (ms1Content != null) File.WriteAllText(outMs1, ms1Content);
            };

        // ── Helper ────────────────────────────────────────────────────────────

        private static MzSpectrum BuildSyntheticSpectrum()
        {
            const double monoMass = 12_221.0, proton = 1.007276, isoDelta = 1.003355;
            var mzs = new List<double>(); var ints = new List<double>();
            for (int z = 8; z <= 16; z++)
            {
                double baseMz = (monoMass + z * proton) / z;
                double baseInt = 1e6 * Math.Exp(-0.5 * Math.Pow((z - 12.0) / 2.5, 2));
                for (int iso = 0; iso <= 2; iso++)
                { mzs.Add(baseMz + iso * isoDelta / z); ints.Add(baseInt * Math.Exp(-iso * 0.4)); }
            }
            var pairs = mzs.Zip(ints, (m, i) => (m, i)).OrderBy(x => x.m).ToList();
            return new MzSpectrum(pairs.Select(x => x.m).ToArray(), pairs.Select(x => x.i).ToArray(), shouldCopy: false);
        }
    }
}