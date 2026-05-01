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
    /// Accepted coverage gaps (intentionally not unit-tested):
    ///   - ResolveExePath null-path fallback (well-known list + PATH search, then
    ///     throw): the well-known list hardcodes installed OpenMS layouts; a test
    ///     that passes a null exe path resolves via well-known on dev machines
    ///     where OpenMS is installed and never reaches the throw, while CI runners
    ///     without OpenMS take the throw path. Same test, different result by
    ///     environment -- flaky. The deterministic failure mode (explicit path
    ///     missing) IS covered by RealFLASH_Factory_BadExePath_ThrowsFileNotFound.
    ///     Production callers set GlobalSettings.FLASHDeconvExecutablePath, so the
    ///     null-path branch is edge-case fallback only.
    ///   - Deconvolute body (process invocation + temp-file orchestration) and
    ///     RunFLASHDeconv (timeout / non-zero-exit handling): require the real
    ///     FLASHDeconv exe. Covered by the RealFLASH_RealExe_* integration tests
    ///     below when FLASHDECONV_EXE is set or OpenMS is installed at KnownExePath.
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
            // We assert a sanity floor rather than the exact count so the test stays
            // green across OpenMS versions.
            Assert.That(total, Is.GreaterThan(1000),
                "Envelope count should be substantial (reference run produced 3149).");
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

            Assume.That(allEnvs.Count, Is.GreaterThan(0));

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

            Assume.That(allEnvs.Count, Is.GreaterThan(0));

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