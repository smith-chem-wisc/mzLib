using MassSpectrometry;
using NUnit.Framework;
using Readers;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Development
{
    /// <summary>
    /// Benchmarks <see cref="MsDataFile.LoadAllStaticData"/> for a local spectra file across a range of
    /// thread counts. It is <c>[Explicit]</c> because it needs a local data file that is not part of the
    /// repository, so it never runs in CI - invoke it by hand when you want to measure reader performance.
    ///
    /// Point it at a file either by setting the <c>MZLIB_BENCHMARK_FILE</c> environment variable or by
    /// editing <see cref="DefaultBenchmarkFile"/> below. Any format the reader supports works
    /// (.raw, .mzML, .mgf, ...), so this doubles as a general file-reading benchmark, not just a Thermo one.
    ///
    /// Run from the CLI, e.g.:
    ///   dotnet test mzLib/Development/Development.csproj --filter "FullyQualifiedName~Benchmark_FileReading"
    /// and read the table it writes to the test output.
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class BenchmarkFileReading
    {
        // Adjust this path to a local file on your machine, or set the MZLIB_BENCHMARK_FILE environment
        // variable (which takes precedence). A ~1 GB .raw file makes the thread-count differences obvious.
        private const string DefaultBenchmarkFile = @"D:\D_Morgen_Glyco\RawFiles\HFX_MB_14751_19.raw";

        // Number of timed repetitions per thread-count configuration; best and median are reported.
        private const int Repetitions = 3;

        [Test]
        [Explicit("Requires a local spectra file; set MZLIB_BENCHMARK_FILE or edit DefaultBenchmarkFile.")]
        public static void Benchmark_FileReading()
        {
            string path = Environment.GetEnvironmentVariable("MZLIB_BENCHMARK_FILE") ?? DefaultBenchmarkFile;

            if (!File.Exists(path))
            {
                Assert.Ignore($"Benchmark file not found: '{path}'. " +
                              "Set the MZLIB_BENCHMARK_FILE environment variable or edit DefaultBenchmarkFile.");
            }

            // Thread counts to sweep, clamped to the machine and de-duplicated (order preserved).
            int cores = Environment.ProcessorCount;
            var threadCounts = new[] { 1, 2, 4, 8, cores }
                .Where(t => t >= 1 && t <= cores)
                .Distinct()
                .ToArray();

            long fileBytes = new FileInfo(path).Length;

            // Warm up: pay the JIT and OS-file-cache costs once, before any timing.
            var warm = MsDataFileReader.GetDataFile(path).LoadAllStaticData(maxThreads: 4);
            int scanCount = warm.NumSpectra;

            TestContext.Out.WriteLine($"File:                {path}");
            TestContext.Out.WriteLine($"Size:                {fileBytes / 1024.0 / 1024.0:F1} MB");
            TestContext.Out.WriteLine($"Scans:               {scanCount}");
            TestContext.Out.WriteLine($"Logical processors:  {cores}");
            TestContext.Out.WriteLine($"Repetitions:         {Repetitions} (best / median reported)");
            TestContext.Out.WriteLine("");
            TestContext.Out.WriteLine("maxThreads=1 is sequential reading; 'speedup' is relative to it.");
            TestContext.Out.WriteLine($"{"maxThreads",-12}{"best (ms)",12}{"median (ms)",14}{"speedup",10}{"peaks",16}");
            TestContext.Out.WriteLine(new string('-', 64));

            long referencePeaks = -1;
            double sequentialBest = 0;

            foreach (int maxThreads in threadCounts)
            {
                var times = new List<double>();
                long peaks = 0;

                for (int rep = 0; rep < Repetitions; rep++)
                {
                    GC.Collect();
                    GC.WaitForPendingFinalizers();
                    GC.Collect();

                    var sw = Stopwatch.StartNew();
                    var file = MsDataFileReader.GetDataFile(path).LoadAllStaticData(maxThreads: maxThreads);
                    sw.Stop();

                    times.Add(sw.Elapsed.TotalMilliseconds);
                    peaks = SumPeaks(file);
                }

                var sorted = times.OrderBy(t => t).ToList();
                double best = sorted[0];
                double median = sorted[sorted.Count / 2];

                // The first configuration is the single-threaded (sequential) baseline.
                if (sequentialBest == 0)
                    sequentialBest = best;
                double speedup = sequentialBest / best;

                TestContext.Out.WriteLine($"{maxThreads,-12}{best,12:F0}{median,14:F0}{speedup,9:F2}x{peaks,16}");

                // Every configuration must read the same data - guards against a threading bug that drops
                // or duplicates scans, and confirms the timings are comparing identical work.
                if (referencePeaks < 0)
                    referencePeaks = peaks;
                else
                    Assert.AreEqual(referencePeaks, peaks, $"maxThreads={maxThreads} read a different peak count.");
            }

            Assert.Greater(scanCount, 0, "No scans were read from the benchmark file.");
        }

        private static long SumPeaks(MsDataFile file)
        {
            long peaks = 0;
            foreach (var scan in file.GetAllScansList())
            {
                if (scan?.MassSpectrum != null)
                    peaks += scan.MassSpectrum.Size;
            }
            return peaks;
        }
    }
}
