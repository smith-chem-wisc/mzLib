using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MassSpectrometry;
using MassSpectrometry.Deconvolution.FeatureTracing;
using NUnit.Framework;
using Readers;

namespace Test.MassSpectrometryTests.Deconvolution
{
    /// <summary>
    /// EXPLICIT, Phase-6 driver: for every .mzML in <see cref="FolderPath"/>, run MetaFlashDecon
    /// per-scan deconvolution and emit a sibling <c>_ms1.feature</c> (RT in seconds) with BOTH
    /// feature tracers — the FLASHDeconv-faithful native tracer and the consensus tracer — into
    /// separate output folders. This produces the head-to-head input for comparing against the
    /// real OpenMS FLASHDeconv <c>_ms1.feature</c> files in
    /// <c>E:\Projects\JurkatTopDown\FlashDecon</c> and for MetaMorpheus FromFile searches.
    ///
    /// Not part of CI; run manually with:
    ///   dotnet test --filter FullyQualifiedName~JurkatMetaFlashDeconGenerator
    /// </summary>
    [TestFixture]
    [Explicit("Generates MetaFlashDecon _ms1.feature files (both tracers) for the local Jurkat TopDown dataset; needs E:\\Projects\\JurkatTopDown.")]
    [ExcludeFromCodeCoverage]
    public class JurkatMetaFlashDeconGenerator
    {
        private const string FolderPath = @"E:\Projects\JurkatTopDown";
        private const string OutputRoot = @"E:\TestData\MetaMorpheus\runs\metaflashdecon_features";

        // Top-down MetaFlashDecon parameters (charge 1..60).
        private const int MinCharge = 1;
        private const int MaxCharge = 60;

        [Test]
        public void GenerateNativeAndConsensusFeatureFiles()
        {
            Assume.That(Directory.Exists(FolderPath), $"Folder not found: {FolderPath}");
            var mzmls = Directory.EnumerateFiles(FolderPath, "*.mzML").OrderBy(p => p).ToList();
            Assume.That(mzmls, Is.Not.Empty, $"No mzML files at {FolderPath}");

            var deconParams = new MetaFlashDeconParameters(minCharge: MinCharge, maxCharge: MaxCharge);

            // Same MetaFlashDecon envelopes, two tracers -> the side-by-side.
            var tracers = new (string Tag, IMassFeatureTracer Tracer)[]
            {
                ("native", new MetaFlashDeconMassFeatureTracer()),
                ("consensus", new ConsensusMassFeatureTracer()),
            };

            var overall = Stopwatch.StartNew();
            int produced = 0, skipped = 0, failed = 0;
            var failures = new List<string>();

            foreach (string mzml in mzmls)
            {
                string name = Path.GetFileNameWithoutExtension(mzml);
                foreach (var (tag, tracer) in tracers)
                {
                    string outDir = Path.Combine(OutputRoot, tag);
                    Directory.CreateDirectory(outDir);
                    string outputPath = Path.Combine(outDir, name + "_ms1.feature");
                    if (File.Exists(outputPath))
                    {
                        skipped++;
                        continue;
                    }

                    try
                    {
                        var sw = Stopwatch.StartNew();
                        var features = Ms1FeatureGenerator.GenerateAndWrite(mzml, deconParams, tracer, outputPath);
                        sw.Stop();
                        TestContext.Progress.WriteLine(
                            $"[{tag}/{name}] features={features.Count} " +
                            $"(multi-z={features.Count(f => f.ChargeCount >= 2)}) in {sw.Elapsed}");
                        produced++;
                    }
                    catch (Exception ex)
                    {
                        failed++;
                        failures.Add($"{tag}/{name}: {ex.GetType().Name}: {ex.Message}");
                        TestContext.Progress.WriteLine($"[{tag}/{name}] FAILED: {ex.Message}");
                    }
                }
            }

            overall.Stop();
            TestContext.Progress.WriteLine(
                $"[MetaFlashDecon] produced={produced}, skipped={skipped}, failed={failed}, elapsed {overall.Elapsed}");
            if (failures.Count > 0)
            {
                foreach (var f in failures) TestContext.Progress.WriteLine("  fail: " + f);
                Assert.Fail($"{failures.Count} generation(s) failed; see output.");
            }
        }
    }
}
