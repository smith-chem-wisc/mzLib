using System;
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
    /// EXPLICIT whole-file evaluation: run MetaFlashDecon on ONE real Jurkat top-down mzML, trace
    /// cross-scan features with the native tracer, and report the WHOLE-FILE feature count against
    /// the real OpenMS FLASHDeconv ground-truth <c>_ms1.feature</c>. Per-scan over-generation can
    /// collapse after tracing because the tracer requires a mass to persist across several scans.
    ///
    /// Two methods so the cheap, high-value GATED result can be obtained without the heavy ungated
    /// tail (ungated traces ~10x more envelopes single-threaded). Each method appends its result to
    /// <see cref="ResultsLog"/> immediately (durable, unbuffered), so progress is visible live.
    ///
    /// Run manually:
    ///   dotnet test --filter FullyQualifiedName~MetaFlashDeconWholeFileEval.Gated
    ///   dotnet test --filter FullyQualifiedName~MetaFlashDeconWholeFileEval.Ungated
    /// </summary>
    [TestFixture]
    [Explicit("Whole-file MetaFlashDecon eval vs real FLASHDeconv ground truth; needs E:\\Projects\\JurkatTopDown.")]
    [ExcludeFromCodeCoverage]
    public class MetaFlashDeconWholeFileEval
    {
        private const string Mzml = @"E:\Projects\JurkatTopDown\02-18-20_jurkat_td_rep2_fract10.mzML";
        private const string GroundTruth = @"E:\Projects\JurkatTopDown\FlashDecon\02-18-20_jurkat_td_rep2_fract10_ms1.feature";
        private const string OutDir = @"E:\TestData\MetaMorpheus\runs\metaflashdecon_eval";
        private const string ResultsLog = @"E:\CodeReview\MetaFlashDecon\deliverables\wholefile_eval_results.txt";

        [Test] public void Gated()   => RunPass(0.5, writeFeatureFile: true);
        [Test] public void Ungated() => RunPass(0.0, writeFeatureFile: false);

        // Recall hunt (2026-05-30): faithful z1-100 default params, with the tracer's per-trace dump enabled
        // so we can split the missed FLASHDeconv features into "trace dropped by 0.85 cosine" vs "no trace".
        private const string TraceDumpPath = @"E:\CodeReview\MetaFlashDecon\deliverables\tracedump_z1_100.txt";
        [Test]
        public void DefaultParamsWithTraceDump()
        {
            Assume.That(File.Exists(Mzml), $"missing mzML: {Mzml}");
            Directory.CreateDirectory(OutDir);
            var sw = Stopwatch.StartNew();
            var dataFile = MsDataFileReader.GetDataFile(Mzml).LoadAllStaticData();
            var ms1 = dataFile.GetAllScansList().Where(s => s.MsnOrder == 1)
                .OrderBy(s => s.OneBasedScanNumber).ToList();
            var p = new MetaFlashDeconParameters(); // FLASHDeconv-faithful defaults: z1-100, cos 0.85
            var perScan = Ms1FeatureGenerator.DeconvolveScans(ms1, p);
            long totalEnv = perScan.Sum(e => (long)e.Count);
            Log($"[tracedump] deconvolved {ms1.Count} scans, envelopes={totalEnv} in {sw.Elapsed}");

            MetaFlashDeconMassFeatureTracer.TraceDump = new System.Collections.Generic.List<string>();
            var features = new MetaFlashDeconMassFeatureTracer().TraceFeatures(ms1, perScan);
            File.WriteAllText(TraceDumpPath,
                "centroid\trtMinSec\trtMaxSec\tsize\tintensity\tcosine\tkept\n" +
                string.Join("\n", MetaFlashDeconMassFeatureTracer.TraceDump));
            int kept = MetaFlashDeconMassFeatureTracer.TraceDump.Count(l => l.EndsWith("\t1"));
            int dropped = MetaFlashDeconMassFeatureTracer.TraceDump.Count - kept;
            MetaFlashDeconMassFeatureTracer.TraceDump = null;
            Log($"[tracedump] traces={MetaFlashDeconMassFeatureTracer.TraceDump?.Count ?? (kept + dropped)} " +
                $"kept(cosine>=0.85)={kept} dropped-by-cosine={dropped}  features={features.Count}  " +
                $"wrote {TraceDumpPath}  (total {sw.Elapsed})");
            string outPath = Path.Combine(OutDir, "rep2_fract10_native_z1_100_ms1.feature");
            Ms1FeatureGenerator.WriteFeatures(features, outPath);
        }

        private static void RunPass(double snrThreshold, bool writeFeatureFile)
        {
            Assume.That(File.Exists(Mzml), $"missing mzML: {Mzml}");
            Directory.CreateDirectory(OutDir);

            Log($"=== pass SNR>={snrThreshold:0.0}  ({DateTime.Now:HH:mm:ss}) ===");

            var sw = Stopwatch.StartNew();
            var dataFile = MsDataFileReader.GetDataFile(Mzml).LoadAllStaticData();
            var ms1 = dataFile.GetAllScansList()
                .Where(s => s.MsnOrder == 1)
                .OrderBy(s => s.OneBasedScanNumber)
                .ToList();
            int gtRows = File.ReadAllLines(GroundTruth).Length - 1;
            Log($"loaded {ms1.Count} MS1 scans in {sw.Elapsed}; ground-truth feature rows = {gtRows}");

            var p = new MetaFlashDeconParameters(minCharge: 1, maxCharge: 60, snrThreshold: snrThreshold);

            var decon = Stopwatch.StartNew();
            var perScan = Ms1FeatureGenerator.DeconvolveScans(ms1, p);
            long totalEnv = perScan.Sum(e => (long)e.Count);
            decon.Stop();
            Log($"deconvolved: per-scan envelopes={totalEnv} (~{totalEnv / (double)ms1.Count:F0}/scan) in {decon.Elapsed}");

            var trace = Stopwatch.StartNew();
            var features = new MetaFlashDeconMassFeatureTracer().TraceFeatures(ms1, perScan);
            trace.Stop();

            Log($"RESULT SNR>={snrThreshold:0.0}: WHOLE-FILE features={features.Count} " +
                $"(multi-z={features.Count(f => f.ChargeCount >= 2)})  " +
                $"mass[{features.Min(f => f.ConsensusMass):F0},{features.Max(f => f.ConsensusMass):F0}]  " +
                $"vs GT {gtRows}   (trace {trace.Elapsed}, total {sw.Elapsed})");

            if (writeFeatureFile)
            {
                string outPath = Path.Combine(OutDir, "rep2_fract10_native_snr05_ms1.feature");
                Ms1FeatureGenerator.WriteFeatures(features, outPath);
                Log($"wrote {outPath}");
            }
        }

        private static void Log(string line)
        {
            TestContext.Progress.WriteLine(line);
            File.AppendAllText(ResultsLog, line + Environment.NewLine);
        }
    }
}
