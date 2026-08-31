using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Readers;

namespace Development.Deconvolution
{
    /// <summary>
    /// Diagnostic tests for target-decoy classic deconvolution on real top-down data.
    ///
    /// Envelopes are scored using <see cref="DeconvolutionScorer"/> (the generic
    /// post-hoc scorer) rather than the raw <see cref="IsotopicEnvelope.Score"/>
    /// from the classic algorithm. The generic scorer uses cosine similarity against
    /// the real Averagine model, so decoy envelopes — whose peak positions are shifted
    /// by ~0.059 Da per isotope — should score much lower than targets even when the
    /// classic algorithm's own score doesn't separate them well.
    ///
    /// These tests are ignored by default. Remove [Ignore] or run explicitly from
    /// the CLI to see the score distribution output.
    /// </summary>
    [TestFixture]
    [Ignore("Diagnostic only — run explicitly to inspect decoy score distributions")]
    public class TestDecoyClassicRealData
    {
        private static readonly ClassicDeconvolutionParameters TopDownParams =
            new ClassicDeconvolutionParameters(1, 60, 4, 3);

        private static readonly AverageResidue ScoringModel = new Averagine();

        // ── Per-file diagnostics ──────────────────────────────────────────────

        [Test]
        public void Ubiquitin_TargetDecoy_ScoreSummary()
            => RunAndPrint("Averaged_221110_UbiqOnly.mzML", TopDownParams);

        [Test]
        public void HGH_TargetDecoy_ScoreSummary()
            => RunAndPrint("Averaged_221110_HGHOnly.mzML", TopDownParams);

        [Test]
        public void Cytochrome_TargetDecoy_ScoreSummary()
            => RunAndPrint("Averaged_221110_CytoOnly.mzML", TopDownParams);

        // ── Assertion-based tests ─────────────────────────────────────────────

        [Test]
        public void Ubiquitin_TargetCount_GreaterThan_DecoyCount()
        {
            string path = FindMzml("Averaged_221110_UbiqOnly.mzML");
            var (targets, decoys) = DeconvoluteFile(path, TopDownParams);

            Console.WriteLine($"Targets: {targets.Count}  Decoys: {decoys.Count}");

            Assert.That(targets.Count, Is.GreaterThan(0),
                "Expected target envelopes from ubiquitin file");
            Assert.That(targets.Count, Is.GreaterThan(decoys.Count),
                "Expected more targets than decoys on real protein data");
        }

        [Test]
        public void Ubiquitin_GenericScore_TargetMean_GreaterThan_DecoyMean()
        {
            string path = FindMzml("Averaged_221110_UbiqOnly.mzML");
            var (targets, decoys) = DeconvoluteFile(path, TopDownParams);

            Assume.That(targets.Count, Is.GreaterThan(0), "No targets found");
            Assume.That(decoys.Count, Is.GreaterThan(0), "No decoys found — cannot compare scores");

            double meanTarget = targets.Average(e => DeconvolutionScorer.ScoreEnvelope(e, ScoringModel));
            double meanDecoy = decoys.Average(e => DeconvolutionScorer.ScoreEnvelope(e, ScoringModel));

            // Diagnostic print only -- Classic AUC ~ 0.50 because Peak2satisfiesRatio
            // pre-filters everything to look Averagine-like, so target and decoy means
            // are not expected to differ. No assertion: the diagnostic value is the
            // printed scores for human inspection when this fixture is run explicitly.
            Console.WriteLine($"Generic score — mean target: {meanTarget:F4}  mean decoy: {meanDecoy:F4}");
        }

        // ── Core helpers ──────────────────────────────────────────────────────

        private static void RunAndPrint(string filename, ClassicDeconvolutionParameters p)
        {
            string path = FindMzml(filename);
            var (targets, decoys) = DeconvoluteFile(path, p);

            // Score all envelopes with the generic scorer
            List<double> targetScores = targets
                .Select(e => DeconvolutionScorer.ScoreEnvelope(e, ScoringModel)).ToList();
            List<double> decoyScores = decoys
                .Select(e => DeconvolutionScorer.ScoreEnvelope(e, ScoringModel)).ToList();

            PrintSummary(filename, targets.Count, decoys.Count, targetScores, decoyScores);
            PrintHistogram("TARGETS (generic score)", targetScores);
            PrintHistogram("DECOYS  (generic score)", decoyScores);

            // Also print the classic raw score for comparison
            List<double> targetRaw = targets.Select(e => e.Score).ToList();
            List<double> decoyRaw = decoys.Select(e => e.Score).ToList();
            PrintHistogram("TARGETS (classic raw score)", targetRaw);
            PrintHistogram("DECOYS  (classic raw score)", decoyRaw);
        }

        private static (List<IsotopicEnvelope> targets, List<IsotopicEnvelope> decoys)
            DeconvoluteFile(string mzmlPath, ClassicDeconvolutionParameters p)
        {
            var allTargets = new List<IsotopicEnvelope>();
            var allDecoys = new List<IsotopicEnvelope>();

            var file = MsDataFileReader.GetDataFile(mzmlPath);
            file.LoadAllStaticData();

            foreach (var scan in file.GetAllScansList()
                                     .Where(s => s.MsnOrder == 1))
            {
                var (t, d) = Deconvoluter.DeconvoluteWithDecoys(scan, p);
                allTargets.AddRange(t);
                allDecoys.AddRange(d);
            }

            return (allTargets, allDecoys);
        }

        private static void PrintSummary(string filename,
            int targetCount, int decoyCount,
            List<double> targetScores, List<double> decoyScores)
        {
            Console.WriteLine($"\n{new string('-', 60)}");
            Console.WriteLine($"FILE: {filename}");
            Console.WriteLine($"{new string('-', 60)}");
            Console.WriteLine($"  Targets : {targetCount,6}");
            Console.WriteLine($"  Decoys  : {decoyCount,6}");

            if (targetScores.Count > 0)
                Console.WriteLine(
                    $"  Target generic scores — " +
                    $"min:{targetScores.Min():F4}  " +
                    $"mean:{targetScores.Average():F4}  " +
                    $"max:{targetScores.Max():F4}");

            if (decoyScores.Count > 0)
                Console.WriteLine(
                    $"  Decoy  generic scores — " +
                    $"min:{decoyScores.Min():F4}  " +
                    $"mean:{decoyScores.Average():F4}  " +
                    $"max:{decoyScores.Max():F4}");

            double ratio = decoyCount > 0
                ? (double)targetCount / decoyCount
                : double.PositiveInfinity;
            Console.WriteLine($"  Target/decoy ratio: {ratio:F2}");

            // AUC-style summary: fraction of targets scoring above median decoy
            if (targetScores.Count > 0 && decoyScores.Count > 0)
            {
                double medianDecoy = Median(decoyScores);
                double fracAbove = targetScores.Count(s => s > medianDecoy)
                                   / (double)targetScores.Count;
                Console.WriteLine($"  Fraction of targets above median decoy score: {fracAbove:P1}");
            }
        }

        private static void PrintHistogram(string label, List<double> scores)
        {
            if (scores.Count == 0)
            {
                Console.WriteLine($"\n  {label}: (none)");
                return;
            }

            const int bins = 10;
            const int barWidth = 30;

            double min = scores.Min();
            double max = scores.Max();
            double range = max - min;
            if (range < 1e-10) range = 1;

            int[] counts = new int[bins];
            foreach (double s in scores)
            {
                int b = (int)Math.Min(bins - 1, Math.Floor((s - min) / range * bins));
                counts[b]++;
            }

            int maxCount = counts.Max();

            Console.WriteLine($"\n  {label}  (n={scores.Count}):");
            Console.WriteLine($"  {"Score range",-22}  {"Count",6}  Bar");
            Console.WriteLine($"  {new string('-', 22)}  {new string('-', 6)}  {new string('-', barWidth)}");

            for (int b = 0; b < bins; b++)
            {
                double lo = min + b * range / bins;
                double hi = min + (b + 1) * range / bins;
                int bar = maxCount > 0 ? counts[b] * barWidth / maxCount : 0;
                Console.WriteLine(
                    $"  [{lo,8:F4} - {hi,8:F4}]  {counts[b],6}  {new string('|', bar)}");
            }
        }

        private static double Median(List<double> values)
        {
            var sorted = values.OrderBy(x => x).ToList();
            int n = sorted.Count;
            return n % 2 == 0
                ? (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
                : sorted[n / 2];
        }

        private static string FindMzml(string filename)
        {
            var candidates = new[]
            {
                Path.Combine(TestContext.CurrentContext.TestDirectory,
                    "DeconvolutionDevelopment", "TestData", filename),
                Path.Combine(TestContext.CurrentContext.TestDirectory,
                    "TestData", filename),
                Path.Combine(TestContext.CurrentContext.TestDirectory, filename),
            };

            string? dir = TestContext.CurrentContext.TestDirectory;
            var extra = new List<string>();
            for (int i = 0; i < 6 && dir != null; i++)
            {
                extra.Add(Path.Combine(dir, "Test", "DeconvolutionDevelopment", "TestData", filename));
                extra.Add(Path.Combine(dir, "Development", "DeconvolutionDevelopment", "TestData", filename));
                dir = Path.GetDirectoryName(dir);
            }

            string? found = candidates.Concat(extra).FirstOrDefault(File.Exists);
            if (found != null) return found;

            Assert.Inconclusive(
                $"{filename} not found. Searched:\n" +
                string.Join("\n", candidates.Concat(extra)));
            return null!;
        }
    }
}