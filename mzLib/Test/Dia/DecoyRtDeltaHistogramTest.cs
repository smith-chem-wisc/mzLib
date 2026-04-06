// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using MassSpectrometry.Dia;
using NUnit.Framework;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

namespace Test.DiaTests
{
    /// <summary>
    /// Diagnostic test: generates decoys from the real HeLa DIA-NN library and
    /// produces an HTML histogram of the RT delta distribution
    /// (decoyRT - targetRT after iRT correction).
    ///
    /// Output: F:\DiaBenchmark\NewLibrary\decoy_rt_delta_histogram.html
    ///
    /// Not run in CI — requires the real library file on disk.
    /// Placement: mzLib/Test/Dia/DecoyRtDeltaHistogramTest.cs
    /// </summary>
    [TestFixture]
    [Explicit("Diagnostic only — requires real library file at F:\\DiaBenchmark\\NewLibrary\\diannFig2helalib.msp")]
    public class DecoyRtDeltaHistogramTest
    {
        private const string LibraryPath = @"F:\DiaBenchmark\NewLibrary\diannFig2helalib.msp";
        private const string OutputPath = @"F:\DiaBenchmark\NewLibrary\decoy_rt_delta_histogram.html";

        [Test]
        public async Task GenerateDecoyRtDeltaHistogram()
        {
            // ── 1. Check library exists ───────────────────────────────────────
            Assert.That(File.Exists(LibraryPath), Is.True,
                $"Library not found: {LibraryPath}");

            // ── 2. Parse targets ──────────────────────────────────────────────
            Console.WriteLine($"Parsing library: {LibraryPath}");
            var targets = KoinaMspParser.Parse(LibraryPath, retentionTimeLookup: null, minIntensity: 0.05f);
            Console.WriteLine($"Parsed {targets.Count:N0} targets");
            Assert.That(targets.Count, Is.GreaterThan(0), "Library must contain at least one entry");

            // ── 3. Generate decoys ────────────────────────────────────────────
            Console.WriteLine("Generating decoys...");
            var decoys = DecoyGenerator.GenerateFromTargets(targets);
            Assert.That(decoys.Count, Is.EqualTo(targets.Count));

            // ── 5. Apply iRT correction ───────────────────────────────────────
            Console.WriteLine("Applying iRT RT correction via Prosit...");
            try
            {
                await ApplyIrtRtCorrection(targets, decoys);
            }
            catch (Exception ex)
            {
                Assert.Inconclusive($"Prosit iRT correction failed — cannot generate post-correction histogram. ({ex.Message})");
                return;
            }

            // ── 6. Capture post-correction RT deltas ──────────────────────────
            var deltasAfter = new List<double>(targets.Count);
            for (int i = 0; i < targets.Count; i++)
            {
                double? tRt = targets[i].RetentionTime ?? targets[i].IrtValue;
                double? dRt = decoys[i].RetentionTime ?? decoys[i].IrtValue;
                if (tRt.HasValue && dRt.HasValue)
                    deltasAfter.Add(dRt.Value - tRt.Value);
            }

            // ── 7. Print summary statistics ───────────────────────────────────
            PrintStats("After iRT correction", deltasAfter);

            // ── 8. Write HTML histogram ───────────────────────────────────────
            string html = BuildHistogramHtml(deltasAfter);
            File.WriteAllText(OutputPath, html);
            Console.WriteLine($"Histogram written to: {OutputPath}");
        }

        // ── iRT correction (mirrors DiaSearchRunner.ApplyIrtRtCorrection) ─────

        private static async Task ApplyIrtRtCorrection(
            IReadOnlyList<LibraryPrecursorInput> targets,
            List<LibraryPrecursorInput> decoys)
        {
            int n = targets.Count;
            var allSequences = new List<string>(n * 2);
            for (int i = 0; i < n; i++)
                allSequences.Add(StripMods(targets[i].Sequence));
            for (int i = 0; i < n; i++)
                allSequences.Add(StripMods(decoys[i].Sequence));

            WarningException warnings = null;
            var model = new Prosit2019iRT(allSequences, out warnings);
            if (warnings != null)
                Console.WriteLine($"  [iRT] {warnings.Message.Split('\n')[0]}");
            await model.RunInferenceAsync();

            var predictions = model.Predictions;
            if (predictions?.Count != n * 2)
                throw new InvalidOperationException(
                    $"Expected {n * 2} predictions, got {predictions?.Count ?? 0}");

            int corrected = 0;
            for (int i = 0; i < n; i++)
            {
                double? targetRt = targets[i].RetentionTime ?? targets[i].IrtValue;
                if (targetRt == null) continue;

                double delta = predictions[i + n].PredictedRetentionTime
                             - predictions[i].PredictedRetentionTime;

                var old = decoys[i];
                decoys[i] = new LibraryPrecursorInput(
                    old.Sequence, old.PrecursorMz, old.ChargeState,
                    retentionTime: targetRt.Value + delta,
                    isDecoy: true,
                    old.FragmentMzs, old.FragmentIntensities,
                    old.IrtValue);
                corrected++;
            }
            Console.WriteLine($"  [iRT] RT correction applied: {corrected:N0} decoys updated.");
        }

        // ── Statistics ────────────────────────────────────────────────────────

        private static void PrintStats(string label, List<double> deltas)
        {
            if (deltas.Count == 0) { Console.WriteLine($"{label}: no data"); return; }
            var sorted = deltas.OrderBy(x => x).ToList();
            double mean = deltas.Average();
            double median = sorted[sorted.Count / 2];
            double std = Math.Sqrt(deltas.Average(x => (x - mean) * (x - mean)));
            double p5 = sorted[(int)(sorted.Count * 0.05)];
            double p95 = sorted[(int)(sorted.Count * 0.95)];
            double absMax = sorted.Select(Math.Abs).Max();
            Console.WriteLine($"\n{label} (n={deltas.Count:N0}):");
            Console.WriteLine($"  mean={mean:F3}  median={median:F3}  std={std:F3}");
            Console.WriteLine($"  P5={p5:F3}  P95={p95:F3}  |max|={absMax:F3}");
            Console.WriteLine($"  |delta|<1: {deltas.Count(d => Math.Abs(d) < 1):N0}  " +
                              $"|delta|<5: {deltas.Count(d => Math.Abs(d) < 5):N0}  " +
                              $"|delta|>=5: {deltas.Count(d => Math.Abs(d) >= 5):N0}");
        }

        // ── HTML histogram builder ────────────────────────────────────────────

        private static string BuildHistogramHtml(List<double> after)
        {
            // Compute histogram bins for the after-correction distribution
            // Use range covering 99th percentile on each side
            var sortedAfter = after.OrderBy(x => x).ToList();
            double lo = sortedAfter[(int)(sortedAfter.Count * 0.005)];
            double hi = sortedAfter[(int)(sortedAfter.Count * 0.995)];
            // Round to nice numbers
            lo = Math.Floor(lo);
            hi = Math.Ceiling(hi);

            int nBins = 80;
            double binWidth = (hi - lo) / nBins;
            if (binWidth <= 0) binWidth = 1;

            int[] binsAfter = ComputeBins(after, lo, hi, nBins, binWidth);

            // Build bin labels
            var labels = new List<string>(nBins);
            var dataAfter = new List<string>(nBins);
            for (int b = 0; b < nBins; b++)
            {
                double center = lo + (b + 0.5) * binWidth;
                labels.Add($"{center:F1}");
                dataAfter.Add(binsAfter[b].ToString());
            }

            // Summary stats for display
            double meanAfter = after.Count > 0 ? after.Average() : 0;
            double stdAfter = after.Count > 0
                ? Math.Sqrt(after.Average(x => (x - meanAfter) * (x - meanAfter))) : 0;
            var sortedAbs = after.Select(Math.Abs).OrderBy(x => x).ToList();
            double p95abs = sortedAbs.Count > 0
                ? sortedAbs[(int)(sortedAbs.Count * 0.95)] : 0;

            return $@"<!DOCTYPE html>
<html lang=""en"">
<head>
<meta charset=""UTF-8"">
<title>Decoy RT Delta Distribution</title>
<script src=""https://cdnjs.cloudflare.com/ajax/libs/Chart.js/4.4.1/chart.umd.min.js""></script>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
          background: #1a1a2e; color: #e0e0e0; margin: 0; padding: 24px; }}
  h1   {{ color: #a78bfa; margin-bottom: 4px; }}
  .subtitle {{ color: #94a3b8; margin-bottom: 24px; font-size: 0.9em; }}
  .stats {{ display: flex; gap: 16px; flex-wrap: wrap; margin-bottom: 24px; }}
  .stat-card {{ background: #16213e; border: 1px solid #2d3748; border-radius: 8px;
                padding: 12px 20px; min-width: 140px; }}
  .stat-label {{ font-size: 0.75em; color: #94a3b8; text-transform: uppercase;
                 letter-spacing: 0.05em; }}
  .stat-value {{ font-size: 1.4em; font-weight: 600; color: #a78bfa; margin-top: 2px; }}
  .chart-container {{ background: #16213e; border: 1px solid #2d3748;
                      border-radius: 12px; padding: 24px; max-width: 1100px; }}
</style>
</head>
<body>
<h1>Decoy RT Delta Distribution</h1>
<p class=""subtitle"">
  decoyRT − targetRT after Prosit iRT correction &nbsp;|&nbsp;
  n = {after.Count:N0} pairs &nbsp;|&nbsp;
  Library: diannFig2helalib.msp
</p>

<div class=""stats"">
  <div class=""stat-card"">
    <div class=""stat-label"">Mean Δ</div>
    <div class=""stat-value"">{meanAfter:+0.00;-0.00} iRT</div>
  </div>
  <div class=""stat-card"">
    <div class=""stat-label"">Std Dev</div>
    <div class=""stat-value"">{stdAfter:F2} iRT</div>
  </div>
  <div class=""stat-card"">
    <div class=""stat-label"">|Δ| P95</div>
    <div class=""stat-value"">{p95abs:F2} iRT</div>
  </div>
  <div class=""stat-card"">
    <div class=""stat-label"">|Δ| &lt; 1</div>
    <div class=""stat-value"">{after.Count(d => Math.Abs(d) < 1):N0}</div>
  </div>
  <div class=""stat-card"">
    <div class=""stat-label"">|Δ| ≥ 5</div>
    <div class=""stat-value"">{after.Count(d => Math.Abs(d) >= 5):N0}</div>
  </div>
  <div class=""stat-card"">
    <div class=""stat-label"">|Δ| ≥ 10</div>
    <div class=""stat-value"">{after.Count(d => Math.Abs(d) >= 10):N0}</div>
  </div>
</div>

<div class=""chart-container"">
  <canvas id=""hist"" height=""90""></canvas>
</div>

<script>
const labels = [{string.Join(",", labels.Select(l => $"\"{l}\""))}];
const after  = [{string.Join(",", dataAfter)}];

new Chart(document.getElementById('hist'), {{
  type: 'bar',
  data: {{
    labels: labels,
    datasets: [
      {{
        label: 'After iRT correction',
        data: after,
        backgroundColor: 'rgba(167,139,250,0.75)',
        borderColor:     'rgba(167,139,250,1)',
        borderWidth: 1
      }}
    ]
  }},
  options: {{
    responsive: true,
    plugins: {{
      legend: {{ labels: {{ color: '#e0e0e0' }} }},
      title: {{
        display: true,
        text: 'Decoy RT Delta (decoyRT − targetRT) after Prosit iRT correction',
        color: '#e0e0e0',
        font: {{ size: 14 }}
      }},
      tooltip: {{
        callbacks: {{
          title: ctx => `Δ ≈ ${{ctx[0].label}} iRT`,
          label: ctx => ` ${{ctx.dataset.label}}: ${{ctx.parsed.y.toLocaleString()}} pairs`
        }}
      }}
    }},
    scales: {{
      x: {{
        ticks: {{ color: '#94a3b8', maxTicksLimit: 20, maxRotation: 45 }},
        grid:  {{ color: '#2d3748' }},
        title: {{ display: true, text: 'RT delta (iRT units)', color: '#94a3b8' }}
      }},
      y: {{
        ticks: {{ color: '#94a3b8' }},
        grid:  {{ color: '#2d3748' }},
        title: {{ display: true, text: 'Precursor count', color: '#94a3b8' }}
      }}
    }},
    barPercentage: 1.0,
    categoryPercentage: 1.0
  }}
}});
</script>
</body>
</html>";
        }

        private static int[] ComputeBins(
            List<double> values, double lo, double hi, int nBins, double binWidth)
        {
            var bins = new int[nBins];
            foreach (double v in values)
            {
                if (v < lo || v > hi) continue;
                int b = (int)((v - lo) / binWidth);
                if (b >= nBins) b = nBins - 1;
                bins[b]++;
            }
            return bins;
        }

        private static string StripMods(string sequence)
        {
            if (string.IsNullOrEmpty(sequence)) return "";
            var sb = new System.Text.StringBuilder(sequence.Length);
            bool inMod = false; char modOpen = ' ';
            foreach (char c in sequence)
            {
                if (!inMod) { if (c == '(' || c == '[') { inMod = true; modOpen = c; } else sb.Append(c); }
                else { char close = modOpen == '(' ? ')' : ']'; if (c == close) inMod = false; }
            }
            return sb.ToString();
        }
    }
}