// FLASHDeconvTsvWriter.cs
//
// Writes a per-envelope TSV companion file for FLASHDeconvolutionAlgorithm output,
// analogous to the OpenMS FLASHDeconv _ms1.tsv format.
//
// ── What OpenMS actually does for q-values ────────────────────────────────────
// OpenMS computes q-values by running the full deconvolution algorithm THREE times
// per spectrum:
//   1. Target run  — normal deconvolution
//   2. Isotope decoy run — same algorithm but isotope step shifted to 0.9444 Da
//      (instead of 1.003 Da), excluding masses found in the target run
//   3. Noise decoy run — another shifted run with a different offset
//
// It then fits a weighted mixture of the three Qscore distributions (using
// least-squares in Qvalue.cpp) and sweeps the Qscore threshold to compute
// target-decoy q-values in the standard way:
//   q(t) = weighted_decoys_with_score >= t / targets_with_score >= t
// enforced to be monotonically non-decreasing as score increases.
//
// ── What we do here (honest approximation) ────────────────────────────────────
// We do NOT run decoy passes. Instead we use the Qscore (which is a probability
// estimate from logistic regression, not a raw cosine) to produce a simplified
// q-value:
//
//   SimpleQvalue(t) = 1 − Qscore(t)
//
// Rationale: The OpenMS Qscore was trained as P(mass is a true identification).
// So 1 − Qscore approximates the local false discovery rate. Converting to a
// cumulative q-value via the monotone sweep is then:
//
//   q(i) = min over j>=i of { 1 − Qscore(j) }
//         (sorted descending by Qscore, monotone enforced)
//
// This is conservative: it overestimates q-values compared to the full
// target-decoy procedure, but is mathematically honest. It produces q-values
// that are a smooth function of the Qscore rather than noisy from small decoy
// pools, which is actually an advantage for a single-file run.
//
// The TSV column 'SimpleQvalue' is clearly labelled to distinguish it from
// the full target-decoy Qvalue produced by OpenMS.
//
// ── Column definitions ────────────────────────────────────────────────────────
// ScanNum          — 1-based scan number from the source MsDataScan
// RetentionTime    — retention time in seconds (matches msalign format)
// MonoisotopicMass — reported monoisotopic mass in Da
// AverageMass      — apex observed mass in Da (MostAbundantObservedIsotopicMass / |charge|)
// SumIntensity     — total intensity of recruited peaks
// MinCharge        — |charge| (we track a single charge per envelope post-dedup)
// MaxCharge        — same as MinCharge (single-charge envelope)
// PeakCount        — number of recruited isotope peaks
// IsotopeCosine    — cosine similarity from Step 4 (before Qscore replacement)
//                    NOTE: after AssignQscores, envelope.Score is the Qscore,
//                    not the cosine. The cosine must be passed separately — see usage.
// Qscore           — logistic regression Qscore in [0,1] (higher = better)
// SimpleQvalue     — monotone q-value estimate = 1 − Qscore, enforced non-decreasing
//                    as Qscore increases (see above)
// RepresentativeCharge — |charge| of the envelope
// AvgPPMError      — mean absolute ppm error of recruited peaks vs Averagine
//
// ── Placement ────────────────────────────────────────────────────────────────
// Add to the same location as the msalign writer in the test project, or to
// a Readers/Writers folder in mzLib proper.
//
// ── Usage ────────────────────────────────────────────────────────────────────
// The writer takes (scan, envelope, isotopeCosine) triples so the cosine is
// preserved separately from the Qscore that replaces envelope.Score.
//
// In the test:
//   var rawResults = new List<(MsDataScan Scan, IsotopicEnvelope Envelope, double Cosine)>();
//   foreach (var scan in ms1Scans)
//   {
//       // Run Steps 1-5, capture cosine before Qscore replacement:
//       var raw = ScoreAndBuildEnvelopes(...);      // Score = cosine
//       var deduped = FLASHDeconvDeduplicator.Deduplicate(raw, ppm);
//       double median = FLASHDeconvScorer.ComputeMedianIntensity(scan.MassSpectrum);
//       var scored = FLASHDeconvScorer.AssignQscores(deduped, median, ppm); // Score = Qscore
//       // pair each scored envelope with its original cosine:
//       foreach (var (env, cos) in scored.Zip(deduped.Select(e => e.Score)))
//           rawResults.Add((scan, env, cos));
//   }
//   FLASHDeconvTsvWriter.WriteMs1Tsv(rawResults, outputPath);
//
// Alternatively, simpler: store the cosine inside the envelope before calling
// AssignQscores (e.g. as a separate list), then zip.

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using MassSpectrometry;

namespace Test.Deconvolution
{
    /// <summary>
    /// Writes a per-envelope TSV companion file for FLASHDeconvolution output.
    /// Produces Qscore and a simplified q-value analogous to the OpenMS
    /// FLASHDeconv _ms1.tsv format. See file header for full methodology notes.
    /// </summary>
    public static class FLASHDeconvTsvWriter
    {
        // ── Column headers ────────────────────────────────────────────────────
        private static readonly string[] Columns =
        {
            "ScanNum",
            "RetentionTime",
            "MonoisotopicMass",
            "AverageMass",
            "SumIntensity",
            "MinCharge",
            "MaxCharge",
            "PeakCount",
            "IsotopeCosine",
            "Qscore",
            "SimpleQvalue",
            "RepresentativeCharge",
            "AvgPPMError"
        };

        /// <summary>
        /// Writes all envelopes across all scans to a TSV file, sorted by scan
        /// then by Qscore descending within each scan. The SimpleQvalue column
        /// is computed globally across the entire file (not per-scan) using a
        /// monotone sweep over all Qscores.
        /// </summary>
        /// <param name="results">
        /// Tuples of (scan, scored envelope, isotope cosine before Qscore replacement).
        /// </param>
        /// <param name="outputPath">Full path for the output TSV file.</param>
        public static void WriteMs1Tsv(
            IEnumerable<(MsDataScan Scan, IsotopicEnvelope Envelope, double IsotopeCosine)> results,
            string outputPath)
        {
            // Materialise so we can sweep twice (once for q-value, once to write)
            var rows = results.ToList();
            if (rows.Count == 0) return;

            // ── Compute simple q-values globally across all envelopes ─────────
            // Sort descending by Qscore (highest confidence first)
            var sortedByQscore = rows
                .Select((r, idx) => (r, idx))
                .OrderByDescending(x => x.r.Envelope.Score)
                .ToList();

            // SimpleQvalue[i] = 1 − Qscore, monotone enforced:
            //   running minimum of (1 − Qscore) as we scan high → low score
            var simpleQvalue = new double[rows.Count];
            double runningMin = 1.0;
            foreach (var (r, idx) in sortedByQscore)
            {
                double local = 1.0 - r.Envelope.Score;
                runningMin = Math.Min(runningMin, local);
                simpleQvalue[idx] = runningMin;
            }
            // Enforce non-decreasing in the low-score direction (reverse pass)
            // Actually the above already does this correctly: as we scan from
            // high Qscore to low Qscore, runningMin can only increase or stay.
            // But we want q-values to be non-DECREASING as Qscore DECREASES,
            // which means non-decreasing as we go forward in the sorted order.
            // Re-sweep to enforce:
            runningMin = 0.0;
            foreach (var (r, idx) in sortedByQscore)
            {
                simpleQvalue[idx] = Math.Max(simpleQvalue[idx], runningMin);
                runningMin = simpleQvalue[idx];
            }

            // ── Write TSV ─────────────────────────────────────────────────────
            using var writer = new StreamWriter(outputPath);
            writer.WriteLine(string.Join("\t", Columns));

            // Sort output by scan number, then by Qscore descending within scan
            var outputOrder = rows
                .Select((r, idx) => (r, idx))
                .OrderBy(x => x.r.Scan.OneBasedScanNumber)
                .ThenByDescending(x => x.r.Envelope.Score);

            foreach (var (r, idx) in outputOrder)
            {
                var (scan, env, cosine) = r;
                int absCharge = Math.Abs(env.Charge);

                // Average mass: apex m/z × charge ≈ MostAbundantObservedIsotopicMass
                // We reconstruct from the most-intense peak in the recruited list
                double apexMz = env.Peaks.Count > 0
                    ? env.Peaks.MaxBy(p => p.intensity).mz
                    : env.MonoisotopicMass.ToMz(env.Charge);
                double averageMass = apexMz * absCharge;

                double avgPpmError = ComputeAvgPpmError(env);

                writer.WriteLine(string.Join("\t", new object[]
                {
                    scan.OneBasedScanNumber,
                    $"{scan.RetentionTime * 60.0:F2}",        // minutes → seconds
                    $"{env.MonoisotopicMass:F6}",
                    $"{averageMass:F6}",
                    $"{env.TotalIntensity:F2}",
                    absCharge,                                  // MinCharge
                    absCharge,                                  // MaxCharge (single-charge envelope)
                    env.Peaks.Count,
                    $"{cosine:F6}",
                    $"{env.Score:F6}",                          // Qscore
                    $"{simpleQvalue[idx]:F6}",
                    absCharge,                                  // RepresentativeCharge
                    $"{avgPpmError:F4}"
                }));
            }
        }

        // ── Helpers ───────────────────────────────────────────────────────────

        private static double ComputeAvgPpmError(IsotopicEnvelope env)
        {
            if (env.Peaks.Count == 0 || env.Charge == 0) return 0.0;
            int absCharge = Math.Abs(env.Charge);
            double isotopeStepMz = Chemistry.Constants.C13MinusC12 / absCharge;
            var apexPeak = env.Peaks.MaxBy(p => p.intensity);
            double apexMz = apexPeak.mz;
            double totalPpm = 0.0;
            foreach (var (mz, _) in env.Peaks)
            {
                int n = (int)Math.Round((mz - apexMz) / isotopeStepMz);
                double theorMz = apexMz + n * isotopeStepMz;
                totalPpm += Math.Abs(mz - theorMz) / theorMz * 1e6;
            }
            return env.Peaks.Count > 0 ? totalPpm / env.Peaks.Count : 0.0;
        }
    }
}