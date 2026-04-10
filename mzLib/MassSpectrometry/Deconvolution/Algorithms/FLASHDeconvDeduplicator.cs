// FLASHDeconvDeduplicator.cs
//
// Deduplicates IsotopicEnvelope results from FLASHDeconvolutionAlgorithm
// within a single spectrum.
//
// MOTIVATION
// ----------
// FLASHDeconv Step 2 generates candidate masses independently for each charge
// state. The same physical protein can produce multiple IsotopicEnvelope objects
// differing only by sub-ppm floating-point jitter (same isotope, same charge,
// slightly different apex selection) or by a small apex shift caused by recruiting
// from a neighbouring isotope peak at the same charge. These are not biologically
// distinct species.
//
// WHAT THIS COLLAPSES (safe)
// --------------------------
//   Population 1 — 0–1 ppm:  identical isotope peak, tiny float variation.
//                             Example: 18786.6084 vs 18786.6099  (0.08 ppm)
//   Population 2 — 1–10 ppm: one-isotope-step apex offset at the given charge.
//                             Example: 7846.0265 vs 7846.0433 at z=6  (2.1 ppm,
//                             = 1 Da / 6 / 7846 × 1e6)
//
// WHAT THIS DOES NOT COLLAPSE (by design)
// ----------------------------------------
//   ~50 ppm gaps: one full C13 step at the monoisotopic mass — missed-monoisotopic
//                 errors. These are real measurement offsets addressed separately.
//   Harmonic species (half-mass, third-mass): handled in Step 2 harmonic rejection.
//   Distinct co-eluting proteins: masses differ by hundreds of ppm or more.
//
// EMPIRICAL BASIS
// ---------------
// The 10 ppm threshold was chosen from analysis of the Filgrastim dataset
// (190226_FIlg_2_FD_500ng.mzML). The within-scan ppm distribution shows:
//   - Two dense populations at 0–1 ppm and 3–4 ppm (true duplicates)
//   - A near-empty gap from 11–49 ppm (0–4 pairs across the entire run)
//   - A third population at ~50 ppm (missed-monoisotopic, intentionally preserved)
// This makes 10 ppm a conservative, data-driven choice.
//
// ALGORITHM
// ---------
// Single-pass greedy clustering on mass-sorted envelopes:
//   1. Sort envelopes by MonoisotopicMass ascending.
//   2. Walk forward; start a new cluster whenever the ppm distance from the
//      cluster's anchor mass exceeds the threshold.
//   3. Within each cluster, keep the envelope with the highest cosine score.
//      Ties broken by TotalIntensity.
//
// Greedy (not transitive-closure) clustering is intentional: it prevents a
// cluster from silently bridging two distinct species through an intermediate.
//
// TWO OVERLOADS
// -------------
// Deduplicate(IEnumerable<EnvelopeScoringData>, tolerancePpm)
//   Primary path. Operates on the full scoring data produced during Step 3 so
//   that per-charge cosine, noise power, and ppm error for the cluster winner
//   are correctly propagated to FLASHDeconvScorer.AssignQscores.
//
// Deduplicate(IEnumerable<IsotopicEnvelope>, tolerancePpm)
//   Legacy / fallback path for callers that do not have EnvelopeScoringData.

using System;
using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry
{
    /// <summary>
    /// Collapses duplicate <see cref="IsotopicEnvelope"/> objects produced by
    /// <see cref="FLASHDeconvolutionAlgorithm"/> for the same physical species.
    /// </summary>
    internal static class FLASHDeconvDeduplicator
    {
        // ── Primary overload: full EnvelopeScoringData ────────────────────────

        /// <summary>
        /// Deduplicates a collection of <see cref="FLASHDeconvScorer.EnvelopeScoringData"/>
        /// from a single spectrum. Within each ppm-tolerance cluster the entry with
        /// the highest global cosine score is kept; ties are broken by TotalIntensity.
        /// The winning entry's full scoring data is preserved so that
        /// <see cref="FLASHDeconvScorer.AssignQscores"/> receives accurate per-charge
        /// cosine and noise values.
        /// </summary>
        /// <param name="scoringData">
        /// Scoring data from Steps 1–5 of <see cref="FLASHDeconvolutionAlgorithm"/>.
        /// </param>
        /// <param name="tolerancePpm">
        /// Maximum ppm difference between two envelopes to be considered the same
        /// species. Recommended: the same ppm tolerance used for peak matching (10 ppm).
        /// </param>
        internal static IEnumerable<FLASHDeconvScorer.EnvelopeScoringData> Deduplicate(
            IEnumerable<FLASHDeconvScorer.EnvelopeScoringData> scoringData,
            double tolerancePpm = 10.0)
        {
            var sorted = scoringData
                .OrderBy(d => d.Envelope.MonoisotopicMass)
                .ToList();

            if (sorted.Count == 0)
                yield break;

            FLASHDeconvScorer.EnvelopeScoringData best = sorted[0];
            double clusterAnchorMass = sorted[0].Envelope.MonoisotopicMass;

            for (int i = 1; i < sorted.Count; i++)
            {
                var current = sorted[i];
                double ppmFromAnchor = PpmDiff(
                    current.Envelope.MonoisotopicMass, clusterAnchorMass);

                if (ppmFromAnchor <= tolerancePpm)
                {
                    if (IsBetter(current.Envelope, best.Envelope))
                        best = current;
                }
                else
                {
                    yield return best;
                    best = current;
                    clusterAnchorMass = current.Envelope.MonoisotopicMass;
                }
            }

            yield return best;
        }

        // ── Legacy overload: envelope-only ───────────────────────────────────

        /// <summary>
        /// Deduplicates a collection of isotopic envelopes from a single spectrum.
        /// Within each ppm-tolerance cluster the envelope with the highest cosine
        /// score is kept; ties are broken by TotalIntensity.
        /// <para>
        /// Use this overload only when <see cref="FLASHDeconvScorer.EnvelopeScoringData"/>
        /// is not available (e.g. when rescoring envelopes read from a file).
        /// </para>
        /// </summary>
        internal static IEnumerable<IsotopicEnvelope> Deduplicate(
            IEnumerable<IsotopicEnvelope> envelopes,
            double tolerancePpm = 10.0)
        {
            var sorted = envelopes
                .OrderBy(e => e.MonoisotopicMass)
                .ToList();

            if (sorted.Count == 0)
                yield break;

            IsotopicEnvelope best = sorted[0];
            double clusterAnchorMass = sorted[0].MonoisotopicMass;

            for (int i = 1; i < sorted.Count; i++)
            {
                IsotopicEnvelope current = sorted[i];
                double ppmFromAnchor = PpmDiff(current.MonoisotopicMass, clusterAnchorMass);

                if (ppmFromAnchor <= tolerancePpm)
                {
                    if (IsBetter(current, best))
                        best = current;
                }
                else
                {
                    yield return best;
                    best = current;
                    clusterAnchorMass = current.MonoisotopicMass;
                }
            }

            yield return best;
        }

        // ── Helpers ──────────────────────────────────────────────────────────

        /// <summary>
        /// Returns true if <paramref name="candidate"/> should replace
        /// <paramref name="incumbent"/> as the cluster representative.
        /// Primary criterion: higher cosine score. Tiebreaker: higher intensity.
        /// </summary>
        private static bool IsBetter(IsotopicEnvelope candidate, IsotopicEnvelope incumbent)
        {
            if (candidate.Score > incumbent.Score) return true;
            if (candidate.Score < incumbent.Score) return false;
            return candidate.TotalIntensity > incumbent.TotalIntensity;
        }

        /// <summary>
        /// Absolute ppm difference between two masses.
        /// Uses the average of the two masses as the denominator (symmetric).
        /// </summary>
        private static double PpmDiff(double massA, double massB)
        {
            double avg = (massA + massB) * 0.5;
            return avg > 0.0
                ? Math.Abs(massA - massB) / avg * 1e6
                : 0.0;
        }
    }
}