// FLASHDeconvDeduplication.cs
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
//                 errors. These are real measurement offsets that will be addressed
//                 separately.
//   Harmonic species (half-mass, third-mass): handled elsewhere by harmonic rejection
//                 in Step 2. This code only sees the Step 5 output.
//   Distinct co-eluting proteins: their masses differ by hundreds of ppm or more.
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
//      cluster's current representative mass exceeds the threshold.
//   3. Within each cluster, keep the envelope with the highest Score (cosine).
//      Ties broken by TotalIntensity.
//
// Greedy (not transitive-closure) clustering is intentional: it means a cluster
// cannot silently expand to bridge two species that happen to have an intermediate
// envelope between them. Each cluster is bounded by consecutive ppm distance, not
// all-pairs distance.
//
// PLACEMENT
// ---------
// Add this static class to MassSpectrometry/Deconvolution/Algorithms/ alongside
// FLASHDeconvolutionAlgorithm.cs. Call it from the tail of Deconvolute() before
// returning, after the existing Step 5 output is assembled.
//
// USAGE IN FLASHDeconvolutionAlgorithm.Deconvolute()
// ---------------------------------------------------
//   // ... existing Steps 1-5 ...
//   var rawEnvelopes = ScoreAndBuildEnvelopes(...);
//   return FLASHDeconvDeduplicator.Deduplicate(rawEnvelopes, p.DeconvolutionTolerancePpm);
//
// The same decon tolerance used for peak matching (default 10 ppm) is a natural
// choice for the deduplication window as well, because a duplicate produced by a
// 1-isotope apex shift will be at most ~(1 Da / mass) ppm — always well under 10
// ppm for proteins above 2 kDa. For smaller peptides (mass ~ 1 kDa), one isotope
// step is ~1 Da / 1000 Da × 1e6 = 1000 ppm, so those will NOT be collapsed at 10
// ppm — which is correct, because at low mass each isotope peak is a genuinely
// distinct deconvolution result.

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
        /// <summary>
        /// Deduplicates a collection of isotopic envelopes from a single spectrum.
        /// Within each ppm-tolerance cluster the envelope with the highest cosine
        /// score is kept; ties are broken by TotalIntensity.
        /// </summary>
        /// <param name="envelopes">
        /// Raw envelopes from Steps 1–5 of <see cref="FLASHDeconvolutionAlgorithm"/>.
        /// </param>
        /// <param name="tolerancePpm">
        /// Maximum ppm difference between two envelopes to be considered the same
        /// species. Default (and recommended) value: the same ppm tolerance used
        /// for peak matching — typically 10 ppm.
        /// </param>
        /// <returns>
        /// One envelope per distinct mass cluster, being the highest-scoring member.
        /// </returns>
        internal static IEnumerable<IsotopicEnvelope> Deduplicate(
            IEnumerable<IsotopicEnvelope> envelopes,
            double tolerancePpm = 10.0)
        {
            // Sort by monoisotopic mass so we can do a single forward pass
            var sorted = envelopes
                .OrderBy(e => e.MonoisotopicMass)
                .ToList();

            if (sorted.Count == 0)
                yield break;

            // Greedy single-pass clustering
            // clusterRep tracks the mass of the FIRST envelope added to the
            // current cluster, which is always the lowest mass in the cluster.
            // All subsequent comparisons are made against this anchor, not
            // against the previously seen envelope. This prevents the cluster
            // from drifting across a large mass range through a chain of small
            // ppm steps.
            IsotopicEnvelope best = sorted[0];
            double clusterAnchorMass = sorted[0].MonoisotopicMass;

            for (int i = 1; i < sorted.Count; i++)
            {
                IsotopicEnvelope current = sorted[i];
                double ppmFromAnchor = PpmDiff(current.MonoisotopicMass, clusterAnchorMass);

                if (ppmFromAnchor <= tolerancePpm)
                {
                    // Same cluster — keep the better-scoring envelope
                    if (IsBetter(current, best))
                        best = current;
                }
                else
                {
                    // New cluster — emit the winner of the previous cluster
                    yield return best;

                    // Start fresh cluster anchored at the current envelope's mass
                    best = current;
                    clusterAnchorMass = current.MonoisotopicMass;
                }
            }

            // Emit the final cluster's winner
            yield return best;
        }

        // ── Helpers ─────────────────────────────────────────────────────────────

        /// <summary>
        /// Returns true if <paramref name="candidate"/> should replace
        /// <paramref name="incumbent"/> as the cluster representative.
        /// Primary criterion: higher cosine score.
        /// Tiebreaker: higher total intensity.
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