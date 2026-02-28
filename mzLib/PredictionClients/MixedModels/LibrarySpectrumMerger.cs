using Omics.SpectrumMatch;
using System.ComponentModel;

namespace PredictionClients.MixedModels
{
    /// <summary>
    /// Stateless, pure merge logic.
    ///
    /// Takes a list of MixedModelResult objects (one per component that ran) and folds
    /// them into a single LibrarySpectrum per peptide. The merger is the only place
    /// in the codebase that understands the rules for combining heterogeneous model outputs.
    ///
    /// MERGE RULES (current)
    /// ----------------------
    /// Fragment ions:
    ///   - All fragment ions from all contributing spectra are unioned into one list.
    ///   - Duplicate ions (same m/z within 0.01 Da AND same annotation) are resolved by
    ///     keeping the one with higher intensity. This should not happen in normal use
    ///     (primary and internal ions occupy different m/z ranges) but is handled defensively.
    ///   - IsInternalFragment on each ion is determined by the ion itself (Product.SecondaryProductType),
    ///     not by which component contributed it.
    ///
    /// Precursor m/z:
    ///   - Taken from the PrimaryFragmentIntensities result if present; otherwise from the
    ///     first available result. The primary model (Prosit) has the most accurate precursor
    ///     m/z since it accounts for modifications precisely.
    ///
    /// Retention time:
    ///   - Taken from the RetentionTime component if present (priority).
    ///   - Falls back to the value in the PrimaryFragmentIntensities spectrum.
    ///   - If neither provides it, remains null.
    ///
    /// ADDING A NEW CONTRIBUTION TYPE
    /// --------------------------------
    /// Add a case to the switch in Merge(). The rest of the infrastructure is unchanged.
    ///
    /// THREAD SAFETY
    /// -------------
    /// All methods are static and stateless. Safe to call from any thread.
    /// </summary>
    public static class LibrarySpectrumMerger
    {
        private const double MzMatchToleranceDa = 0.01;

        /// <summary>
        /// Merges all results into a dictionary of combined LibrarySpectrum objects,
        /// keyed by "Sequence/Charge".
        ///
        /// Results from failed components (Succeeded = false) are skipped with a warning.
        /// Peptides present in some but not all components are still included — they receive
        /// whatever contributions are available.
        /// </summary>
        /// <param name="results">Results from all components that ran.</param>
        /// <param name="warnings">
        /// Collated warnings: per-component warnings + any merge-time conflicts.
        /// Null if everything ran cleanly.
        /// </param>
        public static Dictionary<string, LibrarySpectrum> Merge(
            IReadOnlyList<MixedModelResult> results,
            out WarningException? warnings)
        {
            var warningMessages = new List<string>();

            // Separate failed results and warn about them
            foreach (var failed in results.Where(r => !r.Succeeded))
            {
                warningMessages.Add(
                    $"Component '{failed.ComponentName}' failed and its contributions " +
                    $"were excluded: {failed.Error?.Message ?? "unknown error"}");
            }

            // Collect per-component warnings
            foreach (var w in results.Where(r => r.Succeeded && r.Warning != null))
            {
                warningMessages.Add($"[{w.ComponentName}] {w.Warning!.Message}");
            }

            var succeeded = results.Where(r => r.Succeeded).ToList();

            // Collect all peptide keys across all components
            var allKeys = succeeded
                .SelectMany(r => r.Spectra.Keys.Concat(r.ScalarData.Keys))
                .Distinct()
                .ToHashSet();

            // Separate results by contribution type for clean dispatch
            var primaryResults = succeeded.Where(r => r.ContributionType == ContributionType.PrimaryFragmentIntensities).ToList();
            var internalResults = succeeded.Where(r => r.ContributionType == ContributionType.InternalFragmentIntensities).ToList();
            var rtResults = succeeded.Where(r => r.ContributionType == ContributionType.RetentionTime).ToList();
            // Future: ionMobilityResults, etc.

            var merged = new Dictionary<string, LibrarySpectrum>();

            foreach (var key in allKeys)
            {
                // ── Fragment ions ──────────────────────────────────────────────────────
                var allFragments = new List<Omics.Fragmentation.MatchedFragmentIon>();

                // Primary ions first (defines the anchor spectrum: precursorMz, sequence, charge)
                LibrarySpectrum? anchorSpectrum = null;
                foreach (var r in primaryResults)
                {
                    if (r.Spectra.TryGetValue(key, out var s))
                    {
                        anchorSpectrum = s;
                        allFragments.AddRange(s.MatchedFragmentIons);
                        break; // Only one primary source expected; take the first
                    }
                }

                // Internal fragment ions
                foreach (var r in internalResults)
                {
                    if (r.Spectra.TryGetValue(key, out var s))
                    {
                        // If no primary anchor exists, use the internal spectrum as anchor
                        anchorSpectrum ??= s;
                        allFragments.AddRange(s.MatchedFragmentIons);
                    }
                }

                // Future contribution types: add cases here
                // e.g. foreach (var r in ionMobilityResults) { ... }

                if (anchorSpectrum == null)
                {
                    // This key came only from scalar data (e.g. RT-only component) — skip
                    // until we have fragment ions to anchor the spectrum.
                    warningMessages.Add(
                        $"Key '{key}' had no fragment ion contributions and was excluded " +
                        $"from the merged library.");
                    continue;
                }

                // ── Retention time ─────────────────────────────────────────────────────
                // Priority: dedicated RT component > primary spectrum RT > null
                double? rt = anchorSpectrum.RetentionTime;
                foreach (var r in rtResults)
                {
                    if (r.ScalarData.TryGetValue(key, out var predictedRt))
                    {
                        rt = predictedRt;
                        break; // First RT component wins
                    }
                }

                // ── Deduplicate fragment ions ──────────────────────────────────────────
                // In practice primary and internal ions don't collide, but be defensive.
                var deduped = DeduplicateFragments(allFragments, out var dupCount);
                if (dupCount > 0)
                {
                    warningMessages.Add(
                        $"'{key}': {dupCount} duplicate fragment ion(s) resolved by keeping higher intensity.");
                }

                // ── Build merged spectrum ──────────────────────────────────────────────
                merged[key] = new LibrarySpectrum(
                    sequence: anchorSpectrum.Sequence,
                    precursorMz: anchorSpectrum.PrecursorMz,
                    chargeState: anchorSpectrum.ChargeState,
                    peaks: deduped,
                    rt: rt);
            }

            warnings = warningMessages.Count > 0
                ? new WarningException(string.Join("\n", warningMessages))
                : null;

            return merged;
        }

        /// <summary>
        /// Removes fragment ions that share the same m/z (within MzMatchToleranceDa) AND
        /// the same annotation, keeping the one with higher intensity.
        /// </summary>
        private static List<Omics.Fragmentation.MatchedFragmentIon> DeduplicateFragments(
            List<Omics.Fragmentation.MatchedFragmentIon> ions,
            out int removedCount)
        {
            var result = new List<Omics.Fragmentation.MatchedFragmentIon>(ions.Count);
            var seen = new List<Omics.Fragmentation.MatchedFragmentIon>(ions.Count);
            removedCount = 0;

            foreach (var ion in ions)
            {
                var duplicate = seen.FirstOrDefault(s =>
                    Math.Abs(s.Mz - ion.Mz) < MzMatchToleranceDa &&
                    s.NeutralTheoreticalProduct.Annotation == ion.NeutralTheoreticalProduct.Annotation);

                if (duplicate == null)
                {
                    seen.Add(ion);
                    result.Add(ion);
                }
                else if (ion.Intensity > duplicate.Intensity)
                {
                    // Replace with the higher-intensity version
                    var idx = result.IndexOf(duplicate);
                    result[idx] = ion;
                    seen[seen.IndexOf(duplicate)] = ion;
                    removedCount++;
                }
                else
                {
                    removedCount++;
                }
            }

            return result;
        }
    }
}
