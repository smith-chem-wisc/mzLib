// MetaFlashDeconAveragine.cs
//
// Plan B: faithful port of OpenMS PrecalculatedAveragine's per-mass isotope-distribution
// post-processing (FLASHDeconvHelperStructs.cpp:36-106) — trim the isotope pattern to cover 99.99%
// of its power, L2-normalise it (so FLASHDeconvAlgorithm::getCosine, which divides only by |a|, is a
// true cosine), and record apex / left-count-from-apex / last index (used as recruitment scan bounds
// and the per-charge cosine target b).
//
// The RAW isotope distribution generator (OpenMS CoarseIsotopePatternGenerator vs our
// IsotopicDistribution) is the one differential-testing boundary — it is too entangled to extract
// into a standalone snip — so TrimAndNormalize takes the raw isotope-indexed intensities as input and
// is itself differential-tested (trim_cpp) with a fixed raw distribution.

using System;
using System.Collections.Concurrent;
using Chemistry; // Constants.C13MinusC12

namespace MassSpectrometry
{
    internal sealed class MetaFlashDeconAveragine
    {
        // ── Model: OpenMS-style trimmed/normalised isotope distributions over the mass range ──
        // Wraps an existing AverageResidue (our raw isotope generator) and applies TrimAndNormalize
        // per averagine index, lazily + thread-safely (per-scan deconvolution runs in parallel).
        private readonly AverageResidue _source;
        private readonly double _isoDa;
        private readonly ConcurrentDictionary<int, Entry> _cache = new();

        private readonly record struct Entry(double[] B, int Apex, int LeftFromApex, int RightFromApex, double AverageMassDelta);

        internal MetaFlashDeconAveragine(AverageResidue source, double isoDa)
        {
            _source = source;
            _isoDa = isoDa;
        }

        // TEST SEAM ONLY — a constant averagine (mass-independent apex/left/right/avgDelta), so the
        // candidate-generation logic can be differential-tested vs the C++ on IDENTICAL averagine
        // inputs (isolating logic bugs from the isotope-generator boundary). b is unused by candidate
        // generation, so it is empty here. ⚠ Never use this ctor in production — it returns the same
        // envelope shape for every mass and would destroy scoring; it exists purely for the diff tests.
        private readonly bool _isConstant;
        private readonly Entry _constantEntry;
        internal MetaFlashDeconAveragine(int apex, int left, int right, double avgDelta)
            : this(null, Constants.C13MinusC12)
        {
            _isConstant = true;
            _constantEntry = new Entry(Array.Empty<double>(), apex, left, right, avgDelta);
        }

        // Shared instances so the per-mass trim cache persists across the (per-scan, parallel)
        // deconvolution calls — a fresh algorithm is created per spectrum, but the averagine model is
        // immutable and reusable.
        private static readonly ConcurrentDictionary<(AverageResidue, double), MetaFlashDeconAveragine> _instances = new();
        internal static MetaFlashDeconAveragine For(AverageResidue source, double isoDa)
            => _instances.GetOrAdd((source, isoDa), k => new MetaFlashDeconAveragine(k.Item1, k.Item2));

        private Entry GetEntry(double mass)
        {
            if (_isConstant) return _constantEntry;
            int idx = _source.GetMostIntenseMassIndex(mass);
            return _cache.GetOrAdd(idx, i =>
            {
                double[] masses = _source.GetAllTheoreticalMasses(i);
                double[] intens = _source.GetAllTheoreticalIntensities(i);
                // arrays are intensity-descending: [0] is the most-abundant peak.
                double monoMass = masses[0] - _source.GetDiffToMonoisotopic(i);

                int maxIso = 0;
                for (int k = 0; k < masses.Length; k++)
                {
                    int iso = (int)Math.Round((masses[k] - monoMass) / _isoDa, MidpointRounding.AwayFromZero);
                    if (iso > maxIso) maxIso = iso;
                }
                var raw = new double[maxIso + 1];
                for (int k = 0; k < masses.Length; k++)
                {
                    int iso = (int)Math.Round((masses[k] - monoMass) / _isoDa, MidpointRounding.AwayFromZero);
                    if (iso >= 0 && iso <= maxIso) raw[iso] += intens[k];
                }
                double[] b = TrimAndNormalize(raw, out int apex, out int left, out int right);
                // average_mono_mass_difference_: intensity-weighted mean isotope offset (Da from mono).
                double sw = 0, swi = 0;
                for (int k = 0; k < b.Length; k++) { sw += b[k]; swi += b[k] * k; }
                double avgDelta = sw > 0 ? swi / sw * _isoDa : 0;
                return new Entry(b, apex, left, right, avgDelta);
            });
        }

        /// <summary>Trimmed + L2-normalised isotope distribution (OpenMS <c>avg.get(mass)</c>).</summary>
        internal double[] Get(double mass) => GetEntry(mass).B;
        /// <summary>OpenMS <c>avg.getApexIndex(mass)</c>.</summary>
        internal int GetApexIndex(double mass) => GetEntry(mass).Apex;
        /// <summary>OpenMS <c>avg.getLeftCountFromApex(mass)</c>.</summary>
        internal int GetLeftCountFromApex(double mass) => GetEntry(mass).LeftFromApex;
        /// <summary>OpenMS <c>avg.getRightCountFromApex(mass)</c>.</summary>
        internal int GetRightCountFromApex(double mass) => GetEntry(mass).RightFromApex;
        /// <summary>OpenMS <c>avg.getLastIndex(mass)</c> = apex + right-count-from-apex.</summary>
        internal int GetLastIndex(double mass) { var e = GetEntry(mass); return e.Apex + e.RightFromApex; }
        /// <summary>OpenMS <c>avg.getAverageMassDelta(mass)</c> (intensity-weighted mean isotope offset, Da).</summary>
        internal double GetAverageMassDelta(double mass) => GetEntry(mass).AverageMassDelta;

        /// <summary>
        /// OpenMS <c>PrecalculatedAveragine::getAverageMassDelta</c> at the FLASHDeconv minimum mass (50 Da),
        /// from its <c>CoarseIsotopePatternGenerator</c> (FLASHDeconvHelperStructs.cpp:107). ⚠ This is the ONE
        /// averagine value mzLib cannot reproduce locally: mzLib's <c>IsotopicDistribution</c> + its isotope
        /// abundances give ~0.0251-0.0255 for the mass-50 averagine, vs OpenMS's coarse generator 0.02514582
        /// — the un-snippable isotope-generator boundary. Verified against the OpenMS-linked
        /// <c>avg_snip</c> / <c>FD_SCALARTRACE</c> (full-precision double). Used ONLY to set
        /// <c>mass_bin_min</c>; an exact match removes the ~2-bin axis shift that mis-assigns sparse-scan
        /// candidate charge ranges. See STATUS 2026-05-29.
        /// </summary>
        internal const double OpenMsAvgDeltaAtMass50 = 0.025145816488311823;

        /// <summary>
        /// avgDelta for the bin-axis origin (<c>mass_bin_min = log(min_mass - this)</c>). Returns the
        /// exact-OpenMS value at the FLASHDeconv default min_mass (50); falls back to the model elsewhere.
        /// (min_mass is 50 in every standard top-down run; a non-50 min_mass would need a faithful port of
        /// CoarseIsotopePatternGenerator to be axis-exact.)
        /// </summary>
        internal double GetMassBinMinAvgDelta(double minMass)
            => Math.Abs(minMass - 50.0) < 1e-6 ? OpenMsAvgDeltaAtMass50 : GetAverageMassDelta(minMass);

        /// <summary>
        /// Faithful port of the per-mass loop body in OpenMS
        /// <c>PrecalculatedAveragine::PrecalculatedAveragine</c> (FLASHDeconvHelperStructs.cpp:36-106).
        /// Given the raw isotope-indexed intensities (isotope 0 = monoisotopic), greedily trim the
        /// lower-intensity end peaks while keeping ≥99.99% of the summed power, drop trailing
        /// near-zero peaks (<c>trimRight(1e-10)</c>), and L2-normalise by the trimmed power. Returns
        /// the normalised distribution <c>b</c> and the apex / left-from-apex / right-from-apex counts
        /// (each floored at 2), exactly as OpenMS stores them.
        /// </summary>
        internal static double[] TrimAndNormalize(
            double[] rawIso, out int apexIndex, out int leftCountFromApex, out int rightCountFromApex)
        {
            const double minPwr = 0.9999;
            const int minIsoLength = 2;
            const int minLeftRightCount = 2;

            var iso = (double[])rawIso.Clone();
            double totalPwr = 0.0;
            int mostAbundantIndex = 0;
            double mostAbundantInt = 0.0;
            for (int k = 0; k < iso.Length; k++)
            {
                totalPwr += iso[k] * iso[k];
                if (mostAbundantInt >= iso[k]) continue;   // first-max-wins on ties (OpenMS `>=` continue)
                mostAbundantInt = iso[k];
                mostAbundantIndex = k;
            }

            int leftCount = 0;
            int rightCount = iso.Length - 1;
            int trimCount = 0;
            while (iso.Length - trimCount > minIsoLength && leftCount < rightCount)
            {
                double lint = iso[leftCount];
                double rint = iso[rightCount];
                double pwr;
                bool trimLeft = true;
                if (lint < rint) { pwr = lint * lint; }
                else { pwr = rint * rint; trimLeft = false; }

                if (totalPwr - pwr < totalPwr * minPwr) break;
                totalPwr -= pwr;
                trimCount++;
                if (trimLeft) { iso[leftCount] = 0.0; leftCount++; }
                else { iso[rightCount] = 0.0; rightCount--; }
            }

            leftCount = mostAbundantIndex - leftCount;
            rightCount = rightCount - mostAbundantIndex;

            // trimRight(1e-10): drop trailing near-zero entries (the right-trimmed ones).
            int newLen = iso.Length;
            while (newLen > 0 && iso[newLen - 1] <= 1e-10) newLen--;

            // ⚠ L2 normalisation (divide by sqrt of summed power), NOT L1. OpenMS getCosine divides only
            // by |a| and assumes the averagine target `b` is already unit-L2; an L1 (sum) norm here would
            // silently bias every isotope-cosine. This whole function is differential-tested (trim_cpp).
            double norm = Math.Sqrt(totalPwr);
            var b = new double[newLen];
            for (int k = 0; k < newLen; k++) b[k] = iso[k] / norm;

            leftCount = Math.Max(leftCount, minLeftRightCount);
            rightCount = Math.Max(rightCount, minLeftRightCount);

            apexIndex = mostAbundantIndex;
            leftCountFromApex = leftCount;
            rightCountFromApex = rightCount;
            return b;
        }
    }
}
