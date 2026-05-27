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

namespace MassSpectrometry
{
    internal static class MetaFlashDeconAveragine
    {
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
