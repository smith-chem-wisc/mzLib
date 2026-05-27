// MetaFlashDeconPeakGroup.cs
//
// Plan B: a faithful C# port of the OpenMS FLASHDeconv PeakGroup scoring core, so MetaFlashDecon
// reproduces real FLASHDeconv per-spectrum behaviour (recruitment -> per-charge information ->
// charge-range/charge-fit -> mono-mass -> per-charge cosine -> SNR -> Qscore), enabling a faithful
// removeChargeErrorPeakGroups_ that uses genuine per-charge SNR.
//
// Reference (OpenMS @ E:\GitClones\OpenMS):
//   src/openms/source/ANALYSIS/TOPDOWN/PeakGroup.cpp           (per-charge chain, recruitment, SNR)
//   src/openms/source/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.cpp (getCosine, scoring loop, removeChargeError)
//   src/openms/source/ANALYSIS/TOPDOWN/Qscore.cpp               (logistic Qscore)
//
// Each ported function is pinned to the C++ via the differential-test harness in
// E:\CodeReview\MetaFlashDecon\difftest (MSVC extract vs this C#), see MetaFlashDeconDiffTest.
//
// STATUS: under construction (built incrementally, snip-validated function-by-function).

using System;
using System.Collections.Generic;

namespace MassSpectrometry
{
    internal sealed partial class MetaFlashDeconPeakGroup
    {
        /// <summary>
        /// Faithful port of OpenMS <c>FLASHDeconvAlgorithm::getCosine</c>
        /// (FLASHDeconvAlgorithm.cpp:1315-1386). Cosine of an observed per-isotope intensity vector
        /// <paramref name="a"/> against the (L2-normalised) averagine isotope intensities
        /// <paramref name="b"/>, with <paramref name="b"/> aligned to <paramref name="a"/> by
        /// <paramref name="offset"/> (<c>i = j - offset</c>).
        /// <para>
        /// Faithful subtleties: the denominator <c>a_norm</c> sums a[j]² over the WHOLE window
        /// (including observed peaks that align to no valid b index — so unexplained observed
        /// intensity lowers the score), the numerator only adds aligned terms, and when
        /// <paramref name="minIsoSize"/> &gt; 0 the score is forced to 0 if the maximum-intensity
        /// isotope has no non-zero neighbour. Since OpenMS does NOT divide by |b|, <paramref name="b"/>
        /// is expected L2-normalised (as OpenMS's PrecalculatedAveragine is); callers must normalise.
        /// </para>
        /// </summary>
        internal static double GetCosine(
            IReadOnlyList<double> a, int aStart, int aEnd,
            IReadOnlyList<double> b, int bSize, int offset, int minIsoSize)
        {
            double n = 0.0, aNorm = 0.0;
            aStart = Math.Max(0, aStart);
            aEnd = Math.Min(a.Count, aEnd);

            if (aEnd - aStart < minIsoSize) return 0.0;

            int maxIntensityIndex = 0;
            double maxIntensity = 0.0;

            for (int j = aStart; j < aEnd; j++)
            {
                int i = j - offset;
                aNorm += a[j] * a[j];

                if (maxIntensity < a[j])
                {
                    maxIntensity = a[j];
                    maxIntensityIndex = j;
                }

                if (i < 0 && a[j] > 0)
                {
                    // n -= a[j] * b[0];  // commented out in OpenMS — kept for fidelity
                }
                else if (i >= bSize || i < 0 || b[i] <= 0)
                {
                    continue;
                }
                else
                {
                    n += a[j] * b[i];
                }
            }

            // two consecutive isotopes around the max-intensity isotope
            if (minIsoSize > 0)
            {
                if (maxIntensityIndex == aEnd - 1)
                {
                    if (maxIntensityIndex > 0 && a[maxIntensityIndex - 1] == 0) return 0.0;
                }
                else if (maxIntensityIndex == aStart)
                {
                    if (maxIntensityIndex + 1 < a.Count && a[maxIntensityIndex + 1] == 0) return 0.0;
                }
                else if (maxIntensityIndex > 0 && maxIntensityIndex + 1 < a.Count)
                {
                    if (a[maxIntensityIndex + 1] == 0 && a[maxIntensityIndex - 1] == 0) return 0.0;
                }
            }

            if (aNorm <= 0) return 0.0;
            return n / Math.Sqrt(aNorm);
        }
    }
}
