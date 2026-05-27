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
using Chemistry; // Constants.ProtonMass

namespace MassSpectrometry
{
    using LogMzPeak = MetaFlashDeconAlgorithm.LogMzPeak;

    internal sealed partial class MetaFlashDeconPeakGroup
    {
        // ── PeakGroup state (mirrors OpenMS PeakGroup members) ────────────────
        internal int MinAbsCharge;
        internal int MaxAbsCharge;
        internal bool IsPositive = true;
        internal double IsoDaDistance = Constants.C13MinusC12;     // OpenMS iso_da_distance_ (ISOTOPE_MASSDIFF_55K_U)
        internal int MinNegativeIsotopeIndex = -1;                 // OpenMS min_negative_isotope_index_ default
        internal double MonoisotopicMass;

        /// <summary>Signal peaks with isotopeIndex &gt;= 0 (OpenMS <c>logMzpeaks_</c>).</summary>
        internal readonly List<LogMzPeak> SignalPeaks = new();
        /// <summary>Signal peaks with isotopeIndex &lt; 0 (OpenMS <c>negative_iso_peaks_</c>).</summary>
        internal readonly List<LogMzPeak> NegativeIsoPeaks = new();

        /// <summary>
        /// Faithful port of OpenMS <c>PeakGroup::recruitAllPeaksInSpectrum</c>
        /// (PeakGroup.cpp:532-618). For each charge in [MinAbsCharge, MaxAbsCharge] (scanned high→low),
        /// classify every spectrum peak in the isotope m/z range: on the isotope grid (within ±tol)
        /// → signal (isotopeIndex &gt;= 0 to <see cref="SignalPeaks"/>, &lt; 0 to
        /// <see cref="NegativeIsoPeaks"/>), off-grid but in range (isotopeIndex &gt;= 0) → noise
        /// (returned). Each peak is tagged with <c>abs_charge</c> and <c>isotopeIndex</c>; the same
        /// raw peak can be claimed by several charges (signal for one, noise for another).
        /// <para>
        /// <paramref name="minIsotope"/>/<paramref name="maxIsotope"/> are the averagine-derived scan
        /// bounds (OpenMS <c>avg.getApexIndex/getLeftCountFromApex/getLastIndex</c>); passed in so the
        /// recruitment logic is differential-tested independently of the averagine. OpenMS's
        /// <c>findNearest</c> scan-start is a perf optimisation only — the iso-index break/continue
        /// makes scanning from index 0 yield an identical set.
        /// </para>
        /// </summary>
        internal List<LogMzPeak> RecruitAllPeaksInSpectrum(
            MzSpectrum spectrum, double tol, double monoMass, int minIsotope, int maxIsotope, Polarity polarity)
        {
            var noisyPeaks = new List<LogMzPeak>();
            if (monoMass < 0) return noisyPeaks;

            MonoisotopicMass = monoMass;
            SignalPeaks.Clear();
            NegativeIsoPeaks.Clear();

            int polSign = Math.Sign((int)polarity);
            double chargeMass = polSign * Constants.ProtonMass; // OpenMS getChargeMass(is_positive_)

            for (int c = MaxAbsCharge; c >= MinAbsCharge; c--)
            {
                if (c <= 0) break;
                double cmz = monoMass / c + chargeMass;
                double isoDelta = IsoDaDistance / c;

                for (int index = 0; index < spectrum.Size; index++)
                {
                    double pint = spectrum.YArray[index];
                    if (pint <= 0) continue;
                    double pmz = spectrum.XArray[index];
                    // C++ round() is half-away-from-zero (banker's in C# by default).
                    int isoIndex = (int)Math.Round((pmz - cmz) / isoDelta, MidpointRounding.AwayFromZero);
                    if (isoIndex > maxIsotope) break;
                    if (isoIndex < minIsotope) continue;

                    if (Math.Abs(pmz - cmz - isoIndex * isoDelta) <= pmz * tol)
                    {
                        var p = new LogMzPeak(pmz, pint, MetaFlashDeconAlgorithm.GetLogMz(pmz, polarity), c, isoIndex);
                        if (isoIndex < 0) NegativeIsoPeaks.Add(p);
                        else SignalPeaks.Add(p);
                    }
                    else if (isoIndex >= 0)
                    {
                        noisyPeaks.Add(new LogMzPeak(pmz, pint, MetaFlashDeconAlgorithm.GetLogMz(pmz, polarity), c, isoIndex));
                    }
                }
            }
            return noisyPeaks;
        }

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
