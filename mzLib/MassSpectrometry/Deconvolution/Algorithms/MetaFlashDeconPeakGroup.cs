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

        // ── Per-charge information (OpenMS updatePerChargeInformation_, PeakGroup.cpp:403-445) ──
        internal double[] PerChargeInt;          // per_charge_int_       (Σ signal intensity per charge)
        internal double[] PerChargeSumSignalSq;  // per_charge_sum_signal_squared_
        internal double[] PerChargeNoisePwr;     // per_charge_noise_pwr_ (getNoisePeakPower_ per charge)

        /// <summary>
        /// Faithful port of OpenMS <c>PeakGroup::updatePerChargeInformation_</c>
        /// (PeakGroup.cpp:403-445). Per charge: summed signal intensity and summed squared signal
        /// intensity from <see cref="SignalPeaks"/> (isotopeIndex&gt;=0 only), and the structured
        /// noise power from <see cref="MetaFlashDeconAlgorithm.ComputeNoisePeakPower"/> over that
        /// charge's noisy + signal peaks. Vectors are sized <c>1 + MaxAbsCharge</c>.
        /// </summary>
        internal void UpdatePerChargeInformation(List<LogMzPeak> noisyPeaks)
        {
            int sz = 1 + MaxAbsCharge;
            PerChargeNoisePwr = new double[sz];
            PerChargeSumSignalSq = new double[sz];
            PerChargeInt = new double[sz];

            foreach (var p in SignalPeaks)
            {
                PerChargeInt[p.AbsCharge] += p.Intensity;
                PerChargeSumSignalSq[p.AbsCharge] += p.Intensity * p.Intensity;
            }

            for (int z = MinAbsCharge; z <= MaxAbsCharge; z++)
            {
                var chargeNoisy = new List<(double mz, double intensity)>();
                var chargeSignal = new List<(double mz, double intensity)>();
                foreach (var p in noisyPeaks) if (p.AbsCharge == z) chargeNoisy.Add((p.Mz, p.Intensity));
                foreach (var p in SignalPeaks) if (p.AbsCharge == z) chargeSignal.Add((p.Mz, p.Intensity));
                PerChargeNoisePwr[z] = MetaFlashDeconAlgorithm.ComputeNoisePeakPower(chargeNoisy, chargeSignal, z, IsoDaDistance);
            }
        }

        // ── SNR (OpenMS updateSNR_, PeakGroup.cpp:893-919) ────────────────────
        internal double[] PerChargeCos;   // per_charge_cos_ (set by per-charge cosine step)
        internal double[] PerChargeSnr;   // per_charge_snr_ (output)
        internal double Snr;              // snr_ (overall, output)
        internal double IsotopeCosineScore; // isotope_cosine_score_ (global cosine)

        /// <summary>
        /// Faithful port of OpenMS <c>PeakGroup::updateSNR_</c> (PeakGroup.cpp:893-919). Per charge:
        /// <c>snr_c = cos_c² · int_c² / (1 + noise_c + (1 − cos_c²)·sumSigSq_c)</c>; overall:
        /// <c>snr = globalCos² · Σint_c² / (1 + Σnoise_c + (1 − globalCos²)·ΣsumSigSq_c)</c>.
        /// Consumes <see cref="PerChargeCos"/>, <see cref="PerChargeInt"/>,
        /// <see cref="PerChargeNoisePwr"/>, <see cref="PerChargeSumSignalSq"/> and
        /// <see cref="IsotopeCosineScore"/>; writes <see cref="PerChargeSnr"/> and <see cref="Snr"/>.
        /// </summary>
        internal void UpdateSNR()
        {
            double cosSquared = IsotopeCosineScore * IsotopeCosineScore;
            double signal = 0, noise = 0, sumSignalSquared = 0;
            PerChargeSnr = new double[1 + MaxAbsCharge];

            int upper = Math.Min(PerChargeSumSignalSq.Length, 1 + MaxAbsCharge);
            for (int c = MinAbsCharge; c < upper; c++)
            {
                if (PerChargeCos != null && PerChargeCos.Length > c)
                {
                    double pcc2 = PerChargeCos[c] * PerChargeCos[c];
                    double nom = pcc2 * PerChargeInt[c] * PerChargeInt[c];
                    double denom = 1 + PerChargeNoisePwr[c] + (1 - pcc2) * PerChargeSumSignalSq[c];
                    PerChargeSnr[c] = denom <= 0 ? 0.0 : nom / denom;
                }
                sumSignalSquared += PerChargeSumSignalSq[c];
                signal += PerChargeInt[c] * PerChargeInt[c];
                noise += PerChargeNoisePwr[c];
            }

            double tNom = cosSquared * signal;
            double tDenom = 1 + noise + (1 - cosSquared) * sumSignalSquared;
            Snr = tDenom <= 0 ? 0.0 : tNom / tDenom;
        }

        /// <summary>Per-charge SNR accessor mirroring OpenMS <c>getChargeSNR</c> (0 if out of range).</summary>
        internal double GetChargeSnr(int absCharge)
            => (PerChargeSnr == null || absCharge < 0 || absCharge >= PerChargeSnr.Length) ? 0.0 : PerChargeSnr[absCharge];

        // ── Charge-range refinement (OpenMS updateChargeRange_, PeakGroup.cpp:447-530) ──
        /// <summary>
        /// Faithful port of OpenMS <c>PeakGroup::updateChargeRange_</c>. Finds the max-SNR-proxy
        /// charge, then expands the [min,max] range outward while a charge's signal power exceeds a
        /// threshold (<c>min(maxSig/10, 1)</c>), and re-filters <see cref="SignalPeaks"/> and
        /// <paramref name="noisyPeaks"/> to the refined range. Consumes the per-charge arrays from
        /// <see cref="UpdatePerChargeInformation"/> (must be called first).
        /// </summary>
        internal void UpdateChargeRange(List<LogMzPeak> noisyPeaks)
        {
            int maxSigCharge = 0;
            double maxSig = 0;
            for (int z = MinAbsCharge; z <= MaxAbsCharge; z++)
            {
                double tmpSnr = PerChargeInt[z] * PerChargeInt[z] / (1 + PerChargeNoisePwr[z]);
                if (maxSig < tmpSnr) { maxSig = tmpSnr; maxSigCharge = z; }
            }

            int newMax = maxSigCharge, newMin = maxSigCharge;
            double threshold = Math.Min(maxSig / 10, 1.0);
            for (int z = maxSigCharge; z <= MaxAbsCharge; z++)
            {
                double psp = PerChargeInt[z] * PerChargeInt[z];
                if (psp / (1 + PerChargeNoisePwr[z]) < threshold) break;
                newMax = z;
            }
            for (int z = maxSigCharge; z >= MinAbsCharge; z--)
            {
                double psp = PerChargeInt[z] * PerChargeInt[z];
                if (psp / (1 + PerChargeNoisePwr[z]) < threshold) break;
                newMin = z;
            }

            if (MaxAbsCharge != newMax || MinAbsCharge != newMin)
            {
                SignalPeaks.RemoveAll(p => p.AbsCharge < newMin || p.AbsCharge > newMax);
                noisyPeaks.RemoveAll(p => p.AbsCharge < newMin || p.AbsCharge > newMax);
                MaxAbsCharge = newMax;
                MinAbsCharge = newMin;
            }
            if (MinAbsCharge > MaxAbsCharge)
            {
                SignalPeaks.Clear();
                NegativeIsoPeaks.Clear();
            }
        }

        // ── Charge fit score (OpenMS updateChargeFitScoreAndChargeIntensities_, PeakGroup.cpp:620-686) ──
        internal double ChargeScore;
        /// <summary>
        /// Faithful port of OpenMS <c>PeakGroup::updateChargeFitScoreAndChargeIntensities_</c>.
        /// Measures how unimodal the per-charge intensity profile is: 1 minus the summed positive
        /// "uphill" intensity jumps away from the max charge, over the total (min-subtracted)
        /// intensity. OpenMS gates groups with score &lt; 0.7.
        /// </summary>
        internal void UpdateChargeFitScoreAndChargeIntensities()
        {
            if (MaxAbsCharge == MinAbsCharge) { ChargeScore = 1; return; }

            double maxPerChargeIntensity = 0, summedIntensity = 0;
            int maxIndex = -1, firstIndex = -1, lastIndex = -1;
            double minIntensity = -1;
            for (int c = MinAbsCharge; c <= MaxAbsCharge; c++)
                if (minIntensity < 0 || minIntensity > PerChargeInt[c]) minIntensity = PerChargeInt[c];

            for (int c = MinAbsCharge; c <= MaxAbsCharge; c++)
            {
                summedIntensity += PerChargeInt[c] - minIntensity;
                if (PerChargeInt[c] > 0)
                {
                    if (firstIndex < 0) firstIndex = c;
                    lastIndex = c;
                }
                if (maxPerChargeIntensity > PerChargeInt[c]) continue;
                maxPerChargeIntensity = PerChargeInt[c];
                maxIndex = c;
            }
            if (maxIndex < 0) { ChargeScore = 0; return; }

            firstIndex = firstIndex < 0 ? 0 : firstIndex;
            double p = 0;
            for (int c = maxIndex; c < lastIndex; c++)
            {
                double diff = PerChargeInt[c + 1] - PerChargeInt[c];
                if (diff <= 0) continue;
                p += diff;
            }
            for (int c = maxIndex; c > firstIndex; c--)
            {
                double diff = PerChargeInt[c - 1] - PerChargeInt[c];
                if (diff <= 0) continue;
                p += diff;
            }
            ChargeScore = Math.Max(0.0, 1.0 - p / summedIntensity);
        }

        // ── Mono mass + per-isotope intensities (OpenMS updateMonoMassAndIsotopeIntensities, PeakGroup.cpp:688-726) ──
        internal double[] PerIsotopeInt;  // per_isotope_int_ (indexed by isotopeIndex - MinNegativeIsotopeIndex)
        internal double Intensity;        // intensity_ (Σ signal intensity, iso>=0)
        /// <summary>
        /// Faithful port of OpenMS <c>PeakGroup::updateMonoMassAndIsotopeIntensities</c>. Builds the
        /// per-isotope intensity vector and the intensity-weighted monoisotopic mass
        /// (<c>Σ pi·(unchargedMass − iso·isoDa) / Σ pi</c> over signal peaks, iso&gt;=0). Negative-iso
        /// peaks contribute only to <see cref="PerIsotopeInt"/>.
        /// </summary>
        internal void UpdateMonoMassAndIsotopeIntensities(Polarity polarity)
        {
            if (SignalPeaks.Count == 0) return;

            int maxIsotopeIndex = 0;
            foreach (var p in SignalPeaks) maxIsotopeIndex = Math.Max(maxIsotopeIndex, p.IsotopeIndex);

            PerIsotopeInt = new double[maxIsotopeIndex + 1 - MinNegativeIsotopeIndex];
            Intensity = 0.0;
            double nominator = 0.0;
            double chargeMass = Math.Sign((int)polarity) * Constants.ProtonMass;

            foreach (var p in SignalPeaks)
            {
                double pi = p.Intensity;
                if (p.IsotopeIndex < 0) continue;
                PerIsotopeInt[p.IsotopeIndex - MinNegativeIsotopeIndex] += pi;
                double unchargedMass = (p.Mz - chargeMass) * p.AbsCharge; // OpenMS LogMzPeak::getUnchargedMass
                nominator += pi * (unchargedMass - p.IsotopeIndex * IsoDaDistance);
                Intensity += pi;
            }
            foreach (var p in NegativeIsoPeaks)
            {
                if (p.IsotopeIndex - MinNegativeIsotopeIndex < 0) continue;
                PerIsotopeInt[p.IsotopeIndex - MinNegativeIsotopeIndex] += p.Intensity;
            }

            MonoisotopicMass = nominator / Intensity;
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
