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
        internal double IsoDaDistance = MetaFlashDeconAlgorithm.IsoDaDistance55K; // OpenMS iso_da_distance_ = ISOTOPE_MASSDIFF_55K_U = 1.002371
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
        internal double QscoreValue;      // qscore_
        internal int RepAbsCharge;        // max_snr_abs_charge_ (representative charge)
        internal Polarity Polarity = Polarity.Positive;

        /// <summary>Per-charge isotope cosine accessor mirroring OpenMS <c>getChargeIsotopeCosine</c>.</summary>
        internal double GetChargeIsotopeCosine(int absCharge)
            => (PerChargeCos == null || absCharge < 0 || absCharge >= PerChargeCos.Length) ? 0.0 : PerChargeCos[absCharge];

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

        // ── Per-charge cosine (OpenMS updatePerChargeCos_, PeakGroup.cpp:79-115) ──
        /// <summary>
        /// Faithful port of OpenMS <c>PeakGroup::updatePerChargeCos_</c>. For each charge, build the
        /// observed per-isotope intensity vector (indexed by isotopeIndex) from that charge's signal
        /// peaks and score it against the trimmed/normalised averagine <c>b</c> via
        /// <see cref="GetCosine"/> (offset 0, minIsoSize 0). Requires
        /// <see cref="UpdateMonoMassAndIsotopeIntensities"/> first (sets <see cref="PerIsotopeInt"/>).
        /// </summary>
        internal void UpdatePerChargeCos(MetaFlashDeconAveragine avg)
            => UpdatePerChargeCos(avg.Get(MonoisotopicMass));

        /// <summary>Overload taking the averagine distribution <c>b</c> directly (for differential testing).</summary>
        internal void UpdatePerChargeCos(double[] b)
        {
            int isoSize = b.Length;
            // OpenMS: current_per_isotope_intensities size = getIsotopeIntensities().size() + min_negative_isotope_index_
            int curSize = PerIsotopeInt.Length + MinNegativeIsotopeIndex;
            if (curSize < 0) curSize = 0;
            PerChargeCos = new double[1 + MaxAbsCharge];

            for (int z = MinAbsCharge; z <= MaxAbsCharge; z++)
            {
                var cur = new double[curSize];
                int minIso = curSize, maxIso = -1;
                foreach (var p in SignalPeaks)
                {
                    if (p.AbsCharge != z) continue;
                    if (p.IsotopeIndex >= curSize || p.IsotopeIndex < 0) continue;
                    cur[p.IsotopeIndex] += p.Intensity;
                    if (p.IsotopeIndex < minIso) minIso = p.IsotopeIndex;
                    if (p.IsotopeIndex > maxIso) maxIso = p.IsotopeIndex;
                }
                PerChargeCos[z] = GetCosine(cur, minIso, maxIso + 1, b, isoSize, 0, 0);
            }
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

        /// <summary>
        /// Faithful port of OpenMS <c>FLASHDeconvAlgorithm::getIsotopeCosineAndDetermineIsotopeIndex</c>
        /// (FLASHDeconvAlgorithm.cpp:1226-1313), target path only (the isotope-dummy second-best
        /// branch is omitted — MetaFlashDecon scores targets). Searches isotope offsets in
        /// [-left, right] (where right = apexIndex/4 + 1, shifted by <paramref name="isoIntShift"/>),
        /// returns the best cosine of <paramref name="perIsotopeIntensities"/> vs the normalised
        /// averagine <paramref name="b"/>, and outputs the isotope-index correction
        /// <paramref name="offset"/> (already de-shifted by <paramref name="isoIntShift"/>).
        /// </summary>
        internal static double GetIsotopeCosineAndDetermineIsotopeIndex(
            double[] perIsotopeIntensities, double[] b, int apexIndex,
            int isoIntShift, int windowWidth, int minIsoSize, out int offset)
        {
            offset = 0;
            if (perIsotopeIntensities.Length < minIsoSize + isoIntShift) return 0.0;

            int isoSize = b.Length;
            int right = apexIndex / 4 + 1;
            int left = right;
            if (windowWidth >= 0) { right = Math.Min(right, windowWidth); left = Math.Min(left, windowWidth); }

            double maxCos = -1000.0;
            int maxIsotopeIndex = perIsotopeIntensities.Length; // exclusive
            int minIsotopeIndex = -1;                            // inclusive

            left -= isoIntShift;
            right += isoIntShift;

            for (int i = 0; i < maxIsotopeIndex; i++)
            {
                if (perIsotopeIntensities[i] <= 0) continue;
                if (minIsotopeIndex < 0) minIsotopeIndex = i;
            }
            if (maxIsotopeIndex - minIsotopeIndex < minIsoSize) return 0.0;

            for (int tmpOffset = -left; tmpOffset <= right; tmpOffset++)
            {
                double tmpCos = GetCosine(perIsotopeIntensities, minIsotopeIndex, maxIsotopeIndex, b, isoSize, tmpOffset, minIsoSize);
                if (maxCos < tmpCos) { maxCos = tmpCos; offset = tmpOffset; }
            }

            offset -= isoIntShift;
            return maxCos;
        }

        // ── Qscore (current OpenMS Qscore.cpp:16-61) ──────────────────────────
        /// <summary>
        /// Faithful port of the CURRENT OpenMS <c>Qscore::getQscore</c> (Qscore.cpp). Logistic over
        /// two live features (PPM-error and charge-score weights are 0): weights
        /// <c>{-2.2833, -3.2881, 0, 0, 4.5425(intercept)}</c>, features
        /// <c>f0 = log2(cosine + 1)</c>, <c>f1 = log2(1 + snr/(1+snr))</c>, then
        /// <c>q = 1/(1 + exp(intercept + w0·f0 + w1·f1))</c>. Uses the GLOBAL isotope cosine and
        /// GLOBAL SNR (so it is the same for every charge in OpenMS's updateQscore loop).
        /// </summary>
        internal static double ComputeQscore(double isotopeCosine, double snr)
        {
            double f0 = Math.Log(isotopeCosine + 1.0, 2.0);
            double f1 = Math.Log(1.0 + snr / (1.0 + snr), 2.0);
            double score = 4.5425 + (-2.2833) * f0 + (-3.2881) * f1;
            return 1.0 / (1.0 + Math.Exp(score));
        }

        // ── Orchestrator (OpenMS PeakGroup::updateQscore, PeakGroup.cpp:117-178) ──
        /// <summary>
        /// Faithful port of OpenMS <c>PeakGroup::updateQscore</c>: runs the per-charge chain in order
        /// with the gates (charge-fit ≥0.7, isotope cosine ≥<paramref name="minCos"/>, charge-range
        /// ≥ maxCharge/20), then sets <see cref="QscoreValue"/> (current OpenMS Qscore over the global
        /// cosine + SNR) and <see cref="RepAbsCharge"/> (max-per-charge-SNR charge). Returns the
        /// isotope-index <c>offset</c> for the caller's recruit→score refinement loop (0 = converged).
        /// Requires the signal/negative peaks to be already recruited (<paramref name="noisyPeaks"/>
        /// is the matching noise list from recruitment).
        /// </summary>
        internal int UpdateQscore(List<LogMzPeak> noisyPeaks, MetaFlashDeconAveragine avg, double minCos)
        {
            QscoreValue = 0;
            UpdatePerChargeInformation(noisyPeaks);
            UpdateChargeRange(noisyPeaks);
            if (SignalPeaks.Count == 0) return 0;

            UpdateChargeFitScoreAndChargeIntensities();
            if (ChargeScore < 0.7f) return 0;

            UpdateMonoMassAndIsotopeIntensities(Polarity);
            if (PerIsotopeInt == null || PerIsotopeInt.Length == 0 || MaxAbsCharge < MinAbsCharge) return 0;

            int isoIntShift = -MinNegativeIsotopeIndex; // = 1
            double[] b = avg.Get(MonoisotopicMass);
            int apex = avg.GetApexIndex(MonoisotopicMass);
            IsotopeCosineScore = GetIsotopeCosineAndDetermineIsotopeIndex(
                PerIsotopeInt, b, apex, isoIntShift, -1, 2, out int hOffset);

            if (IsotopeCosineScore < minCos) return 0;
            if (MaxAbsCharge - MinAbsCharge < MaxAbsCharge / 20) return 0;

            UpdatePerChargeCos(avg);
            UpdateSNR();

            RepAbsCharge = 0;
            for (int c = MinAbsCharge; c <= MaxAbsCharge; c++)
            {
                if (GetChargeSnr(c) <= 0 || GetChargeIsotopeCosine(c) <= 0) continue;
                double q = ComputeQscore(IsotopeCosineScore, Snr);
                if (QscoreValue < q) QscoreValue = q;
                if (GetChargeSnr(c) > GetChargeSnr(RepAbsCharge)) RepAbsCharge = c;
            }
            return hOffset;
        }

        /// <summary>Group total signal intensity (OpenMS <c>getIntensity()</c>).</summary>
        internal double GetIntensity() => Intensity;

        internal bool IsTargeted; // we have no targeted masses -> always false

        // ── Charge-error removal (OpenMS removeChargeErrorPeakGroups_, FLASHDeconvAlgorithm.cpp:1388-1464) ──
        /// <summary>
        /// Faithful port of OpenMS <c>FLASHDeconvAlgorithm::removeChargeErrorPeakGroups_</c>. For each
        /// raw signal peak shared by ≥2 groups, a group <c>i</c> is "overlapped" at that peak if some
        /// other group <c>j</c> interprets the peak as a DIFFERENT charge AND is not ≥2× worse by
        /// per-charge SNR (<c>getChargeSNR(repz_i) &gt; getChargeSNR(repz_j)·2</c> ⇒ i wins, skip).
        /// A group whose overlapped intensity reaches ≥50% of its total intensity is dropped as a
        /// charge-error ghost. The charge each shared peak implies is <c>round(mass / (pmz − chargeMass))</c>.
        /// This is the lever for the multi-charge over-generation; it relies on the (now validated)
        /// per-charge SNR — using the Qscore here (the prior attempt) wrongly dropped true masses.
        /// </summary>
        internal static List<MetaFlashDeconPeakGroup> RemoveChargeErrorPeakGroups(
            List<MetaFlashDeconPeakGroup> groups, Polarity polarity)
        {
            double chargeMass = Math.Sign((int)polarity) * Constants.ProtonMass;

            var peakToPgs = new Dictionary<double, HashSet<int>>();
            var mzToIntensity = new Dictionary<double, double>();
            for (int i = 0; i < groups.Count; i++)
            {
                foreach (var p in groups[i].SignalPeaks)
                {
                    if (!peakToPgs.TryGetValue(p.Mz, out var set)) { set = new HashSet<int>(); peakToPgs[p.Mz] = set; }
                    set.Add(i);
                    mzToIntensity[p.Mz] = p.Intensity;
                }
            }

            var overlapIntensity = new double[groups.Count];
            foreach (var kv in peakToPgs)
            {
                var pgIs = kv.Value;
                if (pgIs.Count == 1) continue;
                double pmz = kv.Key;
                double pint = mzToIntensity[pmz];

                foreach (int i in pgIs)
                {
                    bool isOverlap = false;
                    double mass1 = groups[i].MonoisotopicMass;
                    int repz1 = (int)Math.Round(mass1 / (pmz - chargeMass), MidpointRounding.AwayFromZero);
                    foreach (int j in pgIs)
                    {
                        if (i == j) continue;
                        double mass2 = groups[j].MonoisotopicMass;
                        int repz2 = (int)Math.Round(mass2 / (pmz - chargeMass), MidpointRounding.AwayFromZero);
                        if (repz1 == repz2) continue;
                        if (groups[i].GetChargeSnr(repz1) > groups[j].GetChargeSnr(repz2) * 2.0) continue; // i clearly better -> j is the ghost
                        isOverlap = true;
                        break;
                    }
                    if (isOverlap) overlapIntensity[i] += pint;
                }
            }

            var filtered = new List<MetaFlashDeconPeakGroup>(groups.Count);
            for (int i = 0; i < groups.Count; i++)
            {
                if (!groups[i].IsTargeted && overlapIntensity[i] >= groups[i].GetIntensity() * 0.5) continue;
                if (groups[i].RepAbsCharge < groups[i].MinAbsCharge || groups[i].RepAbsCharge > groups[i].MaxAbsCharge) continue;
                filtered.Add(groups[i]);
            }
            return filtered;
        }

        // ── Overlap dedup (OpenMS removeOverlappingPeakGroups_, FLASHDeconvAlgorithm.cpp:1466-1515) ──
        /// <summary>
        /// Faithful port of OpenMS <c>removeOverlappingPeakGroups_</c>: over groups sorted by mono
        /// mass, keep only the highest-SNR group within each <c>mass·windowTol</c> window. (Targeted
        /// handling omitted — MetaFlashDecon has no targeted masses.)
        /// </summary>
        internal static List<MetaFlashDeconPeakGroup> RemoveOverlappingPeakGroups(
            List<MetaFlashDeconPeakGroup> groups, double windowTol)
        {
            if (groups.Count == 0) return groups;
            groups.Sort((a, b) => a.MonoisotopicMass.CompareTo(b.MonoisotopicMass));

            var filtered = new List<MetaFlashDeconPeakGroup>(groups.Count);
            double startMass = groups[0].MonoisotopicMass;
            double localMaxSnr = 0;
            int localMaxIndex = 0;
            for (int i = 0; i < groups.Count; i++)
            {
                double mass = groups[i].MonoisotopicMass;
                if (mass - startMass > mass * windowTol)
                {
                    filtered.Add(groups[localMaxIndex]);
                    startMass = mass;
                    localMaxSnr = 0;
                }
                if (localMaxSnr < groups[i].Snr)
                {
                    localMaxSnr = groups[i].Snr;
                    localMaxIndex = i;
                }
            }
            if (localMaxSnr > 0) filtered.Add(groups[localMaxIndex]);
            return filtered;
        }
    }
}
