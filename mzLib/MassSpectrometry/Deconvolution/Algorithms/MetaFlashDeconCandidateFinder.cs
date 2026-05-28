// MetaFlashDeconCandidateFinder.cs
//
// Plan B candidate generation: faithful port of OpenMS FLASHDeconvAlgorithm's candidate-mass binning
// — updateMzBins_ (:273), updateCandidateMassBins_ (:292-524, continuous-charge support +
// intensity-ratio + the LOW-CHARGE isotope-presence path + two-level harmonic rejection),
// filterMassBins_ (:529-633, top-3 per m/z peak + per-mass charge range). Replaces the old
// FindCandidateMasses, which used continuous-charge support only (so it missed small charge-1/2
// species — the per-spectrum-vs-real low-mass gap) and lacked faithful harmonic rejection.
//
// Output = candidate (neutral mass, min/max abs charge); the Plan B PeakGroup pipeline then recruits
// + scores. Differential-tested vs an MSVC extract (cand_cpp) in MetaFlashDeconDiffTest.

using System;
using System.Collections.Generic;
using Chemistry; // Constants.C13MinusC12

namespace MassSpectrometry
{
    using LogMzPeak = MetaFlashDeconAlgorithm.LogMzPeak;
    using CandidateMass = MetaFlashDeconAlgorithm.CandidateMass;

    internal static class MetaFlashDeconCandidateFinder
    {
        // OpenMS getBinNumber_ : round-half-up ((value-min)*mul + 0.5), 0 if value<min.
        internal static int GetBinNumber(double value, double minValue, double binMulFactor)
            => value < minValue ? 0 : (int)((value - minValue) * binMulFactor + 0.5);

        // OpenMS getBinValue_
        internal static double GetBinValue(int bin, double minValue, double binMulFactor)
            => minValue + bin / binMulFactor;

        // OpenMS setFilters_ (FLASHDeconvAlgorithm.cpp:197-222): universal filter -log(i+1) + the
        // harmonic matrix; charges are 1..maxCharge (abs charge = i+1).
        internal static void SetFilters(int maxCharge, int[] harmonicCharges, out double[] filter, out double[][] harmonicFilter)
        {
            filter = new double[maxCharge];
            for (int i = 0; i < maxCharge; i++) filter[i] = -Math.Log(i + 1);
            harmonicFilter = new double[harmonicCharges.Length][];
            for (int k = 0; k < harmonicCharges.Length; k++)
            {
                int hc = harmonicCharges[k];
                int nn = hc / 2;
                harmonicFilter[k] = new double[maxCharge];
                for (int i = 0; i < maxCharge; i++)
                {
                    double a = i > 0 ? Math.Exp(-filter[i - 1]) : 0;
                    double b = Math.Exp(-filter[i]);
                    harmonicFilter[k][i] = -Math.Log(b - (b - a) * nn / hc);
                }
            }
        }

        /// <summary>
        /// Faithful candidate generation (OpenMS charges 1..<paramref name="maxCharge"/>). Builds the
        /// universal/harmonic filters internally (setFilters_), votes the candidate mass bins, and
        /// returns candidate neutral masses + abs-charge ranges.
        /// </summary>
        internal static List<CandidateMass> FindCandidates(
            List<LogMzPeak> logPeaks, int maxCharge, int[] harmonicCharges,
            double binMulFactor, MetaFlashDeconParameters p, MetaFlashDeconAveragine avg)
        {
            var result = new List<CandidateMass>();
            if (logPeaks.Count == 0) return result;

            SetFilters(maxCharge, harmonicCharges, out double[] filter, out double[][] harmonicFilter);
            int chargeRange = filter.Length;                 // OpenMS current_max_charge_
            int lowCharge = 10;                               // OpenMS low_charge_
            int minSupportPeakCount = 2;                      // OpenMS min_support_peak_count_
            double isoDa = Constants.C13MinusC12;
            int hChargeSize = harmonicCharges.Length;

            double minMass = p.MinMassRange, maxMass = p.MaxMassRange;
            int tmpPeakCntr = Math.Max(0, chargeRange - minSupportPeakCount);

            double mzBinMin = logPeaks[0].LogMz;
            double mzBinMax = logPeaks[^1].LogMz;
            double massBinMin = Math.Log(Math.Max(1.0, minMass - avg.GetAverageMassDelta(minMass)));
            double massBinMax = Math.Min(
                logPeaks[^1].LogMz - filter[tmpPeakCntr],
                Math.Log(maxMass + avg.GetRightCountFromApex(maxMass) + 1.0));

            int massBinNumber = GetBinNumber(massBinMax, massBinMin, binMulFactor) + 1;
            int mzBinNumber = GetBinNumber(mzBinMax, mzBinMin, binMulFactor) + 1;
            if (massBinNumber <= 0 || mzBinNumber <= 0) return result;

            // bin offsets (mz-bin -> mass-bin); C++ round() = away-from-zero
            var binOffsets = new int[chargeRange];
            for (int i = 0; i < chargeRange; i++)
                binOffsets[i] = (int)Math.Round((mzBinMin - filter[i] - massBinMin) * binMulFactor, MidpointRounding.AwayFromZero);
            var harmonicBinOffsets = new int[hChargeSize][];
            for (int k = 0; k < hChargeSize; k++)
            {
                harmonicBinOffsets[k] = new int[chargeRange];
                for (int i = 0; i < chargeRange; i++)
                    harmonicBinOffsets[k][i] = (int)Math.Round((mzBinMin - harmonicFilter[k][i] - massBinMin) * binMulFactor, MidpointRounding.AwayFromZero);
            }

            // updateMzBins_
            var mzBins = new bool[mzBinNumber];
            var mzIntensities = new float[mzBinNumber];
            var mzSetBins = new List<int>();
            foreach (var pk in logPeaks)
            {
                int bi = GetBinNumber(pk.LogMz, mzBinMin, binMulFactor);
                if (bi >= mzBinNumber) continue;
                if (!mzBins[bi]) mzSetBins.Add(bi);
                mzBins[bi] = true;
                mzIntensities[bi] += (float)pk.Intensity;
            }
            mzSetBins.Sort();

            var massBins = new bool[massBinNumber];
            var massIntensities = new float[massBinNumber];

            UpdateCandidateMassBins(mzSetBins, mzBins, mzIntensities, massBins, massIntensities,
                binOffsets, harmonicBinOffsets, harmonicCharges, chargeRange, lowCharge, minSupportPeakCount,
                isoDa, mzBinMin, binMulFactor, massBinNumber);

            var chargeRanges = FilterMassBins(mzSetBins, mzBins, massBins, massIntensities, binOffsets, chargeRange, massBinNumber);

            // getCandidatePeakGroups_ : per surviving mass bin, build the group (local-max guard +
            // harmonic-intensity tracking), apply the group-harmonic gate + >=2-isotope filter, and
            // emit the refined-monoisotopic-mass candidate. (Precision filters.)
            return GetCandidatePeakGroups(logPeaks, massBins, chargeRanges, binOffsets, harmonicCharges,
                mzBinMin, massBinMin, binMulFactor, isoDa, p.DeconvolutionTolerancePpm * 1e-6,
                chargeRange, minMass, maxMass, avg);
        }

        // OpenMS getCandidatePeakGroups_ (FLASHDeconvAlgorithm.cpp:648-959): build + filter candidate
        // groups. Returns refined-mono-mass candidates (CandidateMass.Mass = MONOISOTOPIC mass).
        private static List<CandidateMass> GetCandidatePeakGroups(
            List<LogMzPeak> logPeaks, bool[] massBins, int[][] chargeRanges, int[] binOffsets, int[] harmonicCharges,
            double mzBinMin, double massBinMin, double binMulFactor, double isoDa, double tol,
            int chargeRange, double minMass, double maxMass, MetaFlashDeconAveragine avg)
        {
            var result = new List<CandidateMass>();
            int n = logPeaks.Count;
            int hSize = harmonicCharges.Length;
            const double maxMassDaltonTol = 0.16;

            // per-peak m/z bin number (peaks are logMz-ascending)
            var peakBin = new int[n];
            for (int i = 0; i < n; i++) peakBin[i] = GetBinNumber(logPeaks[i].LogMz, mzBinMin, binMulFactor);

            var currentPeakIndex = new int[chargeRange]; // cpi per charge (monotone as mass bins increase)
            var totalHarmonic = new double[hSize];
            var hPrevIso = new int[hSize];
            var hMaxIsoIntensity = new float[hSize];

            for (int massBinIndex = 0; massBinIndex < massBins.Length; massBinIndex++)
            {
                if (!massBins[massBinIndex]) continue;
                int crLo = chargeRanges[0][massBinIndex], crHi = chargeRanges[1][massBinIndex];
                if (crLo == int.MaxValue || crHi == int.MinValue) continue;

                double mass0 = Math.Exp(GetBinValue(massBinIndex, massBinMin, binMulFactor));
                int rightIndex = avg.GetRightCountFromApex(mass0);
                int leftIndex = avg.GetLeftCountFromApex(mass0);

                var peaks = new List<(double mz, double intensity, int charge, int iso)>();
                double totalSignal = 0;
                Array.Clear(totalHarmonic, 0, hSize);

                for (int j = crLo; j <= crHi; j++)
                {
                    int absCharge = j + 1;
                    int binOffset = binOffsets[j];
                    if (massBinIndex < binOffset) continue;
                    int bIndex = massBinIndex - binOffset;

                    int maxPeakIndex = -1;
                    double maxIntensity = -1;
                    ref int cpi = ref currentPeakIndex[j];
                    while (cpi < n - 1)
                    {
                        if (peakBin[cpi] == bIndex)
                        {
                            double inten = logPeaks[cpi].Intensity;
                            if (inten > maxIntensity) { maxIntensity = inten; maxPeakIndex = cpi; }
                        }
                        else if (peakBin[cpi] > bIndex) break;
                        cpi++;
                    }
                    if (maxPeakIndex < 0) continue;

                    // local-max guard
                    if (maxPeakIndex > 0 && peakBin[maxPeakIndex - 1] == bIndex - 1 && logPeaks[maxPeakIndex - 1].Intensity > maxIntensity) continue;
                    if (maxPeakIndex < n - 1 && peakBin[maxPeakIndex + 1] == bIndex + 1 && logPeaks[maxPeakIndex + 1].Intensity > maxIntensity) continue;

                    double mz = logPeaks[maxPeakIndex].Mz;
                    double isoDelta = isoDa / absCharge;
                    double mzDelta = Math.Min(maxMassDaltonTol / absCharge, 2.0 * tol * mz);
                    double maxMz = mz;
                    float maxPeakIntensity = (float)logPeaks[maxPeakIndex].Intensity;
                    float maxIsoIntensity = 0;
                    int prevIso = -1000;
                    Array.Clear(hPrevIso, 0, hSize);
                    Array.Clear(hMaxIsoIntensity, 0, hSize);

                    // forward
                    for (int pi = maxPeakIndex; pi < n; pi++)
                    {
                        double obsMz = logPeaks[pi].Mz; float inten = (float)logPeaks[pi].Intensity;
                        double mzDiff = obsMz - mz; int ti = (int)Math.Round(mzDiff / isoDelta, MidpointRounding.AwayFromZero);
                        if (obsMz - maxMz > rightIndex * isoDelta + mzDelta) break;
                        if (Math.Abs(mzDiff - ti * isoDelta) < mzDelta)
                        {
                            // OpenMS only counts a signal peak if its mass bin is valid
                            // (FLASHDeconvAlgorithm.cpp:782-783: bin < mass_bin_size && not previously
                            // deconved). Omitting this inflated total_signal -> harmonic gate too lenient.
                            int sbin = peakBin[pi] + binOffset;
                            if (sbin >= 0 && sbin < massBins.Length)
                            {
                                peaks.Add((obsMz, inten, absCharge, ti));
                                if (maxPeakIntensity < inten) maxPeakIntensity = inten;
                                if (prevIso != ti) { totalSignal += maxIsoIntensity; maxIsoIntensity = 0; }
                                maxIsoIntensity = Math.Max(maxIsoIntensity, inten);
                                prevIso = ti;
                            }
                        }
                        else
                        {
                            for (int l = 0; l < hSize; l++)
                            {
                                int hc = harmonicCharges[l]; if (hc * absCharge > chargeRange) break;
                                double hIsoDelta = isoDelta / hc; int thi = (int)Math.Round(mzDiff / hIsoDelta, MidpointRounding.AwayFromZero);
                                if ((double)thi / hc < ti + maxMassDaltonTol) continue;
                                if ((double)thi / hc >= ti + 1 - maxMassDaltonTol) break;
                                if (Math.Abs(mzDiff - thi * hIsoDelta) < mzDelta)
                                {
                                    if (hPrevIso[l] != thi / hc) { totalHarmonic[l] += Math.Min(maxPeakIntensity, hMaxIsoIntensity[l]); hMaxIsoIntensity[l] = 0; }
                                    hMaxIsoIntensity[l] = Math.Max(hMaxIsoIntensity[l], inten); hPrevIso[l] = thi / hc;
                                }
                            }
                        }
                    }
                    totalSignal += maxIsoIntensity;
                    // OpenMS forward-final harmonic add is RAW h_max (FLASHDeconvAlgorithm.cpp:842),
                    // NOT min(max_peak, h_max). (The in-loop add uses min; the final uses raw.)
                    for (int l = 0; l < hSize; l++) totalHarmonic[l] += hMaxIsoIntensity[l];

                    // backward
                    maxIsoIntensity = 0; prevIso = -1000; Array.Clear(hPrevIso, 0, hSize); Array.Clear(hMaxIsoIntensity, 0, hSize);
                    for (int pi = maxPeakIndex - 1; pi >= 0; pi--)
                    {
                        double obsMz = logPeaks[pi].Mz; float inten = (float)logPeaks[pi].Intensity;
                        double mzDiff = obsMz - mz; int ti = (int)Math.Round(mzDiff / isoDelta, MidpointRounding.AwayFromZero);
                        if (maxMz - obsMz > leftIndex * isoDelta + mzDelta) break;
                        if (Math.Abs(mzDiff - ti * isoDelta) < mzDelta)
                        {
                            int sbin = peakBin[pi] + binOffset; // OpenMS bin<mass_bin_size (cpp:863-864)
                            if (sbin >= 0 && sbin < massBins.Length)
                            {
                                peaks.Add((obsMz, inten, absCharge, ti));
                                if (maxPeakIntensity < inten) maxPeakIntensity = inten;
                                if (prevIso != ti) { totalSignal += maxIsoIntensity; maxIsoIntensity = 0; }
                                maxIsoIntensity = Math.Max(maxIsoIntensity, inten);
                                prevIso = ti;
                            }
                        }
                        else
                        {
                            for (int l = 0; l < hSize; l++)
                            {
                                int hc = harmonicCharges[l]; if (hc * absCharge > chargeRange) break;
                                double hIsoDelta = isoDelta / hc; int thi = (int)Math.Round(mzDiff / hIsoDelta, MidpointRounding.AwayFromZero);
                                if ((double)thi / hc > ti - maxMassDaltonTol) continue;
                                if ((double)thi / hc <= ti - 1 + maxMassDaltonTol) break;
                                if (Math.Abs(mzDiff - thi * hIsoDelta) < mzDelta)
                                {
                                    if (hPrevIso[l] != thi / hc) { totalHarmonic[l] += hMaxIsoIntensity[l]; hMaxIsoIntensity[l] = 0; }
                                    hMaxIsoIntensity[l] = Math.Max(hMaxIsoIntensity[l], inten); hPrevIso[l] = thi / hc;
                                }
                            }
                        }
                    }
                    totalSignal += maxIsoIntensity;
                    // OpenMS backward-final harmonic add is min(max_peak, h_max) (FLASHDeconvAlgorithm.cpp:925),
                    // NOT raw. (The backward in-loop add uses raw; the final uses min — mirror of forward.)
                    for (int l = 0; l < hSize; l++) totalHarmonic[l] += Math.Min(maxPeakIntensity, hMaxIsoIntensity[l]);
                }

                // group-harmonic gate
                if (peaks.Count == 0) continue;
                double maxHarm = 0; for (int l = 0; l < hSize; l++) if (totalHarmonic[l] > maxHarm) maxHarm = totalHarmonic[l];
                if (!(totalSignal > maxHarm)) continue;

                // re-anchor isotopes to the most intense peak; >=2-isotope; refined mono mass
                double tMass = 0, mi = -1;
                foreach (var pk in peaks) { double um = (pk.mz - Constants.ProtonMass) * pk.charge; if (pk.intensity > mi) { mi = pk.intensity; tMass = um; } }
                double isoTol = tol * tMass;
                int apexIndex = avg.GetApexIndex(tMass);
                int minOff = 10000, maxOff = -1;
                double nominator = 0, intensity = 0;
                foreach (var pk in peaks)
                {
                    double um = (pk.mz - Constants.ProtonMass) * pk.charge;
                    int iso = (int)Math.Round((um - tMass) / isoDa, MidpointRounding.AwayFromZero);
                    if (Math.Abs(tMass - um + isoDa * iso) > isoTol) continue;
                    iso += apexIndex;
                    if (iso < minOff) minOff = iso;
                    if (iso > maxOff) maxOff = iso;
                    // OpenMS updateMonoMassAndIsotopeIntensities uses the FINAL (apex-anchored) isotope
                    // index here (PeakGroup.cpp:713): mono = Σ pi·(unchargedMass − finalIso·isoDa)/Σpi.
                    nominator += pk.intensity * (um - iso * isoDa);
                    intensity += pk.intensity;
                }
                if (minOff == maxOff || intensity <= 0) continue;
                double monoMass = nominator / intensity;
                if (monoMass < minMass || monoMass > maxMass) continue;

                result.Add(new CandidateMass
                {
                    Mass = monoMass,
                    LogMass = Math.Log(monoMass),
                    // OpenMS getCandidatePeakGroups_ (FLASHDeconvAlgorithm.cpp:678) constructs the
                    // candidate PeakGroup as PeakGroup(1, per_mass_abs_charge_ranges(1,bin)+1, ...):
                    // the MIN abs charge is hardcoded to 1 (NOT crLo+1). The recruit then scans 1..max
                    // and updateChargeRange_ narrows to the real signal range. Using crLo+1 collapsed
                    // candidates to narrow/single-charge -> scoring recruited too few charges, formed
                    // wrong envelopes, dropped the true masses and produced halos.
                    MinAbsCharge = 1,
                    MaxAbsCharge = crHi + 1,
                    SupportIntensity = (float)intensity,
                });
            }
            return result;
        }

        // OpenMS updateCandidateMassBins_ (FLASHDeconvAlgorithm.cpp:292-524) — verbatim logic.
        private static void UpdateCandidateMassBins(
            List<int> mzSetBins, bool[] mzBins, float[] mzIntensities, bool[] massBins, float[] massIntensities,
            int[] binOffsets, int[][] harmonicBinOffsets, int[] harmonicCharges, int chargeRange, int lowCharge,
            int minSupportPeakCount, double isoDa, double mzBinMin, double binMulFactor, int binEnd)
        {
            int n = massBins.Length;
            var supportPeakCount = new ushort[n];
            var prevCharges = new ushort[n];
            for (int i = 0; i < n; i++) prevCharges[i] = (ushort)(chargeRange + 2);
            var prevIntensities = new float[n];
            for (int i = 0; i < n; i++) prevIntensities[i] = 1.0f;
            int hChargeSize = harmonicCharges.Length;
            var subMaxHIntensity = new float[hChargeSize];
            int mzSize = mzBins.Length;

            // iterate set mz bins from high to low (reverse) so charges count small->large
            for (int idx = mzSetBins.Count - 1; idx >= 0; idx--)
            {
                int mzBinIndex = mzSetBins[idx];
                float intensity = mzIntensities[mzBinIndex];
                double logMz = GetBinValue(mzBinIndex, mzBinMin, binMulFactor);
                double mz = Math.Exp(logMz);

                for (int j = 0; j < chargeRange; j++)
                {
                    long massBinIndex = (long)mzBinIndex + binOffsets[j];
                    if (massBinIndex < 0) continue;
                    if (massBinIndex >= binEnd) break;

                    ref ushort spc = ref supportPeakCount[massBinIndex];
                    int absCharge = j + 1;
                    ref float prevIntensity = ref prevIntensities[massBinIndex];
                    ref ushort prevCharge = ref prevCharges[massBinIndex];
                    bool chargeNotContinuous = prevCharge - j != -1 && prevCharge <= chargeRange;

                    float factor = absCharge <= lowCharge ? 10.0f : (5.0f + 5.0f * lowCharge / (float)absCharge);
                    float hfactor = factor / 2.0f;
                    float intensityRatio = intensity / prevIntensity;
                    intensityRatio = intensityRatio < 1 ? 1.0f / intensityRatio : intensityRatio;
                    float supportPeakIntensity = 0;
                    bool passFirstCheck = false;

                    if (chargeNotContinuous || intensityRatio > factor)
                    {
                        spc = 0;
                    }
                    else
                    {
                        passFirstCheck = true;
                        if (spc == 0 && absCharge > lowCharge) supportPeakIntensity = prevIntensity;
                    }

                    float maxHIntensity = 0;
                    if (!passFirstCheck && absCharge <= lowCharge)
                    {
                        Array.Clear(subMaxHIntensity, 0, hChargeSize);
                        for (int d = 1; d >= -1; d -= 2)
                        {
                            bool isoExist = false;
                            double diff = d * isoDa / absCharge / mz;
                            int nextIsoBin = 0;
                            for (int t = -1; t < 2; t++)
                            {
                                int nib = GetBinNumber(logMz + diff, mzBinMin, binMulFactor) + t;
                                if (nib != mzBinIndex && nib > 0 && nib < mzSize && mzBins[nib])
                                {
                                    isoExist = true;
                                    passFirstCheck = true;
                                    if (nextIsoBin == 0 || mzIntensities[nextIsoBin] < mzIntensities[nib]) nextIsoBin = nib;
                                }
                            }
                            if (isoExist)
                            {
                                double hThreshold = intensity + mzIntensities[nextIsoBin];
                                for (int k = 0; k < hChargeSize; k++)
                                {
                                    int hc = harmonicCharges[k];
                                    int harmonicCntr = 0;
                                    if (hc * absCharge > chargeRange) break;
                                    double hdiff = diff / hc;
                                    for (int t = -1; t < 2; t++)
                                    {
                                        int nhib = GetBinNumber(logMz + hdiff, mzBinMin, binMulFactor) + t;
                                        if (nhib != mzBinIndex && nhib >= 0 && nhib < mzSize && mzBins[nhib]
                                            && mzIntensities[nhib] > hThreshold / 2 && mzIntensities[nhib] < hThreshold * 2)
                                        {
                                            harmonicCntr++;
                                            subMaxHIntensity[k] = Math.Max(subMaxHIntensity[k], mzIntensities[nhib]);
                                        }
                                    }
                                    if (harmonicCntr > 0) passFirstCheck = false;
                                }
                                if (passFirstCheck) supportPeakIntensity += mzIntensities[nextIsoBin];
                            }
                        }
                        maxHIntensity = Max(subMaxHIntensity, hChargeSize);
                        passFirstCheck &= maxHIntensity <= 0;
                    }

                    if (passFirstCheck)
                    {
                        if (prevCharge - j == -1) // high-charge continuous path: harmonic artifact check
                        {
                            float maxIntensity = intensity, minIntensity = prevIntensity;
                            if (prevIntensity <= 1.0) { maxIntensity = intensity; minIntensity = intensity; }
                            else if (minIntensity > maxIntensity) { (minIntensity, maxIntensity) = (maxIntensity, minIntensity); }

                            float highThreshold = maxIntensity * hfactor;
                            float lowThreshold = minIntensity / hfactor;
                            bool isHarmonic = false;
                            for (int k = 0; k < hChargeSize; k++)
                            {
                                for (int t = -1; t < 2; t++)
                                {
                                    long hmz = massBinIndex - harmonicBinOffsets[k][j] + t;
                                    if (hmz > 0 && hmz != mzBinIndex && hmz < mzSize && mzBins[hmz])
                                    {
                                        float hi = mzIntensities[hmz];
                                        if (hi > lowThreshold && hi < highThreshold) { maxHIntensity = Math.Max(maxHIntensity, hi); isHarmonic = true; }
                                    }
                                }
                            }
                            if (!isHarmonic)
                            {
                                massIntensities[massBinIndex] += intensity + supportPeakIntensity;
                                spc++;
                                if (spc >= minSupportPeakCount || spc >= absCharge / 2) massBins[massBinIndex] = true;
                            }
                            else
                            {
                                massIntensities[massBinIndex] -= maxHIntensity;
                                if (spc > 0) spc--;
                            }
                        }
                        else if (absCharge <= lowCharge) // low charge: include if isotope present
                        {
                            massIntensities[massBinIndex] += intensity + supportPeakIntensity;
                            spc++;
                            massBins[massBinIndex] = true;
                        }
                    }
                    prevIntensity = intensity;
                    prevCharge = (ushort)j;
                }
            }
        }

        // OpenMS filterMassBins_ (FLASHDeconvAlgorithm.cpp:529-633) — top-3 per m/z peak + charge ranges.
        private static int[][] FilterMassBins(
            List<int> mzSetBins, bool[] mzBins, bool[] massBins, float[] massIntensities,
            int[] binOffsets, int chargeRange, int binSize)
        {
            var ranges = new int[2][];
            ranges[0] = new int[binSize]; ranges[1] = new int[binSize];
            for (int i = 0; i < binSize; i++) { ranges[0][i] = int.MaxValue; ranges[1][i] = int.MinValue; }

            // to_skip = !massBins (mass bins NOT set by candidate stage are skipped); then reset massBins.
            var toSkip = new bool[binSize];
            for (int i = 0; i < binSize; i++) { toSkip[i] = !massBins[i]; massBins[i] = false; }

            const int selectTopN = 3;
            var maxIndices = new long[selectTopN];
            var maxChargeRanges = new int[selectTopN];

            foreach (int mzBinIndex in mzSetBins)
            {
                for (int t = 0; t < selectTopN; t++) { maxIndices[t] = -1; maxChargeRanges[t] = -1; }
                float maxIntensity = -1e11f;
                for (int j = 0; j < chargeRange; j++)
                {
                    long massBinIndex = (long)mzBinIndex + binOffsets[j];
                    if (massBinIndex < 0) continue;
                    if (massBinIndex >= binSize) break;
                    if (toSkip[massBinIndex]) continue;
                    float tval = massIntensities[massBinIndex];
                    if (tval <= 0) continue;
                    if (maxIntensity < tval)
                    {
                        maxIntensity = tval;
                        for (int i = selectTopN - 1; i > 0; i--) { maxIndices[i] = maxIndices[i - 1]; maxChargeRanges[i] = maxChargeRanges[i - 1]; }
                        maxIndices[0] = massBinIndex;
                        maxChargeRanges[0] = j;
                    }
                }
                for (int i = 0; i < selectTopN; i++)
                {
                    long mi = maxIndices[i];
                    int cr = maxChargeRanges[i];
                    if (mi >= 0 && mi < binSize)
                    {
                        ranges[0][mi] = Math.Min(ranges[0][mi], cr);
                        ranges[1][mi] = Math.Max(ranges[1][mi], cr);
                        massBins[mi] = true;
                    }
                }
            }
            return ranges;
        }

        private static float Max(float[] a, int len)
        {
            float m = a[0];
            for (int i = 1; i < len; i++) if (a[i] > m) m = a[i];
            return m;
        }
    }
}
