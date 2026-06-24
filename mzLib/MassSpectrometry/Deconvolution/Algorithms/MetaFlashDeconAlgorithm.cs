// MetaFlashDeconAlgorithm.cs
//
// Status: STEPS 1–5 COMPLETE — full spectral deconvolution implemented
//
// Step 1 — log transform + universal/harmonic pattern construction
// Step 2 — binned candidate mass finding (single-round mass bin computation)
// Step 3 — isotope peak recruitment per candidate per charge
//          + per-charge cosine and noise power collection for Qscoring
// Step 4 — cosine scoring against Averagine; isotope offset optimisation
// Step 5 — IsotopicEnvelope output, deduplication, Qscoring
//
// Reference (paper):
//   Jeong et al. (2020) Cell Systems 10(2):213-218.e6  DOI: 10.1016/j.cels.2020.01.003
//
// Reference (implementation — for faithful reproduction + later fidelity evaluation):
//   OpenMS FLASHDeconv per-spectrum source @ E:\GitClones\OpenMS:
//     src/openms/source/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.cpp  (main pipeline)
//     src/openms/source/ANALYSIS/TOPDOWN/PeakGroup.cpp             (scoring / SNR)
//     src/openms/source/ANALYSIS/TOPDOWN/Qscore.cpp                (logistic Qscore)
//     src/topp/FLASHDeconv.cpp                                     (TOPP defaults)
//   Faithful mapping doc (mzLib step -> OpenMS file:line crosswalk, with the cosine /
//   tol_div_factor / overlap-dedup tolerances we reproduce):
//     E:\CodeReview\MetaFlashDecon\deliverables\reference_flashdecon_faithful_mapping.md
//
// Faithful tolerances (scoped to MetaFlashDecon only; see MetaFlashDeconParameters):
//   - MinCosineScore        = 0.85       (OpenMS min_isotope_cosine MS1; FLASHDeconv.cpp:193)
//   - TolDivFactor          = 2.5        (OpenMS tol_div_factor; FLASHDeconvAlgorithm.cpp:22)
//   - OverlapDedupTolFactor = 1.5        (OpenMS final dedup = 1.5 x ppm; :1223)
//   Candidate voting bins at ppm/2.5; survivors de-duplicated at 1.5 x ppm.

using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;                           // Constants.ProtonMass, ToMz/ToMass

// OpenMS uses Constants::ISOTOPE_MASSDIFF_55K_U = 1.002371 (the average isotope mass difference
// for a 55 kDa peptide, used to bin/anchor isotopes in top-down deconvolution; see
// Constants.h and PeakGroup.h iso_da_distance_). NOT mzLib's Chemistry.Constants.C13MinusC12
// (= 1.0033548, the pure 13C-12C mass diff). Matching this avoids a ~0.001 Da/isotope drift
// that flips borderline isotope assignments in dense spectra.
using MassSpectrometry.MzSpectra;          // SpectralSimilarity.CosineOfAlignedVectors
using MzLibUtil;                           // PpmTolerance

namespace MassSpectrometry
{
    internal class MetaFlashDeconAlgorithm : DeconvolutionAlgorithm
    {
        // ── Algorithm constants ───────────────────────────────────────────────

        // ⚠ CRITICAL CONSTANT — DO NOT change this value or replace it with Chemistry.Constants.C13MinusC12
        // (1.0033548). OpenMS uses Constants::ISOTOPE_MASSDIFF_55K_U (= 1.002371) — the AVERAGE isotope
        // mass difference for a ~55 kDa peptide — as iso_da_distance_ EVERYWHERE in FLASHDeconv
        // (PeakGroup.h:315, FLASHDeconvAlgorithm.cpp:78), NOT the C13−C12 monoisotopic spacing. Using
        // 1.0033548 introduces a ~0.001 Da/isotope drift that flips borderline isotope assignments in
        // dense/high-charge regions and pushes monoisotopic masses off by fractions of a Da — it broke
        // candidate counts and dropped true masses on the densest real scan until this was corrected
        // (commit 3badfcb1: candidates/gates went exact 5840=5840, 32=32). Every iso-offset math in this
        // file, the candidate finder, and the peak group MUST use this single value.
        internal const double IsoDaDistance55K = 1.002371;

        /// <summary>Harmonic denominators (used by the legacy scoring path; candidate gen uses {2,3,5,7,11}).</summary>
        private static readonly int[] HarmonicDenominators = { 2, 3, 5 };

        /// <summary>
        /// Minimum continuous support peaks to mark a mass bin as a candidate.
        /// Note: because the first peak of a run is credited retroactively when
        /// the second arrives, reaching spc==3 means 3 real peaks confirmed.
        /// </summary>
        private const int MinSupportPeakCount = 3;

        /// <summary>
        /// Maximum intensity ratio between adjacent charge states before the
        /// charge run is considered broken and the support count resets.
        /// </summary>
        private const float MaxChargeIntensityRatio = 10.0f;

        // Candidate-mass voting uses bins FINER than the requested tolerance, exactly as
        // OpenMS does: it divides the ppm tolerance by tol_div_factor (=2.5) before computing
        // the bin multiplication factor (FLASHDeconvAlgorithm.cpp:22; updateMembers_:171:
        //   j *= 1e-6; j /= tol_div_factor; bin_mul_factors_.push_back(1.0 / j);).
        // The wider input tolerance is re-applied only at the final overlap dedup
        // (removeOverlappingPeakGroups_, FLASHDeconvAlgorithm.cpp:1223).
        // The single-arg overload (tolDivFactor implicitly 1.0) is retained for callers
        // and tests that compute the bin factor directly at the raw ppm tolerance.
        private static double ComputeBinMulFactor(double tolerancePpm)
            => ComputeBinMulFactor(tolerancePpm, 1.0);

        private static double ComputeBinMulFactor(double tolerancePpm, double tolDivFactor)
        {
            double effectiveTol = tolerancePpm * 1e-6;
            if (tolDivFactor > 0.0) effectiveTol /= tolDivFactor;
            return 1.0 / effectiveTol;
        }

        // ── Constructor ───────────────────────────────────────────────────────
        internal MetaFlashDeconAlgorithm(DeconvolutionParameters deconParameters)
            : base(deconParameters)
        {
        }

        // ══════════════════════════════════════════════════════════════════════
        // ENTRY POINT
        // ══════════════════════════════════════════════════════════════════════

        protected internal override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range)
        {
            var p = DeconvolutionParameters as MetaFlashDeconParameters
                ?? throw new MzLibException("Deconvolution params and algorithm do not match");

            range ??= spectrum.Range;
            if (spectrum.Size == 0) return Enumerable.Empty<IsotopicEnvelope>();

            // ── Step 1 ────────────────────────────────────────────────────────
            var logPeaks = BuildLogMzPeaks(spectrum, range, p.Polarity);
            if (logPeaks.Count == 0) return Enumerable.Empty<IsotopicEnvelope>();

            int maxAbsCharge = Math.Abs(p.MaxAssumedChargeState);
            // Faithful FLASHDeconv: bin the log-m/z axis at ppm / tol_div_factor (finer than the
            // requested tolerance) for candidate voting (FLASHDeconvAlgorithm.cpp:22, :171).
            double binMulFactor = ComputeBinMulFactor(p.DeconvolutionTolerancePpm, p.TolDivFactor);
            var avg = MetaFlashDeconAveragine.For(p.AverageResidueModel, IsoDaDistance55K);

            // ── Step 2 (Plan B): faithful OpenMS candidate generation ─────────
            // updateMzBins_ + updateCandidateMassBins_ (continuous-charge support + the low-charge
            // isotope-presence path + two-level harmonic rejection) + filterMassBins_ (top-3 per m/z
            // peak + charge range). Charges 1..max, harmonic set {2,3,5,7,11}. Replaces the old
            // continuous-charge-only FindCandidateMasses (which missed small charge-1/2 species).
            int[] harmonicCharges = { 2, 3, 5, 7, 11 };
            var candidates = MetaFlashDeconCandidateFinder.FindCandidates(
                logPeaks, maxAbsCharge, harmonicCharges, binMulFactor, p, avg);
            if (candidates.Count == 0) return Enumerable.Empty<IsotopicEnvelope>();

            // ── Steps 3–5 (Plan B): faithful OpenMS PeakGroup pipeline ───────
            // Per candidate, the 10-iteration recruit -> updateQscore offset-refinement loop
            // (MetaFlashDeconPeakGroup) + gates, then removeChargeError -> removeOverlapping, exactly
            // as OpenMS scoreAndFilterPeakGroups_. Each surviving group -> one IsotopicEnvelope whose
            // Score is the OpenMS Qscore. Every per-group function is differential-tested vs the C++.
            var scored = ScoreCandidatesViaPeakGroups(spectrum, candidates, p, avg);

            var afterChargeError = MetaFlashDeconPeakGroup.RemoveChargeErrorPeakGroups(scored, p.Polarity);

            // Final overlap-dedup window = tol × tol_div_factor × 1.5 (FLASHDeconvAlgorithm.cpp:1223):
            // OpenMS passes the RAW tol (=ppm) × tol_div_factor(2.5) × OverlapDedupTolFactor(1.5)
            // = 37.5 ppm at 10 ppm. (Was ppm×1.5 = 15 ppm — 2.5× too narrow, the over-production source.)
            double overlapWindow = p.DeconvolutionTolerancePpm * 1e-6 * p.TolDivFactor * p.OverlapDedupTolFactor;
            var finalGroups = MetaFlashDeconPeakGroup.RemoveOverlappingPeakGroups(afterChargeError, overlapWindow);

            int polSign = Math.Sign((int)p.Polarity);
            var envelopes = new List<IsotopicEnvelope>(finalGroups.Count);
            foreach (var pg in finalGroups)
            {
                var env = new IsotopicEnvelope(
                    id: 0,
                    peaks: pg.SignalPeaks.Select(sp => (sp.Mz, sp.Intensity)).ToList(),
                    monoisotopicmass: pg.MonoisotopicMass,
                    chargestate: polSign * pg.RepAbsCharge,
                    intensity: pg.Intensity,
                    score: pg.QscoreValue);   // Score = OpenMS Qscore
                // Carry the per-isotope vector (OpenMS PeakGroup::getIsotopeIntensities, min-negative
                // isotope-index based) so feature tracing can run the OpenMS MassFeatureTrace
                // isotope-cosine filter. PerIsotopeInt was set by UpdateMonoMassAndIsotopeIntensities.
                env.SetPerIsotopeIntensities(pg.PerIsotopeInt);
                envelopes.Add(env);
            }
            return envelopes;
        }

        // ══════════════════════════════════════════════════════════════════════
        // PLAN B — per-candidate PeakGroup scoring (OpenMS scoreAndFilterPeakGroups_ loop)
        // ══════════════════════════════════════════════════════════════════════
        internal List<MetaFlashDeconPeakGroup> ScoreCandidatesViaPeakGroups(
            MzSpectrum spectrum, List<CandidateMass> candidates, MetaFlashDeconParameters p, MetaFlashDeconAveragine avg)
        {
            var scored = new List<MetaFlashDeconPeakGroup>();
            double isoDa = IsoDaDistance55K;
            // ⚠ CRITICAL: tolerance MUST be narrowed by tol_div_factor here. OpenMS
            // scoreAndFilterPeakGroups_ uses tolerance_[ms_level-1], which updateMembers_ already
            // narrowed (cpp:170-175: j *= 1e-6; j /= tol_div_factor) -> 10 ppm / 2.5 = 4e-6, NOT 1e-5.
            // Do NOT pass the raw `ppm * 1e-6`: the un-narrowed tol made recruitAllPeaksInSpectrum gather
            // ~2.5x too many borderline-isotope peaks, so mis-scored envelopes survived the gates. This
            // single narrowing was the largest fidelity lever on the densest scan (229 -> 16 final masses).
            double tolFraction = p.DeconvolutionTolerancePpm * 1e-6 / p.TolDivFactor;
            const int lowCharge = 10;          // OpenMS low_charge_
            const int minSupportPeakCount = 2; // OpenMS min_support_peak_count_

            foreach (var candidate in candidates)
            {
                double seedMono = candidate.Mass; // getCandidatePeakGroups_ already returns the monoisotopic mass

                var pg = new MetaFlashDeconPeakGroup
                {
                    MinAbsCharge = candidate.MinAbsCharge,
                    MaxAbsCharge = candidate.MaxAbsCharge,
                    IsoDaDistance = isoDa,
                    MinNegativeIsotopeIndex = -1,
                    Polarity = p.Polarity,
                    MonoisotopicMass = seedMono,
                };

                // FD_GATETRACE_CS: which gate (if any) drops each candidate, for a seed-mass window.
                string gtPath = Environment.GetEnvironmentVariable("FD_GATETRACE_CS");
                bool gtOn = gtPath != null && seedMono > 5000 && seedMono < 25000;

                // PRE-LOOP offset seeding (OpenMS scoreAndFilterPeakGroups_, cpp:1090-1101): align
                // the candidate's per-isotope intensity vector against the averagine to find the best
                // initial isotope-index offset, so the 10-iter recruit→updateQscore loop doesn't start
                // anchored to the wrong isotope. Also applies OpenMS's pre-filter `cos < min(.5, min_iso_cos) - .1`.
                int offset = 0;
                if (candidate.PerIsotopeIntensities != null && candidate.PerIsotopeIntensities.Length > 0)
                {
                    double[] bAvg = avg.Get(seedMono);
                    int apexA = avg.GetApexIndex(seedMono);
                    int isoIntShift = -pg.MinNegativeIsotopeIndex; // = 1 (matches OpenMS -peak_group.getMinNegativeIsotopeIndex())
                    // OpenMS min_iso_size_ = 2 (FLASHDeconvAlgorithm.h:140) — applied to both the
                    // input length check and the iso-range check inside getIsotopeCosine.
                    double preCos = MetaFlashDeconPeakGroup.GetIsotopeCosineAndDetermineIsotopeIndex(
                        candidate.PerIsotopeIntensities, bAvg, apexA, isoIntShift, -1, 2, out offset);
                    if (preCos < Math.Min(0.5, p.MinCosineScore) - 0.1)
                    {
                        if (gtOn) System.IO.File.AppendAllText(gtPath, $"GATE\tseed={seedMono:F3}\tz{candidate.MinAbsCharge}-{candidate.MaxAbsCharge}\tdrop=precos\tcos={preCos:F4}\n");
                        continue; // OpenMS cpp:1098
                    }
                    pg.IsotopeCosineScore = preCos;
                }
                double prevMono = seedMono + offset * isoDa; // OpenMS cpp:1093 prev_mono_mass = monoMass + offset * iso_da_distance_
                double iterTraceMono0 = pg.MonoisotopicMass;
                string iterTracePath = Environment.GetEnvironmentVariable("FD_ITERTRACE_CS");
                bool iterTraceOn = iterTracePath != null && iterTraceMono0 > 27880 && iterTraceMono0 < 27900;
                // OpenMS scoreAndFilterPeakGroups_ runs the recruit→updateQscore offset loop up to 10
                // iterations and breaks as soon as the returned offset is 0 (cpp:1103-1124). Keep BOTH
                // the cap of 10 and the `offset == 0` early break — they define the converged mono.
                for (int k = 0; k < 10; k++)
                {
                    double recruitMono = pg.MonoisotopicMass + offset * isoDa;
                    if (recruitMono <= 0) break;
                    int apex = avg.GetApexIndex(recruitMono);
                    int leftCount = avg.GetLeftCountFromApex(recruitMono);
                    int minIso = Math.Max(pg.MinNegativeIsotopeIndex, apex - leftCount + pg.MinNegativeIsotopeIndex);
                    int maxIso = avg.GetLastIndex(recruitMono);
                    var noisy = pg.RecruitAllPeaksInSpectrum(spectrum, tolFraction, recruitMono, minIso, maxIso, p.Polarity);
                    offset = pg.UpdateQscore(noisy, avg, p.MinCosineScore);
                    if (iterTraceOn)
                        System.IO.File.AppendAllText(iterTracePath, $"ITER\t{k}\trecruit={recruitMono:F10}\toffset_ret={offset}\tmono_after={pg.MonoisotopicMass:F10}\n");
                    if (offset == 0) break;
                }

                // FD_PGDUMP_CS: full per-charge recruited-peak dump for a tight seed window.
                if (gtOn && Environment.GetEnvironmentVariable("FD_PGDUMP_CS") is string pgd && seedMono > 8160 && seedMono < 8170)
                {
                    var sbp = new System.Text.StringBuilder();
                    sbp.Append($"PG\tseed={seedMono:F4}\tmono={pg.MonoisotopicMass:F4}\tz{pg.MinAbsCharge}-{pg.MaxAbsCharge}\trep{pg.RepAbsCharge}\tcos={pg.IsotopeCosineScore:F4}\tsnr={pg.Snr:F4}\tq={pg.QscoreValue:F4}\n");
                    foreach (var z in pg.SignalPeaks.GroupBy(sp => sp.AbsCharge).OrderBy(g => g.Key))
                        sbp.Append($"  z{z.Key}: {z.Count()} peaks, iso[{z.Min(x => x.IsotopeIndex)}..{z.Max(x => x.IsotopeIndex)}], sumInt={z.Sum(x => x.Intensity):E2}\n");
                    System.IO.File.AppendAllText(pgd, sbp.ToString());
                }

                // post-gates (OpenMS scoreAndFilterPeakGroups_)
                // ⚠ CRITICAL: this is an EMPTY-ONLY gate. Do NOT change it to
                // `pg.SignalPeaks.Count < p.MinIsotopicPeakCount` (3). OpenMS scoreAndFilterPeakGroups_
                // (cpp:1126) only drops `peak_group.empty()`; a >=3 gate here dropped low-mass /
                // few-isotope truths OpenMS keeps (e.g. 1549.2 z1) — see commit 8af4e4d4. (The faithful
                // floor is 2: a candidate requires min_off != max_off, i.e. >=2 isotopes.)
                string drop = null;
                if (pg.SignalPeaks.Count == 0) drop = "empty";
                else if (pg.MonoisotopicMass < p.MinMassRange || pg.MonoisotopicMass > p.MaxMassRange) drop = "massrange";
                else if (Math.Abs(prevMono - pg.MonoisotopicMass) > 3) drop = "dmono>3";          // moved >3 Da -> different envelope
                else if (pg.MinAbsCharge > lowCharge && (pg.MaxAbsCharge - pg.MinAbsCharge) < minSupportPeakCount) drop = "chargesupport";
                else if (pg.QscoreValue <= 0) drop = "qscore";
                else if (pg.Snr < p.SnrThreshold) drop = "snr";                                   // SNR >= snr_threshold (0.5)
                if (gtOn) System.IO.File.AppendAllText(gtPath, $"GATE\tseed={seedMono:F3}\tmono={pg.MonoisotopicMass:F3}\tz{pg.MinAbsCharge}-{pg.MaxAbsCharge}\tcos={pg.IsotopeCosineScore:F4}\tsnr={pg.Snr:F3}\tq={pg.QscoreValue:F4}\tdrop={drop ?? "KEPT"}\n");
                if (drop != null) continue;

                scored.Add(pg);
            }
            return scored;
        }

        // ══════════════════════════════════════════════════════════════════════
        // STEP 3 — ISOTOPE RECRUITMENT
        // STEP 4 — COSINE SCORING
        // STEP 5 — IsotopicEnvelope OUTPUT
        // ══════════════════════════════════════════════════════════════════════

        /// <summary>
        /// For each candidate mass, recruits isotope peaks from the raw spectrum
        /// for each charge state in the candidate's range, scores the resulting
        /// isotope intensity distribution against Averagine, and returns accepted
        /// envelopes as <see cref="MetaFlashDeconScorer.EnvelopeScoringData"/>.
        ///
        /// During recruitment the following are accumulated simultaneously:
        /// <list type="bullet">
        ///   <item>Global per-isotope intensity vector (all charges combined) for cosine scoring.</item>
        ///   <item>Per-charge isotope intensity vectors for the representative charge.</item>
        ///   <item>Noise peaks — raw spectrum peaks in the same m/z window as each expected
        ///     isotope position but outside the ±ppm tolerance band — for the representative charge.</item>
        /// </list>
        /// </summary>
        private IEnumerable<MetaFlashDeconScorer.EnvelopeScoringData> ScoreAndBuildEnvelopes(
            MzSpectrum spectrum,
            List<CandidateMass> candidates,
            MetaFlashDeconParameters p)
        {
            var results = new List<MetaFlashDeconScorer.EnvelopeScoringData>();
            var tolerance = new PpmTolerance(p.DeconvolutionTolerancePpm);
            int polSign = Math.Sign((int)p.Polarity); // +1 positive, -1 negative

            foreach (var candidate in candidates)
            {
                // ── Step 3a: get Averagine shape for this mass ─────────────────
                // GetAllTheoreticalMasses/Intensities returns arrays sorted
                // INTENSITY-DESCENDING (index 0 = most abundant / apex isotope).
                // We re-sort to MASS-ASCENDING for use in the per-isotope cosine vector.
                //
                // DiffToMonoisotopic[i] = MostIntenseMasses[i] - monoisotopicMass
                //                       = apexMass - monoMass (already computed, exact)
                int avgIdx = p.AverageResidueModel.GetMostIntenseMassIndex(candidate.Mass);
                double apexDaFromMono = p.AverageResidueModel.GetDiffToMonoisotopic(avgIdx);

                double[] avgMasses = p.AverageResidueModel.GetAllTheoreticalMasses(avgIdx);
                double[] avgIntensRaw = p.AverageResidueModel.GetAllTheoreticalIntensities(avgIdx);

                // Sort to mass-ascending order for cosine scoring vector
                var sortedAvg = avgMasses.Zip(avgIntensRaw).OrderBy(pair => pair.First).ToArray();
                double[] avgIntens = sortedAvg.Select(pair => pair.Second).ToArray();
                int nIso = avgIntens.Length;

                // apexIdx in the mass-ascending array — centers the perIsotopeObs vector
                int apexIdx = Array.IndexOf(avgIntens, avgIntens.Max());

                // Representative charge: midpoint of the candidate's charge range.
                // Used for per-charge cosine and noise collection (f[1]–f[3]).
                int repAbsCharge = (candidate.MinAbsCharge + candidate.MaxAbsCharge) / 2;

                // ── Step 3b: recruit peaks per charge state ────────────────────
                int isoWindowSize = nIso + apexIdx;
                var perIsotopeObs = new double[isoWindowSize]; // global (all charges)
                var repChargeObs = new double[isoWindowSize]; // representative charge only

                // recruitedPeaks: (mz, intensity, absCharge, isotopeIndex n from monoisotopic)
                var recruitedPeaks = new List<(double mz, double intensity, int z, int n)>();

                // Faithful global SNR (OpenMS PeakGroup::updateSNR_, PeakGroup.cpp:893-919)
                // aggregates signal / noise / sum-signal-squared across ALL charges; the resulting
                // overall SNR feeds the SNR >= snr_threshold drop gate (getSNR(),
                // FLASHDeconvAlgorithm.cpp:1163). We collect each charge's noisy + signal peaks
                // during recruitment and defer the (relatively expensive) per-charge
                // ComputeNoisePeakPower / getNoisePeakPower_ until AFTER the cosine gate, so noise
                // power is computed only for candidates that survive cosine.
                var perChargeAccum =
                    new List<(List<(double mz, double intensity)> noisy,
                              List<(double mz, double intensity)> signal,
                              double sumSignal, double sumSignalSq, int z, bool isRep)>();

                for (int z = candidate.MinAbsCharge; z <= candidate.MaxAbsCharge; z++)
                {
                    double isoDelta = IsoDaDistance55K / z;
                    double baseMz = candidate.Mass.ToMz(polSign * z);
                    bool isRepZ = (z == repAbsCharge);

                    // Per-charge signal/noise accumulators for the global SNR (reset each charge).
                    var zNoisyPeaks = new List<(double mz, double intensity)>();
                    var zSignalPeaks = new List<(double mz, double intensity)>();
                    double zSumSignal = 0.0;     // Σ signal intensity at z  (-> per_charge_int_[z])
                    double zSumSignalSq = 0.0;   // Σ signal intensity² at z (-> per_charge_sum_signal_squared_[z])

                    // Search forward (n=0 is monoisotopic, n>0 are heavier isotopes)
                    for (int n = 0; n < nIso; n++)
                    {
                        double targetMz = baseMz + n * isoDelta;
                        double windowHalfDa = targetMz * p.DeconvolutionTolerancePpm * 1e-6;
                        double windowLo = targetMz - windowHalfDa;
                        double windowHi = targetMz + windowHalfDa;

                        var signalIndices = spectrum.GetPeakIndicesWithinTolerance(targetMz, tolerance);
                        if (signalIndices.Count == 0)
                        {
                            CollectNoisePeaks(spectrum, windowLo, windowHi,
                                signalIndices, zNoisyPeaks);
                            break;
                        }

                        int bestIdx = signalIndices.OrderByDescending(i => spectrum.YArray[i]).First();
                        double obsIntensity = spectrum.YArray[bestIdx];
                        double obsMz = spectrum.XArray[bestIdx];

                        int isoSlot = apexIdx + n;
                        if (isoSlot < isoWindowSize)
                        {
                            perIsotopeObs[isoSlot] += obsIntensity;
                            if (isRepZ) repChargeObs[isoSlot] += obsIntensity;
                        }

                        recruitedPeaks.Add((obsMz, obsIntensity, z, n));

                        zSumSignal += obsIntensity;
                        zSumSignalSq += obsIntensity * obsIntensity;
                        zSignalPeaks.Add((obsMz, obsIntensity));
                        CollectNoisePeaks(spectrum, windowLo, windowHi,
                            signalIndices, zNoisyPeaks);
                    }

                    // Search backward (lighter than monoisotopic, n=−1,−2,…)
                    for (int n = 1; n <= apexIdx; n++)
                    {
                        double targetMz = baseMz - n * isoDelta;
                        if (targetMz <= 0) break;

                        double windowHalfDa = targetMz * p.DeconvolutionTolerancePpm * 1e-6;
                        double windowLo = targetMz - windowHalfDa;
                        double windowHi = targetMz + windowHalfDa;

                        var signalIndices = spectrum.GetPeakIndicesWithinTolerance(targetMz, tolerance);
                        if (signalIndices.Count == 0)
                        {
                            CollectNoisePeaks(spectrum, windowLo, windowHi,
                                signalIndices, zNoisyPeaks);
                            break;
                        }

                        int bestIdx = signalIndices.OrderByDescending(i => spectrum.YArray[i]).First();
                        double obsIntensity = spectrum.YArray[bestIdx];
                        double obsMz = spectrum.XArray[bestIdx];

                        int isoSlot = apexIdx - n;
                        if (isoSlot >= 0 && isoSlot < isoWindowSize)
                        {
                            perIsotopeObs[isoSlot] += obsIntensity;
                            if (isRepZ) repChargeObs[isoSlot] += obsIntensity;
                        }

                        recruitedPeaks.Add((obsMz, obsIntensity, z, -n));

                        zSumSignal += obsIntensity;
                        zSumSignalSq += obsIntensity * obsIntensity;
                        zSignalPeaks.Add((obsMz, obsIntensity));
                        CollectNoisePeaks(spectrum, windowLo, windowHi,
                            signalIndices, zNoisyPeaks);
                    }

                    // Record this charge's signal/noise for the deferred global-SNR aggregation.
                    perChargeAccum.Add((zNoisyPeaks, zSignalPeaks, zSumSignal, zSumSignalSq, z, isRepZ));
                }

                if (recruitedPeaks.Count < p.MinIsotopicPeakCount)
                    continue;

                // ── Step 4a: find best isotope offset ─────────────────────────
                double bestCosine = -1.0;
                int bestOffset = 0;

                for (int offset = -apexIdx; offset <= apexIdx; offset++)
                {
                    var obsSlice = new double[nIso];
                    for (int i = 0; i < nIso; i++)
                    {
                        int slot = i + offset;
                        if (slot >= 0 && slot < isoWindowSize)
                            obsSlice[i] = perIsotopeObs[slot];
                    }

                    double cos = SpectralSimilarity.CosineOfAlignedVectors(obsSlice, avgIntens);
                    if (cos > bestCosine)
                    {
                        bestCosine = cos;
                        bestOffset = offset;
                    }
                }

                // ── Step 4b: per-charge cosine for representative charge ────────
                // Apply the same offset alignment to the representative-charge-only
                // isotope vector to get chargeIsotopeCosine (used in f[1]).
                var repObsSlice = new double[nIso];
                for (int i = 0; i < nIso; i++)
                {
                    int slot = i + bestOffset;
                    if (slot >= 0 && slot < isoWindowSize)
                        repObsSlice[i] = repChargeObs[slot];
                }
                double repChargeCosine = SpectralSimilarity.CosineOfAlignedVectors(repObsSlice, avgIntens);

                // ── Step 4c: cosine gate ───────────────────────────────────────
                if (bestCosine < p.MinCosineScore) continue;

                // ── Step 4d: faithful all-charge SNR gate ──────────────────────
                // The candidate passed cosine, so now compute each charge's structured-noise
                // power (ComputeNoisePeakPower / getNoisePeakPower_) and aggregate the OpenMS
                // overall SNR across all charges (updateSNR_, PeakGroup.cpp:893-919):
                //   snr = cos² · Σ_z(sigInt_z)²  /  (1 + Σ_z noisePwr_z + (1-cos²)·Σ_z sumSigSq_z)
                // The representative charge's noise / sum-signal-squared are captured here for the
                // Qscore (f[2]/f[3]) -- identical values to before, so Qscores are unchanged. Drop
                // the group when SNR < snr_threshold (getSNR() gate, FLASHDeconvAlgorithm.cpp:1163);
                // this is what discriminates harmonics / noise now that noise power is realistic.
                double globalSignalPower = 0.0;   // Σ_z (per-charge summed signal intensity)²
                double globalNoisePower = 0.0;    // Σ_z noisePwr(z)
                double globalSumSignalSq = 0.0;   // Σ_z Σ(signal intensity²)
                double repChargeNoisePower = 0.0; // representative charge only (Qscore)
                double repChargeSumSignalSq = 0.0;
                foreach (var pc in perChargeAccum)
                {
                    double zNoisePwr = ComputeNoisePeakPower(
                        pc.noisy, pc.signal, pc.z, IsoDaDistance55K);
                    globalSignalPower += pc.sumSignal * pc.sumSignal;
                    globalNoisePower += zNoisePwr;
                    globalSumSignalSq += pc.sumSignalSq;
                    if (pc.isRep)
                    {
                        repChargeNoisePower = zNoisePwr;
                        repChargeSumSignalSq = pc.sumSignalSq;
                    }
                }
                double globalSnr = MetaFlashDeconScorer.ComputeGlobalSnr(
                    bestCosine, globalSignalPower, globalNoisePower, globalSumSignalSq);
                if (globalSnr < p.SnrThreshold) continue;

                // Deduplicate: same raw spectrum peak recruited from multiple charges.
                var uniquePeaks = recruitedPeaks
                    .GroupBy(pk => spectrum.GetClosestPeakIndex(pk.mz))
                    .Select(g => g.OrderByDescending(pk => pk.intensity).First())
                    .ToList();

                if (uniquePeaks.Count < p.MinIsotopicPeakCount) continue;

                // monoMass = candidate.Mass − apexDaFromMono + bestOffset × C13MinusC12
                double monoMass = candidate.Mass - apexDaFromMono + bestOffset * IsoDaDistance55K;
                if (monoMass < p.MinMassRange || monoMass > p.MaxMassRange) continue;

                int signedCharge = polSign * repAbsCharge;
                double totalIntensity = uniquePeaks.Sum(pk => pk.intensity);
                var outputPeaks = uniquePeaks.Select(pk => (pk.mz, pk.intensity)).ToList();

                // ── Compute avg ppm error (isotope index n still known) ────────
                double totalPpmError = 0.0;
                foreach (var (obsMz, _, z, n) in uniquePeaks)
                {
                    double isoDelta = IsoDaDistance55K / z;
                    double theorMz = candidate.Mass.ToMz(polSign * z) + n * isoDelta;
                    totalPpmError += Math.Abs(obsMz - theorMz) / theorMz * 1e6;
                }
                double avgPpmError = uniquePeaks.Count > 0
                    ? totalPpmError / uniquePeaks.Count
                    : 0.0;

                var envelope = new IsotopicEnvelope(
                    id: 0,
                    peaks: outputPeaks,
                    monoisotopicmass: monoMass,
                    chargestate: signedCharge,
                    intensity: totalIntensity,
                    score: bestCosine);   // Score = global cosine; replaced by Qscore after scoring

                results.Add(new MetaFlashDeconScorer.EnvelopeScoringData(
                    envelope: envelope,
                    avgPpmError: avgPpmError,
                    repChargeIsotopeCosine: repChargeCosine,
                    repChargeNoisePower: repChargeNoisePower,
                    repChargeSumSignalSquared: repChargeSumSignalSq));
            }

            return results;
        }

        // -- Noise collection helper --
        /// <summary>
        /// Appends to <paramref name="noisyPeaks"/> every raw spectrum peak within
        /// [windowLo, windowHi] that is NOT in <paramref name="signalIndices"/> (i.e. peaks in
        /// the isotope window but outside the +/-ppm band around the expected isotope position).
        /// These noisy peaks feed <see cref="ComputeNoisePeakPower"/>.
        /// </summary>
        private static void CollectNoisePeaks(
            MzSpectrum spectrum,
            double windowLo,
            double windowHi,
            IReadOnlyCollection<int> signalIndices,
            List<(double mz, double intensity)> noisyPeaks)
        {
            int lo = spectrum.GetClosestPeakIndex(windowLo);
            for (int i = lo; i < spectrum.Size; i++)
            {
                double mz = spectrum.XArray[i];
                if (mz > windowHi) break;
                if (mz < windowLo) continue;
                if (signalIndices.Contains(i)) continue; // signal peak, skip
                noisyPeaks.Add((mz, spectrum.YArray[i]));
            }
        }

        /// <summary>
        /// Faithful port of OpenMS <c>PeakGroup::getNoisePeakPower_</c> (PeakGroup.cpp:180-401).
        /// Estimates the noise power for one charge state by detecting STRUCTURED noise: chains of
        /// peaks (signal + noisy) separated by a CONSISTENT, non-isotopic m/z spacing - the
        /// signature of a co-eluting contaminant ladder rather than the target envelope. Each such
        /// path's noisy intensity is summed and squared; leftover unpaired noisy peaks contribute
        /// their individual power. This is far larger and more realistic than "sum of squared
        /// out-of-band intensity", which is why the prior estimate was near-zero on dense spectra.
        /// <para>
        /// Mirrors the C++ exactly: 50-peak intensity cap, 29 distance bins (24 + 5), the
        /// 0.9-1.1 isotope-distance exclusion, the 0.75 both-signal exclusion, the >=2-edge path
        /// requirement, and the used/skipped power accounting.
        /// </para>
        /// </summary>
        /// <param name="noisyPeaks">Noisy peaks (out-of-band window peaks) for this charge.</param>
        /// <param name="signalPeaks">Recruited signal peaks for this charge.</param>
        /// <param name="absCharge">Absolute charge z (the OpenMS <c>z</c>).</param>
        /// <param name="isoDaDistance">Isotope spacing in Da (C13 - C12 for targets).</param>
        internal static double ComputeNoisePeakPower(
            List<(double mz, double intensity)> noisyPeaks,
            List<(double mz, double intensity)> signalPeaks,
            int absCharge,
            double isoDaDistance)
        {
            const int maxNoisyPeaks = 50; // too many noise peaks slow the process
            const int maxBinNumber = 29;  // 24 bins + 5 extra
            if (absCharge <= 0) return 0.0;
            // OpenMS derives the charge z from the signal peaks; with no signal peaks z stays 0 and
            // getNoisePeakPower_ returns 0 (PeakGroup.cpp:226-229) regardless of noisy peaks.
            if (signalPeaks.Count == 0) return 0.0;

            double threshold = 0.0;
            int total = noisyPeaks.Count + signalPeaks.Count;
            if (total > maxNoisyPeaks)
            {
                // Faithful to OpenMS getNoisePeakPower_ (PeakGroup.cpp:195-202), INCLUDING the quirk:
                // the intensities vector is pre-sized with `total` ZEROS, then ONLY the noisy-peak
                // intensities are appended (signal intensities are NOT included). threshold is then
                // the element at index [size - 50]. (Differential-tested vs the C++ extract.)
                var intens = new List<double>(total + noisyPeaks.Count);
                for (int i = 0; i < total; i++) intens.Add(0.0);
                foreach (var pk in noisyPeaks) intens.Add(pk.intensity);
                intens.Sort();
                threshold = intens[intens.Count - maxNoisyPeaks];
            }

            var allPeaks = new List<(double mz, double intensity, bool isSignal)>(total);
            var signalMzs = new HashSet<double>();
            foreach (var pk in noisyPeaks)
                if (pk.intensity >= threshold) allPeaks.Add((pk.mz, pk.intensity, false));
            foreach (var pk in signalPeaks)
                if (pk.intensity >= threshold) { allPeaks.Add((pk.mz, pk.intensity, true)); signalMzs.Add(pk.mz); }

            allPeaks.Sort((a, b) => a.mz.CompareTo(b.mz));
            int nAll = allPeaks.Count;
            if (nAll == 0) return 0.0;

            var isSignal = new bool[nAll];
            for (int i = 0; i < nAll; i++)
                isSignal[i] = signalMzs.Contains(allPeaks[i].mz);

            var perBinEdges = new int[maxBinNumber][];
            var perBinStartIndex = new int[maxBinNumber];
            for (int k = 0; k < maxBinNumber; k++)
            {
                perBinEdges[k] = new int[nAll];           // 0 = no edge (matches OpenMS sentinel)
                perBinStartIndex[k] = -2;                 // -2 empty, -1 used, >=0 path start
            }

            for (int i = 0; i < nAll; i++)
            {
                bool p1Signal = isSignal[i];
                for (int j = i + 1; j < nAll; j++)
                {
                    double normalizedDist = (allPeaks[j].mz - allPeaks[i].mz) * absCharge / isoDaDistance;
                    if (normalizedDist > 0.9 && normalizedDist < 1.1) continue; // ~ isotope distance: skip
                    // ⚠ MidpointRounding.AwayFromZero REQUIRED: C++ round() is half-away-from-zero; C#'s
                    // default Math.Round is banker's (half-to-even). Dropping this flag was one of the two
                    // real bugs the differential harness caught in the original noise-power port (commit
                    // 786b56d0) — it mis-bins noise peaks and shifts the SNR. Keep the explicit flag.
                    int bin = (int)Math.Round(normalizedDist * (maxBinNumber - 5), MidpointRounding.AwayFromZero);
                    if (bin == 0) continue;
                    if (bin >= maxBinNumber) break;
                    if (p1Signal && isSignal[j] && normalizedDist >= 0.75) continue; // two signals ~1 iso apart
                    perBinEdges[bin][i] = j;
                    perBinStartIndex[bin] = -1;
                }
            }

            var maxIntensitySumToBin = new SortedDictionary<double, int>();
            for (int k = 0; k < maxBinNumber; k++)
            {
                if (perBinStartIndex[k] == -2) continue;
                var edges = perBinEdges[k];
                double maxSumIntensity = 0.0;
                for (int i = 0; i < edges.Length; i++)
                {
                    if (edges[i] == 0) continue;
                    double intensity = isSignal[i] ? 0.0 : allPeaks[i].intensity;
                    double sumIntensity = intensity;
                    int j = edges[i];
                    int cntr = 0;
                    while (j < edges.Length)
                    {
                        cntr++;
                        j = edges[j];
                        if (j <= 0) break;
                        sumIntensity += intensity;
                        if (!isSignal[j]) intensity = allPeaks[j].intensity;
                    }
                    if (cntr > 2 && maxSumIntensity < sumIntensity)
                    {
                        maxSumIntensity = sumIntensity;
                        perBinStartIndex[k] = i;
                    }
                }
                maxIntensitySumToBin[maxSumIntensity] = k; // later duplicate sums overwrite, as in OpenMS
            }

            double chargeNoisePwr = 0.0;
            var used = new bool[nAll];

            foreach (var kv in maxIntensitySumToBin.Reverse())
            {
                int bin = kv.Value;
                int index = perBinStartIndex[bin];
                if (index < 0) continue;

                var edges = perBinEdges[bin];
                double intensity = isSignal[index] ? 0.0 : allPeaks[index].intensity;
                double sumIntensity = 0.0, sumSquaredIntensity = 0.0;
                int skippedPeakCntr = 0;

                if (!used[index]) sumIntensity += intensity;
                else { sumSquaredIntensity += intensity * intensity; skippedPeakCntr++; }
                used[index] = true;

                int j = edges[index];
                while (j < edges.Length)
                {
                    j = edges[j];
                    if (j <= 0) break;
                    if (!used[j]) sumIntensity += intensity;
                    else { sumSquaredIntensity += intensity * intensity; skippedPeakCntr++; }
                    used[j] = true;
                    if (!isSignal[j]) intensity = allPeaks[j].intensity;
                }

                chargeNoisePwr += skippedPeakCntr < 2
                    ? sumIntensity * sumIntensity
                    : sumSquaredIntensity;
            }

            for (int i = 0; i < nAll; i++)
            {
                if (used[i] || isSignal[i]) continue;
                chargeNoisePwr += allPeaks[i].intensity * allPeaks[i].intensity;
            }
            return chargeNoisePwr;
        }

        // ══════════════════════════════════════════════════════════════════════
        // STEP 2 — CANDIDATE MASS FINDING (unchanged)
        // ══════════════════════════════════════════════════════════════════════

        internal sealed class CandidateMass
        {
            public double Mass { get; init; }
            public double LogMass { get; init; }
            public int MinAbsCharge { get; init; }
            public int MaxAbsCharge { get; init; }
            public float SupportIntensity { get; init; }
            /// <summary>
            /// Per-isotope intensity vector from candidate-gen, indexed by `isotopeIndex -
            /// MinNegativeIsotopeIndex` (same shape as <see cref="MetaFlashDeconPeakGroup.PerIsotopeInt"/>,
            /// i.e. index 0 = first allowed negative isotope, isoIntShift = 1). Used by the pre-loop
            /// `GetIsotopeCosineAndDetermineIsotopeIndex` call in scoreAndFilterPeakGroups_
            /// (FLASHDeconvAlgorithm.cpp:1091-1093) to seed the initial isotope-index offset before
            /// the 10-iter recruit→updateQscore loop.
            /// </summary>
            public double[] PerIsotopeIntensities { get; init; }
        }

        internal static List<CandidateMass> FindCandidateMasses(
            List<LogMzPeak> logPeaks,
            double[] universalPattern,
            double[][] harmonicPatterns,
            double binMulFactor,
            MetaFlashDeconParameters deconParams)
        {
            if (logPeaks.Count == 0) return new List<CandidateMass>();

            int chargeRange = universalPattern.Length;
            int minAbsCharge = Math.Abs(deconParams.MinAssumedChargeState);

            double mzBinMinValue = logPeaks[0].LogMz;
            double mzBinMaxValue = logPeaks[^1].LogMz;
            double massBinMinValue = Math.Log(Math.Max(1.0, deconParams.MinMassRange - 50.0));
            double massBinMaxValue = Math.Log(deconParams.MaxMassRange + 200.0);

            int mzBinCount = BinNumber(mzBinMaxValue, mzBinMinValue, binMulFactor) + 1;
            int massBinCount = BinNumber(massBinMaxValue, massBinMinValue, binMulFactor) + 1;

            // Pre-compute integer harmonic bin offsets in mz-bin space.
            // The harmonic check compares mz-bin positions of the current peak against
            // positions that would be occupied by harmonic species at the same mass bin.
            var binHarmonicPatterns = new int[harmonicPatterns.Length][];
            for (int k = 0; k < harmonicPatterns.Length; k++)
            {
                binHarmonicPatterns[k] = new int[chargeRange];
                for (int j = 0; j < chargeRange; j++)
                    binHarmonicPatterns[k][j] = (int)Math.Round(
                        (mzBinMinValue - harmonicPatterns[k][j] - massBinMinValue) * binMulFactor);
            }

            // mz-bin occupancy map — used by the harmonic check
            var mzBinOccupied = new bool[mzBinCount];
            var mzBinIntensity = new float[mzBinCount];
            var peakBinIndex = new int[logPeaks.Count];

            for (int i = 0; i < logPeaks.Count; i++)
            {
                int b = BinNumber(logPeaks[i].LogMz, mzBinMinValue, binMulFactor);
                if (b < 0 || b >= mzBinCount) continue;
                peakBinIndex[i] = b;
                mzBinOccupied[b] = true;
                mzBinIntensity[b] += (float)logPeaks[i].Intensity;
            }

            var massBinOccupied = new bool[massBinCount];
            var massBinIntensity = new float[massBinCount];
            var massBinSupportCount = new int[massBinCount];
            var massBinMinChargeIdx = new int[massBinCount];
            var massBinMaxChargeIdx = new int[massBinCount];
            for (int i = 0; i < massBinCount; i++)
            {
                massBinMinChargeIdx[i] = int.MaxValue;
                massBinMaxChargeIdx[i] = int.MinValue;
            }

            var prevChargeIdx = new int[massBinCount];
            var prevIntensity = new float[massBinCount];
            for (int i = 0; i < massBinCount; i++) prevChargeIdx[i] = chargeRange + 2;

            // Main scan: right-to-left (high logMz → low logMz).
            // Scanning in this direction means charge index j increases as we move
            // to higher logMz, so consecutive charges for the same mass appear as
            // j, j-1 — i.e. prevChargeIdx[massBin] - j == -1 when continuous.
            for (int i = logPeaks.Count - 1; i >= 0; i--)
            {
                int mzBin = peakBinIndex[i];
                double logMz = logPeaks[i].LogMz;
                float intensity = (float)logPeaks[i].Intensity;

                for (int j = 0; j < chargeRange; j++)
                {
                    int massBin = BinNumber(logMz - universalPattern[j], massBinMinValue, binMulFactor);
                    if (massBin < 0 || massBin >= massBinCount) continue;

                    bool chargeIsContinuous = (prevChargeIdx[massBin] - j) == -1
                                             && prevChargeIdx[massBin] <= chargeRange;
                    float ratio = prevIntensity[massBin] <= 0f
                        ? MaxChargeIntensityRatio + 1f
                        : intensity / prevIntensity[massBin];
                    if (ratio < 1f) ratio = 1f / ratio;

                    bool passCheck = chargeIsContinuous && ratio <= MaxChargeIntensityRatio;

                    if (!passCheck)
                    {
                        massBinSupportCount[massBin] = 0;
                    }
                    else
                    {
                        bool isHarmonic = false;
                        float maxHarmonicIntensity = 0f;

                        for (int k = 0; k < binHarmonicPatterns.Length && !isHarmonic; k++)
                        {
                            int hBin = mzBin - binHarmonicPatterns[k][j];
                            if (hBin > 0 && hBin < mzBinCount && hBin != mzBin
                                && mzBinOccupied[hBin])
                            {
                                float hIntensity = mzBinIntensity[hBin];
                                float hRatio = hIntensity > intensity
                                    ? hIntensity / intensity
                                    : intensity / hIntensity;
                                if (hRatio <= MaxChargeIntensityRatio)
                                {
                                    isHarmonic = true;
                                    maxHarmonicIntensity = Math.Max(maxHarmonicIntensity, hIntensity);
                                }
                            }
                        }

                        if (!isHarmonic)
                        {
                            if (massBinSupportCount[massBin] == 0)
                            {
                                // Retroactive credit: this is the second confirming peak,
                                // so we count both this and the previous peak → support = 2.
                                massBinSupportCount[massBin] = 2;
                                massBinIntensity[massBin] += prevIntensity[massBin] + intensity;
                                int prevJ = prevChargeIdx[massBin];
                                if (prevJ < chargeRange)
                                {
                                    if (prevJ < massBinMinChargeIdx[massBin]) massBinMinChargeIdx[massBin] = prevJ;
                                    if (prevJ > massBinMaxChargeIdx[massBin]) massBinMaxChargeIdx[massBin] = prevJ;
                                }
                            }
                            else
                            {
                                massBinSupportCount[massBin]++;
                                massBinIntensity[massBin] += intensity;
                            }

                            if (!massBinOccupied[massBin]
                                && massBinSupportCount[massBin] >= MinSupportPeakCount)
                                massBinOccupied[massBin] = true;

                            if (j < massBinMinChargeIdx[massBin]) massBinMinChargeIdx[massBin] = j;
                            if (j > massBinMaxChargeIdx[massBin]) massBinMaxChargeIdx[massBin] = j;
                        }
                        else
                        {
                            massBinIntensity[massBin] -= maxHarmonicIntensity;
                            if (massBinSupportCount[massBin] > 0) massBinSupportCount[massBin]--;
                        }
                    }

                    prevIntensity[massBin] = intensity;
                    prevChargeIdx[massBin] = j;
                }
            }

            // Conflict resolution — faithful to OpenMS FLASHDeconv filterMassBins_
            // (FLASHDeconvAlgorithm.cpp:543, select_top_N = 3): each m/z peak matches one
            // candidate mass per charge, so it is allowed to "claim" only its TOP-3
            // most-intense occupied mass bins (top-N=3 "to consider frequent coelution"),
            // not just its single best. This bounds how many masses one peak can seed.
            const int SelectTopN = 3;
            var finalMassBins = new bool[massBinCount];
            Span<long> topBins = stackalloc long[SelectTopN];
            Span<float> topInts = stackalloc float[SelectTopN];
            for (int i = 0; i < logPeaks.Count; i++)
            {
                double logMz = logPeaks[i].LogMz;
                for (int t = 0; t < SelectTopN; t++) { topBins[t] = -1; topInts[t] = -1f; }

                for (int j = 0; j < chargeRange; j++)
                {
                    int massBin = BinNumber(logMz - universalPattern[j], massBinMinValue, binMulFactor);
                    if (massBin < 0 || massBin >= massBinCount) continue;
                    if (!massBinOccupied[massBin]) continue;

                    float mbi = massBinIntensity[massBin];
                    if (mbi <= topInts[SelectTopN - 1]) continue;

                    // Insert into the descending-by-intensity top-N (rolling shift).
                    int pos = SelectTopN - 1;
                    while (pos > 0 && topInts[pos - 1] < mbi)
                    {
                        topInts[pos] = topInts[pos - 1];
                        topBins[pos] = topBins[pos - 1];
                        pos--;
                    }
                    topInts[pos] = mbi;
                    topBins[pos] = massBin;
                }

                for (int t = 0; t < SelectTopN; t++)
                    if (topBins[t] >= 0) finalMassBins[topBins[t]] = true;
            }

            var candidates = new List<CandidateMass>();
            for (int b = 0; b < massBinCount; b++)
            {
                if (!finalMassBins[b]) continue;
                if (massBinMinChargeIdx[b] == int.MaxValue || massBinMaxChargeIdx[b] == int.MinValue) continue;

                double logMass = BinValue(b, massBinMinValue, binMulFactor);
                double mass = Math.Exp(logMass);

                if (mass < deconParams.MinMassRange || mass > deconParams.MaxMassRange) continue;

                candidates.Add(new CandidateMass
                {
                    Mass = mass,
                    LogMass = logMass,
                    MinAbsCharge = minAbsCharge + massBinMinChargeIdx[b],
                    MaxAbsCharge = minAbsCharge + massBinMaxChargeIdx[b],
                    SupportIntensity = massBinIntensity[b]
                });
            }

            return candidates;
        }

        // ══════════════════════════════════════════════════════════════════════
        // BIN HELPERS
        // ══════════════════════════════════════════════════════════════════════

        internal static int BinNumber(double value, double minValue, double binMulFactor)
        {
            double raw = (value - minValue) * binMulFactor;
            return raw < 0.0 ? 0 : (int)Math.Round(raw);
        }

        internal static double BinValue(int bin, double minValue, double binMulFactor)
            => minValue + bin / binMulFactor;

        // ══════════════════════════════════════════════════════════════════════
        // STEP 1 HELPERS
        // ══════════════════════════════════════════════════════════════════════

        internal static double GetLogMz(double mz, Polarity polarity)
        {
            double uncharged = mz - Math.Sign((int)polarity) * Constants.ProtonMass;
            return uncharged > 0.0 ? Math.Log(uncharged) : double.NegativeInfinity;
        }

        internal static List<LogMzPeak> BuildLogMzPeaks(
            MzSpectrum spectrum, MzRange range, Polarity polarity)
        {
            var peaks = new List<LogMzPeak>(spectrum.Size);
            for (int i = 0; i < spectrum.Size; i++)
            {
                double mz = spectrum.XArray[i];
                double intensity = spectrum.YArray[i];
                if (mz < range.Minimum || mz > range.Maximum) continue;
                if (intensity <= 0.0) continue;
                double logMz = GetLogMz(mz, polarity);
                if (!double.IsFinite(logMz)) continue;
                peaks.Add(new LogMzPeak(mz, intensity, logMz));
            }
            peaks.Sort((a, b) => a.LogMz.CompareTo(b.LogMz));
            return peaks;
        }

        internal static double[] BuildUniversalPattern(int minAbsCharge, int chargeRange)
        {
            if (chargeRange <= 0) throw new ArgumentOutOfRangeException(nameof(chargeRange), "Must be > 0");
            if (minAbsCharge <= 0) throw new ArgumentOutOfRangeException(nameof(minAbsCharge), "Must be > 0");
            var p = new double[chargeRange];
            for (int j = 0; j < chargeRange; j++) p[j] = -Math.Log(minAbsCharge + j);
            return p;
        }

        internal static double[][] BuildHarmonicPatterns(int minAbsCharge, int chargeRange)
        {
            double[] up = BuildUniversalPattern(minAbsCharge, chargeRange);
            var h = new double[HarmonicDenominators.Length][];
            for (int k = 0; k < HarmonicDenominators.Length; k++)
            {
                int d = HarmonicDenominators[k];
                int n = d / 2;
                h[k] = new double[chargeRange];
                for (int j = 0; j < chargeRange; j++)
                {
                    double b = Math.Exp(-up[j]);
                    double a = j > 0 ? Math.Exp(-up[j - 1]) : 0.0;
                    double hMz = b - (b - a) * (double)n / d;
                    h[k][j] = hMz > 0.0 ? -Math.Log(hMz) : up[j];
                }
            }
            return h;
        }

        // ══════════════════════════════════════════════════════════════════════
        // LogMzPeak STRUCT
        // ══════════════════════════════════════════════════════════════════════

        // Mirror of OpenMS FLASHDeconvHelperStructs::LogMzPeak (double mz, FLOAT intensity, double logMz).
        // ⚠ Intensity is stored as double for interop, but the values are float-origin (mzML 32-bit
        // intensities widened to double), and OpenMS treats intensity as float throughout. The
        // intensity-weighted mono computations therefore cast to (float) and sum the denominator in float
        // ON PURPOSE — see MetaFlashDeconPeakGroup.UpdateMonoMassAndIsotopeIntensities. Do not assume the
        // extra double precision here is meaningful; matching OpenMS's float arithmetic is what's correct.
        internal readonly struct LogMzPeak
        {
            public readonly double Mz;
            public readonly double Intensity;
            public readonly double LogMz;
            public readonly int AbsCharge;
            public readonly int IsotopeIndex;

            internal LogMzPeak(double mz, double intensity, double logMz,
                               int absCharge = 0, int isotopeIndex = 0)
            {
                Mz = mz;
                Intensity = intensity;
                LogMz = logMz;
                AbsCharge = absCharge;
                IsotopeIndex = isotopeIndex;
            }

            public LogMzPeak WithAssignment(int absCharge, int isotopeIndex)
                => new LogMzPeak(Mz, Intensity, LogMz, absCharge, isotopeIndex);

            public override string ToString()
                => $"mz={Mz:F4} logMz={LogMz:F6} I={Intensity:G4} z={AbsCharge} iso={IsotopeIndex}";
        }
    }
}