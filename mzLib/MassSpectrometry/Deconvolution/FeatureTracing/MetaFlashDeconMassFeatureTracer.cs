using System;
using System.Collections.Generic;
using System.Linq;
using MassSpectrometry.Deconvolution.Consensus;

namespace MassSpectrometry.Deconvolution.FeatureTracing
{
    /// <summary>
    /// FLASHDeconv-faithful <see cref="IMassFeatureTracer"/>. Traces features over
    /// <b>neutral mass</b> (charge collapsed per scan), as in OpenMS FLASHDeconv's
    /// <c>MassFeatureTrace</c> step — fundamentally different from
    /// <see cref="ConsensusMassFeatureTracer"/>, which builds charge-locked traces first and
    /// stitches charges afterward. Steps (see deliverables/research_flashdeconv_feature_tracing.md):
    /// <list type="number">
    ///   <item>Collapse each MS1 scan's per-charge envelopes into neutral-mass peaks
    ///   (charge intensities summed) — the per-scan PeakGroup analog.</item>
    ///   <item>Detect mass traces across retention time with an apex-sorted, bidirectional,
    ///   intensity-weighted-mean-centroid extension (the Kenar elution-profile algorithm OpenMS
    ///   reuses), with outlier-run termination and a minimum length / sample-rate filter.</item>
    ///   <item>Assemble each surviving trace into a <see cref="MassFeature"/>, decomposed into
    ///   one <see cref="CorrectedTrace"/> per charge so the existing feature/writer model is reused.</item>
    /// </list>
    /// <para>
    /// Parameter defaults mirror the FLASHDeconv TOPP-tool overrides (10 ppm mass tolerance,
    /// 10 s min trace length, 0.05 min sample rate, 20-scan outlier termination). NOTE:
    /// <paramref name="minTraceLengthRt"/> is in the SAME units as
    /// <see cref="MsDataScan.RetentionTime"/> — minutes for mzML — so the FLASHDeconv 10 s default
    /// is expressed as 10/60 min.
    /// </para>
    /// <para>
    /// Deferred refinements (Phase 3b, validated in Phase 6 vs real FLASHDeconv): the
    /// integer isotope-offset mass correction (averagine cosine on the trace-summed isotope
    /// envelope) and the feature-level isotope-cosine ≥ 0.85 filter. The input envelopes are
    /// already cosine- and Qscore-filtered per scan by <see cref="MetaFlashDeconAlgorithm"/>.
    /// </para>
    /// </summary>
    public sealed class MetaFlashDeconMassFeatureTracer : IMassFeatureTracer
    {
        /// <summary>Mass tolerance (ppm) for collapsing charges into one neutral-mass peak per scan.</summary>
        public const double DefaultMassTolerancePpm = 10.0;
        /// <summary>OpenMS MassTraceDetection mass_error_ppm (NOT overridden by MassFeatureTrace, so the
        /// default 20). The trace-extension window is ±3·(centroid/1e6·mass_error_ppm) = ±60 ppm.</summary>
        public const double DefaultMassErrorPpm = 20.0;
        /// <summary>Minimum trace RT span to keep, in RetentionTime units (minutes); OpenMS min_trace_length
        /// = 10 s (MassFeatureTrace.cpp:19) = 10/60 min.</summary>
        public const double DefaultMinTraceLengthRt = 10.0 / 60.0;
        /// <summary>OpenMS MassFeatureTrace min_sample_rate (MassFeatureTrace.cpp:18) = 0.1.</summary>
        public const double DefaultMinSampleRate = 0.1;
        /// <summary>OpenMS MassTraceDetection trace_termination_outliers (NOT overridden, default 5):
        /// a direction stops after this many CONSECUTIVE non-empty scans with no acceptable peak.</summary>
        public const int DefaultMaxOutlierScans = 5;

        /// <summary>OpenMS MassFeatureTrace min_isotope_cosine for MS1 (MassFeatureTrace.cpp:32). Traces
        /// whose summed isotope pattern has cosine &lt; this vs the averagine are dropped.</summary>
        public const double DefaultMinIsotopeCosine = 0.75;
        // Per-isotope accumulation vector size. >= OpenMS avg.getMaxIsotopeIndex() (=299 at max_mass
        // 100k); the cosine ignores trailing zeros, so any size >= the highest occupied isotope is exact
        // (verified vs massfeature_snip: cos within 2e-3, offset/filter identical).
        private const int IsotopeVectorSize = 400;

        private readonly double _massTolerancePpm;
        private readonly double _massErrorPpm = DefaultMassErrorPpm;
        private readonly double _minTraceLengthRt;
        private readonly double _minSampleRate;
        private readonly int _maxOutlierScans;
        private readonly double _minIsotopeCosine = DefaultMinIsotopeCosine;
        private readonly double _isoDaDistance;
        private readonly MetaFlashDeconAveragine _averagine;

        /// <param name="massTolerancePpm">Mass tolerance (ppm) for charge collapse and trace extension.</param>
        /// <param name="minTraceLengthRt">Minimum trace RT span (RetentionTime units, i.e. minutes for mzML).</param>
        /// <param name="minSampleRate">Minimum fraction of spanned scans carrying a peak.</param>
        /// <param name="maxOutlierScans">Consecutive empty scans that terminate a trace direction.</param>
        public MetaFlashDeconMassFeatureTracer(
            double massTolerancePpm = DefaultMassTolerancePpm,
            double minTraceLengthRt = DefaultMinTraceLengthRt,
            double minSampleRate = DefaultMinSampleRate,
            int maxOutlierScans = DefaultMaxOutlierScans)
        {
            _massTolerancePpm = massTolerancePpm;
            _minTraceLengthRt = minTraceLengthRt;
            _minSampleRate = minSampleRate;
            _maxOutlierScans = maxOutlierScans;
            // Averagine for the OpenMS MassFeatureTrace feature-level isotope-cosine filter. Same model +
            // iso_da_distance the decon uses (FLASHDeconvAlgorithm passes avg_ to MassFeatureTrace).
            _isoDaDistance = MetaFlashDeconAlgorithm.IsoDaDistance55K;
            _averagine = MetaFlashDeconAveragine.For(new Averagine(), _isoDaDistance);
        }

        /// <inheritdoc />
        public List<MassFeature> TraceFeatures(
            IReadOnlyList<MsDataScan> ms1Scans,
            IReadOnlyList<IReadOnlyList<IsotopicEnvelope>> perScanEnvelopes)
        {
            if (perScanEnvelopes.Count != ms1Scans.Count)
                throw new ArgumentException(
                    $"perScanEnvelopes ({perScanEnvelopes.Count}) must align 1:1 with ms1Scans ({ms1Scans.Count}).",
                    nameof(perScanEnvelopes));

            var peaksByScan = CollapseAllScans(ms1Scans, perScanEnvelopes);

            // Step 2 — neutral-mass trace detection.
            var traces = DetectTraces(peaksByScan);

            // Step 3 — assemble surviving traces into MassFeatures (dropping those that fail the
            // OpenMS MassFeatureTrace feature-level isotope-cosine filter).
            var features = new List<MassFeature>(traces.Count);
            int nextId = 1;
            foreach (var trace in traces)
            {
                var feature = BuildFeature(trace, nextId);
                if (feature != null) { features.Add(feature); nextId++; }
            }
            return features;
        }

        private List<MassPeak>[] CollapseAllScans(
            IReadOnlyList<MsDataScan> ms1Scans, IReadOnlyList<IReadOnlyList<IsotopicEnvelope>> perScanEnvelopes)
        {
            var peaksByScan = new List<MassPeak>[ms1Scans.Count];
            for (int s = 0; s < ms1Scans.Count; s++)
                peaksByScan[s] = CollapseToNeutralMassPeaks(
                    s, ms1Scans[s].OneBasedScanNumber, ms1Scans[s].RetentionTime, perScanEnvelopes[s]);
            return peaksByScan;
        }

        // Test hook for the MassTraceDetection differential snip: returns the per-scan collapsed input
        // peaks (mass, intensity, RT in seconds) AND the raw detected traces (pre-cosine-filter) summarised
        // as (size, rtMin/rtMax in seconds, intensity-weighted centroid mass) — so masstrace_snip can be
        // fed identical input and compared trace-for-trace.
        internal (List<(int scanIndex, double rtSeconds, List<(double mass, double intensity)> peaks)> input,
                  List<(int size, double rtMinSec, double rtMaxSec, double centroid)> traces)
            DetectTracesDiagnostic(IReadOnlyList<MsDataScan> ms1Scans, IReadOnlyList<IReadOnlyList<IsotopicEnvelope>> perScanEnvelopes)
        {
            var peaksByScan = CollapseAllScans(ms1Scans, perScanEnvelopes);
            var rawTraces = DetectTraces(peaksByScan);

            var input = new List<(int, double, List<(double, double)>)>();
            for (int s = 0; s < peaksByScan.Length; s++)
            {
                var pks = peaksByScan[s].Select(p => (p.Mass, p.Intensity)).ToList();
                double rtSec = peaksByScan[s].Count > 0 ? peaksByScan[s][0].Rt * 60.0 : ms1Scans[s].RetentionTime * 60.0;
                input.Add((s, rtSec, pks));
            }
            var traceSummaries = new List<(int, double, double, double)>();
            foreach (var t in rawTraces)
            {
                double ti = t.Sum(p => p.Intensity);
                double centroid = ti > 0 ? t.Sum(p => p.Mass * p.Intensity) / ti : t[0].Mass;
                traceSummaries.Add((t.Count, t.Min(p => p.Rt) * 60.0, t.Max(p => p.Rt) * 60.0, centroid));
            }
            return (input, traceSummaries);
        }

        // ── Step 1: per-scan neutral-mass collapse ────────────────────────────
        // Builds the per-RT neutral-mass peak list — OpenMS MassFeatureTrace's
        // storeInformationFromDeconvolvedSpectrum keys one Peak1D(monoMass, intensity) per PeakGroup
        // (MassFeatureTrace.cpp:171-187), i.e. NO ppm merging. A MetaFlashDecon decon emits exactly one
        // envelope per neutral mass per scan (charges already aggregated in the peak group, RepAbsCharge
        // labelled), and RemoveOverlappingPeakGroups keeps survivors apart — so this is effectively 1:1
        // and FAITHFUL to OpenMS on the real path (VERIFIED: 0/912 envelopes merged over 60 real scans,
        // Collapse_MergeStats). The within-ppm merge below is dormant for MetaFlashDecon input; it only
        // fires for non-faithful sources that hand the IMassFeatureTracer multiple same-mass per-charge
        // envelopes (e.g. the synthetic NativeTracer_ tests). Do not rely on it for production fidelity.
        private List<MassPeak> CollapseToNeutralMassPeaks(
            int scanIndex, int scanNumber, double rt, IReadOnlyList<IsotopicEnvelope> envelopes)
        {
            var peaks = new List<MassPeak>();
            if (envelopes.Count == 0) return peaks;

            foreach (var env in envelopes.OrderBy(e => e.MonoisotopicMass))
            {
                var last = peaks.Count > 0 ? peaks[^1] : null;
                double tolDa = last != null ? last.Mass * _massTolerancePpm * 1e-6 : 0;
                if (last != null && Math.Abs(env.MonoisotopicMass - last.Mass) <= tolDa)
                {
                    // Merge into the open group: update the intensity-weighted mean mass.
                    double newInt = last.Intensity + env.TotalIntensity;
                    last.Mass = newInt > 0
                        ? (last.Mass * last.Intensity + env.MonoisotopicMass * env.TotalIntensity) / newInt
                        : last.Mass;
                    last.Intensity = newInt;
                    last.Envelopes.Add(env);
                }
                else
                {
                    peaks.Add(new MassPeak
                    {
                        ScanIndex = scanIndex,
                        ScanNumber = scanNumber,
                        Rt = rt,
                        Mass = env.MonoisotopicMass,
                        Intensity = env.TotalIntensity,
                        Envelopes = new List<IsotopicEnvelope> { env },
                    });
                }
            }
            return peaks;
        }

        // ── Step 2: neutral-mass trace detection — faithful port of OpenMS MassTraceDetection::run_ ──
        // (MassTraceDetection.cpp:449-663, the algorithm MassFeatureTrace runs via mtdet.run). Each
        // collapsed neutral-mass peak is an MSExperiment Peak1D(mono, intensity); we trace over neutral
        // mass exactly as the real tool. ⚠ Faithful details (copy, do not "improve"):
        //   - apices sorted by intensity ASCENDING then iterated in reverse (stable), so the most intense
        //     unused peak seeds first.
        //   - tolerance window is ±3·ftl_sd where ftl_sd = centroid/1e6·mass_error_ppm; with the OpenMS
        //     default mass_error_ppm=20 that is ±60 ppm. ftl_sd is fixed (reestimate_mt_sd=false).
        //   - centroid is the ITERATIVE weighted mean (UpdateIterativeWeightedMean), apex counted twice
        //     at init exactly as OpenMS — the multiplicative form is not algebraically identical to a
        //     naive Σwv/Σw in floating point, so it is replicated verbatim.
        //   - a candidate is the single findNearest peak, accepted only if inside the window AND not yet
        //     visited; outlier (consecutive-missed > trace_termination_outliers=5) terminates a direction;
        //     empty scans advance the index but do NOT increment the missed counter.
        //   - peaks are marked visited ONLY after the trace passes isTraceValid_ (min_trace_length +
        //     min_sample_rate over total_scans − missed_down − missed_up); failed traces leave their peaks
        //     available to other seeds.
        private List<List<MassPeak>> DetectTraces(List<MassPeak>[] peaksByScan)
        {
            int nScans = peaksByScan.Length;
            var visited = new bool[nScans][];
            for (int s = 0; s < nScans; s++) visited[s] = new bool[peaksByScan[s].Count];

            // apices: every peak, stable-sorted by (float) intensity ascending; iterate descending.
            var apices = new List<(int scan, int idx)>();
            for (int s = 0; s < nScans; s++)
                for (int i = 0; i < peaksByScan[s].Count; i++) apices.Add((s, i));
            apices = apices.OrderBy(a => (float)peaksByScan[a.scan][a.idx].Intensity).ToList(); // OrderBy is stable

            var traces = new List<List<MassPeak>>();
            for (int ai = apices.Count - 1; ai >= 0; ai--)
            {
                var (apexScan, apexIdx) = apices[ai];
                if (visited[apexScan][apexIdx]) continue;
                var apex = peaksByScan[apexScan][apexIdx];

                var members = new List<MassPeak> { apex };
                var gathered = new List<(int s, int i)> { (apexScan, apexIdx) };

                double centroid = apex.Mass;
                double counter = (float)apex.Intensity * apex.Mass;
                double denom = (float)apex.Intensity;
                UpdateIterativeWeightedMean(apex.Mass, (float)apex.Intensity, ref centroid, ref counter, ref denom);

                double ftlSd = centroid / 1e6 * _massErrorPpm; // fixed (reestimate_mt_sd = false)
                int downIdx = apexScan, upIdx = apexScan;
                bool downActive = true, upActive = true;
                int missedDown = 0, missedUp = 0, scanDown = 0, scanUp = 0;

                while ((downIdx > 0 && downActive) || (upIdx < nScans - 1 && upActive))
                {
                    if (downIdx > 0 && downActive)
                    {
                        var spec = peaksByScan[downIdx - 1];
                        if (spec.Count > 0)
                        {
                            int ci = FindNearest(spec, centroid);
                            double cmz = spec[ci].Mass;
                            if (cmz <= centroid + 3 * ftlSd && cmz >= centroid - 3 * ftlSd && !visited[downIdx - 1][ci])
                            {
                                members.Insert(0, spec[ci]); // push_front (RT-ascending)
                                UpdateIterativeWeightedMean(cmz, (float)spec[ci].Intensity, ref centroid, ref counter, ref denom);
                                gathered.Add((downIdx - 1, ci));
                                missedDown = 0;
                            }
                            else missedDown++;
                        }
                        downIdx--; scanDown++;
                        if (missedDown > _maxOutlierScans) downActive = false;
                    }
                    if (upIdx < nScans - 1 && upActive)
                    {
                        var spec = peaksByScan[upIdx + 1];
                        if (spec.Count > 0)
                        {
                            int ci = FindNearest(spec, centroid);
                            double cmz = spec[ci].Mass;
                            if (cmz <= centroid + 3 * ftlSd && cmz >= centroid - 3 * ftlSd && !visited[upIdx + 1][ci])
                            {
                                members.Add(spec[ci]); // push_back
                                UpdateIterativeWeightedMean(cmz, (float)spec[ci].Intensity, ref centroid, ref counter, ref denom);
                                gathered.Add((upIdx + 1, ci));
                                missedUp = 0;
                            }
                            else missedUp++;
                        }
                        upIdx++; scanUp++;
                        if (missedUp > _maxOutlierScans) upActive = false;
                    }
                }

                // isTraceValid_ (MassTraceDetection.cpp:423-447)
                int totalScans = scanDown + scanUp + 1;
                double rtRange = Math.Abs(members[^1].Rt - members[0].Rt);
                if (rtRange < _minTraceLengthRt) continue;
                int adjusted = totalScans - missedDown - missedUp;
                double quality = (double)members.Count / adjusted;
                if (quality < _minSampleRate) continue;

                foreach (var (s, i) in gathered) visited[s][i] = true;
                traces.Add(members);
            }
            return traces;
        }

        // OpenMS MassTraceDetection::updateIterativeWeightedMean_ (MassTraceDetection.cpp:58-73) — the
        // exact multiplicative iterative form (NOT a naive Σwv/Σw; replicated for FP fidelity).
        private static void UpdateIterativeWeightedMean(double addedValue, double addedIntensity,
            ref double centroid, ref double prevCounter, ref double prevDenom)
        {
            double counterTmp = 1.0 + addedIntensity * addedValue / prevCounter;
            double denomTmp = 1.0 + addedIntensity / prevDenom;
            centroid *= counterTmp / denomTmp;
            prevCounter *= counterTmp;
            prevDenom *= denomTmp;
        }

        // OpenMS MSSpectrum::findNearest: index of the peak whose mass is closest to target (mass-sorted).
        private static int FindNearest(List<MassPeak> spec, double target)
        {
            int lo = 0, hi = spec.Count - 1;
            while (lo < hi)
            {
                int mid = (lo + hi) / 2;
                if (spec[mid].Mass < target) lo = mid + 1; else hi = mid;
            }
            int best = lo;
            if (lo > 0 && Math.Abs(spec[lo - 1].Mass - target) <= Math.Abs(spec[lo].Mass - target)) best = lo - 1;
            return best;
        }

        // ── Step 3: trace → MassFeature (decomposed into per-charge traces) ───
        // Returns null if the trace fails the OpenMS MassFeatureTrace feature-level isotope-cosine
        // filter (MassFeatureTrace.cpp:109-143): accumulate each member peak-group's per-isotope
        // intensities into a trace-level pattern (shifted by iso_off relative to the trace centroid),
        // then drop the trace if its averagine cosine < min_isotope_cosine (0.75). This is the dominant
        // over-production reducer that our independent tracer was missing (deferred "Phase 3b").
        // Differential-tested vs the real OpenMS static getIsotopeCosineAndDetermineIsotopeIndex
        // (massfeature_snip): offset & filter decisions identical, cosine within 2e-3 (averagine boundary).
        private MassFeature BuildFeature(List<MassPeak> tracePeaks, int id)
        {
            // ── Feature-level isotope-cosine filter (OpenMS MassFeatureTrace.cpp:84-143) ──
            // Trace centroid = intensity-weighted mean neutral mass (mt.getCentroidMZ()).
            double totalForCentroid = tracePeaks.Sum(p => p.Intensity);
            double centroid = totalForCentroid > 0
                ? tracePeaks.Sum(p => p.Mass * p.Intensity) / totalForCentroid
                : tracePeaks[0].Mass;

            // per_isotope_intensity accumulated in FLOAT (OpenMS std::vector<float>). Each member peak
            // group contributes its per-isotope vector shifted by iso_off = (int)(0.5 + dMass/isoDa)
            // (C# (int) truncates toward zero == C++ int()).
            var perIso = new float[IsotopeVectorSize];
            bool anyIsotopeData = false;
            foreach (var pk in tracePeaks)
                foreach (var env in pk.Envelopes)
                {
                    var iso = env.PerIsotopeIntensities; // double[], min-negative-isotope based (== OpenMS getIsotopeIntensities)
                    if (iso == null) continue;
                    anyIsotopeData = true;
                    int isoOff = (int)(0.5 + (env.MonoisotopicMass - centroid) / _isoDaDistance);
                    for (int i = 0; i < perIso.Length - isoOff; i++)
                    {
                        if (i + isoOff < 0 || i >= iso.Count) continue;
                        perIso[i + isoOff] += (float)iso[i];
                    }
                }

            // Real MetaFlashDecon envelopes always carry per-isotope intensities (set in Deconvolute), so
            // the filter always applies in production — matching OpenMS, which computes this cosine for
            // every trace. Only skip it when NO envelope supplies per-isotope data (synthetic envelopes or
            // a non-MetaFlashDecon source), where the cosine is unevaluable rather than genuinely zero.
            if (anyIsotopeData)
            {
            var perIsoD = new double[IsotopeVectorSize];
            for (int i = 0; i < IsotopeVectorSize; i++) perIsoD[i] = perIso[i];
            // window_width = 0 (so the trace offset is always 0) + iso_int_shift = 0 + min_iso_size = 2,
            // exactly the args OpenMS MassFeatureTrace passes (MassFeatureTrace.cpp:138).
            double isotopeCosine = MetaFlashDeconPeakGroup.GetIsotopeCosineAndDetermineIsotopeIndex(
                perIsoD, _averagine.Get(centroid), _averagine.GetApexIndex(centroid), 0, 0, 2, out _);
            if (isotopeCosine < _minIsotopeCosine) return null;
            }

            // Group the trace's constituent envelopes by charge; each charge becomes one
            // CorrectedTrace so MassFeature.Finalise() / ToMs1Feature() work unchanged.
            var byCharge = new Dictionary<int, List<CorrectedEnvelope>>();
            foreach (var pk in tracePeaks)
                foreach (var env in pk.Envelopes)
                {
                    if (!byCharge.TryGetValue(env.Charge, out var list))
                        byCharge[env.Charge] = list = new List<CorrectedEnvelope>();
                    list.Add(new CorrectedEnvelope
                    {
                        ScanIndex = pk.ScanIndex,
                        ScanNumber = pk.ScanNumber,
                        RT = pk.Rt,
                        OriginalMass = env.MonoisotopicMass,
                        CorrectedMass = env.MonoisotopicMass,
                        Charge = env.Charge,
                        Intensity = env.TotalIntensity,
                        WasCorrected = false,
                    });
                }

            var feature = new MassFeature { Id = id };
            foreach (var (charge, envs) in byCharge)
            {
                envs.Sort((a, b) => a.ScanIndex.CompareTo(b.ScanIndex));
                double intensity = envs.Sum(e => e.Intensity);
                double consensusMass = intensity > 0
                    ? envs.Sum(e => e.CorrectedMass * e.Intensity) / intensity
                    : envs.Average(e => e.CorrectedMass);
                feature.Traces.Add(new CorrectedTrace
                {
                    Id = charge,
                    Charge = charge,
                    ConsensusMass = consensusMass,
                    Envelopes = envs,
                });
            }

            feature.Finalise();
            return feature;
        }

        // Per-scan collapsed neutral-mass peak: one neutral mass with charge intensities summed.
        private sealed class MassPeak
        {
            public int ScanIndex;
            public int ScanNumber;
            public double Rt;
            public double Mass;
            public double Intensity;
            public List<IsotopicEnvelope> Envelopes;
            public bool Visited;
        }
    }
}
