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
        /// <summary>Mass tolerance (ppm) for collapsing charges per scan and for trace extension.</summary>
        public const double DefaultMassTolerancePpm = 10.0;
        /// <summary>Minimum trace RT span to keep, in RetentionTime units (minutes); FLASHDeconv 10 s = 10/60.</summary>
        public const double DefaultMinTraceLengthRt = 10.0 / 60.0;
        /// <summary>Minimum fraction of spanned scans that must carry a peak (FLASHDeconv min_sample_rate).</summary>
        public const double DefaultMinSampleRate = 0.05;
        /// <summary>Consecutive empty scans that terminate a trace direction (FLASHDeconv trace_termination_outliers).</summary>
        public const int DefaultMaxOutlierScans = 20;

        /// <summary>OpenMS MassFeatureTrace min_isotope_cosine for MS1 (MassFeatureTrace.cpp:32). Traces
        /// whose summed isotope pattern has cosine &lt; this vs the averagine are dropped.</summary>
        public const double DefaultMinIsotopeCosine = 0.75;
        // Per-isotope accumulation vector size. >= OpenMS avg.getMaxIsotopeIndex() (=299 at max_mass
        // 100k); the cosine ignores trailing zeros, so any size >= the highest occupied isotope is exact
        // (verified vs massfeature_snip: cos within 2e-3, offset/filter identical).
        private const int IsotopeVectorSize = 400;

        private readonly double _massTolerancePpm;
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

            // Step 1 — collapse each scan to neutral-mass peaks (charges summed).
            var peaksByScan = new List<MassPeak>[ms1Scans.Count];
            for (int s = 0; s < ms1Scans.Count; s++)
                peaksByScan[s] = CollapseToNeutralMassPeaks(
                    s, ms1Scans[s].OneBasedScanNumber, ms1Scans[s].RetentionTime, perScanEnvelopes[s]);

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

        // ── Step 1: per-scan neutral-mass collapse ────────────────────────────
        // Envelopes of the same proteoform at different charges share a neutral mass; group
        // them (within ppm), sum intensities, and keep the constituent per-charge envelopes.
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

        // ── Step 2: Kenar-style neutral-mass trace detection ──────────────────
        private List<List<MassPeak>> DetectTraces(List<MassPeak>[] peaksByScan)
        {
            int nScans = peaksByScan.Length;

            // Apex-sorted seeding: most intense unvisited peak first.
            var allPeaks = peaksByScan.SelectMany(p => p).OrderByDescending(p => p.Intensity).ToList();

            var traces = new List<List<MassPeak>>();
            foreach (var seed in allPeaks)
            {
                if (seed.Visited) continue;
                seed.Visited = true;

                var members = new List<MassPeak> { seed };
                double wMass = seed.Mass * seed.Intensity;
                double wInt = seed.Intensity;

                // Extend forward then backward from the seed scan, tracking the
                // intensity-weighted-mean centroid mass; terminate a direction after
                // _maxOutlierScans consecutive scans with no acceptable peak.
                Extend(peaksByScan, seed.ScanIndex, +1, nScans, members, ref wMass, ref wInt);
                Extend(peaksByScan, seed.ScanIndex, -1, nScans, members, ref wMass, ref wInt);

                members.Sort((a, b) => a.ScanIndex.CompareTo(b.ScanIndex));
                if (PassesLengthFilter(members)) traces.Add(members);
            }
            return traces;
        }

        private void Extend(List<MassPeak>[] peaksByScan, int seedScan, int step, int nScans,
            List<MassPeak> members, ref double wMass, ref double wInt)
        {
            int outliers = 0;
            for (int s = seedScan + step; s >= 0 && s < nScans; s += step)
            {
                double centroid = wInt > 0 ? wMass / wInt : members[0].Mass;
                double tolDa = centroid * _massTolerancePpm * 1e-6;

                MassPeak best = null;
                double bestDelta = double.MaxValue;
                foreach (var p in peaksByScan[s])
                {
                    if (p.Visited) continue;
                    double delta = Math.Abs(p.Mass - centroid);
                    if (delta <= tolDa && delta < bestDelta) { best = p; bestDelta = delta; }
                }

                if (best != null)
                {
                    best.Visited = true;
                    members.Add(best);
                    wMass += best.Mass * best.Intensity;
                    wInt += best.Intensity;
                    outliers = 0;
                }
                else if (++outliers > _maxOutlierScans)
                {
                    break;
                }
            }
        }

        private bool PassesLengthFilter(List<MassPeak> members)
        {
            double rtSpan = members[^1].Rt - members[0].Rt;
            if (rtSpan < _minTraceLengthRt) return false;

            int scansSpanned = members[^1].ScanIndex - members[0].ScanIndex + 1;
            double sampleRate = (double)members.Count / scansSpanned;
            return sampleRate >= _minSampleRate;
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
