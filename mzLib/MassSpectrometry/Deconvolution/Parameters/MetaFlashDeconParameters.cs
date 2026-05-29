#nullable enable
using System;

namespace MassSpectrometry
{
    /// <summary>
    /// Parameters for the FLASHDeconv deconvolution algorithm.
    /// FLASHDeconv uses a harmonic charge-mass transformation to rapidly identify
    /// isotopic envelopes across charge states without brute-force charge enumeration.
    /// See: Kim et al. (2020) DOI: 10.1021/acs.jproteome.9b00738
    ///
    /// ⚠ The defaults below (minCosineScore 0.85, snrThreshold 0.5, tolDivFactor 2.5,
    /// overlapDedupTolFactor 1.5, charge 1–60, tol 10 ppm, mass 50–100000) are the LITERAL OpenMS
    /// FLASHDeconv values, each cited to its source line. They were tuned by differential-testing the
    /// whole pipeline against the real OpenMS library to exact per-scan agreement — they are NOT
    /// free knobs. Changing a default re-opens the fidelity gap; treat them as load-bearing constants.
    /// </summary>
    public class MetaFlashDeconParameters : DeconvolutionParameters
    {
        /// <inheritdoc />
        public override DeconvolutionType DeconvolutionType { get; protected set; }
            = DeconvolutionType.MetaFlashDecon;

        /// <summary>
        /// Number of harmonic charge relationships considered when scoring an isotopic envelope.
        /// Higher values increase sensitivity at the cost of speed.
        /// </summary>
        public int PrecursorHarmonicCount { get; set; }

        /// <summary>
        /// Maximum number of isotopic peaks to consider per candidate envelope.
        /// </summary>
        public int MaxIsotopicPeakCount { get; set; }

        /// <summary>
        /// Minimum number of isotopic peaks required to report an envelope.
        /// Envelopes with fewer observed peaks than this threshold are discarded.
        /// </summary>
        public int MinIsotopicPeakCount { get; set; }

        /// <summary>
        /// Peak matching tolerance in parts per million.
        /// Used when assigning observed peaks to theoretical isotope positions.
        /// </summary>
        public double DeconvolutionTolerancePpm { get; set; }

        /// <summary>
        /// Minimum cosine similarity score between the observed and Averagine-predicted
        /// isotope distributions required to report an envelope.
        /// Range: 0.0 (no similarity) to 1.0 (perfect match).
        /// </summary>
        public double MinCosineScore { get; set; }

        /// <summary>
        /// Minimum fraction of total theoretical isotope intensity that must be
        /// accounted for by observed peaks. Envelopes below this threshold are discarded.
        /// </summary>
        public double IsotopeIntensityRatioThreshold { get; set; }

        /// <summary>
        /// Maximum neutral monoisotopic mass (Da) to report.
        /// Candidate envelopes with inferred masses above this value are discarded.
        /// </summary>
        public double MaxMassRange { get; set; }

        /// <summary>
        /// Minimum neutral monoisotopic mass (Da) to report.
        /// Candidate envelopes with inferred masses below this value are discarded.
        /// </summary>
        public double MinMassRange { get; set; }

        // ── Faithful OpenMS FLASHDeconv tolerance mechanism ───────────────────
        // OpenMS does NOT bin the log-m/z axis at the requested ppm tolerance. It bins
        // FINER — at ppm / TolDivFactor — during candidate-mass voting (for sensitivity),
        // then re-applies the (wider) input tolerance only at the final overlap dedup.
        // This narrow-then-wide trick is load-bearing for reproducing FLASHDeconv's
        // candidate set. These properties are consumed ONLY by MetaFlashDeconAlgorithm;
        // they do not affect any other deconvolution algorithm. Defaults are the literal
        // OpenMS values. See FLASHDeconvAlgorithm.cpp:22 and :1223, and the citations in
        // MetaFlashDeconAlgorithm.cs.

        /// <summary>
        /// Bin-width divisor used during candidate-mass voting: the log-m/z axis is binned
        /// at <c>DeconvolutionTolerancePpm / TolDivFactor</c> (finer than the requested
        /// tolerance). OpenMS: <c>tol_div_factor = 2.5</c> (FLASHDeconvAlgorithm.cpp:22,
        /// applied in updateMembers_:171).
        /// </summary>
        public double TolDivFactor { get; set; }

        /// <summary>
        /// Tolerance multiplier for the final overlap-dedup window. OpenMS calls
        /// <c>removeOverlappingPeakGroups_</c> with a window of
        /// <c>(ppm / TolDivFactor) × TolDivFactor × 1.5</c> = 1.5 × the input ppm
        /// (FLASHDeconvAlgorithm.cpp:1223). Setting <see cref="OverlapDedupTolFactor"/>
        /// to 1.5 reproduces that.
        /// </summary>
        public double OverlapDedupTolFactor { get; set; }

        /// <summary>
        /// Minimum overall signal-to-noise ratio required to report an envelope. The SNR is the
        /// faithful OpenMS all-charge value (<c>getSNR()</c>, aggregated by
        /// <c>updateSNR_</c>); envelopes below this are dropped as harmonics / noise. OpenMS:
        /// <c>snr_threshold = 0.5</c> (FLASHDeconvAlgorithm.cpp:1163). Consumed only by
        /// <see cref="MetaFlashDeconAlgorithm"/>.
        /// </summary>
        public double SnrThreshold { get; set; }

        /// <summary>
        /// Constructs FLASHDeconv parameters with sensible defaults suitable for
        /// top-down proteomics. All parameters have defaults so that
        /// <c>new MetaFlashDeconParameters()</c> is a valid zero-argument call.
        /// </summary>
        /// <param name="minCharge">Minimum assumed charge state (inclusive).</param>
        /// <param name="maxCharge">Maximum assumed charge state (inclusive).</param>
        /// <param name="deconvolutionTolerancePpm">Peak matching tolerance in ppm.</param>
        /// <param name="minIsotopicPeakCount">Minimum isotopic peaks required per envelope.</param>
        /// <param name="maxIsotopicPeakCount">Maximum isotopic peaks considered per envelope.</param>
        /// <param name="precursorHarmonicCount">Number of harmonic charge relationships to consider.</param>
        /// <param name="minCosineScore">Minimum cosine similarity score (0–1).</param>
        /// <param name="isotopeIntensityRatioThreshold">Minimum fraction of theoretical intensity that must be observed.</param>
        /// <param name="minMassRange">Minimum neutral mass (Da) to report.</param>
        /// <param name="maxMassRange">Maximum neutral mass (Da) to report.</param>
        /// <param name="polarity">Ionisation polarity.</param>
        /// <param name="averageResidueModel">Averagine model used for theoretical isotope distributions. Defaults to <see cref="Averagine"/>.</param>
        /// <param name="tolDivFactor">Bin-width divisor for candidate voting (OpenMS <c>tol_div_factor</c>, default 2.5).</param>
        /// <param name="overlapDedupTolFactor">Overlap-dedup window multiplier (OpenMS final dedup = 1.5 × input ppm, default 1.5).</param>
        /// <param name="snrThreshold">Minimum all-charge SNR to report an envelope (OpenMS <c>snr_threshold</c>, default 0.5).</param>
        public MetaFlashDeconParameters(
            int minCharge = 1,
            int maxCharge = 60,
            double deconvolutionTolerancePpm = 10.0,
            int minIsotopicPeakCount = 3,
            int maxIsotopicPeakCount = 50,
            int precursorHarmonicCount = 2,
            double minCosineScore = 0.85, // FLASHDeconv min_isotope_cosine default (MS1); FLASHDeconvAlgorithm.cpp / FLASHDeconv.cpp:193
            double isotopeIntensityRatioThreshold = 0.2,
            double minMassRange = 50.0,
            double maxMassRange = 100000.0,
            Polarity polarity = Polarity.Positive,
            AverageResidue? averageResidueModel = null,
            double tolDivFactor = 2.5,          // OpenMS tol_div_factor, FLASHDeconvAlgorithm.cpp:22
            double overlapDedupTolFactor = 1.5, // OpenMS final dedup = 1.5 × input ppm, FLASHDeconvAlgorithm.cpp:1223
            double snrThreshold = 0.5)          // OpenMS snr_threshold, FLASHDeconvAlgorithm.cpp:1163
            : base(minCharge, maxCharge, polarity, averageResidueModel)
        {
            DeconvolutionTolerancePpm = deconvolutionTolerancePpm;
            MinIsotopicPeakCount = minIsotopicPeakCount;
            MaxIsotopicPeakCount = maxIsotopicPeakCount;
            PrecursorHarmonicCount = precursorHarmonicCount;
            MinCosineScore = minCosineScore;
            IsotopeIntensityRatioThreshold = isotopeIntensityRatioThreshold;
            MinMassRange = minMassRange;
            MaxMassRange = maxMassRange;
            TolDivFactor = tolDivFactor;
            OverlapDedupTolFactor = overlapDedupTolFactor;
            SnrThreshold = snrThreshold;
        }

        #region IEquatable<DeconvolutionParameters>

        protected override bool EqualProperties(DeconvolutionParameters other)
        {
            var o = (MetaFlashDeconParameters)other;
            return PrecursorHarmonicCount == o.PrecursorHarmonicCount
                && MaxIsotopicPeakCount == o.MaxIsotopicPeakCount
                && MinIsotopicPeakCount == o.MinIsotopicPeakCount
                && DeconvolutionTolerancePpm.Equals(o.DeconvolutionTolerancePpm)
                && MinCosineScore.Equals(o.MinCosineScore)
                && IsotopeIntensityRatioThreshold.Equals(o.IsotopeIntensityRatioThreshold)
                && MinMassRange.Equals(o.MinMassRange)
                && MaxMassRange.Equals(o.MaxMassRange)
                && TolDivFactor.Equals(o.TolDivFactor)
                && OverlapDedupTolFactor.Equals(o.OverlapDedupTolFactor)
                && SnrThreshold.Equals(o.SnrThreshold);
        }

        protected override void AddHashCodes(HashCode hash)
        {
            hash.Add(PrecursorHarmonicCount);
            hash.Add(MaxIsotopicPeakCount);
            hash.Add(MinIsotopicPeakCount);
            hash.Add(DeconvolutionTolerancePpm);
            hash.Add(MinCosineScore);
            hash.Add(IsotopeIntensityRatioThreshold);
            hash.Add(MinMassRange);
            hash.Add(MaxMassRange);
            hash.Add(TolDivFactor);
            hash.Add(OverlapDedupTolFactor);
            hash.Add(SnrThreshold);
        }

        public override MetaFlashDeconParameters Clone()
        {
            return new MetaFlashDeconParameters(
                MinAssumedChargeState, MaxAssumedChargeState,
                DeconvolutionTolerancePpm, MinIsotopicPeakCount, MaxIsotopicPeakCount,
                PrecursorHarmonicCount, MinCosineScore, IsotopeIntensityRatioThreshold,
                MinMassRange, MaxMassRange, Polarity, AverageResidueModel,
                TolDivFactor, OverlapDedupTolFactor, SnrThreshold)
            {
                ExpectedIsotopeSpacing = ExpectedIsotopeSpacing,
                UseGenericScore = UseGenericScore
            };
        }

        // Thread-safe lazy caching of decoy parameters (mirrors ClassicDeconvolutionParameters).
        private volatile DeconvolutionParameters? _decoyParams = null;
        private readonly object _decoyParamsLock = new();
        public override DeconvolutionParameters ToDecoyParameters()
        {
            if (_decoyParams is not null)
                return _decoyParams;

            lock (_decoyParamsLock)
            {
                return _decoyParams ??= new MetaFlashDeconParameters(
                    MinAssumedChargeState, MaxAssumedChargeState,
                    DeconvolutionTolerancePpm, MinIsotopicPeakCount, MaxIsotopicPeakCount,
                    PrecursorHarmonicCount, MinCosineScore, IsotopeIntensityRatioThreshold,
                    MinMassRange, MaxMassRange, Polarity,
                    averageResidueModel: new DecoyAveragine(AverageResidueModel,
                        DecoyAveragine.DefaultDecoyIsotopeSpacing, ExpectedIsotopeSpacing),
                    tolDivFactor: TolDivFactor,
                    overlapDedupTolFactor: OverlapDedupTolFactor,
                    snrThreshold: SnrThreshold)
                {
                    ExpectedIsotopeSpacing = DecoyAveragine.DefaultDecoyIsotopeSpacing,
                    UseGenericScore = UseGenericScore
                };
            }
        }

        #endregion
    }
}