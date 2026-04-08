#nullable enable
namespace MassSpectrometry
{
    /// <summary>
    /// Parameters for the FLASHDeconv deconvolution algorithm.
    /// FLASHDeconv uses a harmonic charge-mass transformation to rapidly identify
    /// isotopic envelopes across charge states without brute-force charge enumeration.
    /// See: Kim et al. (2020) DOI: 10.1021/acs.jproteome.9b00738
    /// </summary>
    public class FLASHDeconvolutionParameters : DeconvolutionParameters
    {
        /// <inheritdoc />
        public override DeconvolutionType DeconvolutionType { get; protected set; }
            = DeconvolutionType.FLASHDeconvolution;

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

        /// <summary>
        /// Constructs FLASHDeconv parameters with sensible defaults suitable for
        /// top-down proteomics. All parameters have defaults so that
        /// <c>new FLASHDeconvolutionParameters()</c> is a valid zero-argument call.
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
        public FLASHDeconvolutionParameters(
            int minCharge = 1,
            int maxCharge = 60,
            double deconvolutionTolerancePpm = 10.0,
            int minIsotopicPeakCount = 3,
            int maxIsotopicPeakCount = 50,
            int precursorHarmonicCount = 2,
            double minCosineScore = 0.6,
            double isotopeIntensityRatioThreshold = 0.2,
            double minMassRange = 50.0,
            double maxMassRange = 100000.0,
            Polarity polarity = Polarity.Positive,
            AverageResidue? averageResidueModel = null)
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
        }
    }
}