#nullable enable
using System.Runtime.InteropServices;
using Chemistry;

namespace MassSpectrometry;

public class IsoDecDeconvolutionParameters : DeconvolutionParameters
{
    public override DeconvolutionType DeconvolutionType { get; protected set; } = DeconvolutionType.IsoDecDeconvolution;

    #region Interop Parameters 

    /// <summary>
    /// The struct that is passed into the isodec.dll
    /// </summary>
    public struct IsoSettings
    {
        public int phaseres; // Precision of encoding matrix
        public int verbose; // Verbose output
        public int peakwindow; // Peak Detection Window
        public float peakthresh; // Peak Detection Threshold
        public int minpeaks; // Minimum Peaks for an allowed peak
        public float css_thresh; // Minimum cosine similarity score for isotope distribution
        public float matchtol; // Match Tolerance for peak detection in ppm
        public int maxshift; // Maximum shift allowed for isotope distribution
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 2)]
        public float[] mzwindow; // MZ Window for isotope distribution
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 2)]
        public float[] plusoneintwindow; // Plus One Intensity range. Will be used for charge state 1
        public int knockdown_rounds; // Number of knockdown rounds
        public float min_score_diff; // Minimum score difference for isotope distribution to allow missed monoisotopic peaks
        public float minareacovered; // Minimum area covered by isotope distribution. Use in or with css_thresh
        public int isolength; // Isotope Distribution Length
        public double mass_diff_c; // Mass difference between isotopes
        public float adductmass; // Adduct Mass
        public int minusoneaszero; // Use set the -1 isotope as 0 to help force better alignments
        public float isotopethreshold; // Threshold for isotope distribution. Will remove relative intensities below this.
        public float datathreshold; // Threshold for data. Will remove relative intensities below this relative to max intensity in each cluster
        public float zscore_threshold; //Ratio above which a secondary charge state prediction will be returned.
    }

    private IsoSettings? _isoSettings;

    internal IsoSettings ToIsoSettings()
    {
        if (_isoSettings != null)
            return _isoSettings.Value;

        IsoSettings result = new IsoSettings
        {
            phaseres = PhaseRes,
            verbose = Verbose,
            peakwindow = PeakWindow,
            peakthresh = PeakThreshold,
            minpeaks = MinPeaks,
            css_thresh = CssThreshold,
            matchtol = MatchTolerance,
            maxshift = MaxShift,
            mzwindow = MzWindow,
            plusoneintwindow = PlusOneIntWindow,
            knockdown_rounds = KnockdownRounds,
            min_score_diff = MinScoreDiff,
            minareacovered = MinAreaCovered,
            isolength = IsoLength,
            // NOTE: We use ExpectedIsotopeSpacing (from the base class) rather than the
            // legacy MassDiffC field so that decoy runs (ExpectedIsotopeSpacing = 0.9444 Da,
            // set via ToDecoyParameters()) are correctly forwarded to the DLL.
            // For normal runs ExpectedIsotopeSpacing == Constants.C13MinusC12 == MassDiffC.
            mass_diff_c = ExpectedIsotopeSpacing,
            adductmass = AdductMass,
            minusoneaszero = MinusOneAreasZero,
            isotopethreshold = IsotopeThreshold,
            datathreshold = DataThreshold,
            zscore_threshold = ZScoreThreshold
        };

        _isoSettings = result;
        return result;
    }

    #endregion

    #region User-Accessible Parameters

    /// <summary>
    /// Precision of encoding matrix
    /// </summary>
    public int PhaseRes { get; set; }

    /// <summary>
    /// Minimum cosine similarity score for isotope distribution
    /// </summary>
    public float CssThreshold { get; set; }

    /// <summary>
    /// Match Tolerance for peak detection in ppm
    /// </summary>
    public float MatchTolerance { get; set; }

    /// <summary>
    /// Maximum shift allowed for isotope distribution
    /// </summary>
    public int MaxShift { get; set; }

    /// <summary>
    /// MZ Window for isotope distribution
    /// </summary>
    public float[] MzWindow { get; set; }

    /// <summary>
    /// Number of knockdown rounds
    /// </summary>
    public int KnockdownRounds { get; set; }

    /// <summary>
    /// Minimum area covered by isotope distribution. Use in or with css_thresh
    /// </summary>
    public float MinAreaCovered { get; set; }

    /// <summary>
    /// Threshold for data. Will remove relative intensities below this relative to max intensity in each cluster
    /// </summary>
    public float DataThreshold { get; set; }

    /// <summary>
    /// Report multiple monoisotopic peaks
    /// </summary>
    public bool ReportMulitpleMonoisos { get; set; }

    #endregion User-Accessible Parameters

    #region Hard-Coded Parameters

    /// <summary>
    /// Verbose output
    /// </summary>
    public int Verbose { get; protected set; } = 0;

    /// <summary>
    /// Peak Detection Window
    /// </summary>
    public int PeakWindow { get; protected set; } = 80;

    /// <summary>
    /// Peak Detection Threshold
    /// </summary>
    public float PeakThreshold { get; protected set; } = (float)0.0001;

    /// <summary>
    /// Minimum Peaks for an allowed peak
    /// </summary>
    public int MinPeaks { get; protected set; } = 3;

    /// <summary>
    /// Plus One Intensity range. Will be used for charge state 1
    /// </summary>
    public float[] PlusOneIntWindow { get; protected set; } = new float[] { (float)0.1, (float)0.6 };

    /// <summary>
    /// Minimum score difference for isotope distribution to allow missed monoisotopic peaks
    /// </summary>
    public float MinScoreDiff { get; protected set; } = (float)0.1;

    /// <summary>
    /// Isotope Distribution Length
    /// </summary>
    public int IsoLength { get; protected set; } = 64;

    /// <summary>
    /// Legacy field kept for reference. The value actually passed to the DLL is
    /// <see cref="DeconvolutionParameters.ExpectedIsotopeSpacing"/>, which defaults
    /// to <see cref="Constants.C13MinusC12"/> and is set to
    /// <see cref="DecoyAveragine.DefaultDecoyIsotopeSpacing"/> by
    /// <see cref="ToDecoyParameters"/> for decoy runs.
    /// </summary>
    public double MassDiffC { get; protected set; } = 1.0033;

    /// <summary>
    /// Adduct Mass
    /// </summary>
    public float AdductMass { get; protected set; } = (float)1.00727276467;

    /// <summary>
    /// Use set the -1 isotope as 0 to help force better alignments
    /// </summary>
    public int MinusOneAreasZero { get; protected set; } = 1;

    /// <summary>
    /// Threshold for isotope distribution. Will remove relative intensities below this.
    /// </summary>
    public float IsotopeThreshold { get; protected set; } = (float)0.01;

    /// <summary>
    /// Ratio above which a secondary charge state prediction will be returned.
    /// </summary>
    public float ZScoreThreshold { get; protected set; } = (float)0.95;

    #endregion Hard-Coded Parameters

    public IsoDecDeconvolutionParameters(
        Polarity polarity = Polarity.Positive,
        int phaseRes = 8,
        bool reportMultipleMonoisos = true,
        float cssThreshold = (float)0.7,
        float matchTolerance = (float)5,
        int maxShift = 3,
        float[] mzWindow = null,
        int knockdownRounds = 5,
        float minAreaCovered = (float)0.20,
        float relativeDataThreshold = (float)0.05,
        AverageResidue? averageResidueModel = null,
        double expectedIsotopeSpacing = Constants.C13MinusC12)
        : base(1, 50, polarity, averageResidueModel, expectedIsotopeSpacing)
    {
        PhaseRes = phaseRes;
        ReportMulitpleMonoisos = reportMultipleMonoisos;
        CssThreshold = cssThreshold;
        MatchTolerance = matchTolerance;
        MaxShift = maxShift;
        MzWindow = mzWindow ?? new float[] { (float)-1.05, (float)2.05 };
        KnockdownRounds = knockdownRounds;
        MinAreaCovered = minAreaCovered;
        DataThreshold = relativeDataThreshold;
        if (Polarity == Polarity.Negative)
            AdductMass = (float)-1.00727276467;
    }

    private DeconvolutionParameters? _decoyParams = null;

    /// <summary>
    /// Returns a version of these parameters configured for decoy deconvolution,
    /// with <see cref="DeconvolutionParameters.ExpectedIsotopeSpacing"/> set to
    /// <see cref="DecoyAveragine.DefaultDecoyIsotopeSpacing"/> (0.9444 Da).
    /// The cached IsoSettings struct is invalidated on the clone so that the next
    /// call to <see cref="ToIsoSettings"/> rebuilds it with the decoy spacing.
    /// </summary>
    public override DeconvolutionParameters ToDecoyParameters()
    {
        if (_decoyParams != null) return _decoyParams;

        var decoy = new IsoDecDeconvolutionParameters(
            polarity: Polarity,
            phaseRes: PhaseRes,
            reportMultipleMonoisos: ReportMulitpleMonoisos,
            cssThreshold: CssThreshold,
            matchTolerance: MatchTolerance,
            maxShift: MaxShift,
            mzWindow: MzWindow,
            knockdownRounds: KnockdownRounds,
            minAreaCovered: MinAreaCovered,
            relativeDataThreshold: DataThreshold,
            averageResidueModel: new DecoyAveragine(AverageResidueModel),
            expectedIsotopeSpacing: DecoyAveragine.DefaultDecoyIsotopeSpacing);

        return _decoyParams = decoy;
    }
}