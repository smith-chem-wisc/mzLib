namespace MzLibSpectralAveraging
{
    public class MzLibSpectralAveragingOptions
    {
        #region Public Properties
        public SpectralAveragingOptions SpectralAveragingOptions { get; set; }
        public SpectraFileProcessingType SpectraFileProcessingType { get; set; }
        public int NumberOfScansToAverage { get; set; }
        public int ScanOverlap { get; set; }
        public OutputType OutputType { get; set; }

        #region SpectralAveragingOptions Wrapped Properties

        public RejectionType RejectionType
        {
            get => SpectralAveragingOptions.RejectionType;
            set => SpectralAveragingOptions.RejectionType = value;
        }

        public WeightingType WeightingType
        {
            get => SpectralAveragingOptions.WeightingType;
            set => SpectralAveragingOptions.WeightingType = value;
        }

        public SpectrumMergingType SpectrumMergingType
        {
            get => SpectralAveragingOptions.SpectrumMergingType;
            set => SpectralAveragingOptions.SpectrumMergingType = value;
        }

        public bool PerformNormalization
        {
            get => SpectralAveragingOptions.PerformNormalization;
            set => SpectralAveragingOptions.PerformNormalization = value;
        }

        public double Percentile
        {
            get => SpectralAveragingOptions.Percentile;
            set => SpectralAveragingOptions.Percentile = value;
        }

        public double MinSigmaValue
        {
            get => SpectralAveragingOptions.MinSigmaValue;
            set => SpectralAveragingOptions.MinSigmaValue = value;
        }

        public double MaxSigmaValue
        {
            get => SpectralAveragingOptions.MaxSigmaValue;
            set => SpectralAveragingOptions.MaxSigmaValue = value;  
        }

        public double BinSize
        {
            get => SpectralAveragingOptions.BinSize;
            set => SpectralAveragingOptions.BinSize = value;
        }

        #endregion

        #endregion

        public MzLibSpectralAveragingOptions(SpectralAveragingOptions spectralAveragingOptions)
        {
            SpectralAveragingOptions = spectralAveragingOptions;
            SetDefaultValues(true);
        }

        public MzLibSpectralAveragingOptions() => this.SetDefaultValues();

        /// <summary>
        /// Can be used to set the values of the options in a single call
        /// </summary>
        /// <param name="spectraFileProcessingType"></param>
        /// <param name="numberOfScansToAverage"></param>
        /// <param name="scanOverlap"></param>
        /// <param name="outputType"></param>
        public void SetValues(SpectraFileProcessingType spectraFileProcessingType = SpectraFileProcessingType.AverageAll, 
            int numberOfScansToAverage = 5, int scanOverlap = 2, OutputType outputType = OutputType.mzML)
        {
            SpectraFileProcessingType = spectraFileProcessingType;
            NumberOfScansToAverage = numberOfScansToAverage;
            ScanOverlap = scanOverlap;
            OutputType = outputType;
        }

        /// <summary>
        /// Sets the values of the options to their defaults
        /// </summary>
        public void SetDefaultValues(bool resetOnlyMzLibValues = false)
        {
            SpectraFileProcessingType = SpectraFileProcessingType.AverageAll;
            NumberOfScansToAverage = 5;
            ScanOverlap = 2;
            OutputType = OutputType.mzML;
            if (SpectralAveragingOptions == null)
                SpectralAveragingOptions = new();
            else if (!resetOnlyMzLibValues)
                SpectralAveragingOptions.SetDefaultValues();
        }
    }
}
