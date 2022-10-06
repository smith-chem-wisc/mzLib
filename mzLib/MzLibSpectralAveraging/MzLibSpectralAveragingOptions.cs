using SpectralAveraging;

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

        #endregion

        public MzLibSpectralAveragingOptions(SpectralAveragingOptions? spectralAveragingOptions = null)
        {
            SpectralAveragingOptions = spectralAveragingOptions ?? new SpectralAveragingOptions();
            SetDefaultValues();
        }

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
        public void SetDefaultValues()
        {
            SpectraFileProcessingType = SpectraFileProcessingType.AverageAll;
            NumberOfScansToAverage = 5;
            ScanOverlap = 2;
            OutputType = OutputType.mzML;
        }
    }
}
