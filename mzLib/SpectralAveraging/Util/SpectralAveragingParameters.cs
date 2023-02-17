using System.Text;

namespace SpectralAveraging
{
    public class SpectralAveragingParameters
    {
        #region Averaging Options
        public OutlierRejectionType OutlierRejectionType { get; set; }
        public SpectraWeightingType SpectralWeightingType { get; set; }
        public NormalizationType NormalizationType { get; set; }
        public SpectralAveragingType SpectralAveragingType { get; set; }
        public SpectraFileAveragingType SpectraFileAveragingType { get; set; }
        public OutputType OutputType { get; set; }
        public double Percentile { get; set; }
        public double MinSigmaValue { get; set; }
        public double MaxSigmaValue { get; set; }
        public double BinSize { get; set; }
        public int NumberOfScansToAverage { get; set; }
        public int ScanOverlap { get; set; }
        public int MaxThreadsToUsePerFile { get; set; } = 1;

        #endregion

        public SpectralAveragingParameters()
        {
            SetDefaultValues();
        }

        /// <summary>
        ///     Can be used to set the values of the options class in one method call
        /// </summary>
        public void SetValues(OutlierRejectionType outlierRejectionType = OutlierRejectionType.NoRejection,
            SpectraWeightingType spectraWeighingType = SpectraWeightingType.WeightEvenly,
            SpectralAveragingType spectralAveragingType = SpectralAveragingType.MzBinning,
            NormalizationType normalizationType = NormalizationType.RelativeToTics,
            SpectraFileAveragingType specAveragingType = SpectraFileAveragingType.AverageAll,
            OutputType outputType = OutputType.MzML, int numToAverage = 5, int overlap = 2,
            double percentile = 0.1, double minSigma = 1.5, double maxSigma = 1.5, double binSize = 0.01,
            int maxThreads = 1)
        {
            OutlierRejectionType = outlierRejectionType;
            SpectralWeightingType = spectraWeighingType;
            SpectralAveragingType = spectralAveragingType;
            NormalizationType = normalizationType;
            SpectraFileAveragingType = specAveragingType;
            OutputType = outputType;
            NumberOfScansToAverage = numToAverage;
            ScanOverlap = overlap;
            Percentile = percentile;
            MinSigmaValue = minSigma;
            MaxSigmaValue = maxSigma;
            BinSize = binSize;
            MaxThreadsToUsePerFile = maxThreads;
        }

        /// <summary>
        ///     Sets the values of the options to their defaults
        /// </summary>
        public void SetDefaultValues()
        {
            OutlierRejectionType = OutlierRejectionType.NoRejection;
            SpectralWeightingType = SpectraWeightingType.WeightEvenly;
            SpectralAveragingType = SpectralAveragingType.MzBinning;
            SpectraFileAveragingType = SpectraFileAveragingType.AverageAll;
            NormalizationType = NormalizationType.RelativeToTics;
            OutputType = OutputType.MzML;
            ScanOverlap = 2;
            NumberOfScansToAverage = 5;
            Percentile = 0.1;
            MinSigmaValue = 1.5;
            MaxSigmaValue = 1.5;
            BinSize = 0.01;
            MaxThreadsToUsePerFile = 1;
        }

        /// <summary>
        ///     Override for the ToString method that can be used for file output naming
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            var stringBuilder = new StringBuilder();
            stringBuilder.Append(OutlierRejectionType.ToString() + '_');
            stringBuilder.Append(SpectralWeightingType.ToString() + '_');
            stringBuilder.Append(NormalizationType + "_");

            if (SpectralAveragingType == SpectralAveragingType.MzBinning)
                stringBuilder.Append("BinSize-" + BinSize + "_");

            // rejection type specific 
            if (OutlierRejectionType == OutlierRejectionType.PercentileClipping)
                stringBuilder.Append("Percentile-" + Percentile + '_');
            if (OutlierRejectionType is OutlierRejectionType.WinsorizedSigmaClipping
                or OutlierRejectionType.AveragedSigmaClipping
                or OutlierRejectionType.SigmaClipping)
            {
                stringBuilder.Append("MinSigma-" + MinSigmaValue + '_');
                stringBuilder.Append("MaxSigma-" + MaxSigmaValue + '_');
            }

            stringBuilder.Remove(stringBuilder.ToString().LastIndexOf('_'), 1);

            return stringBuilder.ToString();
        }
    }
}