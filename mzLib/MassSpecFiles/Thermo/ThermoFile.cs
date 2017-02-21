using MassSpectrometry;
using MSFileReaderLib;
using MzLibUtil;
using System;
using System.Globalization;
using System.Text.RegularExpressions;

namespace IO.Thermo
{
    public abstract class ThermoFile : MsDataFile<IThermoScan>
    {

        #region Private Fields

        private static readonly Regex PolarityRegex = new Regex(@"\+ ", RegexOptions.Compiled);

        private static bool?[] couldBePrecursor;

        #endregion Private Fields

        #region Public Constructors

        public ThermoFile(IThermoScan[] scans, ThermoGlobalParams thermoGlobalParams) : base(scans)
        {
            this.thermoGlobalParams = thermoGlobalParams;
        }

        public ThermoFile(IXRawfile5 _rawConnection, int numSpectra) : base(numSpectra)
        {
            this.thermoGlobalParams = GetAllGlobalStuff(_rawConnection);
        }

        #endregion Public Constructors

        #region Private Enums

        private enum ThermoMzAnalyzer
        {
            None = -1,
            ITMS = 0,
            TQMS = 1,
            SQMS = 2,
            TOFMS = 3,
            FTMS = 4,
            Sector = 5
        }

        #endregion Private Enums

        #region Public Properties

        public ThermoGlobalParams thermoGlobalParams { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public static ThermoGlobalParams GetAllGlobalStuff(IXRawfile5 _rawConnection)
        {
            int pnNumInstMethods = 0;
            _rawConnection.GetNumInstMethods(ref pnNumInstMethods);

            string[] instrumentMethods = new string[pnNumInstMethods];
            for (int nInstMethodItem = 0; nInstMethodItem < pnNumInstMethods; nInstMethodItem++)
            {
                string pbstrInstMethod = null;
                _rawConnection.GetInstMethod(nInstMethodItem, ref pbstrInstMethod);
                instrumentMethods[nInstMethodItem] = pbstrInstMethod;
            }

            string pbstrInstSoftwareVersion = null;
            _rawConnection.GetInstSoftwareVersion(ref pbstrInstSoftwareVersion);

            string pbstrInstName = null;
            _rawConnection.GetInstName(ref pbstrInstName);

            string pbstrInstModel = null;
            _rawConnection.GetInstModel(ref pbstrInstModel);

            int pnControllerNumber = 0;
            int pnControllerType = 0;
            _rawConnection.GetCurrentController(ref pnControllerType, ref pnControllerNumber);

            return new ThermoGlobalParams(pnNumInstMethods, instrumentMethods, pbstrInstSoftwareVersion, pbstrInstName, pbstrInstModel, pnControllerType, pnControllerNumber);
        }

        #endregion Public Methods

        #region Protected Methods

        protected static IThermoScan GetMsDataOneBasedScanFromThermoFile(int nScanNumber, IXRawfile5 theConnection)
        {
            if (couldBePrecursor == null)
            {
                int pnFirstSpectrum = 0;
                theConnection.GetFirstSpectrumNumber(ref pnFirstSpectrum);
                int pnLastSpectrum = 0;
                theConnection.GetLastSpectrumNumber(ref pnLastSpectrum);
                couldBePrecursor = new bool?[pnLastSpectrum - pnFirstSpectrum + 1];
            }

            int pnNumPackets = 0;
            double pdLowMass = 0;
            double pdHighMass = 0;
            double pdTIC = 0;
            double pdBasePeakMass = 0;
            double pdBasePeakIntensity = 0;
            int pnNumChannels = 0;
            int pbUniformTime = 0;
            double pdFrequency = 0;
            double pdStartTime = 0;
            theConnection.GetScanHeaderInfoForScanNum(nScanNumber, ref pnNumPackets, ref pdStartTime, ref pdLowMass, ref pdHighMass, ref pdTIC, ref pdBasePeakMass, ref pdBasePeakIntensity, ref pnNumChannels, ref pbUniformTime, ref pdFrequency);

            double? ms2isolationWidth = null;
            double? precursorMonoisotopicMZfromTrailierExtra = null;
            int? chargeState = null;
            double? injectionTime = null;

            object pvarValues = null;
            object pvarLables = null;
            int pnArraySize = 0;
            theConnection.GetTrailerExtraForScanNum(nScanNumber, ref pvarLables, ref pvarValues, ref pnArraySize);

            string[] labels = (string[])pvarLables;
            string[] values = (string[])pvarValues;
            for (int i = labels.GetLowerBound(0); i <= labels.GetUpperBound(0); i++)
            {
                if (labels[i].StartsWith("MS2 Isolation Width", StringComparison.Ordinal))
                {
                    ms2isolationWidth = double.Parse(values[i], CultureInfo.InvariantCulture) == 0 ?
                        (double?)null :
                        double.Parse(values[i], CultureInfo.InvariantCulture);
                }
                if (labels[i].StartsWith("Monoisotopic M/Z", StringComparison.Ordinal))
                {
                    precursorMonoisotopicMZfromTrailierExtra = double.Parse(values[i], CultureInfo.InvariantCulture) == 0 ?
                        (double?)null :
                        double.Parse(values[i], CultureInfo.InvariantCulture);
                }
                if (labels[i].StartsWith("Charge State", StringComparison.Ordinal))
                {
                    chargeState = int.Parse(values[i], CultureInfo.InvariantCulture) == 0 ?
                        (int?)null :
                        int.Parse(values[i], CultureInfo.InvariantCulture);
                }
                if (labels[i].StartsWith("Ion Injection Time (ms)", StringComparison.Ordinal))
                {
                    injectionTime = double.Parse(values[i], CultureInfo.InvariantCulture) == 0 ?
                        (double?)null :
                        double.Parse(values[i], CultureInfo.InvariantCulture);
                }
            }

            string pbstrFilter = null;
            theConnection.GetFilterForScanNum(nScanNumber, ref pbstrFilter);

            int pbIsCentroidScan = 0;
            theConnection.IsCentroidScanForScanNum(nScanNumber, ref pbIsCentroidScan);

            int pnMSOrder = 0;
            theConnection.GetMSOrderForScanNum(nScanNumber, ref pnMSOrder);

            int pnMassAnalyzerType = 0;
            theConnection.GetMassAnalyzerTypeForScanNum(nScanNumber, ref pnMassAnalyzerType);

            double[,] data;
            try
            {
                object pvarFlags = null;
                object pvarLabels = null;
                theConnection.GetLabelData(ref pvarLabels, ref pvarFlags, ref nScanNumber);
                data = pvarLabels as double[,];
                if (data == null || data.Length == 0)
                    throw new ArgumentException("For spectrum number " + nScanNumber + " the data is null!");
            }
            catch (ArgumentException)
            {
                string bstrFilter = null;
                int nIntensityCutoffType = 0;
                int nIntensityCutoffValue = 0;
                int nMaxNumberOfPeaks = 0;
                int bCentroidResult = 1;
                double pdCentroidPeakWidth = 0;
                object pvarnMassList = null;
                object pvarPeakFlags = null;
                theConnection.GetMassListFromScanNum(ref nScanNumber,
                                bstrFilter,
                                nIntensityCutoffType,
                                nIntensityCutoffValue,
                                nMaxNumberOfPeaks,
                                bCentroidResult,
                                ref pdCentroidPeakWidth,
                                ref pvarnMassList,
                                ref pvarPeakFlags,
                                ref pnArraySize);
                data = (double[,])pvarnMassList;
            }

            MZAnalyzerType mzAnalyzerType;
            switch ((ThermoMzAnalyzer)pnMassAnalyzerType)
            {
                case ThermoMzAnalyzer.FTMS:
                    mzAnalyzerType = MZAnalyzerType.Orbitrap;
                    break;

                default:
                    mzAnalyzerType = MZAnalyzerType.Unknown;
                    break;
            }

            if (pnMSOrder > 1)
            {
                int pnActivationType = 0;
                theConnection.GetActivationTypeForScanNum(nScanNumber, pnMSOrder, ref pnActivationType);

                double pdMass = 0;
                theConnection.GetPrecursorMassForScanNum(nScanNumber, pnMSOrder, ref pdMass);

                int pnChargeState = 0;
                double pdHeaderMass = 0;
                double pdFoundMass = 0;
                int pnMasterScan = 0;
                theConnection.FindPrecursorMassInFullScan(nScanNumber, ref pnMasterScan, ref pdFoundMass, ref pdHeaderMass, ref pnChargeState);

                if (pnMasterScan == 0)
                    pnMasterScan = GetLastScanEventThatIs1(theConnection, nScanNumber);

                return new ThermoScanWithPrecursor(
                    nScanNumber,
                    new ThermoSpectrum(data),
                    pnMSOrder,
                    PolarityRegex.IsMatch(pbstrFilter) ? Polarity.Positive : Polarity.Negative,
                    pdStartTime,
                    new MzRange(pdLowMass, pdHighMass),
                    pbstrFilter,
                    mzAnalyzerType,
                    pdTIC,
                    pdFoundMass == 0 ? (double?)null : pdFoundMass,
                    pnChargeState == 0 ? chargeState : pnChargeState,
                    pdMass,
                    ms2isolationWidth,
                    (DissociationType)pnActivationType,
                    pnMasterScan,
                    precursorMonoisotopicMZfromTrailierExtra,
                    injectionTime);
            }
            else
            {
                return new ThermoScan(nScanNumber,
                    new ThermoSpectrum(data),
                    1,
                    PolarityRegex.IsMatch(pbstrFilter) ? Polarity.Positive : Polarity.Negative,
                    pdStartTime,
                    new MzRange(pdLowMass, pdHighMass),
                    pbstrFilter,
                    mzAnalyzerType,
                    pdTIC,
                    injectionTime);
            }
        }

        #endregion Protected Methods

        #region Private Methods

        private static int GetLastScanEventThatIs1(IXRawfile5 _rawConnection, int scanNumber)
        {
            int oneBasedPrecursorNumber = scanNumber;
            while (oneBasedPrecursorNumber > 0)
            {
                if (!couldBePrecursor[oneBasedPrecursorNumber - 1].HasValue)
                {
                    object pvarValue = null;
                    _rawConnection.GetTrailerExtraValueForScanNum(oneBasedPrecursorNumber, "Scan Event:", ref pvarValue);
                    couldBePrecursor[oneBasedPrecursorNumber - 1] = pvarValue.ToString().Equals("1");
                }
                if (couldBePrecursor[oneBasedPrecursorNumber - 1] == true)
                    return oneBasedPrecursorNumber;
                oneBasedPrecursorNumber--;
            }
            throw new ArgumentException("Could not find precursor for scan number " + scanNumber);
        }

        #endregion Private Methods

        //protected static ThermoSpectrum GetSpectrumFromRawFile(int spectrumNumber)
        //{
        //    double[,] data;
        //    try
        //    {
        //        data = GetLabeledData(spectrumNumber);
        //    }
        //    catch (ArgumentException)
        //    {
        //        data = GetUnlabeledData(spectrumNumber, true);
        //    }

        //    //int arrayLength = data.GetLength(1);

        //    //if (arrayLength > maxPeaksPerScan)
        //    //{
        //    //    double[] intensityArray = new double[arrayLength];
        //    //    for (int i = 0; i < arrayLength; i++)
        //    //        intensityArray[i] = data[1, i];
        //    //    var cutoffIntensity = intensityArray.Quantile(1.0 - (double)maxPeaksPerScan / arrayLength);

        //    //    int thiscOUNT = 0;
        //    //    for (int j = 0; j < arrayLength; j++)
        //    //        if (data[1, j] >= cutoffIntensity)
        //    //            thiscOUNT++;

        //    //    double[,] newData = new double[data.GetLength(0), thiscOUNT];
        //    //    int okIndex = 0;
        //    //    for (int j = 0; j < arrayLength; j++)
        //    //        if (data[1, j] >= cutoffIntensity)
        //    //        {
        //    //            for (int i = 0; i < data.GetLength(0); i++)
        //    //                newData[i, okIndex] = data[i, j];
        //    //            okIndex++;
        //    //        }
        //    //    data = newData;
        //    //}

        //    return new ThermoSpectrum(data);
        //}
    }
}