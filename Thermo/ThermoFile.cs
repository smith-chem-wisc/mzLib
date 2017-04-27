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
        private static readonly Regex mFindParentIonOnlyNonMsx = new Regex(@"[Mm][Ss]\d*[^\[\r\n]* (?<ParentMZ>[0-9.]+)@?[A-Za-z]*\d*\.?\d*(\[[^\]\r\n]\])?", RegexOptions.IgnoreCase | RegexOptions.Compiled);
        private static readonly Regex mFindParentIonOnlyMsx = new Regex(@"[Mm][Ss]\d* (?<ParentMZ>[0-9.]+)@?[A-Za-z]*\d*\.?\d*[^\[\r\n]*(\[[^\]\r\n]+\])?", RegexOptions.IgnoreCase | RegexOptions.Compiled);

        #endregion Private Fields

        #region Public Constructors

        public ThermoFile(IThermoScan[] scans, ThermoGlobalParams thermoGlobalParams) : base(scans)
        {
            this.ThermoGlobalParams = thermoGlobalParams;
        }

        public ThermoFile(IXRawfile5 _rawConnection, int numSpectra, ManagedThermoHelperLayer.PrecursorInfo[] couldBePrecursor) : base(numSpectra)
        {
            this.ThermoGlobalParams = GetAllGlobalStuff(_rawConnection, couldBePrecursor);
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

        public ThermoGlobalParams ThermoGlobalParams { get; }

        #endregion Public Properties

        #region Public Methods

        public static ThermoGlobalParams GetAllGlobalStuff(IXRawfile5 _rawConnection, ManagedThermoHelperLayer.PrecursorInfo[] couldBePrecursor)
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

            return new ThermoGlobalParams(pnNumInstMethods, instrumentMethods, pbstrInstSoftwareVersion, pbstrInstName, pbstrInstModel, pnControllerType, pnControllerNumber, couldBePrecursor);
        }

        #endregion Public Methods

        #region Protected Methods

        protected static IThermoScan GetMsDataOneBasedScanFromThermoFile(int nScanNumber, IXRawfile5 theConnection, ThermoGlobalParams globalParams)
        {
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

            double? ms2isolationWidthFromTrailerExtra = null;
            double? injectionTimeFromTrailerExtra = null;
            double? precursorMonoisotopicMZfromTrailierExtra = null;
            int? chargeStatefromTrailierExtra = null;

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
                    ms2isolationWidthFromTrailerExtra = double.Parse(values[i], CultureInfo.InvariantCulture) == 0 ?
                        (double?)null :
                        double.Parse(values[i], CultureInfo.InvariantCulture);
                }
                if (labels[i].StartsWith("Ion Injection Time (ms)", StringComparison.Ordinal))
                {
                    injectionTimeFromTrailerExtra = double.Parse(values[i], CultureInfo.InvariantCulture) == 0 ?
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
                    chargeStatefromTrailierExtra = int.Parse(values[i], CultureInfo.InvariantCulture) == 0 ?
                        (int?)null :
                        int.Parse(values[i], CultureInfo.InvariantCulture);
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

            object pvarNoisePacket = null;
            theConnection.GetNoiseData(ref pvarNoisePacket, nScanNumber);
            double[,] noiseData = pvarNoisePacket as double[,];

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

                // Trust this first
                var precursorInfo = globalParams.couldBePrecursor[nScanNumber - 1];

                // THIS METHOD IS BUGGY!!! DO NOT USE
                //theConnection.FindPrecursorMassInFullScan(nScanNumber, ref pnMasterScan, ref pdFoundMass, ref pdHeaderMass, ref pnChargeState);

                // TODO: Get it from more locations
                int oneBasedPrecursorScanNumber;
                if (precursorInfo.nScanNumber > 0)
                    oneBasedPrecursorScanNumber = precursorInfo.nScanNumber;
                else
                {
                    oneBasedPrecursorScanNumber = nScanNumber - 1;
                    // Use info from ScanEvent! Loop back until scan event is equal to 1
                    while (true)
                    {
                        if (globalParams.scanEvent[oneBasedPrecursorScanNumber - 1] == 0)
                        {
                            object pvarValuesHere = null;
                            object pvarLablesHere = null;
                            int pnArraySizeHere = 0;
                            theConnection.GetTrailerExtraForScanNum(nScanNumber, ref pvarLablesHere, ref pvarValuesHere, ref pnArraySizeHere);
                            string[] labelsHere = (string[])pvarLablesHere;
                            string[] valuesHere = (string[])pvarValuesHere;
                            for (int i = labelsHere.GetLowerBound(0); i <= labelsHere.GetUpperBound(0); i++)
                                if (labelsHere[i].StartsWith("Scan Event", StringComparison.Ordinal))
                                    globalParams.scanEvent[oneBasedPrecursorScanNumber - 1] = int.Parse(valuesHere[i], CultureInfo.InvariantCulture);
                        }
                        if (globalParams.scanEvent[oneBasedPrecursorScanNumber - 1] == 1)
                            break;
                        oneBasedPrecursorScanNumber--;
                    }
                }

                int? selectedIonGuessChargeStateGuess = null;
                if (precursorInfo.nChargeState > 0)
                    selectedIonGuessChargeStateGuess = precursorInfo.nChargeState;
                else if (chargeStatefromTrailierExtra.HasValue)
                    selectedIonGuessChargeStateGuess = chargeStatefromTrailierExtra;

                double? selectedIonGuessMonoisotopicMz = null;
                if (precursorMonoisotopicMZfromTrailierExtra.HasValue && precursorMonoisotopicMZfromTrailierExtra.Value > 0)
                    selectedIonGuessMonoisotopicMz = precursorMonoisotopicMZfromTrailierExtra;
                if (precursorInfo.dMonoIsoMass > 0 && !selectedIonGuessMonoisotopicMz.HasValue)
                    selectedIonGuessMonoisotopicMz = precursorInfo.dMonoIsoMass;

                // TODO: Get it from more locations
                double selectedIonGuessMZ;
                if (precursorInfo.dIsolationMass > 0)
                    selectedIonGuessMZ = precursorInfo.dIsolationMass;
                else
                {
                    Regex matcher;
                    if (pbstrFilter.ToLower().Contains("msx"))
                        matcher = mFindParentIonOnlyMsx;
                    else
                        matcher = mFindParentIonOnlyNonMsx;
                    selectedIonGuessMZ = double.Parse(matcher.Match(pbstrFilter).Groups["ParentMZ"].Value);
                }

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
                    selectedIonGuessMZ,
                    selectedIonGuessChargeStateGuess,
                    ms2isolationWidthFromTrailerExtra,
                    (DissociationType)pnActivationType,
                    oneBasedPrecursorScanNumber,
                    selectedIonGuessMonoisotopicMz,
                    injectionTimeFromTrailerExtra,
                    noiseData);
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
                    injectionTimeFromTrailerExtra,
                    noiseData);
            }
        }

        #endregion Protected Methods

    }
}