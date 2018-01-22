using MassSpectrometry;
using MSFileReaderLib;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text.RegularExpressions;

namespace IO.Thermo
{
    public abstract class ThermoFile : MsDataFile<IThermoScan>
    {
        #region Private Fields
        #region test data
        private static int? windowSize;
        private static bool windowMode;
        #endregion
        private static readonly Regex PolarityRegex = new Regex(@"\+ ", RegexOptions.Compiled);
        private static readonly Regex mFindParentIonOnlyNonMsx = new Regex(@"[Mm][Ss]\d*[^\[\r\n]* (?<ParentMZ>[0-9.]+)@?[A-Za-z]*\d*\.?\d*(\[[^\]\r\n]\])?", RegexOptions.IgnoreCase | RegexOptions.Compiled);
        private static readonly Regex mFindParentIonOnlyMsx = new Regex(@"[Mm][Ss]\d* (?<ParentMZ>[0-9.]+)@?[A-Za-z]*\d*\.?\d*[^\[\r\n]*(\[[^\]\r\n]+\])?", RegexOptions.IgnoreCase | RegexOptions.Compiled);

        #endregion Private Fields

        #region Protected Constructors

        protected ThermoFile(IThermoScan[] scans, ThermoGlobalParams thermoGlobalParams, SourceFile sourceFile) : base(scans, sourceFile)
        {
            this.ThermoGlobalParams = thermoGlobalParams;
        }

        protected ThermoFile(IXRawfile5 _rawConnection, int numSpectra, ManagedThermoHelperLayer.PrecursorInfo[] couldBePrecursor, SourceFile sourceFile, ThermoGlobalParams thermoGlobalParams) : base(numSpectra, sourceFile)
        {
            this.ThermoGlobalParams = thermoGlobalParams;
        }

        #endregion Protected Constructors

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

        public static ThermoGlobalParams GetAllGlobalStuff(IXRawfile5 _rawConnection, ManagedThermoHelperLayer.PrecursorInfo[] couldBePrecursor, string filePath)
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

            int[] msOrderByScan = new int[couldBePrecursor.Length];
            for (int i = 0; i < couldBePrecursor.Length; i++)
                _rawConnection.GetMSOrderForScanNum((i + 1), ref msOrderByScan[i]);

            return new ThermoGlobalParams(pnNumInstMethods, instrumentMethods, pbstrInstSoftwareVersion, pbstrInstName, pbstrInstModel, pnControllerType, pnControllerNumber, couldBePrecursor, filePath, msOrderByScan);
        }

        public static bool CheckForMsFileReader()
        {
            const string THERMO_READER_CLSID = "{1d23188d-53fe-4c25-b032-dc70acdbdc02}";
            //Check if Thermo File Reader Exists
            try
            {
                var thermoReader = Type.GetTypeFromCLSID(Guid.Parse(THERMO_READER_CLSID));
                Activator.CreateInstance(thermoReader);
            }
            catch (COMException ex)
            {
                return false;
            }
            return true;
        }

        #endregion Public Methods

        #region Protected Methods
        protected static IThermoScan GetMsDataOneBasedScanFromThermoFile(
            int nScanNumber, IXRawfile5 theConnection, ThermoGlobalParams globalParams,
            FilteringParams ThermoParams, bool trimMs1Peaks, bool trimMsMsPeaks)
        {
            FilteringParams parameters = new FilteringParams(1, 2, false, 10);
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
            int? masterScanfromTrailierExtra = null;

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
                if (labels[i].StartsWith("Master Scan Number", StringComparison.Ordinal))
                {
                    masterScanfromTrailierExtra = int.Parse(values[i], CultureInfo.InvariantCulture) == 0 ?
                        (int?)null :
                        int.Parse(values[i], CultureInfo.InvariantCulture);
                }
            }

            string pbstrFilter = null;
            theConnection.GetFilterForScanNum(nScanNumber, ref pbstrFilter);

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
                    throw new MzLibException("For spectrum number " + nScanNumber + " the data is null!");
            }
            catch (MzLibException)
            {
                // Warning: the masses reported by GetMassListFromScanNum when centroiding are not properly calibrated and thus could be off by 0.3 m/z or more

                double pdCentroidPeakWidth = 0;
                object pvarnMassList = null;
                object pvarPeakFlags = null;
                theConnection.GetMassListFromScanNum(ref nScanNumber,
                                null,
                                0,
                                0,
                                0,
                                1,
                                ref pdCentroidPeakWidth,
                                ref pvarnMassList,
                                ref pvarPeakFlags,
                                ref pnArraySize);
                data = (double[,])pvarnMassList;
            }

            ThermoSpectrum thermoSpectrum;
            //modification Mark
            if ((ThermoParams.minRatio.HasValue || ThermoParams.topNpeaks.HasValue) && ((trimMs1Peaks && pnMSOrder == 1) || (trimMsMsPeaks && pnMSOrder > 1)))
            {
                var count = data.GetLength(1);

                var mzArray = new double[count];
                var intensityArray = new double[count];
                Buffer.BlockCopy(data, 0, mzArray, 0, sizeof(double) * count);
                Buffer.BlockCopy(data, sizeof(double) * count, intensityArray, 0, sizeof(double) * count);
                if (!windowMode)
                {
                    int numPeaks = ThermoParams.TopNpeakHelper(intensityArray, mzArray);
                    //the following arrays are modified after TopN helper
                    Array.Resize(ref intensityArray, numPeaks);
                    Array.Resize(ref mzArray, numPeaks);


                }
                //Array reference passed by value, array calues will be modified after calling
                else
                {
                    int temp = count / windowSize.Value;
                    var mzTemp = new double[temp];
                    var intensityTemp = new double[temp];
                    List<double> mzResults = new List<double>();
                    List<double> intensityResults = new List<double>();

                    for (int i = 0; i < windowSize; i++)
                    {
                        Buffer.BlockCopy(mzArray, sizeof(double) * temp * i, mzTemp, 0, sizeof(double) * temp);
                        Buffer.BlockCopy(intensityArray, sizeof(double) * temp * i, intensityTemp, 0, sizeof(double) * temp);
                        int numPeaks = ThermoParams.TopNpeakHelper(intensityTemp, mzTemp);
                        Array.Resize(ref intensityTemp, numPeaks);
                        Array.Resize(ref mzTemp, numPeaks);
                        mzResults.AddRange(mzTemp);
                        intensityResults.AddRange(intensityTemp);
                    }
                    mzArray = mzResults.ToArray();
                }
                Array.Sort(mzArray, intensityArray);
                thermoSpectrum = new ThermoSpectrum(mzArray, intensityArray, false);
            }
            else
                thermoSpectrum = new ThermoSpectrum(data);
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
            string nativeId = "controllerType=0 controllerNumber=1 scan=" + nScanNumber;

            if (pnMSOrder > 1)
            {
                int pnActivationType = 0;
                theConnection.GetActivationTypeForScanNum(nScanNumber, pnMSOrder, ref pnActivationType);

                // INITIALIZE globalParams.couldBePrecursor[nScanNumber - 1] (for dynamic connections that don't have it initialized yet)
                if (globalParams.couldBePrecursor[nScanNumber - 1].Equals(default(ManagedThermoHelperLayer.PrecursorInfo)))
                {
                    var ok = new ManagedThermoHelperLayer.HelperClass();
                    globalParams.couldBePrecursor[nScanNumber - 1] = ok.GetSingleScanPrecursorInfo(nScanNumber, globalParams.filePath);
                }

                var precursorInfo = globalParams.couldBePrecursor[nScanNumber - 1];

                // THIS METHOD IS BUGGY!!! DO NOT USE
                //theConnection.FindPrecursorMassInFullScan(nScanNumber, ref pnMasterScan, ref pdFoundMass, ref pdHeaderMass, ref pnChargeState);

                int oneBasedPrecursorScanNumber;
                if (precursorInfo.nScanNumber > 0)
                    oneBasedPrecursorScanNumber = precursorInfo.nScanNumber;
                else if (masterScanfromTrailierExtra.HasValue)
                    oneBasedPrecursorScanNumber = masterScanfromTrailierExtra.Value;
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
                            theConnection.GetTrailerExtraForScanNum(oneBasedPrecursorScanNumber, ref pvarLablesHere, ref pvarValuesHere, ref pnArraySizeHere);
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

                Regex matcher;
                if (pbstrFilter.ToLower().Contains("msx"))
                    matcher = mFindParentIonOnlyMsx;
                else
                    matcher = mFindParentIonOnlyNonMsx;
                double selectedIonGuessMZ = double.Parse(matcher.Match(pbstrFilter).Groups["ParentMZ"].Value);

                return new ThermoScanWithPrecursor(
                    nScanNumber,
                    thermoSpectrum,
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
                    noiseData,
                    nativeId);
            }
            else
            {
                return new ThermoScan(nScanNumber,
                    thermoSpectrum,
                    1,
                    PolarityRegex.IsMatch(pbstrFilter) ? Polarity.Positive : Polarity.Negative,
                    pdStartTime,
                    new MzRange(pdLowMass, pdHighMass),
                    pbstrFilter,
                    mzAnalyzerType,
                    pdTIC,
                    injectionTimeFromTrailerExtra,
                    noiseData,
                    nativeId);
            }
        }

        #endregion Protected Methods
    }
}