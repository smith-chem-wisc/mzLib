// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work Copyright 2016 Stefan Solntsev
//
// This file (ThermoRawFile.cs) is part of MassSpecFiles.
//
// MassSpecFiles is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MassSpecFiles is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with MassSpecFiles. If not, see <http://www.gnu.org/licenses/>.

using MassSpectrometry;
using MSFileReaderLib;
using MzLibUtil;
using System;
using System.Globalization;
using System.IO;
using System.Runtime.InteropServices;
using System.Security.Cryptography;
using System.Text.RegularExpressions;

namespace IO.Thermo
{
    public class ThermoStaticData : ThermoDataFile
    {
        private static readonly Regex PolarityRegex = new Regex(@"\+ ", RegexOptions.Compiled);
        private static readonly Regex mFindParentIonOnlyNonMsx = new Regex(@"[Mm][Ss]\d*[^\[\r\n]* (?<ParentMZ>[0-9.]+)@?[A-Za-z]*\d*\.?\d*(\[[^\]\r\n]\])?", RegexOptions.IgnoreCase | RegexOptions.Compiled);
        private static readonly Regex mFindParentIonOnlyMsx = new Regex(@"[Mm][Ss]\d* (?<ParentMZ>[0-9.]+)@?[A-Za-z]*\d*\.?\d*[^\[\r\n]*(\[[^\]\r\n]+\])?", RegexOptions.IgnoreCase | RegexOptions.Compiled);

        private ThermoStaticData(MsDataScan[] scans, ThermoGlobalParams thermoGlobalParams, SourceFile sourceFile) : base(scans, sourceFile)
        {
            this.ThermoGlobalParams = thermoGlobalParams;
        }

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

        public static ThermoStaticData LoadAllStaticData(string filePath, IFilteringParams filterParams = null)
        {
            if (!File.Exists(filePath))
            {
                throw new FileNotFoundException();
            }

            if (CheckForMsFileReader() == false)
            {
                throw new MzLibException("MsFileReader Not Installed");
            }

            var ok = new ManagedThermoHelperLayer.HelperClass();
            IXRawfile5 theConnection = (IXRawfile5)new MSFileReader_XRawfile();
            theConnection.Open(filePath);
            int pbSMData = 0;
            theConnection.IsThereMSData(ref pbSMData);
            if (pbSMData == 0)
            {
                throw new MzLibException("Could not read data from file " + Path.GetFileName(filePath));
            }

            theConnection.SetCurrentController(0, 1);

            var precursorInfosArray = ok.GetAllPrecursorInfos(filePath);
            for (int i = 0; i < precursorInfosArray.Length; i++)
            {
                if (precursorInfosArray[i].nScanNumber == 0)
                {
                    precursorInfosArray[i].nScanNumber = -1;
                }
            }

            int pnFirstSpectrum = 0;
            theConnection.GetFirstSpectrumNumber(ref pnFirstSpectrum);
            int pnLastSpectrum = 0;
            theConnection.GetLastSpectrumNumber(ref pnLastSpectrum);

            ThermoGlobalParams p = GetAllGlobalStuff(theConnection, precursorInfosArray, filePath);

            MsDataScan[] scans = new MsDataScan[pnLastSpectrum - pnFirstSpectrum + 1];
            for (int nScanNumber = pnFirstSpectrum; nScanNumber <= pnLastSpectrum; nScanNumber++)
            {
                scans[nScanNumber - pnFirstSpectrum] = GetMsDataOneBasedScanFromThermoFile(theConnection, nScanNumber, p, filterParams);
            }

            theConnection.Close();

            string sendCheckSum;
            using (FileStream stream = File.OpenRead(filePath))
            {
                using (SHA1Managed sha = new SHA1Managed())
                {
                    byte[] checksum = sha.ComputeHash(stream);
                    sendCheckSum = BitConverter.ToString(checksum)
                        .Replace("-", string.Empty);
                }
            }
            SourceFile sourceFile = new SourceFile(
                @"Thermo nativeID format",
                @"Thermo RAW format",
                sendCheckSum,
                @"SHA-1",
                filePath,
                Path.GetFileNameWithoutExtension(filePath));

            return new ThermoStaticData(scans, p, sourceFile);
        }

        public static MsDataScan GetMsDataOneBasedScanFromThermoFile(IXRawfile5 theConnection, int nScanNumber, ThermoGlobalParams globalParams, IFilteringParams filterParams = null)
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
            try //if there is no noise data
            {
                theConnection.GetNoiseData(ref pvarNoisePacket, nScanNumber);
            }
            catch
            {
                //pvarNoisePAcket is already null
            }
            double[,] noiseData = pvarNoisePacket as double[,];

            double[,] data;
            try
            {
                object pvarFlags = null;
                object pvarLabels = null;
                theConnection.GetLabelData(ref pvarLabels, ref pvarFlags, ref nScanNumber);
                data = pvarLabels as double[,];
                if (data == null || data.Length == 0)
                {
                    throw new MzLibException("For spectrum number " + nScanNumber + " the data is null!");
                }
            }
            catch (MzLibException)
            {
                // Warning: the masses reported by GetMassListFromScanNum when centroiding are not properly calibrated and thus could be off by 0.3 m/z or more

                double pdCentroidPeakWidth = 0;
                object pvarnMassList = null;
                object pvarPeakFlags = null;
                theConnection.GetMassListFromScanNum(ref nScanNumber, null, 0, 0, 0, 1, ref pdCentroidPeakWidth, ref pvarnMassList, ref pvarPeakFlags, ref pnArraySize);
                data = (double[,])pvarnMassList;
            }

            MzSpectrum thermoSpectrum;
            if (filterParams != null && data.GetLength(1) > 0 && (filterParams.MinimumAllowedIntensityRatioToBasePeakM.HasValue || filterParams.NumberOfPeaksToKeepPerWindow.HasValue) && ((filterParams.ApplyTrimmingToMs1 && pnMSOrder == 1) || (filterParams.ApplyTrimmingToMsMs && pnMSOrder > 1)))
            {
                var count = data.GetLength(1);

                var mzArray = new double[count];
                var intensityArray = new double[count];
                Buffer.BlockCopy(data, 0, mzArray, 0, sizeof(double) * count);
                Buffer.BlockCopy(data, sizeof(double) * count, intensityArray, 0, sizeof(double) * count);
                if (filterParams.NumberOfWindows == null)
                {
                    int numPeaks = TopNpeakHelper(ref intensityArray, ref mzArray, filterParams);
                    //the following arrays are modified after TopN helper
                    Array.Resize(ref intensityArray, numPeaks);
                    Array.Resize(ref mzArray, numPeaks);
                }
                //Array reference passed by value, array calues will be modified after calling
                else
                {
                    WindowModeHelper(ref intensityArray, ref mzArray, filterParams);
                }
                Array.Sort(mzArray, intensityArray);
                thermoSpectrum = new MzSpectrum(mzArray, intensityArray, false);
            }
            else
            {
                thermoSpectrum = new MzSpectrum(data);
            }
            MZAnalyzerType mzAnalyzerType;
            if ((ThermoMzAnalyzer)pnMassAnalyzerType == ThermoMzAnalyzer.FTMS)
            {
                mzAnalyzerType = MZAnalyzerType.Orbitrap;
            }
            else
            {
                mzAnalyzerType = MZAnalyzerType.Unknown;
            }
            string nativeId = "controllerType=0 controllerNumber=1 scan=" + nScanNumber;

            if (pnMSOrder > 1)
            {
                int pnActivationType = 0;
                theConnection.GetActivationTypeForScanNum(nScanNumber, pnMSOrder, ref pnActivationType);

                // INITIALIZE globalParams.couldBePrecursor[nScanNumber - 1] (for dynamic connections that don't have it initialized yet)
                if (globalParams.CouldBePrecursor[nScanNumber - 1].Equals(default(ManagedThermoHelperLayer.PrecursorInfo)))
                {
                    var ok = new ManagedThermoHelperLayer.HelperClass();
                    globalParams.CouldBePrecursor[nScanNumber - 1] = ok.GetSingleScanPrecursorInfo(nScanNumber, globalParams.FilePath);
                }

                var precursorInfo = globalParams.CouldBePrecursor[nScanNumber - 1];

                // THIS METHOD IS BUGGY!!! DO NOT USE
                //theConnection.FindPrecursorMassInFullScan(nScanNumber, ref pnMasterScan, ref pdFoundMass, ref pdHeaderMass, ref pnChargeState);

                int oneBasedPrecursorScanNumber = -1;
                if (precursorInfo.nScanNumber > 0)
                {
                    oneBasedPrecursorScanNumber = precursorInfo.nScanNumber;
                }
                else if (masterScanfromTrailierExtra.HasValue && masterScanfromTrailierExtra > 0)
                {
                    oneBasedPrecursorScanNumber = masterScanfromTrailierExtra.Value;
                }
                else
                {
                    // we weren't able to get the precursor scan number, so we'll have to guess;
                    // loop back to find precursor scan
                    // (assumed to be the first scan before this scan with an MS order of this scan's MS order - 1)
                    // e.g., if this is an MS2 scan, find the first MS1 scan before this and assume that's the precursor scan
                    int scanOrder = globalParams.MsOrderByScan[nScanNumber - 1];
                    int precursorScanOrder = scanOrder - 1;

                    for (int i = nScanNumber - 1; i >= 0; i--)
                    {
                        int msOrder = globalParams.MsOrderByScan[i];

                        if (msOrder == precursorScanOrder)
                        {
                            oneBasedPrecursorScanNumber = i + 1;
                            break;
                        }
                    }
                }
                if (oneBasedPrecursorScanNumber == -1)
                {
                    throw new MzLibException("Could not find precursor info for scan #" + nScanNumber);
                }

                int? selectedIonGuessChargeStateGuess = null;
                if (precursorInfo.nChargeState > 0)
                {
                    selectedIonGuessChargeStateGuess = precursorInfo.nChargeState;
                }
                else if (chargeStatefromTrailierExtra.HasValue)
                {
                    selectedIonGuessChargeStateGuess = chargeStatefromTrailierExtra;
                }

                double? selectedIonGuessMonoisotopicMz = null;
                if (precursorMonoisotopicMZfromTrailierExtra.HasValue && precursorMonoisotopicMZfromTrailierExtra.Value > 0)
                {
                    selectedIonGuessMonoisotopicMz = precursorMonoisotopicMZfromTrailierExtra;
                }
                if (precursorInfo.dMonoIsoMass > 0 && !selectedIonGuessMonoisotopicMz.HasValue)
                {
                    selectedIonGuessMonoisotopicMz = precursorInfo.dMonoIsoMass;
                }

                Regex matcher;
                if (pbstrFilter.ToLower().Contains("msx"))
                {
                    matcher = mFindParentIonOnlyMsx;
                }
                else
                {
                    matcher = mFindParentIonOnlyNonMsx;
                }
                double selectedIonGuessMZ = double.Parse(matcher.Match(pbstrFilter).Groups["ParentMZ"].Value);

                //   int? selectedIonChargeStateGuess, double? selectedIonIntensity, double? isolationMZ, double? isolationWidth, DissociationType dissociationType, int? oneBasedPrecursorScanNumber, double? selectedIonMonoisotopicGuessMz, double? injectionTime, double[,] noiseData, string nativeId)
                // double TotalIonCurrent, double selectedIonMZ, int? selectedIonChargeStateGuess, double? selectedIonIntensity, double? isolationMZ, double? isolationWidth, DissociationType dissociationType, int? oneBasedPrecursorScanNumber, double? selectedIonMonoisotopicGuessMz, double? injectionTime, double[,] noiseData, string nativeId)
                return new MsDataScan(
                    thermoSpectrum,
                    nScanNumber,
                    pnMSOrder,
                    true,
                    PolarityRegex.IsMatch(pbstrFilter) ? Polarity.Positive : Polarity.Negative,
                    pdStartTime,
                    new MzRange(pdLowMass, pdHighMass),
                    pbstrFilter,
                    mzAnalyzerType,
                    pdTIC,
                    injectionTimeFromTrailerExtra,
                    noiseData,
                    nativeId,
                    selectedIonGuessMZ,
                    selectedIonGuessChargeStateGuess,
                    null,
                    selectedIonGuessMZ,
                    ms2isolationWidthFromTrailerExtra,
                    (DissociationType)pnActivationType,
                    oneBasedPrecursorScanNumber,
                    selectedIonGuessMonoisotopicMz
                    );
            }
            else
            {
                return new MsDataScan(
                    thermoSpectrum,
                    nScanNumber,
                    1,
                    true,
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
            {
                _rawConnection.GetMSOrderForScanNum((i + 1), ref msOrderByScan[i]);
            }

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
            catch (COMException)
            {
                return false;
            }
            return true;
        }
    }
}