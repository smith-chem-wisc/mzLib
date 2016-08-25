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
using Spectra;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace IO.Thermo
{
    public class ThermoRawFile : MsDataFile<ThermoSpectrum>
    {
        internal enum RawLabelDataColumn
        {
            MZ = 0,
            Intensity = 1,
            Resolution = 2,
            NoiseBaseline = 3,
            NoiseLevel = 4,
            Charge = 5
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

        public enum Smoothing
        {
            None = 0,
            Boxcar = 1,
            Gauusian = 2
        }

        public enum IntensityCutoffType
        {
            None = 0,
            Absolute = 1,
            Relative = 2
        };

        private IXRawfile5 _rawConnection;

        public ThermoRawFile(string filePath)
            : base(filePath, true, MsDataFileType.ThermoRawFile)
        {
        }

        public override void Open()
        {
            if (_rawConnection != null)
                return;
            if (!File.Exists(FilePath) && !Directory.Exists(FilePath))
            {
                throw new IOException(string.Format("The MS data file {0} does not currently exist", FilePath));
            }

            _rawConnection = (IXRawfile5)new MSFileReader_XRawfile();
            _rawConnection.Open(FilePath);
            _rawConnection.SetCurrentController(0, 1); // first 0 is for mass spectrometer
        }

        protected override int GetFirstSpectrumNumber()
        {
            int spectrumNumber = 0;
            _rawConnection.GetFirstSpectrumNumber(ref spectrumNumber);
            return spectrumNumber;
        }

        protected override int GetLastSpectrumNumber()
        {
            int spectrumNumber = 0;
            _rawConnection.GetLastSpectrumNumber(ref spectrumNumber);
            return spectrumNumber;
        }

        private int GetParentSpectrumNumber(int spectrumNumber)
        {
            return Convert.ToInt32(Regex.Match(GetPrecursorID(spectrumNumber), @"\d+$").Value);
        }

        public string GetSofwareVersion()
        {
            string softwareVersion = null;
            _rawConnection.GetInstSoftwareVersion(ref softwareVersion);
            return softwareVersion;
        }

        private object GetExtraValue(int spectrumNumber, string filter)
        {
            object value = null;
            _rawConnection.GetTrailerExtraValueForScanNum(spectrumNumber, filter, ref value);
            return value;
        }

        private string GetScanFilter(int spectrumNumber)
        {
            string filter = null;
            _rawConnection.GetFilterForScanNum(spectrumNumber, ref filter);
            return filter;
        }


        private string GetSpectrumID(int spectrumNumber)
        {
            int pnControllerType = 0;
            int pnControllerNumber = 0;
            _rawConnection.GetCurrentController(ref pnControllerType, ref pnControllerNumber);
            return "controllerType=" + pnControllerType + " controllerNumber=" + pnControllerNumber + " scan=" + spectrumNumber;
        }


        private static readonly Regex PolarityRegex = new Regex(@"\+ ", RegexOptions.Compiled);

        private Polarity GetPolarity(int spectrumNumber)
        {
            string filter = GetScanFilter(spectrumNumber);
            return PolarityRegex.IsMatch(filter) ? Polarity.Positive : Polarity.Negative;
        }

        protected ThermoSpectrum GetSpectrumFromRawFile(int spectrumNumber)
        {
            try
            {
                return new ThermoSpectrum(GetLabeledData(spectrumNumber));
            }
            catch (ArgumentNullException)
            {
                return new ThermoSpectrum(GetUnlabeledData(spectrumNumber, true));
            }
        }

        private double[,] GetUnlabeledData(int spectrumNumber, bool useCentroid)
        {
            object massList = null;
            object peakFlags = null;
            int arraySize = -1;
            double centroidPeakWidth = 0.001;
            _rawConnection.GetMassListFromScanNum(ref spectrumNumber, null, 0, 0, 0, Convert.ToInt32(useCentroid), ref centroidPeakWidth, ref massList, ref peakFlags, ref arraySize);
            return (double[,])massList;
        }

        private double[,] GetLabeledData(int spectrumNumber)
        {
            object labels = null;
            object flags = null;
            _rawConnection.GetLabelData(ref labels, ref flags, ref spectrumNumber);
            double[,] data = labels as double[,];
            if (data == null || data.Length == 0)
                throw new ArgumentNullException("For spectrum number " + spectrumNumber + " the data is null!");
            return data;
        }

        private MZAnalyzerType GetMzAnalyzer(int spectrumNumber)
        {
            int mzanalyzer = 0;
            _rawConnection.GetMassAnalyzerTypeForScanNum(spectrumNumber, ref mzanalyzer);

            switch ((ThermoMzAnalyzer)mzanalyzer)
            {
                case ThermoMzAnalyzer.FTMS:
                    return MZAnalyzerType.Orbitrap;
                default:
                    return MZAnalyzerType.Unknown;
            }
        }

        private double GetPrecursorMonoisotopicMZ(int spectrumNumber)
        {
            int parentScanNumber = GetParentSpectrumNumber(spectrumNumber);
            var ms1Spectrum = GetScan(parentScanNumber).MassSpectrum;
            double trailerMZ = GetPrecursorMonoisotopicMZfromTrailierExtra(spectrumNumber);
            if (trailerMZ == -1)
                return GetSelectedIonMZ(spectrumNumber);
            else
            {
                return ms1Spectrum.GetClosestPeak(trailerMZ).MZ;
            }
        }

        private double GetPrecursorMonoisotopicMZfromTrailierExtra(int scanNumber)
        {
            object labels_obj = null;
            object values_obj = null;
            int array_size = -1;
            _rawConnection.GetTrailerExtraForScanNum(scanNumber, ref labels_obj, ref values_obj, ref array_size);
            string[] labels = (string[])labels_obj;
            string[] values = (string[])values_obj;
            for (int i = labels.GetLowerBound(0); i <= labels.GetUpperBound(0); i++)
            {
                if (labels[i].StartsWith("Monoisotopic M/Z"))
                {
                    double monoisotopic_mz = double.Parse(values[i], CultureInfo.InvariantCulture);
                    if (monoisotopic_mz > 0.0)
                        return monoisotopic_mz;
                    else
                        return -1;
                }
            }
            return -1;
        }

        private double GetIsolationWidth(int spectrumNumber)
        {
            object width = GetExtraValue(spectrumNumber, "MS2 Isolation Width:");
            return Convert.ToDouble(width);
        }

        public double GetElapsedScanTime(int spectrumNumber)
        {
            object elapsedScanTime = GetExtraValue(spectrumNumber, "Elapsed Scan Time (sec):");
            return Convert.ToDouble(elapsedScanTime);
        }

        public double GetTIC(int spectrumNumber)
        {
            int numberOfPackets = -1;
            double startTime = double.NaN;
            double lowMass = double.NaN;
            double highMass = double.NaN;
            double totalIonCurrent = double.NaN;
            double basePeakMass = double.NaN;
            double basePeakIntensity = double.NaN;
            int numberOfChannels = -1;
            int uniformTime = -1;
            double frequency = double.NaN;
            _rawConnection.GetScanHeaderInfoForScanNum(spectrumNumber, ref numberOfPackets, ref startTime, ref lowMass,
                ref highMass,
                ref totalIonCurrent, ref basePeakMass, ref basePeakIntensity,
                ref numberOfChannels, ref uniformTime, ref frequency);

            return totalIonCurrent;
        }

        private DissociationType GetDissociationType(int spectrumNumber, int msnOrder = 2)
        {
            int type = 0;
            _rawConnection.GetActivationTypeForScanNum(spectrumNumber, msnOrder, ref type);
            return (DissociationType)type;
        }

        private int GetPrecusorCharge(int spectrumNumber)
        {
            short charge = Convert.ToInt16(GetExtraValue(spectrumNumber, "Charge State:"));
            return charge * (int)GetPolarity(spectrumNumber);
        }

        public override int GetSpectrumNumber(double retentionTime)
        {
            int spectrumNumber = 0;
            _rawConnection.ScanNumFromRT(retentionTime, ref spectrumNumber);
            return spectrumNumber;
        }

        private double GetInjectionTime(int spectrumNumber)
        {
            object time = GetExtraValue(spectrumNumber, "Ion Injection Time (ms):");
            return Convert.ToDouble(time);
        }

        public string GetInstrumentName()
        {
            string name = null;
            _rawConnection.GetInstName(ref name);
            return name;
        }

        public string GetInstrumentModel()
        {
            string model = null;
            _rawConnection.GetInstModel(ref model);
            return model;
        }

        public Chromatogram GetTICChroma()
        {
            int nChroType1 = 1; //1=TIC 0=MassRange
            int nChroOperator = 0;
            int nChroType2 = 0;
            string bstrFilter = null;
            string bstrMassRanges1 = null;
            string bstrMassRanges2 = null;
            double dDelay = 0.0;
            double dStartTime = 0.0;
            double dEndTime = 0.0;
            int nSmoothingType = 1; //0=None 1=Boxcar 2=Gaussian
            int nSmoothingValue = 7;
            object pvarChroData = null;
            object pvarPeakFlags = null;
            int pnArraySize = 0;

            //(int nChroType1, int nChroOperator, int nChroType2, string bstrFilter, string bstrMassRanges1, string bstrMassRanges2, double dDelay, ref double pdStartTime, 
            //ref double pdEndTime, int nSmoothingType, int nSmoothingValue, ref object pvarChroData, ref object pvarPeakFlags, ref int pnArraySize);
            _rawConnection.GetChroData(nChroType1, nChroOperator, nChroType2, bstrFilter, bstrMassRanges1, bstrMassRanges2, dDelay, dStartTime, dEndTime, nSmoothingType, nSmoothingValue, ref pvarChroData, ref pvarPeakFlags, ref pnArraySize);

            double[,] pvarArray = (double[,])pvarChroData;

            return new Chromatogram(pvarArray);
        }

        private readonly static Regex _msxRegex = new Regex(@"([\d.]+)@", RegexOptions.Compiled);

        public List<double> GetMSXPrecursors(int spectrumNumber)
        {
            string scanheader = GetScanFilter(spectrumNumber);

            int msxnumber = -1;
            _rawConnection.GetMSXMultiplexValueFromScanNum(spectrumNumber, ref msxnumber);

            var matches = _msxRegex.Matches(scanheader);

            return (from Match match in matches select double.Parse(match.Groups[1].Value)).ToList();
        }

        private bool GetIsCentroid(int spectrumNumber)
        {
            int isCentroid = -1;
            _rawConnection.IsCentroidScanForScanNum(spectrumNumber, ref isCentroid);
            return isCentroid > 0;
        }

        private string GetPrecursorID(int spectrumNumber)
        {
            return GetSpectrumID(GetPrecursor(spectrumNumber));
        }

        private double GetSelectedIonIntensity(int scanNumber)
        {
            double mz = -1;
            _rawConnection.GetPrecursorMassForScanNum(scanNumber, 2, ref mz);
            return GetScan(GetPrecursor(scanNumber)).MassSpectrum.GetClosestPeak(mz).Intensity;
        }



        private int GetPrecursor(int spectrumNumber)
        {
            int ms_order = -1;
            while (spectrumNumber >= 1 && ms_order != 1)
            {
                _rawConnection.GetMSOrderForScanNum(spectrumNumber, ref ms_order);
                spectrumNumber--;
            }
            return spectrumNumber + 1;
        }

        protected override MsDataScan<ThermoSpectrum> GetMsDataScanFromFile(int spectrumNumber)
        {
            var precursorID = GetPrecursorID(spectrumNumber);

            int numberOfPackets = -1;
            double startTime = double.NaN;
            double lowMass = double.NaN;
            double highMass = double.NaN;
            double totalIonCurrent = double.NaN;
            double basePeakMass = double.NaN;
            double basePeakIntensity = double.NaN;
            int numberOfChannels = -1;
            int uniformTime = -1;
            double frequency = double.NaN;
            _rawConnection.GetScanHeaderInfoForScanNum(spectrumNumber, ref numberOfPackets, ref startTime, ref lowMass,
                ref highMass, ref totalIonCurrent, ref basePeakMass, ref basePeakIntensity,
                ref numberOfChannels, ref uniformTime, ref frequency);


            MzRange ScanWindowRange = new MzRange(lowMass, highMass);



            double retentionTime = 0;
            _rawConnection.RTFromScanNum(spectrumNumber, ref retentionTime);
            int msnOrder = 0;
            _rawConnection.GetMSOrderForScanNum(spectrumNumber, ref msnOrder);


            if (precursorID.Equals(GetSpectrumID(spectrumNumber)))
                return new MsDataScan<ThermoSpectrum>(spectrumNumber, GetSpectrumFromRawFile(spectrumNumber), GetSpectrumID(spectrumNumber), msnOrder, GetIsCentroid(spectrumNumber), GetPolarity(spectrumNumber), retentionTime, ScanWindowRange, GetScanFilter(spectrumNumber), GetMzAnalyzer(spectrumNumber), GetInjectionTime(spectrumNumber), totalIonCurrent);
            else
                return new MsDataScan<ThermoSpectrum>(spectrumNumber, GetSpectrumFromRawFile(spectrumNumber), GetSpectrumID(spectrumNumber), msnOrder, GetIsCentroid(spectrumNumber), GetPolarity(spectrumNumber), retentionTime, ScanWindowRange, GetScanFilter(spectrumNumber), GetMzAnalyzer(spectrumNumber), GetInjectionTime(spectrumNumber), totalIonCurrent, precursorID, GetSelectedIonMZ(spectrumNumber), GetPrecusorCharge(spectrumNumber), GetSelectedIonIntensity(spectrumNumber), GetIsolationMZ(spectrumNumber), GetIsolationWidth(spectrumNumber), GetDissociationType(spectrumNumber), GetParentSpectrumNumber(spectrumNumber), GetPrecursorMonoisotopicIntensity(spectrumNumber), GetPrecursorMonoisotopicMZ(spectrumNumber));
        }

        private double GetPrecursorMonoisotopicIntensity(int spectrumNumber)
        {
            int parentScanNumber = GetParentSpectrumNumber(spectrumNumber);
            var ms1Spectrum = GetScan(parentScanNumber).MassSpectrum;
            double trailerMZ = GetPrecursorMonoisotopicMZfromTrailierExtra(spectrumNumber);
            if (trailerMZ == -1)
                return GetSelectedIonIntensity(spectrumNumber);
            else
            {
                return ms1Spectrum.GetClosestPeak(trailerMZ).Intensity;
            }
        }

        private double GetSelectedIonMZ(int spectrumNumber)
        {
            double mz = double.NaN;
            _rawConnection.GetPrecursorMassForScanNum(spectrumNumber, 2, ref mz);

            int parentScanNumber = GetParentSpectrumNumber(spectrumNumber);
            var ms1Spectrum = GetScan(parentScanNumber).MassSpectrum;
            MzPeak peak = ms1Spectrum.GetClosestPeak(mz);
            return peak.MZ;
        }

        private double GetIsolationMZ(int spectrumNumber)
        {
            double mz = double.NaN;
            _rawConnection.GetPrecursorMassForScanNum(spectrumNumber, 2, ref mz);
            return mz;
        }

        public bool monoisotopicPrecursorSelectionEnabled
        {
            get
            {
                int n = 0;
                _rawConnection.GetNumInstMethods(ref n);
                Console.WriteLine("number of instrument methods:" + n);

                string s;
                for (int i = 0; i < n; i++)
                {
                    s = null;
                    _rawConnection.GetInstMethod(i, ref s);
                    if (Regex.IsMatch(s, "Monoisotopic precursor selection enabled"))
                        return true;
                }
                return false;
            }
        }
    }
}