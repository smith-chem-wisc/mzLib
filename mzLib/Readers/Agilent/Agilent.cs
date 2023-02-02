using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using UsefulProteomicsDatabases;
using Agilent.MassSpectrometry.DataAnalysis;
using MzLibUtil;
using System.Text.RegularExpressions;

namespace Readers
{
    public class Agilent : MsDataFile
    {
        private IMsdrDataReader _reader;

        public Agilent()
        {
            InitiateDynamicConnection();
        }

        public Agilent(string filepath) : base(filepath)
        {
            InitiateDynamicConnection();
        }

        public override MsDataFile LoadAllStaticData(FilteringParams filteringParams = null, int maxThreads = 1)
        {
            if (!File.Exists(FilePath))
            {
                throw new FileNotFoundException();
            }

            Loaders.LoadElements();


            throw new NotImplementedException();
        }

        public override SourceFile GetSourceFile()
        {
            throw new NotImplementedException();
        }

        public override MsDataScan GetOneBasedScanFromDynamicConnection(int oneBasedScanNumber, IFilteringParams filterParams = null)
        {
            if (_reader == null)
            {
                throw new MzLibException("Cannot get scan; the dynamic data connection to " + FilePath + " has been closed!");
            }

            // MS1
            MzSpectrum mzSpectrum = GetSpectrum(oneBasedScanNumber);
            int oneBasedIndex;
            int msnOrder = GetMsnOrder(oneBasedScanNumber);
            bool isCentroid;
            Polarity polarity = GetPolarity(oneBasedScanNumber);
            double retentionTime = GetRetentionTime(oneBasedScanNumber);
            MzRange range = GetMzRange(oneBasedScanNumber);
            string scanFilter;
            MZAnalyzerType analyzer = GetMzAnalyzer(oneBasedScanNumber);
            double tic;
            double injectionTime = GetInjectionTime(oneBasedScanNumber);
            string nativeId;



            // MS2
            double selectedIonMz = GetPrecursorMz(oneBasedScanNumber);
            double selectedIonCharge = GetPrecusorCharge(oneBasedScanNumber);
            double selectedIonIntensity;
            double isolationMz;
            MzRange isolationWidth;
            DissociationType dissociationType = GetDissociationType(oneBasedScanNumber);
            double precursorScanNumber;
            double monoisotopicMz;

            throw new NotImplementedException();
        }

        public override void CloseDynamicConnection()
        {
            if (_reader == null) return;
            _reader.CloseDataFile();

        }

        public sealed override void InitiateDynamicConnection()
        {
            if (!File.Exists(FilePath))
                throw new FileNotFoundException();
            _reader = new MassSpecDataReader();
            _reader.OpenDataFile(FilePath);
            
        }

        #region Private Helpers

        public MzSpectrum GetSpectrum(int spectrumNumber)
        {
            IBDASpecData spectrum = _reader.GetSpectrum(spectrumNumber - 1);

            double[] doubleArray = new double[spectrum.YArray.Length];
            for (int i = 0; i < doubleArray.Length; i++)
            {
                doubleArray[i] = spectrum.YArray[i];
            }

            return new MzSpectrum(spectrum.XArray, doubleArray, true);
        }

        public double GetRetentionTime(int spectrumNumber)
        {
            IMSScanRecord scan_record = _reader.GetScanRecord(spectrumNumber - 1);
            return scan_record.RetentionTime;
        }

        public int GetMsnOrder(int spectrumNumber)
        {
            IMSScanRecord scan_record = _reader.GetScanRecord(spectrumNumber - 1);
            return scan_record.MSLevel == MSLevel.MSMS ? 2 : 1;
        }

        public Polarity GetPolarity(int spectrumNumber)
        {
            IMSScanRecord scan_record = _reader.GetScanRecord(spectrumNumber - 1);
            switch (scan_record.IonPolarity)
            {
                case IonPolarity.Positive:
                    return Polarity.Positive;
                case IonPolarity.Negative:
                    return Polarity.Negative;
                default:
                    return Polarity.Unknown;
            }
        }

        public MZAnalyzerType GetMzAnalyzer(int spectrumNumber)
        {
            IBDASpecData spectrum = _reader.GetSpectrum(spectrumNumber - 1);
            switch (spectrum.DeviceType)
            {
                case DeviceType.IonTrap:
                    return MZAnalyzerType.IonTrap3D;
                case DeviceType.Quadrupole:
                    return MZAnalyzerType.Quadrupole;
                case DeviceType.TandemQuadrupole:
                    return MZAnalyzerType.Quadrupole;
                case DeviceType.TimeOfFlight:
                    return MZAnalyzerType.TOF;
                default:
                    return MZAnalyzerType.Unknown;
            }
        }

        public double GetPrecursorMz(int spectrumNumber, int msnOrder = 2)
        {
            IMSScanRecord scan_record = _reader.GetScanRecord(spectrumNumber - 1);
            return scan_record.MZOfInterest;
        }

        private static Regex ISOLATION_WIDTH_REGEX;
        public double GetIsolationWidth(int spectrumNumber, int msnOrder = 2)
        {
            string acquisition_method;
            using (StreamReader acquisition_method_sr = new StreamReader(Path.Combine(_reader.FileInformation.DataFileName, @"/AcqData/AcqMethod.xml")))
            {
                acquisition_method = acquisition_method_sr.ReadToEnd();
            }
            if (ISOLATION_WIDTH_REGEX == null)
            {
                ISOLATION_WIDTH_REGEX = new Regex(@"\s*(?:&lt;|<)ID(?:&gt;|>)TargetIsolationWidth(?:&lt;|<)/ID(?:&gt;|>)\s*(?:&lt;|<)Value(?:&gt;|>).*\(~([0-9.])+ amu\)(?:&lt;|<)/Value(?:&gt;|>)");
            }
            Match match = ISOLATION_WIDTH_REGEX.Match(acquisition_method);
            return double.Parse(match.Groups[1].Value);
        }

        public DissociationType GetDissociationType(int spectrumNumber, int msnOrder = 2)
        {
            return DissociationType.CID;
        }

        public MzRange GetMzRange(int spectrumNumber)
        {
            IBDASpecData spectrum = _reader.GetSpectrum(spectrumNumber - 1);
            return new MzRange(spectrum.MeasuredMassRange.Start, spectrum.MeasuredMassRange.End);
        }

        public int GetPrecusorCharge(int spectrumNumber, int msnOrder = 2)
        {
            IBDASpecData spectrum = _reader.GetSpectrum(spectrumNumber - 1);
            int precursor_charge;
            spectrum.GetPrecursorCharge(out precursor_charge);
            return (short)precursor_charge;
        }

        public int GetSpectrumNumber(double retentionTime)
        {
            IBDAChromData tic = _reader.GetTIC();
            int index = -1;
            for (int i = 0; i < tic.TotalDataPoints; i++)
            {
                if (index < 0 || Math.Abs(tic.XArray[i] - retentionTime) < Math.Abs(tic.XArray[index] - retentionTime))
                {
                    index = i;
                }
            }
            return index + 1;
        }

        public double GetInjectionTime(int spectrumNumber)
        {
            IMSScanRecord scan_record = _reader.GetScanRecord(spectrumNumber - 1);
            int num_transients = 0;
            double length_transient = double.NaN;
            IBDAActualData[] actuals = _reader.ActualsInformation.GetActualCollection(scan_record.RetentionTime);
            foreach (IBDAActualData actual in actuals)
            {
                if (actual.DisplayName == "Number of Transients")
                {
                    num_transients = int.Parse(actual.DisplayValue);
                }
                else if (actual.DisplayName == "Length of Transients")
                {
                    length_transient = double.Parse(actual.DisplayValue);
                }
            }
            return num_transients * length_transient; // may be off by a factor of two for extended dynamic range mode
        }

        #endregion
    }
}
