using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace IO.Mgf
{
    public class Mgf : MsDataFile
    {
        private MsDataScan[] indexedScans;

        private Mgf(MsDataScan[] scans, SourceFile sourceFile) : base(scans, sourceFile)
        {
            indexedScans = new MsDataScan[scans[scans.Length - 1].OneBasedScanNumber];
            foreach (MsDataScan scan in scans)
            {
                indexedScans[scan.OneBasedScanNumber - 1] = scan;
            }
        }

        public static Mgf LoadAllStaticData(string filePath, FilteringParams filterParams = null)
        {
            if (!File.Exists(filePath))
            {
                throw new FileNotFoundException();
            }

            Loaders.LoadElements();

            List<MsDataScan> scans = new List<MsDataScan>();
            HashSet<int> checkForDuplicateScans = new HashSet<int>();

            using (FileStream fs = new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                using (BufferedStream bs = new BufferedStream(fs))
                {
                    using (StreamReader sr = new StreamReader(bs))
                    {
                        while (sr.Peek() > 0)
                        {
                            string line = sr.ReadLine();
                            if (line != "BEGIN IONS")
                            {
                                continue;
                            }

                            var scan = GetNextMsDataOneBasedScanFromConnection(sr, checkForDuplicateScans, filterParams);

                            scans.Add(scan);
                        }
                    }
                }
            }

            SourceFile sourceFile = new SourceFile("no nativeID format", "mgf format", null, null, null);

            return new Mgf(scans.OrderBy(x => x.OneBasedScanNumber).ToArray(), sourceFile);
        }

        public override MsDataScan GetOneBasedScan(int scanNumber)
        {
            return indexedScans[scanNumber - 1];
        }

        public static MsDataScan GetNextMsDataOneBasedScanFromConnection(StreamReader sr, HashSet<int> scanNumbersAlreadyObserved, 
            IFilteringParams filterParams = null, int? alreadyKnownScanNumber = null)
        {
            List<double> mzs = new List<double>();
            List<double> intensities = new List<double>();
            int charge = 2; //default when unknown
            double precursorMz = 0;
            double rtInMinutes = double.NaN; //default when unknown

            int oldScanNumber = scanNumbersAlreadyObserved.Count > 0 ? scanNumbersAlreadyObserved.Max() : 0;
            int scanNumber = alreadyKnownScanNumber.HasValue ? alreadyKnownScanNumber.Value : 0;

            // read the scan data
            while (sr.Peek() > 0)
            {
                string line = sr.ReadLine();
                string[] sArray = line.Split('=');

                if (String.IsNullOrWhiteSpace(line))
                {
                    continue;
                }

                if (char.IsDigit(line[0]) && sArray.Length == 1)
                {
                    ParsePeakLine(line, mzs, intensities);
                }
                else if (line.StartsWith("PEPMASS"))
                {
                    sArray = sArray[1].Split(' ');
                    precursorMz = Convert.ToDouble(sArray[0], CultureInfo.InvariantCulture);
                }
                else if (line.StartsWith("CHARGE"))
                {
                    string entry = sArray[1];
                    charge = Convert.ToInt32(entry.Substring(0, entry.Length - 1));
                    if (entry[entry.Length - 1].Equals("-"))
                    {
                        charge *= -1;
                    }
                }
                else if (line.StartsWith("SCANS"))
                {
                    scanNumber = Convert.ToInt32(sArray[1]);
                }
                else if (line.StartsWith("RTINSECONDS"))
                {
                    rtInMinutes = Convert.ToDouble(sArray[sArray.Length - 1], CultureInfo.InvariantCulture) / 60.0;
                }
                else if (line.StartsWith("END IONS"))
                {
                    break;
                }
            }

            double[] mzArray = mzs.ToArray();
            double[] intensityArray = intensities.ToArray();

            Array.Sort(mzArray, intensityArray);
            MzRange scanRange = new MzRange(mzArray[0], mzArray[mzArray.Length - 1]);

            // peak filtering
            if (filterParams != null && intensityArray.Length > 0 && filterParams.ApplyTrimmingToMsMs)
            {
                MsDataFile.WindowModeHelper(ref intensityArray, ref mzArray, filterParams, scanRange.Minimum, scanRange.Maximum);
            }

            MzSpectrum spectrum = new MzSpectrum(mzArray, intensityArray, false);

            if (scanNumber == 0)
            {
                scanNumber = oldScanNumber + 1;
            }

            scanNumbersAlreadyObserved.Add(scanNumber);

            return new MsDataScan(spectrum, scanNumber, 2, true, charge > 0 ? Polarity.Positive : Polarity.Negative,
                rtInMinutes, scanRange, null, MZAnalyzerType.Unknown,
                intensities.Sum(), 0, null, null, precursorMz, charge, null, precursorMz, null,
                DissociationType.Unknown, null, precursorMz);
        }

        private static void ParsePeakLine(string line, List<double> mzs, List<double> intensities)
        {
            var sArray = line.Split(' ', StringSplitOptions.RemoveEmptyEntries);

            mzs.Add(Convert.ToDouble(sArray[0], CultureInfo.InvariantCulture));
            intensities.Add(Convert.ToDouble(sArray[1], CultureInfo.InvariantCulture));
        }
    }
}
