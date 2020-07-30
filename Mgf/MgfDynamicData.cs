using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using UsefulProteomicsDatabases;

namespace IO.Mgf
{
    public class MgfDynamicData : DynamicDataConnection
    {
        private StreamReader reader;
        private Dictionary<int, long> scanToByteOffset;
        private static Regex scanNumberParser = new Regex(@"(^|\s)SCANS=(.*?)($|\D)");

        public MgfDynamicData(string filepath) : base(filepath)
        {
            InitiateDynamicConnection();
        }

        public override MsDataScan GetOneBasedScanFromDynamicConnection(int scanNumber, IFilteringParams filterParams = null)
        {
            MsDataScan scan = null;

            if (reader == null)
            {
                throw new MzLibException("Cannot get scan; the dynamic data connection to " + FilePath + " has been closed!");
            }

            if (scanToByteOffset.TryGetValue(scanNumber, out long byteOffset))
            {
                // seek to the byte of the scan
                reader.BaseStream.Position = byteOffset;
                reader.DiscardBufferedData();

                char[] peakSplitter = new char[] { ' ' };
                List<(double mz, double intensity)> peaks = new List<(double, double)>();
                int charge = 2; //default when unknown
                double precursorMz = 0;
                double rtInMinutes = double.NaN; //default when unknown

                // read the scan data
                while (reader.Peek() > 0)
                {
                    string line = reader.ReadLine();
                    string[] sArray = line.Split('=');

                    if (char.IsDigit(line[0]) && sArray.Length == 1)
                    {
                        string[] split = line.Split(peakSplitter, StringSplitOptions.RemoveEmptyEntries);
                        (double mz, double intensity) peak = (double.Parse(split[0], CultureInfo.InvariantCulture), double.Parse(split[1], CultureInfo.InvariantCulture));
                        peaks.Add(peak);
                    }
                    else if (line.StartsWith("PEPMASS="))
                    {
                        sArray = sArray[1].Split(' ');
                        precursorMz = double.Parse(sArray[0]);
                    }
                    else if (line.StartsWith("CHARGE="))
                    {
                        string entry = sArray[1];
                        charge = int.Parse(entry.Substring(0, entry.Length - 1));
                        if (entry[entry.Length - 1].Equals("-"))
                        {
                            charge *= -1;
                        }
                    }
                    else if (line.StartsWith("RTINSECONDS="))
                    {
                        rtInMinutes = double.Parse(sArray[sArray.Length - 1]) / 60.0;
                    }
                    else if (line.StartsWith("END IONS"))
                    {
                        break;
                    }
                }

                double[] mzs = peaks.Select(p => p.mz).ToArray();
                double[] intensities = peaks.Select(p => p.intensity).ToArray();

                Array.Sort(mzs, intensities);
                MzRange scanRange = new MzRange(mzs[0], mzs[mzs.Length - 1]);

                // peak filtering
                if (filterParams != null && intensities.Length > 0 && filterParams.ApplyTrimmingToMsMs)
                {
                    MsDataFile.WindowModeHelper(ref intensities, ref mzs, filterParams, scanRange.Minimum, scanRange.Maximum);
                }

                MzSpectrum spectrum = new MzSpectrum(mzs, intensities, false);

                scan = new MsDataScan(spectrum, scanNumber, 2, true, charge > 0 ? Polarity.Positive : Polarity.Negative,
                    rtInMinutes, scanRange, null, MZAnalyzerType.Unknown,
                    intensities.Sum(), 0, null, null, precursorMz, charge, null, precursorMz, null,
                    DissociationType.Unknown, null, precursorMz);
            }

            return scan;
        }

        public override void CloseDynamicConnection()
        {
            if (reader != null)
            {
                reader.Dispose();
            }
        }

        protected override void InitiateDynamicConnection()
        {
            if (!File.Exists(FilePath))
            {
                throw new FileNotFoundException();
            }

            if (Path.GetExtension(FilePath).ToUpper() != ".MGF")
            {
                throw new InvalidDataException();
            }

            Loaders.LoadElements();
            reader = new StreamReader(FilePath);

            BuildIndex();
        }

        /// <summary>
        /// Mgf files are not indexed, so the index will need to be built.
        /// This method builds the index.
        /// </summary>
        private void BuildIndex()
        {
            scanToByteOffset = new Dictionary<int, long>();
            int oneBasedScanNumber = 0;
            long currentPositionByteOffset = 0;
            long oneBasedScanByteOffset = 0;
            bool scanHasAScanNumber = false;

            while (reader.Peek() > 0)
            {
                currentPositionByteOffset = TextFileReading.GetByteOffsetAtCurrentPosition(reader);

                string line = reader.ReadLine();

                if (line.StartsWith("BEGIN IONS", StringComparison.InvariantCultureIgnoreCase))
                {
                    oneBasedScanByteOffset = currentPositionByteOffset;
                    scanHasAScanNumber = false;
                }
                else if (line.StartsWith("SCANS=", StringComparison.InvariantCultureIgnoreCase))
                {
                    scanHasAScanNumber = true;

                    Match result = scanNumberParser.Match(line);
                    var scanString = result.Groups[2].Value;
                    oneBasedScanNumber = int.Parse(scanString);
                }
                else if (line.StartsWith("END IONS", StringComparison.InvariantCultureIgnoreCase))
                {
                    if (!scanHasAScanNumber)
                    {
                        oneBasedScanNumber++;
                    }

                    if (scanToByteOffset.ContainsKey(oneBasedScanNumber))
                    {
                        throw new MzLibException("Scan number " + oneBasedScanNumber.ToString() +
                            " appeared multiple times in " + FilePath + ", which is not allowed because we assume all scan numbers are unique.");
                    }

                    scanToByteOffset.Add(oneBasedScanNumber, oneBasedScanByteOffset);
                }
            }
        }
    }
}
