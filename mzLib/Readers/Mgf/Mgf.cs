using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using UsefulProteomicsDatabases;

// old namespace to ensure backwards compatibility
namespace IO.Mgf
{
    public class Mgf : Readers.Mgf
    {
        public Mgf(string filePath) : base(filePath) { }
    }
}

namespace Readers
{
    public class Mgf : MsDataFile
    {

        protected MsDataScan[] IndexedScans { get; set; }
        public Mgf(string filePath) : base(filePath)
        {
            
        }

        public override MsDataFile LoadAllStaticData(FilteringParams filterParams = null, int maxThreads = 1)
        {
            if (!File.Exists(FilePath))
            {
                throw new FileNotFoundException();
            }

            Loaders.LoadElements();

            List<MsDataScan> scans = new List<MsDataScan>();
            HashSet<int> checkForDuplicateScans = new HashSet<int>();

            using (FileStream fs = new FileStream(FilePath, FileMode.Open, FileAccess.Read, FileShare.Read))
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

            SourceFile = GetSourceFile();

            // ensures that if a scan (OneBasedScanNumber) does not exist,
            // the final scans array will contain a null value  
            // this unique case is due to the nature of loading MGF files
            var orderedScans = scans.OrderBy(x => x.OneBasedScanNumber).ToArray();
            var indexedScans = new MsDataScan[orderedScans[^1].OneBasedScanNumber];
            foreach (var scan in orderedScans)
                indexedScans[scan.OneBasedScanNumber - 1] = scan;

            IndexedScans = indexedScans;
            Scans = orderedScans;
            return this;
        }

        public override SourceFile GetSourceFile()
        {
            return new SourceFile("no nativeID format", "mgf format", null, null, null);
        }

        public override MsDataScan GetOneBasedScan(int scanNumber)
        {
            return IndexedScans[scanNumber - 1];
        }

        public override MsDataScan GetOneBasedScanFromDynamicConnection(int scanNumber, IFilteringParams filterParams = null)
        {
            if (_streamReader == null)
            {
                throw new MzLibException("Cannot get scan; the dynamic data connection to " + FilePath + " has been closed!");
            }

            if (_scanByteOffset.TryGetValue(scanNumber, out long byteOffset))
            {
                // seek to the byte of the scan
                _streamReader.BaseStream.Position = byteOffset;
                _streamReader.DiscardBufferedData();

                return Mgf.GetNextMsDataOneBasedScanFromConnection(_streamReader, new HashSet<int>(), filterParams, scanNumber);
            }
            else
            {
                throw new MzLibException("The specified scan number: " + scanNumber + " does not exist in " + FilePath);
            }
        }

        public override void CloseDynamicConnection()
        {
            if (_streamReader != null)
            {
                _streamReader.Dispose();
            }
        }

        public override void InitiateDynamicConnection()
        {
            if (!File.Exists(FilePath))
            {
                throw new FileNotFoundException();
            }
            Loaders.LoadElements();
            _streamReader = new StreamReader(FilePath);

            BuildIndex();
        }

        /// <summary>
        /// This method ensures backwards compatibility with previous mzLib implementations
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="filteringParams"></param>
        /// <param name="maxThreads"></param>
        /// <returns></returns>
        public static MsDataFile LoadAllStaticData(string filePath, FilteringParams filteringParams = null,
            int maxThreads = 1) => MsDataFileReader.GetDataFile(filePath).LoadAllStaticData(filteringParams, maxThreads);

        private static MsDataScan GetNextMsDataOneBasedScanFromConnection(StreamReader sr, HashSet<int> scanNumbersAlreadyObserved, 
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

            //Remove Zero Intensity Peaks
            double zeroEquivalentIntensity = 0.01;
            int zeroIntensityCount = intensityArray.Count(i => i < zeroEquivalentIntensity);
            int intensityValueCount = intensityArray.Count();
            if (zeroIntensityCount > 0 && zeroIntensityCount < intensityValueCount)
            {
                Array.Sort(intensityArray, mzArray);
                double[] nonZeroIntensities = new double[intensityValueCount - zeroIntensityCount];
                double[] nonZeroMzs = new double[intensityValueCount - zeroIntensityCount];
                intensityArray = intensityArray.SubArray(zeroIntensityCount, intensityValueCount - zeroIntensityCount);
                mzArray = mzArray.SubArray(zeroIntensityCount, intensityValueCount - zeroIntensityCount);
                Array.Sort(mzArray, intensityArray);
            }


            MzRange scanRange = new MzRange(mzArray[0], mzArray[mzArray.Length - 1]);

            // peak filtering
            if (filterParams != null && intensityArray.Length > 0 && filterParams.ApplyTrimmingToMsMs)
            {
                WindowModeHelper.Run(ref intensityArray, ref mzArray, 
                    filterParams, scanRange.Minimum, scanRange.Maximum);
            }

            MzSpectrum spectrum = new MzSpectrum(mzArray, intensityArray, false);

            if (scanNumber == 0)
            {
                scanNumber = oldScanNumber + 1;
            }

            scanNumbersAlreadyObserved.Add(scanNumber);

            return new MsDataScan(spectrum, scanNumber, 2, true,
                charge > 0 ? Polarity.Positive : Polarity.Negative,
                rtInMinutes, scanRange, null, MZAnalyzerType.Unknown,
                intensities.Sum(), 0, null, null, precursorMz, charge, 
                null, precursorMz, null,  DissociationType.Unknown, 
                null, precursorMz);
        }

        private static void ParsePeakLine(string line, List<double> mzs, List<double> intensities)
        {
            string[] sArray = line.Split(new Char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);

            mzs.Add(Convert.ToDouble(sArray[0], CultureInfo.InvariantCulture));
            intensities.Add(Convert.ToDouble(sArray[1], CultureInfo.InvariantCulture));
        }


        private StreamReader _streamReader;
        private Dictionary<int, long> _scanByteOffset; 
        private static Regex _scanNumberparser = new Regex(@"(^|\s)SCANS=(.*?)($|\D)");

        private void BuildIndex()
        {
            _scanByteOffset = new Dictionary<int, long>();
            int oneBasedScanNumber = 0;
            long currentPositionByteOffset = 0;
            long oneBasedScanByteOffset = 0;
            bool scanHasAScanNumber = false;

            while (_streamReader.Peek() > 0)
            {
                currentPositionByteOffset = TextFileReading.GetByteOffsetAtCurrentPosition(_streamReader);

                string line = _streamReader.ReadLine();

                if (line.StartsWith("BEGIN IONS", StringComparison.InvariantCultureIgnoreCase))
                {
                    oneBasedScanByteOffset = currentPositionByteOffset;
                    scanHasAScanNumber = false;
                }
                else if (line.StartsWith("SCANS=", StringComparison.InvariantCultureIgnoreCase))
                {
                    scanHasAScanNumber = true;

                    Match result = _scanNumberparser.Match(line);
                    var scanString = result.Groups[2].Value;
                    oneBasedScanNumber = int.Parse(scanString);
                }
                else if (line.StartsWith("END IONS", StringComparison.InvariantCultureIgnoreCase))
                {
                    if (!scanHasAScanNumber)
                    {
                        oneBasedScanNumber++;
                    }

                    if (_scanByteOffset.ContainsKey(oneBasedScanNumber))
                    {
                        throw new MzLibException("Scan number " + oneBasedScanNumber.ToString() +
                                                 " appeared multiple times in " + FilePath + ", which is not allowed because we assume all scan numbers are unique.");
                    }

                    _scanByteOffset.Add(oneBasedScanNumber, oneBasedScanByteOffset);
                }
            }
        }
    }
}
