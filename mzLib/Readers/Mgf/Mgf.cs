using MassSpectrometry;
using MzLibUtil;
using System.Globalization;
using System.Reflection.Metadata.Ecma335;
using System.Text.RegularExpressions;

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

                            if (scan is not null)
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
            if (CheckIfScansLoaded() && scanNumber <= IndexedScans.Length)
                return GetOneBasedScan(scanNumber);

            lock (DynamicReadingLock)
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

        private static MsDataScan? GetNextMsDataOneBasedScanFromConnection(StreamReader sr, HashSet<int> scanNumbersAlreadyObserved, 
            IFilteringParams filterParams = null, int? alreadyKnownScanNumber = null)
        {
            List<double> mzs = new List<double>();
            List<double> intensities = new List<double>();
            int charge = 2; //default when unknown
            double precursorMz = 0;
            double? precursorIntensity = null; //default when unknown
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

                // A header whose value is missing entirely -- "CHARGE" with no "=" -- is skipped rather
                // than read, the same as any other unrecognised line. Indexing sArray[1] regardless is
                // what made a stray header take down the whole file.
                bool hasValue = sArray.Length > 1;

                if (char.IsDigit(line[0]) && sArray.Length == 1)
                {
                    TryParsePeakLine(line, mzs, intensities);
                }
                else if (line.StartsWith("PEPMASS") && hasValue)
                {
                    sArray = sArray[1].Split(' ');
                    precursorMz = Convert.ToDouble(sArray[0], CultureInfo.InvariantCulture);
                    if (sArray.Length > 1)
                        precursorIntensity = Convert.ToDouble(sArray[1], CultureInfo.InvariantCulture);
                }
                else if (line.StartsWith("CHARGE") && hasValue && sArray[1].Length > 0)
                {
                    // The sign suffix is optional: "2+", "2-" and "2" are all valid charge states.
                    // Strip it only when present -- stripping unconditionally drops a digit, so "12"
                    // read as 1. Polarity is derived from the sign of this charge.
                    string entry = sArray[1];
                    char sign = entry[entry.Length - 1];
                    bool signed = sign == '+' || sign == '-';

                    charge = Convert.ToInt32(signed ? entry.Substring(0, entry.Length - 1) : entry);
                    if (sign == '-')
                    {
                        charge *= -1;
                    }
                }
                else if (line.StartsWith("SCANS") && hasValue)
                {
                    scanNumber = Convert.ToInt32(sArray[1]);
                }
                else if (line.StartsWith("RTINSECONDS") && hasValue)
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
            if (mzArray.IsNullOrEmpty() || intensityArray.IsNullOrEmpty())
                return null;

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
                precursorIntensity, precursorMz, null, DissociationType.Unknown,
                null, precursorMz);
        }

        /// <summary>
        /// Adds a peak to <paramref name="mzs"/> and <paramref name="intensities"/>, and reports whether
        /// the line was a well-formed peak at all.
        ///
        /// "Starts with a digit" is the only test that gets a line here, which is a weak classifier: a
        /// stray line carrying one number, or a number followed by something that is not one, reaches
        /// this method too. Those are skipped like any other unrecognised line instead of indexing past
        /// the end of the split, which surfaced as an IndexOutOfRangeException out of the reader and
        /// took down the entire file.
        /// </summary>
        private static bool TryParsePeakLine(string line, List<double> mzs, List<double> intensities)
        {
            string[] sArray = line.Split(new Char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);

            if (sArray.Length < 2
                || !double.TryParse(sArray[0], NumberStyles.Any, CultureInfo.InvariantCulture, out double mz)
                || !double.TryParse(sArray[1], NumberStyles.Any, CultureInfo.InvariantCulture, out double intensity))
            {
                return false;
            }

            mzs.Add(mz);
            intensities.Add(intensity);
            return true;
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
