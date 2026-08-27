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
            int? msLevel = null; //from the MSLEVEL line when present
            DissociationType? dissociationType = null; //from ACTIVATIONMETHOD, an mzLib extension header
            int? precursorScanNumber = null;           //from PRECURSORSCAN, an mzLib extension header
            double? isolationWidth = null;             //from ISOLATIONWIDTH, an mzLib extension header
            double? isolationMz = null;                //from ISOLATIONMZ, an mzLib extension header
            double? totalIonCurrent = null;            //from TIC, an mzLib extension header
            bool sawPrecursorMz = false;
            bool sawCharge = false;

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
                    if (sArray.Length > 1)
                        precursorIntensity = Convert.ToDouble(sArray[1], CultureInfo.InvariantCulture);
                    sawPrecursorMz = true;
                }
                else if (line.StartsWith("MSLEVEL"))
                {
                    msLevel = Convert.ToInt32(sArray[1], CultureInfo.InvariantCulture);
                }
                else if (line.StartsWith("CHARGE"))
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
                    sawCharge = true;
                }
                else if (line.StartsWith("SCANS"))
                {
                    scanNumber = Convert.ToInt32(sArray[1]);
                }
                else if (line.StartsWith("RTINSECONDS"))
                {
                    rtInMinutes = Convert.ToDouble(sArray[sArray.Length - 1], CultureInfo.InvariantCulture) / 60.0;
                }
                // mzLib extension headers -- see MgfMethods.WriteScan for why these exist. TryParse
                // throughout: these come from files we did not necessarily write, and one unparseable
                // value must not take the whole file down.
                else if (line.StartsWith("TIC") && sArray.Length > 1)
                {
                    if (double.TryParse(sArray[1], NumberStyles.Float | NumberStyles.AllowThousands, CultureInfo.InvariantCulture, out double parsedTic))
                    {
                        totalIonCurrent = parsedTic;
                    }
                }
                else if (line.StartsWith("ACTIVATIONMETHOD") && sArray.Length > 1)
                {
                    if (Enum.TryParse(sArray[1].Trim(), ignoreCase: true, out DissociationType parsedDissociation))
                    {
                        dissociationType = parsedDissociation;
                    }
                }
                else if (line.StartsWith("PRECURSORSCAN") && sArray.Length > 1)
                {
                    if (int.TryParse(sArray[1], NumberStyles.Integer, CultureInfo.InvariantCulture, out int parsedPrecursorScan))
                    {
                        precursorScanNumber = parsedPrecursorScan;
                    }
                }
                else if (line.StartsWith("ISOLATIONMZ") && sArray.Length > 1)
                {
                    if (double.TryParse(sArray[1], NumberStyles.Float | NumberStyles.AllowThousands, CultureInfo.InvariantCulture, out double parsedIsolationMz))
                    {
                        isolationMz = parsedIsolationMz;
                    }
                }
                else if (line.StartsWith("ISOLATIONWIDTH") && sArray.Length > 1)
                {
                    if (double.TryParse(sArray[1], NumberStyles.Float | NumberStyles.AllowThousands, CultureInfo.InvariantCulture, out double parsedIsolationWidth))
                    {
                        isolationWidth = parsedIsolationWidth;
                    }
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

            // MSLEVEL when the writer supplied one; otherwise a block with a precursor is MS2 and a
            // block without one is MS1. Files predating MSLEVEL all carry PEPMASS, so they read as before.
            int msnOrder = msLevel ?? (sawPrecursorMz ? 2 : 1);

            // peak filtering
            if (filterParams != null && intensityArray.Length > 0 && ShouldTrim(filterParams, msnOrder))
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

            // MGF has no polarity field, so it is only inferrable from the sign of CHARGE. A block with
            // no CHARGE line reads as positive on both the static and dynamic paths -- inheriting it from
            // a neighbouring block would make the two disagree, since a random-access read cannot see
            // what preceded it.
            Polarity polarity = charge > 0 ? Polarity.Positive : Polarity.Negative;

            if (msnOrder == 1)
            {
                return new MsDataScan(spectrum, scanNumber, msnOrder, true, polarity,
                    rtInMinutes, scanRange, null, MZAnalyzerType.Unknown,
                    totalIonCurrent ?? intensities.Sum(), 0, null, null);
            }

            return new MsDataScan(spectrum, scanNumber, msnOrder, true, polarity,
                rtInMinutes, scanRange, null, MZAnalyzerType.Unknown,
                totalIonCurrent ?? intensities.Sum(), 0, null, null, precursorMz, charge,
                precursorIntensity, isolationMz ?? precursorMz, isolationWidth,
                dissociationType ?? DissociationType.Unknown,
                precursorScanNumber, precursorMz);
        }

        // Peak trimming must honor the per-MS-level FilteringParams flags, exactly as the mzML/Thermo/Bruker
        // readers do. Trimming an MS1 precursor scan strips its isotope envelopes and precursor
        // deconvolution then finds nothing.
        private static bool ShouldTrim(IFilteringParams filterParams, int msnOrder)
        {
            return (filterParams.ApplyTrimmingToMs1 && msnOrder == 1)
                   || (filterParams.ApplyTrimmingToMsMs && msnOrder == 2)
                   || (filterParams.ApplyTrimmingToMsN && msnOrder > 2);
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
