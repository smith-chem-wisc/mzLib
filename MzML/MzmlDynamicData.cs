using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Xml;
using UsefulProteomicsDatabases;

namespace IO.MzML
{
    public class MzmlDynamicData
    {
        public readonly string FilePath;

        private Dictionary<int, long> ScanNumberToByteOffset;
        private StreamReader reader;

        //private XmlSerializer serializer;
        private static Regex nativeIdScanNumberParser = new Regex(@"(^|\s)scan=(.*?)($|\D)");

        public MzmlDynamicData(string filepath)
        {
            //XmlRootAttribute xRoot = new XmlRootAttribute();
            //xRoot.ElementName = "spectrum";
            //xRoot.IsNullable = false;
            //serializer = new XmlSerializer(typeof(IO.MzML.Generated.SpectrumType), xRoot);

            this.FilePath = filepath;
            InitiateDynamicConnection();
        }

        /// <summary>
        /// Closes the dynamic connection with the .mzML file.
        /// </summary>
        public void CloseDynamicConnection()
        {
            reader.Dispose();
        }

        /// <summary>
        /// Gets the scan with the specified one-based scan number.
        /// </summary>
        public MsDataScan GetOneBasedScanFromDynamicConnection(int oneBasedScanNumber, IFilteringParams filterParams = null)
        {
            MsDataScan scan = null;

            if (ScanNumberToByteOffset.TryGetValue(oneBasedScanNumber, out long byteOffset))
            {
                // seek to the byte of the scan
                reader.BaseStream.Position = byteOffset;
                reader.DiscardBufferedData();

                // DO NOT USE THIS METHOD! it does not seek reliably
                //stream.BaseStream.Seek(byteOffset, SeekOrigin.Begin);

                // read the scan
                using (XmlReader xmlReader = XmlReader.Create(reader))
                {
                    string nativeId = null;
                    while (xmlReader.Read())
                    {
                        // this skips whitespace
                        string upperName = xmlReader.Name.ToUpper();
                        if (upperName == "SPECTRUM" && xmlReader.IsStartElement())
                        {
                            nativeId = xmlReader["id"];
                            break;
                        }
                    }

                    // deserializing the scan's data doesn't work well. the spectrum type is deserialized
                    // but sub-elements aren't. this is probably because we're trying to deserialize only
                    // a part of the XML file... deserialization would probably be cleaner code than
                    // using the current implementation but I couldn't get it to work
                    //var deserializedSpectrum = (IO.MzML.Generated.SpectrumType)serializer.Deserialize(xmlReader.ReadSubtree());

                    MzSpectrum spectrum = null;
                    int? msOrder = 0;
                    bool? isCentroid = false;
                    Polarity polarity = Polarity.Unknown;
                    double retentionTime = 0;
                    MzRange range = null;
                    string scanFilter = null;
                    MZAnalyzerType mzAnalyzerType = MZAnalyzerType.Unknown;
                    double tic = 0;
                    double? injTime = null;
                    double[,] noiseData = null; // TODO: read this
                    double? selectedIonMz = null;
                    int? selectedCharge = null;
                    double? selectedIonIntensity = null;
                    double? isolationMz = null; // TODO: should this be refined? or taken from the scan header?
                    double? isolationWidth = null;
                    DissociationType? dissociationType = null;
                    int? oneBasedPrecursorScanNumber = null;
                    double? selectedIonMonoisotopicGuessMz = null; // TODO: not sure what this is

                    double scanLowerLimit = 0;
                    double scanUpperLimit = 0;
                    double isolationWindowLowerOffset = 0;
                    double isolationWindowUpperOffset = 0;

                    bool compressed = false;
                    bool readingMzs = false;
                    bool readingIntensities = false;
                    bool is32bit = true;
                    double[] mzs = null;
                    double[] intensities = null;

                    while (xmlReader.Read())
                    {
                        switch (xmlReader.Name.ToUpper())
                        {
                            // controlled vocabulary parameter
                            case "CVPARAM":
                                switch (xmlReader["accession"])
                                {
                                    // MS order
                                    case "MS:1000511":
                                        msOrder = int.Parse(xmlReader["value"]);
                                        break;

                                    // centroid mode
                                    case "MS:1000127":
                                        isCentroid = true;
                                        break;

                                    // profile mode
                                    case "MS:1000128":
                                        isCentroid = false;
                                        throw new Exception("Reading profile mode mzmls not supported");
                                        break;

                                    // positive scan mode
                                    case "MS:1000130":
                                        polarity = Polarity.Positive;
                                        break;

                                    // negative scan mode
                                    case "MS:1000129":
                                        polarity = Polarity.Negative;
                                        break;

                                    // total ion current
                                    case "MS:1000285":
                                        tic = double.Parse(xmlReader["value"]);
                                        break;

                                    // retention time
                                    case "MS:1000016":
                                        retentionTime = double.Parse(xmlReader["value"]);
                                        break;

                                    // filter string
                                    case "MS:1000512":
                                        scanFilter = xmlReader["value"];
                                        break;

                                    // ion injection time
                                    case "MS:1000927":
                                        injTime = double.Parse(xmlReader["value"]);
                                        break;

                                    // scan lower limit
                                    case "MS:1000501":
                                        scanLowerLimit = double.Parse(xmlReader["value"]);
                                        break;

                                    // scan upper limit
                                    case "MS:1000500":
                                        scanUpperLimit = double.Parse(xmlReader["value"]);
                                        break;

                                    // isolation window lower offset
                                    case "MS:1000828":
                                        isolationWindowLowerOffset = double.Parse(xmlReader["value"]);
                                        break;

                                    // isolation window upper offset
                                    case "MS:1000829":
                                        isolationWindowUpperOffset = double.Parse(xmlReader["value"]);
                                        break;

                                    // isolated m/z
                                    case "MS:1000827":
                                        isolationMz = double.Parse(xmlReader["value"]);
                                        break;

                                    // selected ion m/z
                                    case "MS:1000744":
                                        selectedIonMz = double.Parse(xmlReader["value"]);
                                        break;

                                    // selected charge state
                                    case "MS:1000041":
                                        selectedCharge = int.Parse(xmlReader["value"]);
                                        break;

                                    // activation types
                                    case "MS:1000133":
                                        dissociationType = DissociationType.CID;
                                        break;

                                    case "MS:1001880":
                                        dissociationType = DissociationType.ISCID;
                                        break;

                                    case "MS:1000422":
                                        dissociationType = DissociationType.HCD;
                                        break;

                                    case "MS:1000598":
                                        dissociationType = DissociationType.ETD;
                                        break;

                                    case "MS:1000435":
                                        dissociationType = DissociationType.IRMPD;
                                        break;

                                    case "MS:1000599":
                                        dissociationType = DissociationType.PQD;
                                        break;

                                    // mass analyzer types
                                    case "MS:1000081":
                                        mzAnalyzerType = MZAnalyzerType.Quadrupole;
                                        break;

                                    case "MS:1000291":
                                        mzAnalyzerType = MZAnalyzerType.IonTrap2D;
                                        break;

                                    case "MS:1000082":
                                        mzAnalyzerType = MZAnalyzerType.IonTrap3D;
                                        break;

                                    case "MS:1000484":
                                        mzAnalyzerType = MZAnalyzerType.Orbitrap;
                                        break;

                                    case "MS:1000084":
                                        mzAnalyzerType = MZAnalyzerType.TOF;
                                        break;

                                    case "MS:1000079":
                                        mzAnalyzerType = MZAnalyzerType.FTICR;
                                        break;

                                    case "MS:1000080":
                                        mzAnalyzerType = MZAnalyzerType.Sector;
                                        break;

                                    case "MS:1000523":
                                        is32bit = false;
                                        break;

                                    case "MS:1000574":
                                        compressed = true;
                                        break;

                                    case "MS:1000514":
                                        readingMzs = true;
                                        break;

                                    case "MS:1000515":
                                        readingIntensities = true;
                                        break;
                                }
                                break;

                            // binary data array (e.g., m/z or intensity array)
                            case "BINARY":
                                if (!readingMzs && !readingIntensities)
                                {
                                    break;
                                }

                                while (string.IsNullOrWhiteSpace(xmlReader.Value))
                                {
                                    xmlReader.Read();
                                }

                                string binaryString = xmlReader.Value;

                                byte[] binaryData = null;

                                if (!is32bit)
                                {
                                    binaryData = Convert.FromBase64String(binaryString);
                                }
                                else
                                {
                                    // todo: check. not sure if this is right
                                    binaryData = Encoding.UTF8.GetBytes(binaryString);
                                }

                                double[] data = Mzml.ConvertBase64ToDoubles(binaryData, compressed, is32bit);

                                if (readingMzs)
                                {
                                    mzs = data;
                                    readingMzs = false;
                                }
                                else if (readingIntensities)
                                {
                                    intensities = data;
                                    readingIntensities = false;
                                }

                                break;

                            case "PRECURSOR":
                                if (xmlReader.IsStartElement())
                                {
                                    string precursorScanInfo = xmlReader["spectrumRef"];

                                    Match result = nativeIdScanNumberParser.Match(precursorScanInfo);
                                    var scanString = result.Groups[2].Value;
                                    oneBasedPrecursorScanNumber = int.Parse(scanString);
                                }
                                break;

                            // done reading spectrum
                            case "SPECTRUM":
                                if (!xmlReader.IsStartElement())
                                {
                                    if (msOrder > 1)
                                    {
                                        isolationWidth = isolationWindowUpperOffset + isolationWindowLowerOffset;

                                        // look for precursor ion info in the precursor scan
                                        if (oneBasedPrecursorScanNumber != null && selectedIonMz.HasValue)
                                        {
                                            MsDataScan precursorScan = GetOneBasedScanFromDynamicConnection(oneBasedPrecursorScanNumber.Value);
                                            int? bestIndex = precursorScan.MassSpectrum.GetClosestPeakIndex(selectedIonMz.Value);

                                            if (bestIndex.HasValue &&
                                                Math.Abs(precursorScan.MassSpectrum.XArray[bestIndex.Value] - selectedIonMz.Value) < 0.1)
                                            {
                                                selectedIonIntensity = precursorScan.MassSpectrum.YArray[bestIndex.Value];
                                                selectedIonMz = precursorScan.MassSpectrum.XArray[bestIndex.Value];
                                            }
                                        }
                                    }

                                    if (!msOrder.HasValue || !isCentroid.HasValue)
                                    {
                                        throw new MzLibException("Could not determine the MS order or centroid/profile status");
                                    }

                                    // peak filtering
                                    if (filterParams != null && intensities.Length > 0 &&
                                        ((filterParams.ApplyTrimmingToMs1 && msOrder.Value == 1) || (filterParams.ApplyTrimmingToMsMs && msOrder.Value > 1)))
                                    {
                                        MsDataFile.WindowModeHelper(ref intensities, ref mzs, filterParams, scanLowerLimit, scanUpperLimit);
                                    }

                                    Array.Sort(mzs, intensities);

                                    range = new MzRange(scanLowerLimit, scanUpperLimit);
                                    spectrum = new MzSpectrum(mzs, intensities, false);

                                    scan = new MsDataScan(spectrum, oneBasedScanNumber, msOrder.Value, isCentroid.Value, polarity,
                                        retentionTime, range, scanFilter, mzAnalyzerType, tic, injTime, noiseData,
                                        nativeId, selectedIonMz, selectedCharge, selectedIonIntensity, isolationMz, isolationWidth,
                                        dissociationType, oneBasedPrecursorScanNumber, selectedIonMonoisotopicGuessMz);

                                    return scan;
                                }
                                else
                                {
                                    throw new Exception("Spectrum data is malformed");
                                }
                        }
                    }
                }
            }

            return scan;
        }

        /// <summary>
        /// Initiates the dynamic connection with the .mzML file.
        /// </summary>
        private void InitiateDynamicConnection()
        {
            if (!File.Exists(FilePath))
            {
                throw new FileNotFoundException();
            }

            if (Path.GetExtension(FilePath).ToUpper() != ".MZML")
            {
                throw new InvalidDataException();
            }

            Loaders.LoadElements();
            reader = new StreamReader(FilePath);

            ScanNumberToByteOffset = new Dictionary<int, long>();
            FindOrCreateIndex(FilePath);
        }

        /// <summary>
        /// Finds the index in the .mzML file. If the index doesn't exist or can't be found,
        /// then an index is created by the method BuildIndexFromUnindexedMzml().
        /// </summary>
        private void FindOrCreateIndex(string filepath)
        {
            // look for the index in the mzML file
            LookForIndex();

            // the index could not be found. we will need to build it
            if (!ScanNumberToByteOffset.Any())
            {
                CreateIndexFromUnindexedMzml();
            }
        }

        /// <summary>
        /// Looks for the index in the .mzML file.
        /// </summary>
        private void LookForIndex()
        {
            long? indexByteOffset = null;

            // check the top of the file for the index
            while (reader.Peek() > 0)
            {
                var line = reader.ReadLine();

                string trimmedline = line.Trim();
                if (trimmedline.StartsWith("<indexListOffset>", StringComparison.InvariantCultureIgnoreCase))
                {
                    trimmedline = trimmedline.Replace("<indexListOffset>", "");
                    trimmedline = trimmedline.Replace("</indexListOffset>", "");

                    indexByteOffset = long.Parse(trimmedline);
                    break;
                }
                else if (trimmedline.StartsWith("<mzML", StringComparison.InvariantCultureIgnoreCase))
                {
                    break;
                }
            }

            if (!indexByteOffset.HasValue)
            {
                // check the bottom of the file for the index
                // this is super annoying... we need to read the file backwards starting from the end 
                // and then parse the xml...
                ReverseLineReader rlr = new ReverseLineReader(FilePath);

                foreach (string line in rlr)
                {
                    string trimmedline = line.Trim();

                    if (trimmedline.StartsWith("<indexListOffset>", StringComparison.InvariantCultureIgnoreCase))
                    {
                        trimmedline = trimmedline.Replace("<indexListOffset>", "");
                        trimmedline = trimmedline.Replace("</indexListOffset>", "");

                        indexByteOffset = long.Parse(trimmedline);
                        break;
                    }
                    else if (trimmedline.StartsWith("</mzML", StringComparison.InvariantCultureIgnoreCase))
                    {
                        break;
                    }
                }
            }

            if (indexByteOffset.HasValue)
            {
                ReadIndex(indexByteOffset.Value);
            }
        }

        /// <summary>
        /// Reads the index, once its position in the file is known.
        /// </summary>
        private void ReadIndex(long byteOffsetOfIndex)
        {
            // seek to the byte of the scan
            reader.BaseStream.Position = byteOffsetOfIndex;
            reader.DiscardBufferedData();

            bool readingScanIndex = false;

            using (XmlReader xmlReader = XmlReader.Create(reader))
            {
                while (xmlReader.Read())
                {
                    if (xmlReader.IsStartElement())
                    {
                        string upperName = xmlReader.Name.ToUpper();

                        // found the scan index
                        if (upperName == "INDEX" && xmlReader["name"] != null && xmlReader["name"].ToUpper() == "SPECTRUM")
                        {
                            readingScanIndex = true;
                        }

                        // read element of the scan index
                        else if (readingScanIndex && upperName == "OFFSET")
                        {
                            var spectrumInfo = xmlReader["idRef"];
                            long byteOffset = 0;

                            Match result = nativeIdScanNumberParser.Match(spectrumInfo);
                            var scanString = result.Groups[2].Value;

                            if (xmlReader.Read())
                            {
                                var textNode = xmlReader.Value.Trim();

                                // TODO: throw some kind of exception here if it doesn't parse
                                byteOffset = long.Parse(textNode);
                            }
                            else
                            {
                                // TODO: throw exception
                            }

                            if (int.TryParse(scanString, out int scanNumber))
                            {
                                ScanNumberToByteOffset.Add(scanNumber, byteOffset);
                            }
                            else
                            {
                                //TODO: throw an exception or something.. this index is messed up
                            }
                        }
                    }
                    else if (readingScanIndex && xmlReader.NodeType == XmlNodeType.EndElement)
                    {
                        string upperName = xmlReader.Name.ToUpper();

                        // reached the end of the scan index
                        if (upperName == "INDEX")
                        {
                            readingScanIndex = false;
                            break;
                        }
                    }
                }
            }
        }

        /// <summary>
        /// This creates a scan number-to-byte index for .mzML files that have no index.
        /// This means that a dynamic connection can be created with unindexed .mzML files,
        /// it just takes longer because we have to read the entire file one line at a time 
        /// to create the index via code.
        /// 
        /// This does NOT write an index to the .mzML file. This is intentional. This .mzML reader
        /// does not modify the .mzML data file at all. It just creates an index in memory 
        /// via code if one is not provided in the .mzML file.
        /// </summary>
        private void CreateIndexFromUnindexedMzml()
        {
            while (reader.Peek() > 0)
            {
                // this byte offset might be a little different than what it technically 
                // should be because of white space but it will work out ok
                long byteOffset = GetByteOffsetAtCurrentPosition();
                var line = reader.ReadLine();

                line = line.Trim();

                if (line.StartsWith("<spectrum ", StringComparison.InvariantCultureIgnoreCase))
                {
                    Match result = nativeIdScanNumberParser.Match(line);
                    string scanString = result.Groups[2].Value;

                    if (int.TryParse(scanString, out int scanNumber))
                    {
                        var tryLine = reader.ReadLine();

                        ScanNumberToByteOffset.Add(scanNumber, byteOffset);
                    }
                    else
                    {
                        throw new Exception("Could not resolve this information into a scan number: " + line);
                    }
                }
            }
        }

        /// <summary>
        /// The StreamReader object does not have a reliable method to get a byte position in
        /// a file. This method calculates it.
        /// 
        /// Code from: https://stackoverflow.com/questions/5404267/streamreader-and-seeking
        /// </summary>
        private long GetByteOffsetAtCurrentPosition()
        {
            System.Reflection.BindingFlags flags = System.Reflection.BindingFlags.DeclaredOnly | System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance | System.Reflection.BindingFlags.GetField;

            // The current buffer of decoded characters
            char[] charBuffer = (char[])reader.GetType().InvokeMember("charBuffer", flags, null, reader, null);

            // The index of the next char to be read from charBuffer
            int charPos = (int)reader.GetType().InvokeMember("charPos", flags, null, reader, null);

            // The number of decoded chars presently used in charBuffer
            int charLen = (int)reader.GetType().InvokeMember("charLen", flags, null, reader, null);

            // The current buffer of read bytes (byteBuffer.Length = 1024; this is critical).
            byte[] byteBuffer = (byte[])reader.GetType().InvokeMember("byteBuffer", flags, null, reader, null);

            // The number of bytes read while advancing reader.BaseStream.Position to (re)fill charBuffer
            int byteLen = (int)reader.GetType().InvokeMember("byteLen", flags, null, reader, null);

            // The number of bytes the remaining chars use in the original encoding.
            int numBytesLeft = reader.CurrentEncoding.GetByteCount(charBuffer, charPos, charLen - charPos);

            // For variable-byte encodings, deal with partial chars at the end of the buffer
            int numFragments = 0;
            if (byteLen > 0 && !reader.CurrentEncoding.IsSingleByte)
            {
                if (reader.CurrentEncoding.CodePage == 65001) // UTF-8
                {
                    byte byteCountMask = 0;
                    while ((byteBuffer[byteLen - numFragments - 1] >> 6) == 2) // if the byte is "10xx xxxx", it's a continuation-byte
                        byteCountMask |= (byte)(1 << ++numFragments); // count bytes & build the "complete char" mask
                    if ((byteBuffer[byteLen - numFragments - 1] >> 6) == 3) // if the byte is "11xx xxxx", it starts a multi-byte char.
                        byteCountMask |= (byte)(1 << ++numFragments); // count bytes & build the "complete char" mask
                                                                      // see if we found as many bytes as the leading-byte says to expect
                    if (numFragments > 1 && ((byteBuffer[byteLen - numFragments] >> 7 - numFragments) == byteCountMask))
                        numFragments = 0; // no partial-char in the byte-buffer to account for
                }
                else if (reader.CurrentEncoding.CodePage == 1200) // UTF-16LE
                {
                    if (byteBuffer[byteLen - 1] >= 0xd8) // high-surrogate
                        numFragments = 2; // account for the partial character
                }
                else if (reader.CurrentEncoding.CodePage == 1201) // UTF-16BE
                {
                    if (byteBuffer[byteLen - 2] >= 0xd8) // high-surrogate
                        numFragments = 2; // account for the partial character
                }
            }
            return reader.BaseStream.Position - numBytesLeft - numFragments;
        }
    }

    /// <summary>
    /// Takes an encoding (defaulting to UTF-8) and a function which produces a seekable stream
    /// (or a filename for convenience) and yields lines from the end of the stream backwards.
    /// Only single byte encodings, and UTF-8 and Unicode, are supported. The stream
    /// returned by the function must be seekable.
    /// </summary>
    public sealed class ReverseLineReader : IEnumerable<string>
    {
        /// <summary>
        /// Buffer size to use by default. Classes with internal access can specify
        /// a different buffer size - this is useful for testing.
        /// </summary>
        private const int DefaultBufferSize = 4096;

        /// <summary>
        /// Means of creating a Stream to read from.
        /// </summary>
        private readonly Func<Stream> streamSource;

        /// <summary>
        /// Encoding to use when converting bytes to text
        /// </summary>
        private readonly Encoding encoding;

        /// <summary>
        /// Size of buffer (in bytes) to read each time we read from the
        /// stream. This must be at least as big as the maximum number of
        /// bytes for a single character.
        /// </summary>
        private readonly int bufferSize;

        /// <summary>
        /// Function which, when given a position within a file and a byte, states whether
        /// or not the byte represents the start of a character.
        /// </summary>
        private Func<long, byte, bool> characterStartDetector;

        /// <summary>
        /// Creates a LineReader from a stream source. The delegate is only
        /// called when the enumerator is fetched. UTF-8 is used to decode
        /// the stream into text.
        /// </summary>
        /// <param name="streamSource">Data source</param>
        public ReverseLineReader(Func<Stream> streamSource)
            : this(streamSource, Encoding.UTF8)
        {
        }

        /// <summary>
        /// Creates a LineReader from a filename. The file is only opened
        /// (or even checked for existence) when the enumerator is fetched.
        /// UTF8 is used to decode the file into text.
        /// </summary>
        /// <param name="filename">File to read from</param>
        public ReverseLineReader(string filename)
            : this(filename, Encoding.UTF8)
        {
        }

        /// <summary>
        /// Creates a LineReader from a filename. The file is only opened
        /// (or even checked for existence) when the enumerator is fetched.
        /// </summary>
        /// <param name="filename">File to read from</param>
        /// <param name="encoding">Encoding to use to decode the file into text</param>
        public ReverseLineReader(string filename, Encoding encoding)
            : this(() => File.OpenRead(filename), encoding)
        {
        }

        /// <summary>
        /// Creates a LineReader from a stream source. The delegate is only
        /// called when the enumerator is fetched.
        /// </summary>
        /// <param name="streamSource">Data source</param>
        /// <param name="encoding">Encoding to use to decode the stream into text</param>
        public ReverseLineReader(Func<Stream> streamSource, Encoding encoding)
            : this(streamSource, encoding, DefaultBufferSize)
        {
        }

        internal ReverseLineReader(Func<Stream> streamSource, Encoding encoding, int bufferSize)
        {
            this.streamSource = streamSource;
            this.encoding = encoding;
            this.bufferSize = bufferSize;
            if (encoding.IsSingleByte)
            {
                // For a single byte encoding, every byte is the start (and end) of a character
                characterStartDetector = (pos, data) => true;
            }
            else if (encoding is UnicodeEncoding)
            {
                // For UTF-16, even-numbered positions are the start of a character.
                // TODO: This assumes no surrogate pairs. More work required
                // to handle that.
                characterStartDetector = (pos, data) => (pos & 1) == 0;
            }
            else if (encoding is UTF8Encoding)
            {
                // For UTF-8, bytes with the top bit clear or the second bit set are the start of a character
                // See http://www.cl.cam.ac.uk/~mgk25/unicode.html
                characterStartDetector = (pos, data) => (data & 0x80) == 0 || (data & 0x40) != 0;
            }
            else
            {
                throw new ArgumentException("Only single byte, UTF-8 and Unicode encodings are permitted");
            }
        }

        /// <summary>
        /// Returns the enumerator reading strings backwards. If this method discovers that
        /// the returned stream is either unreadable or unseekable, a NotSupportedException is thrown.
        /// </summary>
        public IEnumerator<string> GetEnumerator()
        {
            Stream stream = streamSource();
            if (!stream.CanSeek)
            {
                stream.Dispose();
                throw new NotSupportedException("Unable to seek within stream");
            }
            if (!stream.CanRead)
            {
                stream.Dispose();
                throw new NotSupportedException("Unable to read within stream");
            }
            return GetEnumeratorImpl(stream);
        }

        private IEnumerator<string> GetEnumeratorImpl(Stream stream)
        {
            try
            {
                long position = stream.Length;

                if (encoding is UnicodeEncoding && (position & 1) != 0)
                {
                    throw new InvalidDataException("UTF-16 encoding provided, but stream has odd length.");
                }

                // Allow up to two bytes for data from the start of the previous
                // read which didn't quite make it as full characters
                byte[] buffer = new byte[bufferSize + 2];
                char[] charBuffer = new char[encoding.GetMaxCharCount(buffer.Length)];
                int leftOverData = 0;
                String previousEnd = null;
                // TextReader doesn't return an empty string if there's line break at the end
                // of the data. Therefore we don't return an empty string if it's our *first*
                // return.
                bool firstYield = true;

                // A line-feed at the start of the previous buffer means we need to swallow
                // the carriage-return at the end of this buffer - hence this needs declaring
                // way up here!
                bool swallowCarriageReturn = false;

                while (position > 0)
                {
                    int bytesToRead = Math.Min(position > int.MaxValue ? bufferSize : (int)position, bufferSize);

                    position -= bytesToRead;
                    stream.Position = position;
                    StreamUtil.ReadExactly(stream, buffer, bytesToRead);
                    // If we haven't read a full buffer, but we had bytes left
                    // over from before, copy them to the end of the buffer
                    if (leftOverData > 0 && bytesToRead != bufferSize)
                    {
                        // Buffer.BlockCopy doesn't document its behaviour with respect
                        // to overlapping data: we *might* just have read 7 bytes instead of
                        // 8, and have two bytes to copy...
                        Array.Copy(buffer, bufferSize, buffer, bytesToRead, leftOverData);
                    }
                    // We've now *effectively* read this much data.
                    bytesToRead += leftOverData;

                    int firstCharPosition = 0;
                    while (!characterStartDetector(position + firstCharPosition, buffer[firstCharPosition]))
                    {
                        firstCharPosition++;
                        // Bad UTF-8 sequences could trigger this. For UTF-8 we should always
                        // see a valid character start in every 3 bytes, and if this is the start of the file
                        // so we've done a short read, we should have the character start
                        // somewhere in the usable buffer.
                        if (firstCharPosition == 3 || firstCharPosition == bytesToRead)
                        {
                            throw new InvalidDataException("Invalid UTF-8 data");
                        }
                    }
                    leftOverData = firstCharPosition;

                    int charsRead = encoding.GetChars(buffer, firstCharPosition, bytesToRead - firstCharPosition, charBuffer, 0);
                    int endExclusive = charsRead;

                    for (int i = charsRead - 1; i >= 0; i--)
                    {
                        char lookingAt = charBuffer[i];
                        if (swallowCarriageReturn)
                        {
                            swallowCarriageReturn = false;
                            if (lookingAt == '\r')
                            {
                                endExclusive--;
                                continue;
                            }
                        }
                        // Anything non-line-breaking, just keep looking backwards
                        if (lookingAt != '\n' && lookingAt != '\r')
                        {
                            continue;
                        }
                        // End of CRLF? Swallow the preceding CR
                        if (lookingAt == '\n')
                        {
                            swallowCarriageReturn = true;
                        }
                        int start = i + 1;
                        string bufferContents = new string(charBuffer, start, endExclusive - start);
                        endExclusive = i;
                        string stringToYield = previousEnd == null ? bufferContents : bufferContents + previousEnd;
                        if (!firstYield || stringToYield.Length != 0)
                        {
                            yield return stringToYield;
                        }
                        firstYield = false;
                        previousEnd = null;
                    }

                    previousEnd = endExclusive == 0 ? null : (new string(charBuffer, 0, endExclusive) + previousEnd);

                    // If we didn't decode the start of the array, put it at the end for next time
                    if (leftOverData != 0)
                    {
                        Buffer.BlockCopy(buffer, 0, buffer, bufferSize, leftOverData);
                    }
                }
                if (leftOverData != 0)
                {
                    // At the start of the final buffer, we had the end of another character.
                    throw new InvalidDataException("Invalid UTF-8 data at start of stream");
                }
                if (firstYield && string.IsNullOrEmpty(previousEnd))
                {
                    yield break;
                }
                yield return previousEnd ?? "";
            }
            finally
            {
                stream.Dispose();
            }
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }
    }

    public static class StreamUtil
    {
        public static void ReadExactly(Stream input, byte[] buffer, int bytesToRead)
        {
            int index = 0;
            while (index < bytesToRead)
            {
                int read = input.Read(buffer, index, bytesToRead - index);
                if (read == 0)
                {
                    throw new EndOfStreamException
                        (String.Format("End of stream reached with {0} byte{1} left to read.",
                                       bytesToRead - index,
                                       bytesToRead - index == 1 ? "s" : ""));
                }
                index += read;
            }
        }

    }
}