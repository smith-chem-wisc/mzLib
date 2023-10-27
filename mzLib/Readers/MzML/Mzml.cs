// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work Copyright 2016, 2017 Stefan Solntsev
//
// This file (Mzml.cs) is part of MassSpecFiles.
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
using MzLibUtil;
using System.Collections.Concurrent;
using System.Globalization;
using System.IO.Compression;
using System.Security.Cryptography;
using System.Text;
using System.Text.RegularExpressions;
using System.Xml;
using UsefulProteomicsDatabases;

// old namespace to ensure backwards compatibility
namespace IO.MzML
{
    public class Mzml : Readers.Mzml
    {
        public Mzml(string filePath) : base(filePath) { }
    }
}

namespace Readers
{
    public class Mzml : MsDataFile
    {
        #region Constants

        private const string _zlibCompression = "MS:1000574";
        private const string _64bit = "MS:1000523";
        private const string _32bit = "MS:1000521";
        private const string _filterString = "MS:1000512";
        private const string _centroidSpectrum = "MS:1000127";
        private const string _profileSpectrum = "MS:1000128";
        private const string _peakIntensity = "MS:1000042";
        private const string _totalIonCurrent = "MS:1000285";
        private const string _scanWindowLowerLimit = "MS:1000501";
        private const string _scanWindowUpperLimit = "MS:1000500";
        private const string _msnOrderAccession = "MS:1000511";
        private const string _precursorCharge = "MS:1000041";
        private const string _selectedIonMz = "MS:1000744";
        private const string _isolationWindowTargetMZ = "MS:1000827";
        private const string _isolationWindowLowerOffset = "MS:1000828";
        private const string _isolationWindowUpperOffset = "MS:1000829";
        private const string _oneBasedScanNumber = "MS:1000797";
        private const string _retentionTime = "MS:1000016";
        private const string _ionInjectionTime = "MS:1000927";
        private const string _mzArray = "MS:1000514";
        private const string _intensityArray = "MS:1000515";
        private static readonly Regex MZAnalyzerTypeRegex = new Regex(@"^[a-zA-Z]*", RegexOptions.Compiled);

        public static readonly Dictionary<string, Polarity> PolarityDictionary = new Dictionary<string, Polarity>
        {
            {"MS:1000129", Polarity.Negative},
            {"MS:1000130", Polarity.Positive}
        };

        public static readonly Dictionary<string, MZAnalyzerType> AnalyzerDictionary = new Dictionary<string, MZAnalyzerType>
        {
            { "MS:1000443", MZAnalyzerType.Unknown},
            { "MS:1000081", MZAnalyzerType.Quadrupole},
            { "MS:1000291", MZAnalyzerType.IonTrap2D},
            { "MS:1000082", MZAnalyzerType.IonTrap3D},
            { "MS:1000484", MZAnalyzerType.Orbitrap},
            { "MS:1000084", MZAnalyzerType.TOF},
            { "MS:1000079", MZAnalyzerType.FTICR},
            { "MS:1000080", MZAnalyzerType.Sector}
        };

        public static readonly Dictionary<string, DissociationType> DissociationDictionary = new Dictionary<string, DissociationType>
        {
            { "MS:1000133", DissociationType.CID},
            { "MS:1000134", DissociationType.PD},
            { "MS:1000135", DissociationType.PSD},
            { "MS:1000136", DissociationType.SID},
            { "MS:1000242", DissociationType.BIRD},
            { "MS:1000250", DissociationType.ECD},
            { "MS:1000262", DissociationType.IRMPD},
            { "MS:1000282", DissociationType.SORI},
            { "MS:1000435", DissociationType.MPD},
            { "MS:1000598", DissociationType.ETD},
            { "MS:1000599", DissociationType.PQD},
            { "MS:1001880", DissociationType.ISCID},
            { "MS:1000422", DissociationType.HCD},
            { "MS:1002631", DissociationType.EThcD},
            { "MS:1000044", DissociationType.Unknown},
        };

        public static readonly Dictionary<string, DissociationType> DissociationTypeNames = new Dictionary<string, DissociationType>
        {
            { "collision-induced dissociation", DissociationType.CID},
            { "plasma desorption", DissociationType.PD},
            { "post-source decay", DissociationType.PSD},
            { "surface-induced dissociation", DissociationType.SID},
            { "blackbody infrared radiative dissociation", DissociationType.BIRD},
            { "electron capture dissociation", DissociationType.ECD},
            { "infrared multiphoton dissociation", DissociationType.IRMPD},
            { "sustained off-resonance irradiation", DissociationType.SORI},
            { "photodissociation", DissociationType.MPD},
            { "electron transfer dissociation", DissociationType.ETD},
            { "pulsed q dissociation", DissociationType.PQD},
            { "in-source collision-induced dissociation", DissociationType.ISCID},
            { "higher energy beam-type collision-induced dissociation", DissociationType.HCD},
            { "Electron-Transfer/Higher-Energy Collision Dissociation (EThcD)", DissociationType.EThcD},
            { "unknown dissociation type", DissociationType.Unknown},
        };

        #endregion


        Generated.mzMLType _mzMLConnection;

        public Mzml(string filePath) : base(filePath)
        {
            InitializeConnection();
        }

        private void InitializeConnection()
        {
            try
            {
                using (FileStream fs = new FileStream(FilePath, FileMode.Open, FileAccess.Read, FileShare.Read))
                {
                    var _indexedmzMLConnection = (Generated.indexedmzML)MzmlMethods.indexedSerializer.Deserialize(fs);
                    _mzMLConnection = _indexedmzMLConnection.mzML;
                }
            }
            catch
            {
                using (FileStream fs = new FileStream(FilePath, FileMode.Open, FileAccess.Read, FileShare.Read))
                    _mzMLConnection = (Generated.mzMLType)MzmlMethods.mzmlSerializer.Deserialize(fs);
            }
        }

        public override MsDataFile LoadAllStaticData(FilteringParams filterParams = null, int maxThreads = 1)
        {
            if (!File.Exists(FilePath))
            {
                throw new FileNotFoundException();
            }

            Loaders.LoadElements();

            SourceFile = GetSourceFile();

            var numSpecta = _mzMLConnection.run.spectrumList.spectrum.Length;
            MsDataScan[] scans = new MsDataScan[numSpecta];

            Parallel.ForEach(Partitioner.Create(0, numSpecta), new ParallelOptions
            { MaxDegreeOfParallelism = maxThreads }, fff =>
        {
            for (int i = fff.Item1; i < fff.Item2; i++)
            {
                scans[i] = GetMsDataOneBasedScanFromConnection(_mzMLConnection, i + 1, filterParams);
            }
        });

            scans = scans.Where(s => s.MassSpectrum != null).ToArray();

            //Mzml sometimes have scan numbers specified, but usually not.
            //In the event that they do, the iterator above unintentionally assigned them to an incorrect index.
            //Check to make sure that the scans are in order and that there are no duplicate scans
            HashSet<int> checkForDuplicateScans = new HashSet<int>();
            bool ordered = true;
            int previousScanNumber = -1;
            foreach (MsDataScan scan in scans)
            {
                //check if no duplicates
                if (!checkForDuplicateScans.Add(scan.OneBasedScanNumber)) //returns false if the scan already exists
                {
                    throw new MzLibException("Scan number " + scan.OneBasedScanNumber.ToString() + " appeared multiple times in " + FilePath);
                }
                //check if scans are in order
                if (previousScanNumber > scan.OneBasedScanNumber)
                {
                    ordered = false;
                }
                previousScanNumber = scan.OneBasedScanNumber;
            }

            if (!ordered) //reassign indexes if not ordered
            {
                MsDataScan[] indexedScans = new MsDataScan[checkForDuplicateScans.Max()];
                foreach (MsDataScan scan in scans)
                {
                    indexedScans[scan.OneBasedScanNumber - 1] = scan;
                }
                scans = indexedScans;
            }

            //make reference pervious ms1 scan
            // we weren't able to get the precursor scan number, so we'll have to guess;
            // loop back to find precursor scan
            // (assumed to be the first scan before this scan with an MS order of this scan's MS order - 1)
            // e.g., if this is an MS2 scan, find the first MS1 scan before this and assume that's the precursor scan
            for (int i = 0; i < scans.Length; i++)
            {
                if (scans[i].MsnOrder > 1 && scans[i].OneBasedPrecursorScanNumber == null)
                {
                    for (int j = i; j >= 0; j--)
                    {
                        if (scans[i].MsnOrder - scans[j].MsnOrder == 1)
                        {
                            scans[i].SetOneBasedPrecursorScanNumber(scans[j].OneBasedScanNumber);
                            break;
                        }
                    }
                }
            }

            Scans = scans;
            return this;
        }

        public override SourceFile GetSourceFile()
        {
            SourceFile sourceFile;
            if (_mzMLConnection.fileDescription.sourceFileList != null && _mzMLConnection.fileDescription.sourceFileList.sourceFile
                != null && _mzMLConnection.fileDescription.sourceFileList.sourceFile[0] != null
                && _mzMLConnection.fileDescription.sourceFileList.sourceFile[0].cvParam != null)
            {
                var simpler = _mzMLConnection.fileDescription.sourceFileList.sourceFile[0];
                string nativeIdFormat = null;
                string fileFormat = null;
                string checkSum = null;
                string checkSumType = null;
                foreach (var cv in simpler.cvParam)
                {
                    if (cv.accession.Equals(@"MS:1000563"))
                    {
                        fileFormat = "Thermo RAW format";
                    }
                    if (cv.accession.Equals(@"MS:1000584"))
                    {
                        fileFormat = "mzML format";
                    }

                    if (cv.accession.Equals(@"MS:1000768"))
                    {
                        nativeIdFormat = "Thermo nativeID format";
                    }
                    if (cv.accession.Equals(@"MS:1000776"))
                    {
                        nativeIdFormat = "scan number only nativeID format";
                    }
                    if (cv.accession.Equals(@"MS:1000824"))
                    {
                        nativeIdFormat = "no nativeID format";
                    }

                    if (cv.accession.Equals(@"MS:1000568"))
                    {
                        checkSum = cv.value;
                        checkSumType = "MD5";
                    }
                    if (cv.accession.Equals(@"MS:1000569"))
                    {
                        checkSum = cv.value;
                        checkSumType = "SHA-1";
                    }
                }

                sourceFile = new SourceFile(
                    nativeIdFormat,
                    fileFormat,
                    checkSum,
                    checkSumType,
                    new Uri(simpler.location),
                    simpler.id,
                    simpler.name);
            }
            else
            {
                string sendCheckSum;
                using (FileStream stream = File.OpenRead(FilePath))
                {
                    SHA1 sha = SHA1.Create();
                    byte[] checksum = sha.ComputeHash(stream);
                    sendCheckSum = BitConverter.ToString(checksum)
                        .Replace("-", string.Empty);
                }
                sourceFile = new SourceFile(
                    @"no nativeID format",
                    @"mzML format",
                    sendCheckSum,
                    @"SHA-1",
                    Path.GetFullPath(FilePath),
                    Path.GetFileNameWithoutExtension(FilePath));
            }
            return sourceFile;
        }

        /// <summary>
        /// Gets the scan with the specified one-based scan number.
        /// </summary>
        public override MsDataScan GetOneBasedScanFromDynamicConnection(int oneBasedScanNumber, IFilteringParams filterParams = null)
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
                    double retentionTime = double.NaN;
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
                    double? selectedIonMonoisotopicGuessMz = null;

                    double scanLowerLimit = double.NaN;
                    double scanUpperLimit = double.NaN;
                    double isolationWindowLowerOffset = double.NaN;
                    double isolationWindowUpperOffset = double.NaN;

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
                                string cvParamAccession = xmlReader["accession"];

                                if (Mzml.DissociationDictionary.ContainsKey(cvParamAccession))
                                {
                                    dissociationType = Mzml.DissociationDictionary[cvParamAccession];
                                    break;
                                }

                                if (Mzml.PolarityDictionary.ContainsKey(cvParamAccession))
                                {
                                    polarity = Mzml.PolarityDictionary[cvParamAccession];
                                    break;
                                }

                                switch (cvParamAccession)
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
                                        throw new MzLibException("Reading profile mode mzmls not supported");
                                    //break;

                                    // total ion current
                                    case "MS:1000285":
                                        tic = double.Parse(xmlReader["value"]);
                                        break;

                                    // retention time
                                    case "MS:1000016":
                                        retentionTime = double.Parse(xmlReader["value"]);

                                        // determine units (e.g., minutes or seconds)
                                        string units = xmlReader["unitAccession"];

                                        if (units != null && units == "UO:0000010")
                                        {
                                            // convert from seconds to minutes
                                            retentionTime /= 60;
                                        }
                                        else if (units != null && units == "UO:0000031")
                                        {
                                            // do nothing; the RT is already in minutes
                                        }
                                        else
                                        {
                                            throw new MzLibException("The retention time for scan " + oneBasedScanNumber + " could not be interpreted because there was " +
                                                "no value for units (e.g., minutes or seconds)");
                                        }

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

                                    // selected intensity
                                    case "MS:1000042":
                                        selectedIonIntensity = double.Parse(xmlReader["value"]);
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

                                    case "MS:1000521":
                                        is32bit = true;
                                        break;

                                    case "MS:1000576":
                                        compressed = false;
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

                                byte[] binaryData = Convert.FromBase64String(binaryString);

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
                                    // TODO: note that the precursor scan info may not be available in the .mzML. in this case the precursor
                                    // scan number will incorrectly be null. one fix would be to go backwards through the scans to find
                                    // the precursor scan and then set the scan num here, which would be very time consuming.
                                    string precursorScanInfo = xmlReader["spectrumRef"];

                                    if (precursorScanInfo != null)
                                    {
                                        oneBasedPrecursorScanNumber = NativeIdToScanNumber[precursorScanInfo];
                                    }
                                }
                                break;

                            case "USERPARAM":
                                if (xmlReader.IsStartElement() && xmlReader["name"] != null && xmlReader["name"] == "[mzLib]Monoisotopic M/Z:")
                                {
                                    selectedIonMonoisotopicGuessMz = double.Parse(xmlReader["value"]);
                                }
                                break;

                            // done reading spectrum
                            case "SPECTRUM":
                                if (!xmlReader.IsStartElement())
                                {
                                    if (msOrder > 1)
                                    {
                                        isolationWidth = isolationWindowUpperOffset + isolationWindowLowerOffset;

                                        if (dissociationType == null)
                                        {
                                            dissociationType = DissociationType.Unknown;
                                        }
                                    }

                                    if (!msOrder.HasValue || !isCentroid.HasValue)
                                    {
                                        throw new MzLibException("Could not determine the MS order or centroid/profile status");
                                    }

                                    //Remove Zero Intensity Peaks
                                    double zeroEquivalentIntensity = 0.01;
                                    int zeroIntensityCount = intensities.Count(i => i < zeroEquivalentIntensity);
                                    int intensityValueCount = intensities.Count();
                                    if (zeroIntensityCount > 0 && zeroIntensityCount < intensityValueCount)
                                    {
                                        Array.Sort(intensities, mzs);
                                        double[] nonZeroIntensities = new double[intensityValueCount - zeroIntensityCount];
                                        double[] nonZeroMzs = new double[intensityValueCount - zeroIntensityCount];
                                        intensities = intensities.SubArray(zeroIntensityCount, intensityValueCount - zeroIntensityCount);
                                        mzs = mzs.SubArray(zeroIntensityCount, intensityValueCount - zeroIntensityCount);
                                        Array.Sort(mzs, intensities);
                                    }


                                    // peak filtering
                                    if (filterParams != null && intensities.Length > 0 &&
                                        ((filterParams.ApplyTrimmingToMs1 && msOrder.Value == 1) || (filterParams.ApplyTrimmingToMsMs && msOrder.Value > 1)))
                                    {
                                        WindowModeHelper.Run(ref intensities, ref mzs, filterParams, scanLowerLimit, scanUpperLimit);
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
                                    throw new MzLibException("Spectrum data is malformed");
                                }
                        }
                    }
                }
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

        public override void InitiateDynamicConnection()
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
            NativeIdToScanNumber = new Dictionary<string, int>();

            FindOrCreateIndex();
        }

        private static MsDataScan GetMsDataOneBasedScanFromConnection(Generated.mzMLType _mzMLConnection, int oneBasedIndex, IFilteringParams filterParams)
        {
            // Read in the instrument configuration types from connection (in mzml it's at the start)

            Generated.InstrumentConfigurationType[] configs = new Generated.InstrumentConfigurationType[_mzMLConnection.instrumentConfigurationList.instrumentConfiguration.Length];
            for (int i = 0; i < _mzMLConnection.instrumentConfigurationList.instrumentConfiguration.Length; i++)
            {
                configs[i] = _mzMLConnection.instrumentConfigurationList.instrumentConfiguration[i];
            }

            var defaultInstrumentConfig = _mzMLConnection.run.defaultInstrumentConfigurationRef;
            // May be null!
            var scanSpecificInsturmentConfig = _mzMLConnection.run.spectrumList.spectrum[oneBasedIndex - 1].scanList.scan[0].instrumentConfigurationRef;

            MZAnalyzerType analyzer = default(MZAnalyzerType);
            // use default
            if (scanSpecificInsturmentConfig == null || scanSpecificInsturmentConfig == defaultInstrumentConfig)
            {
                if (configs[0].componentList == null)
                {
                    analyzer = default(MZAnalyzerType);
                }
                else if (AnalyzerDictionary.TryGetValue(configs[0].componentList.analyzer[0].cvParam[0].accession, out MZAnalyzerType returnVal))
                {
                    analyzer = returnVal;
                }
            }
            // use scan-specific
            else
            {
                for (int i = 0; i < _mzMLConnection.instrumentConfigurationList.instrumentConfiguration.Length; i++)
                {
                    if (configs[i].id.Equals(scanSpecificInsturmentConfig))
                    {
                        AnalyzerDictionary.TryGetValue(configs[i].componentList.analyzer[0].cvParam[0].accession, out MZAnalyzerType returnVal);
                        analyzer = returnVal;
                    }
                }
            }

            string nativeId = _mzMLConnection.run.spectrumList.spectrum[oneBasedIndex - 1].id;

            int? msOrder = null;
            bool? isCentroid = null;
            Polarity polarity = Polarity.Unknown;
            double tic = double.NaN;

            foreach (Generated.CVParamType cv in _mzMLConnection.run.spectrumList.spectrum[oneBasedIndex - 1].cvParam)
            {
                if (cv.accession.Equals(_msnOrderAccession))
                {
                    msOrder = int.Parse(cv.value);
                }
                if (cv.accession.Equals(_centroidSpectrum))
                {
                    isCentroid = true;
                }
                if (cv.accession.Equals(_profileSpectrum))
                {
                    throw new MzLibException("Reading profile mode mzmls not supported");
                }
                if (cv.accession.Equals(_totalIonCurrent))
                {
                    tic = double.Parse(cv.value, CultureInfo.InvariantCulture);
                }
                if (polarity.Equals(Polarity.Unknown))
                {
                    PolarityDictionary.TryGetValue(cv.accession, out polarity);
                }
            }

            double rtInMinutes = double.NaN;
            string scanFilter = null;
            double? injectionTime = null;
            int oneBasedScanNumber = oneBasedIndex;
            if (_mzMLConnection.run.spectrumList.spectrum[oneBasedIndex - 1].scanList.scan[0].cvParam != null)
            {
                foreach (Generated.CVParamType cv in _mzMLConnection.run.spectrumList.spectrum[oneBasedIndex - 1].scanList.scan[0].cvParam)
                {
                    if (cv.accession.Equals(_retentionTime))
                    {
                        rtInMinutes = double.Parse(cv.value, CultureInfo.InvariantCulture);
                        if (cv.unitName == "second")
                        {
                            rtInMinutes /= 60;
                        }
                    }
                    if (cv.accession.Equals(_filterString))
                    {
                        scanFilter = cv.value;
                    }
                    if (cv.accession.Equals(_ionInjectionTime))
                    {
                        injectionTime = double.Parse(cv.value, CultureInfo.InvariantCulture);
                    }
                    if (cv.accession.Equals(_oneBasedScanNumber)) //get the real one based spectrum number (if available), the other assumes they are in order. This is present in .mgf->.mzml conversions from MSConvert
                    {
                        oneBasedScanNumber = int.Parse(cv.value);
                    }
                }
            }

            if (!msOrder.HasValue || !isCentroid.HasValue)
                //one instance when this if statment is true (i.e. not false) is when there is no mz/intensity data
                //so, we return the MsDataScan object with a null spectrum
                //scans w/ null spectra are checked later and the scan numbers associated w those scans are returned to the reader.
                return new MsDataScan(
                    null,
                    oneBasedScanNumber,
                    msOrder.Value,
                    false, //have to return a value here b/c it is not nullable
                    polarity,
                    rtInMinutes,
                    null,
                    scanFilter,
                    analyzer,
                    tic,
                    injectionTime,
                    null,
                    nativeId);

            double[] masses = new double[0];
            double[] intensities = new double[0];

            foreach (Generated.BinaryDataArrayType binaryData in _mzMLConnection.run.spectrumList.spectrum[oneBasedIndex - 1].binaryDataArrayList.binaryDataArray)
            {
                bool compressed = false;
                bool mzArray = false;
                bool intensityArray = false;
                bool is32bit = true;
                foreach (Generated.CVParamType cv in binaryData.cvParam)
                {
                    compressed |= cv.accession.Equals(_zlibCompression);
                    is32bit &= !cv.accession.Equals(_64bit);
                    is32bit |= cv.accession.Equals(_32bit);
                    mzArray |= cv.accession.Equals(_mzArray);
                    intensityArray |= cv.accession.Equals(_intensityArray);
                }

                //in the futurem we may see scass w/ no data and there will be a crash here. if that happens, you can retrun an MsDataScan with null as the mzSpectrum
                //the scans with no spectra will be reported to the reader and left out of the scan list.
                double[] data = ConvertBase64ToDoubles(binaryData.binary, compressed, is32bit);
                if (mzArray)
                {
                    masses = data;
                }

                if (intensityArray)
                {
                    intensities = data;
                }
            }

            double high = double.NaN;
            double low = double.NaN;

            var aScanWindowList = _mzMLConnection.run.spectrumList.spectrum[oneBasedIndex - 1].scanList.scan[0].scanWindowList;

            if (aScanWindowList != null)
            {
                foreach (Generated.CVParamType cv in _mzMLConnection.run.spectrumList.spectrum[oneBasedIndex - 1].scanList.scan[0].scanWindowList.scanWindow[0].cvParam)
                {
                    if (cv.accession.Equals(_scanWindowLowerLimit))
                    {
                        low = double.Parse(cv.value, CultureInfo.InvariantCulture);
                    }
                    else if (cv.accession.Equals(_scanWindowUpperLimit))
                    {
                        high = double.Parse(cv.value, CultureInfo.InvariantCulture);
                    }
                }
            }

            //Remove Zero Intensity Peaks
            double zeroEquivalentIntensity = 0.01;
            int zeroIntensityCount = intensities.Count(i => i < zeroEquivalentIntensity);
            int intensityValueCount = intensities.Count();
            if (zeroIntensityCount > 0 && zeroIntensityCount < intensityValueCount)
            {
                Array.Sort(intensities, masses);
                double[] nonZeroIntensities = new double[intensityValueCount - zeroIntensityCount];
                double[] nonZeroMzs = new double[intensityValueCount - zeroIntensityCount];
                intensities = intensities.SubArray(zeroIntensityCount, intensityValueCount - zeroIntensityCount);
                masses = masses.SubArray(zeroIntensityCount, intensityValueCount - zeroIntensityCount);
                Array.Sort(masses, intensities);
            }

            if (filterParams != null && intensities.Length > 0 && ((filterParams.ApplyTrimmingToMs1 && msOrder.Value == 1) || (filterParams.ApplyTrimmingToMsMs && msOrder.Value > 1)))
            {
                WindowModeHelper.Run(ref intensities, ref masses, filterParams, low, high);
            }

            Array.Sort(masses, intensities);
            var mzmlMzSpectrum = new MzSpectrum(masses, intensities, false);

            if (msOrder.Value == 1)
            {
                return new MsDataScan(
                    mzmlMzSpectrum,
                    oneBasedScanNumber,
                    msOrder.Value,
                    isCentroid.Value,
                    polarity,
                    rtInMinutes,
                    new MzRange(low, high),
                    scanFilter,
                    analyzer,
                    tic,
                    injectionTime,
                    null,
                    nativeId);
            }

            double selectedIonMz = double.NaN;
            int? selectedIonCharge = null;
            double? selectedIonIntensity = null;
            foreach (Generated.CVParamType cv in _mzMLConnection.run.spectrumList.spectrum[oneBasedIndex - 1].precursorList.precursor[0].selectedIonList.selectedIon[0].cvParam)
            {
                if (cv.accession.Equals(_selectedIonMz))
                {
                    selectedIonMz = double.Parse(cv.value, CultureInfo.InvariantCulture);
                }
                if (cv.accession.Equals(_precursorCharge))
                {
                    selectedIonCharge = int.Parse(cv.value, CultureInfo.InvariantCulture);
                }
                if (cv.accession.Equals(_peakIntensity))
                {
                    selectedIonIntensity = double.Parse(cv.value, CultureInfo.InvariantCulture);
                }
            }

            double? isolationMz = null;
            double lowIsolation = double.NaN;
            double highIsolation = double.NaN;
            if (_mzMLConnection.run.spectrumList.spectrum[oneBasedIndex - 1].precursorList.precursor[0].isolationWindow != null)
            {
                foreach (Generated.CVParamType cv in _mzMLConnection.run.spectrumList.spectrum[oneBasedIndex - 1].precursorList.precursor[0].isolationWindow.cvParam)
                {
                    if (cv.accession.Equals(_isolationWindowTargetMZ))
                    {
                        isolationMz = double.Parse(cv.value, CultureInfo.InvariantCulture);
                    }
                    if (cv.accession.Equals(_isolationWindowLowerOffset))
                    {
                        lowIsolation = double.Parse(cv.value, CultureInfo.InvariantCulture);
                    }
                    if (cv.accession.Equals(_isolationWindowUpperOffset))
                    {
                        highIsolation = double.Parse(cv.value, CultureInfo.InvariantCulture);
                    }
                }
            }

            DissociationType dissociationType = DissociationType.Unknown;

            if (_mzMLConnection.run.spectrumList.spectrum[oneBasedIndex - 1].precursorList.precursor[0].activation.cvParam != null)
            {
                // for EThcD scans, the dissociation type will not be listed as EThcD. it will be 2 different dissociation types
                // in the list, one as ETD and one with HCD. so we need to check for that case and interpret it as EThcD.
                List<DissociationType> scanDissociationTypes = new List<DissociationType>();

                foreach (Generated.CVParamType cv in _mzMLConnection.run.spectrumList.spectrum[oneBasedIndex - 1].precursorList.precursor[0].activation.cvParam)
                {
                    if (DissociationDictionary.TryGetValue(cv.accession, out var scanDissociationType))
                    {
                        scanDissociationTypes.Add(scanDissociationType);
                    }
                }

                if (scanDissociationTypes.Contains(DissociationType.ETD) && scanDissociationTypes.Contains(DissociationType.HCD))
                {
                    dissociationType = DissociationType.EThcD;
                }
                else if (scanDissociationTypes.Any())
                {
                    dissociationType = scanDissociationTypes.First();
                }
                else
                {
                    dissociationType = DissociationType.Unknown;
                }
            }
            double? monoisotopicMz = null;
            if (_mzMLConnection.run.spectrumList.spectrum[oneBasedIndex - 1].scanList.scan[0].userParam != null)
            {
                foreach (var userParam in _mzMLConnection.run.spectrumList.spectrum[oneBasedIndex - 1].scanList.scan[0].userParam)
                {
                    if (userParam.name.EndsWith("Monoisotopic M/Z:"))
                    {
                        monoisotopicMz = double.Parse(userParam.value, CultureInfo.InvariantCulture);
                    }
                }
            }

            int? precursorScanNumber;
            if (_mzMLConnection.run.spectrumList.spectrum[oneBasedIndex - 1].precursorList.precursor[0].spectrumRef == null)
            {
                precursorScanNumber = null;
            }
            else
            {
                precursorScanNumber = GetOneBasedPrecursorScanNumber(_mzMLConnection, oneBasedIndex);
            }

            return new MsDataScan(
                mzmlMzSpectrum,
                oneBasedIndex,
                msOrder.Value,
                isCentroid.Value,
                polarity,
                rtInMinutes,
                new MzRange(low, high),
                scanFilter,
                analyzer,
                tic,
                injectionTime,
                null,
                nativeId,
                selectedIonMz,
                selectedIonCharge,
                selectedIonIntensity,
                isolationMz,
                lowIsolation + highIsolation,
                dissociationType,
                precursorScanNumber,
                monoisotopicMz
                );
        }

        /// <summary>
        /// Converts a 64-based encoded byte array into an double[]
        /// </summary>
        /// <param name="bytes">the 64-bit encoded byte array</param>
        /// <param name="zlibCompressed">Specifies if the byte array is zlib compressed</param>
        /// <returns>a decompressed, de-encoded double[]</returns>
        public static double[] ConvertBase64ToDoubles(byte[] bytes, bool zlibCompressed = false, bool is32bit = true)
        {
            // Add capability of compressed data

            if (zlibCompressed)
            {
                var output = new MemoryStream();
                using (var compressStream = new MemoryStream(bytes))
                {
                    compressStream.ReadByte();
                    compressStream.ReadByte();
                    using (var decompressor = new DeflateStream(compressStream, CompressionMode.Decompress))
                    {
                        decompressor.CopyTo(output);
                        decompressor.Close();
                        output.Position = 0;
                        bytes = output.ToArray();
                    }
                }
            }

            int size = is32bit ? sizeof(float) : sizeof(double);

            int length = bytes.Length / size;
            double[] convertedArray = new double[length];

            for (int i = 0; i < length; i++)
            {
                if (is32bit)
                {
                    convertedArray[i] = BitConverter.ToSingle(bytes, i * size);
                }
                else
                {
                    convertedArray[i] = BitConverter.ToDouble(bytes, i * size);
                }
            }
            return convertedArray;
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

        private static int GetOneBasedPrecursorScanNumber(Generated.mzMLType _mzMLConnection, int oneBasedSpectrumNumber)
        {
            string precursorID = _mzMLConnection.run.spectrumList.spectrum[oneBasedSpectrumNumber - 1].precursorList.precursor[0].spectrumRef;
            do
            {
                oneBasedSpectrumNumber--;
            } while (!precursorID.Equals(_mzMLConnection.run.spectrumList.spectrum[oneBasedSpectrumNumber - 1].id));
            return oneBasedSpectrumNumber;
        }


        private Dictionary<int, long> ScanNumberToByteOffset;
        private Dictionary<string, int> NativeIdToScanNumber;
        private StreamReader reader;
        public static readonly Regex nativeIdScanNumberParser = new Regex(@"(^|\s)scan=(.*?)($|\D)");

        /// <summary>
        /// Finds the index in the .mzML file. If the index doesn't exist or can't be found,
        /// then an index is created by the method CreateIndexFromUnindexedMzml().
        /// </summary>
        private void FindOrCreateIndex()
        {
            // look for the index in the mzML file
            try
            {
                LookForIndex();
            }
            catch (Exception)
            {
                // something went wrong reading the index
                // the index will need to be created instead
                ScanNumberToByteOffset.Clear();
                NativeIdToScanNumber.Clear();
            }

            // the index could not be found or could not be read. we will need to build it
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
            // seek to the position of the index in the file
            reader.BaseStream.Position = byteOffsetOfIndex;
            reader.DiscardBufferedData();

            bool readingScanIndex = false;

            // some nativeID formats don't have a scan number specified. 
            // in this case, use a counter to determine the scan number.
            int scanNumber = 0;

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

                            if (result.Groups[2].Success)
                            {
                                // TODO: throw some kind of exception here if it doesn't parse
                                scanNumber = int.Parse(result.Groups[2].Value);
                            }
                            else
                            {
                                scanNumber++;
                            }

                            NativeIdToScanNumber.Add(spectrumInfo, scanNumber);

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

                            ScanNumberToByteOffset.Add(scanNumber, byteOffset);
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
            reader.BaseStream.Position = 0;
            reader.DiscardBufferedData();

            int scanNumber = 0;

            while (reader.Peek() > 0)
            {
                // this byte offset might be a little different than what it technically 
                // should be because of white space but it will work out ok
                long byteOffset = TextFileReading.GetByteOffsetAtCurrentPosition(reader);
                var line = reader.ReadLine();

                line = line.Trim();

                if (line.StartsWith("<spectrum ", StringComparison.InvariantCultureIgnoreCase))
                {
                    Match result = nativeIdScanNumberParser.Match(line);
                    int ind = line.IndexOf("id=\"");
                    string nativeId;

                    if (ind >= 0)
                    {
                        StringBuilder nativeIdBuilder = new StringBuilder();

                        for (int r = ind + 4; r < line.Length; r++)
                        {
                            if (line[r] == '"')
                            {
                                break;
                            }

                            nativeIdBuilder.Append(line[r]);
                        }

                        nativeId = nativeIdBuilder.ToString();
                    }
                    else
                    {
                        throw new MzLibException("Could not get nativeID from line: " + line);
                    }

                    if (result.Groups[2].Success)
                    {
                        scanNumber = int.Parse(result.Groups[2].Value);
                    }
                    else
                    {
                        scanNumber++;
                    }

                    NativeIdToScanNumber.Add(nativeId, scanNumber);
                    ScanNumberToByteOffset.Add(scanNumber, byteOffset);
                }
            }
        }
    }
}