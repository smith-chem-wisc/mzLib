// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work Copyright 2016 Stefan Solntsev
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

using Ionic.Zlib;
using MassSpectrometry;
using Spectra;
using System;
using System.IO;
using System.Text.RegularExpressions;
using System.Xml.Serialization;

namespace IO.MzML
{
    public class Mzml : MsDataFile<DefaultMzSpectrum>
    {
        private static string _msnOrderAccession = "MS:1000511";
        private static string _precursorCharge = "MS:1000041";
        private static string _precursorMass = "MS:1000744";
        private static string _isolationWindowTargetMZ = "MS:1000827";
        private static string _isolationWindowLowerOffset = "MS:1000828";
        private static string _isolationWindowUpperOffset = "MS:1000829";
        private static string _retentionTime = "MS:1000016";
        private static string _ionInjectionTime = "MS:1000927";
        private static string _mzArray = "MS:1000514";
        private static string _intensityArray = "MS:1000515";
        private const string _CID = "MS:1000133";
        private const string _ISCID = "MS:1001880";
        private const string _HCD = "MS:1000422";
        private const string _ETD = "MS:1000598";
        private const string _MPD = "MS:1000435";
        private const string _ECD = "MS:1000250";
        private const string _PQD = "MS:1000599";
        private const string _DefaultDissociation = "MS:1000044";
        private const string _quadrupole = "MS:1000081";
        private const string _linearIonTrap = "MS:1000291";
        private const string _IonTrap2DAxialEject = "MS:1000078";
        private const string _IonTrap2DRadialEject = "MS:1000083";
        private const string _IonTrap3D = "MS:1000082";
        private const string _orbitrap = "MS:1000484";
        private const string _TOF = "MS:1000084";
        private const string _FTICR = "MS:1000079";
        private const string _magneticSector = "MS:1000080";
        private const string _nozlibCompress = "MS:1000576";
        private const string _zlibCompression = "MS:1000574";
        private const string _64bit = "MS:1000523";
        private const string _32bit = "MS:1000521";
        private const string _negativePolarity = "MS:1000129";
        private const string _positivePolarity = "MS:1000130";
        private const string _filterString = "MS:1000512";
        private const string _centroidSpectrum = "MS:1000127";
        private const string _profileSpectrum = "MS:1000128";
        private const string _peakIntensity = "MS:1000042";
        private const string _totalIonCurrent = "MS:1000285";
        private const string _scanWindowLowerLimit = "MS:1000501";
        private const string _scanWindowUpperLimit = "MS:1000500";


        private static XmlSerializer _indexedSerializer = new XmlSerializer(typeof(Generated.indexedmzML));
        private static XmlSerializer _mzMLSerializer = new XmlSerializer(typeof(Generated.mzMLType));

        private Generated.indexedmzML _indexedmzMLConnection;
        private Generated.mzMLType _mzMLConnection;

        public Mzml(string filePath)
            : base(filePath, true, MsDataFileType.Mzml)
        {
        }
        public override void Open()
        {
            if (_mzMLConnection == null)
            {
                Stream stream = new FileStream(FilePath, FileMode.Open);
                try
                {
                    _indexedmzMLConnection = _indexedSerializer.Deserialize(stream) as Generated.indexedmzML;
                    _mzMLConnection = _indexedmzMLConnection.mzML;
                }
                catch (Exception)
                {
                    try
                    {
                        _mzMLConnection = _mzMLSerializer.Deserialize(stream) as Generated.mzMLType;
                    }
                    catch (Exception)
                    {
                        throw new InvalidDataException("Unable to parse " + FilePath + " as an mzML file!");
                    }
                }
            }
        }

        public static void Write(string filePath, Generated.indexedmzML _indexedmzMLConnection)
        {
            TextWriter writer = new StreamWriter(filePath);
            _indexedSerializer.Serialize(writer, _indexedmzMLConnection);
            writer.Close();
        }

        public bool IsIndexedMzML
        {
            get { return _indexedmzMLConnection != null; }
        }

        private DissociationType GetDissociationType(int spectrumNumber)
        {
            spectrumNumber--;
            foreach (Generated.CVParamType cv in _mzMLConnection.run.spectrumList.spectrum[spectrumNumber].precursorList.precursor[0].activation.cvParam)
            {
                switch (cv.accession)
                {
                    case _CID:
                        return DissociationType.CID;
                    case _ISCID:
                        return DissociationType.ISCID;
                    case _HCD:
                        return DissociationType.HCD;
                    case _ETD:
                        return DissociationType.ETD;
                    case _MPD:
                        return DissociationType.MPD;
                    case _PQD:
                        return DissociationType.PQD;
                    case _DefaultDissociation:
                        return DissociationType.Unknown;
                }
            }
            throw new ArgumentNullException("Could not find dissociation type for spectrum number " + spectrumNumber + 1);
        }

        private int GetMsnOrder(int spectrumNumber)
        {
            spectrumNumber--;
            foreach (Generated.CVParamType cv in _mzMLConnection.run.spectrumList.spectrum[spectrumNumber].cvParam)
            {
                if (cv.accession.Equals(_msnOrderAccession))
                {
                    return int.Parse(cv.value);
                }
            }
            throw new ArgumentNullException("Could not find MSn level for spectrum number " + spectrumNumber + 1);
        }

        private int GetPrecusorCharge(int spectrumNumber)
        {
            // PRECURSOR ARE ONLY IN MS2 SPECTRA!!!
            spectrumNumber--;

            if (_mzMLConnection.run.spectrumList.spectrum[spectrumNumber].precursorList == null)
                throw new ArgumentNullException("Couldn't find precursor charge in spectrum number " + spectrumNumber + 1 + ", possibly an MS1 Spectrum!");

            foreach (Generated.CVParamType cv in _mzMLConnection.run.spectrumList.spectrum[spectrumNumber].precursorList.precursor[0].selectedIonList.selectedIon[0].cvParam)
            {
                if (cv.accession.Equals(_precursorCharge))
                {
                    return short.Parse(cv.value);
                }
            }
            throw new ArgumentNullException("Couldn't find precursor charge in spectrum number " + spectrumNumber + 1);
        }

        private MzRange GetScanWindowMzRange(int spectrumNumber)
        {
            spectrumNumber--;
            double high = double.NaN;
            double low = double.NaN;

            foreach (Generated.CVParamType cv in _mzMLConnection.run.spectrumList.spectrum[spectrumNumber].scanList.scan[0].scanWindowList.scanWindow[0].cvParam)
            {
                if (cv.accession.Equals(_scanWindowLowerLimit))
                {
                    low = double.Parse(cv.value);
                }
                if (cv.accession.Equals(_scanWindowUpperLimit))
                {
                    high = double.Parse(cv.value);
                }
            }
            if (double.IsNaN(low) || double.IsNaN(high))
            {
                throw new ArgumentNullException("Could not find scan window MZ range for " + spectrumNumber + 1);
            }
            return new MzRange(low, high);
        }


        private double GetIsolationWidth(int spectrumNumber)
        {
            spectrumNumber--;
            double low = double.NaN;
            double high = double.NaN;

            foreach (Generated.CVParamType cv in _mzMLConnection.run.spectrumList.spectrum[spectrumNumber].precursorList.precursor[0].isolationWindow.cvParam)
            {
                if (cv.accession.Equals(_isolationWindowLowerOffset))
                {
                    low = double.Parse(cv.value);
                }
                if (cv.accession.Equals(_isolationWindowUpperOffset))
                {
                    high = double.Parse(cv.value);
                }
            }
            if (double.IsNaN(low) || double.IsNaN(high))
            {
                throw new ArgumentNullException("Could not determine isolation width for " + spectrumNumber + 1);
            }
            return high - low;
        }

        private string GetScanFilter(int spectrumNumber)
        {
            spectrumNumber--;

            foreach (Generated.CVParamType cv in _mzMLConnection.run.spectrumList.spectrum[spectrumNumber].scanList.scan[0].cvParam)
            {
                if (cv.accession.Equals(_filterString))
                {
                    return cv.value;
                }
            }
            // Not a problem if null, scan filter is optional!
            return null;
        }

        private static readonly Regex MZAnalyzerTypeRegex = new Regex(@"^[a-zA-Z]*", RegexOptions.Compiled);

        private MZAnalyzerType GetMzAnalyzer(int spectrumNumber)
        {

            string filter = GetScanFilter(spectrumNumber);

            if (filter == null)
                throw new ArgumentNullException("Cannot get analyzer for spectrum number " + spectrumNumber + " because scan filter is not present!");

            string type = MZAnalyzerTypeRegex.Match(filter).Captures[0].Value;

            switch (type)
            {
                case "ITMS":
                    return MZAnalyzerType.IonTrap2D;
                case "TQMS":
                    throw new InvalidDataException("Not sure what TQMS is");
                case "SQMS":
                    throw new InvalidDataException("Not sure what SQMS is");
                case "TOFMS":
                    return MZAnalyzerType.TOF;
                case "FTMS":
                    return MZAnalyzerType.Orbitrap;
                case "Sector":
                    return MZAnalyzerType.Sector;
            }

            // Maybe in the beginning of the file, there is a single analyzer?
            // Gets the first analyzer used.        
            string analyzer = _mzMLConnection.instrumentConfigurationList.instrumentConfiguration[0].cvParam[0].accession;

            switch (analyzer)
            {
                case _quadrupole:
                    return MZAnalyzerType.Quadrupole;
                case _linearIonTrap:
                    return MZAnalyzerType.IonTrap2D;
                case _IonTrap3D:
                    return MZAnalyzerType.IonTrap3D;
                case _orbitrap:
                    return MZAnalyzerType.Orbitrap;
                case _TOF:
                    return MZAnalyzerType.TOF;
                case _FTICR:
                    return MZAnalyzerType.FTICR;
                case _magneticSector:
                    return MZAnalyzerType.Sector;
                default:
                    return MZAnalyzerType.Unknown;
            }
        }

        private Polarity GetPolarity(int spectrumNumber)
        {

            spectrumNumber--;
            foreach (Generated.CVParamType cv in _mzMLConnection.run.spectrumList.spectrum[spectrumNumber].cvParam)
            {
                if (cv.accession.Equals(_negativePolarity))
                {
                    return Polarity.Negative;
                }
                if (cv.accession.Equals(_positivePolarity))
                {
                    return Polarity.Positive;
                }
            }
            //return Polarity.Neutral;
            throw new ArgumentNullException("Could not find polarity for spectrum number " + spectrumNumber + 1);
        }

        private double GetRetentionTime(int spectrumNumber)
        {
            spectrumNumber--;
            if (_mzMLConnection.run.spectrumList.spectrum[spectrumNumber].scanList.scan[0].cvParam == null)
            {
                return double.NaN;
            }
            double rt = -1;
            foreach (Generated.CVParamType cv in _mzMLConnection.run.spectrumList.spectrum[spectrumNumber].scanList.scan[0].cvParam)
            {
                if (cv.accession.Equals(_retentionTime))
                {
                    rt = double.Parse(cv.value);
                }

                if (cv.unitName == "second")
                    rt /= 60;
            }

            if (rt >= 0)
                return rt;

            throw new ArgumentNullException("Could not determine retention time for " + spectrumNumber + 1);
        }

        private double GetInjectionTime(int spectrumNumber)
        {
            spectrumNumber--;
            foreach (Generated.CVParamType cv in _mzMLConnection.run.spectrumList.spectrum[spectrumNumber].scanList.scan[0].cvParam)
            {
                if (cv.accession.Equals(_ionInjectionTime))
                {
                    return double.Parse(cv.value);
                }
            }
            // HACK
            return -1;
        }

        protected override int GetFirstSpectrumNumber()
        {
            return 1;
        }

        protected override int GetLastSpectrumNumber()
        {
            return Convert.ToInt32(_mzMLConnection.run.spectrumList.count);
        }

        public override int GetSpectrumNumber(double retentionTime)
        {
            // TODO need to convert this to a binary search of some sort. Or if the data is indexedMZML see if the indices work better.
            int totalSpectra = Convert.ToInt32(_mzMLConnection.run.spectrumList.count);
            for (int i = 0; i < totalSpectra; i++)
            {
                foreach (Generated.CVParamType cv in _mzMLConnection.run.spectrumList.spectrum[i].scanList.scan[0].cvParam)
                {
                    if (cv.accession.Equals(_retentionTime))
                    {
                        if (double.Parse(cv.value) == retentionTime)
                        {
                            return i + 1;
                        }
                    }
                }
            }
            throw new ArgumentNullException("Could not determine spectrum number");
        }

        public static byte[] ConvertDoublestoBase64(double[] toConvert, bool zlibCompressed)
        {

            var mem = new MemoryStream();
            for (int i = 0; i < toConvert.Length; i++)
            {
                byte[] ok = BitConverter.GetBytes(toConvert[i]);
                mem.Write(ok, 0, ok.Length);
            }
            mem.Position = 0;

            byte[] bytes = mem.ToArray();
            if (zlibCompressed)
                bytes = ZlibStream.CompressBuffer(bytes);

            return bytes;
        }



        /// <summary>
        /// Converts a 64-based encoded byte array into an double[]
        /// </summary>
        /// <param name="bytes">the 64-bit encoded byte array</param>
        /// <param name="zlibCompressed">Specifies if the byte array is zlib compressed</param>
        /// <returns>a decompressed, de-encoded double[]</returns>
        private static double[] ConvertBase64ToDoubles(byte[] bytes, bool zlibCompressed = false, bool is32bit = true)
        {

            // Add capability of compressed data
            if (zlibCompressed)
                bytes = ZlibStream.UncompressBuffer(bytes);

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

        private bool GetIsCentroid(int spectrumNumber)
        {
            spectrumNumber--;
            foreach (Generated.CVParamType cv in _mzMLConnection.run.spectrumList.spectrum[spectrumNumber].cvParam)
            {
                if (cv.accession.Equals(_centroidSpectrum))
                {
                    return true;
                }
                if (cv.accession.Equals(_profileSpectrum))
                {
                    return false;
                }
            }
            throw new ArgumentNullException("Could not determine if spectrum " + spectrumNumber + 1 + " is centroid or profile");
        }

        private string GetSpectrumID(int spectrumNumber)
        {
            spectrumNumber--;
            return _mzMLConnection.run.spectrumList.spectrum[spectrumNumber].id;
        }

        private string GetPrecursorID(int spectrumNumber)
        {
            spectrumNumber--;
            return _mzMLConnection.run.spectrumList.spectrum[spectrumNumber].precursorList.precursor[0].spectrumRef;
        }

        private double GetPrecursorIntensity(int spectrumNumber)
        {
            spectrumNumber--;
            foreach (Generated.CVParamType cv in _mzMLConnection.run.spectrumList.spectrum[spectrumNumber].precursorList.precursor[0].selectedIonList.selectedIon[0].cvParam)
            {
                if (cv.accession.Equals(_peakIntensity))
                {
                    return Convert.ToDouble(cv.value);
                }
            }
            throw new ArgumentNullException("Could not determine precursor intensity of spectrum " + spectrumNumber + 1);
        }

        private double GetPrecursorMz(int spectrumNumber)
        {
            spectrumNumber--;
            foreach (Generated.CVParamType cv in _mzMLConnection.run.spectrumList.spectrum[spectrumNumber].precursorList.precursor[0].selectedIonList.selectedIon[0].cvParam)
            {
                if (cv.accession.Equals(_precursorMass))
                {
                    return double.Parse(cv.value);
                }
            }
            throw new ArgumentNullException("Could not determine precursor monoisotopic mass for " + spectrumNumber + 1);
        }

        private double GetIsolationMz(int spectrumNumber)
        {
            spectrumNumber--;
            foreach (Generated.CVParamType cv in _mzMLConnection.run.spectrumList.spectrum[spectrumNumber].precursorList.precursor[0].isolationWindow.cvParam)
            {
                if (cv.accession.Equals(_isolationWindowTargetMZ))
                {
                    return double.Parse(cv.value);
                }
            }
            throw new ArgumentNullException("Could not determine isolation mz for " + spectrumNumber + 1);

        }

        protected override MsDataScan<DefaultMzSpectrum> GetMsDataScanFromFile(int spectrumNumber)
        {

            spectrumNumber--; // 0-based indexing

            double[] masses = null;
            double[] intensities = null;

            foreach (Generated.BinaryDataArrayType binaryData in _mzMLConnection.run.spectrumList.spectrum[spectrumNumber].binaryDataArrayList.binaryDataArray)
            {
                bool compressed = false;
                bool mzArray = false;
                bool intensityArray = false;
                bool is32bit = true;
                foreach (Generated.CVParamType cv in binaryData.cvParam)
                {
                    if (cv.accession.Equals(_zlibCompression))
                    {
                        compressed = true;
                    }
                    if (cv.accession.Equals(_64bit))
                    {
                        is32bit = false;
                    }
                    if (cv.accession.Equals(_32bit))
                    {
                        is32bit = true;
                    }
                    if (cv.accession.Equals(_mzArray))
                    {
                        mzArray = true;
                    }
                    if (cv.accession.Equals(_intensityArray))
                    {
                        intensityArray = true;
                    }
                }

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

            if (masses == null || intensities == null)
            {
                throw new InvalidDataException("Unable to find spectral data for spectrum number " + spectrumNumber + 1);
            }

            var ok = new DefaultMzSpectrum(masses, intensities, false);

            if (GetMsnOrder(spectrumNumber + 1) == 1)
                return new MsDataScan<DefaultMzSpectrum>(spectrumNumber, ok, GetSpectrumID(spectrumNumber + 1), GetMsnOrder(spectrumNumber + 1), GetIsCentroid(spectrumNumber + 1), GetPolarity(spectrumNumber + 1), GetRetentionTime(spectrumNumber + 1), GetScanWindowMzRange(spectrumNumber + 1), GetScanFilter(spectrumNumber + 1), GetMzAnalyzer(spectrumNumber + 1), GetInjectionTime(spectrumNumber + 1), GetTotalIonCurrent(spectrumNumber + 1));
            else
                return new MsDataScan<DefaultMzSpectrum>(spectrumNumber, ok, GetSpectrumID(spectrumNumber + 1), GetMsnOrder(spectrumNumber + 1), GetIsCentroid(spectrumNumber + 1), GetPolarity(spectrumNumber + 1), GetRetentionTime(spectrumNumber + 1), GetScanWindowMzRange(spectrumNumber + 1), GetScanFilter(spectrumNumber + 1), GetMzAnalyzer(spectrumNumber + 1), GetInjectionTime(spectrumNumber + 1), GetTotalIonCurrent(spectrumNumber + 1), GetPrecursorID(spectrumNumber + 1), GetPrecursorMz(spectrumNumber + 1), GetPrecusorCharge(spectrumNumber + 1), GetPrecursorIntensity(spectrumNumber + 1), GetIsolationMz(spectrumNumber + 1), GetIsolationWidth(spectrumNumber + 1), GetDissociationType(spectrumNumber + 1), GetPrecursorScanNumber(spectrumNumber + 1), GetPrecursorMonoisotopicIntensity(spectrumNumber + 1), GetPrecursorMonoisotopicMZ(spectrumNumber + 1));

        }

        private double GetPrecursorMonoisotopicIntensity(int spectrumNumber)
        {
            spectrumNumber--;
            foreach (Generated.CVParamType cv in _mzMLConnection.run.spectrumList.spectrum[spectrumNumber].precursorList.precursor[0].selectedIonList.selectedIon[0].cvParam)
            {
                if (cv.accession.Equals(_peakIntensity))
                {
                    return Convert.ToDouble(cv.value);
                }
            }
            throw new ArgumentNullException("Could not determine precursor intensity of spectrum " + spectrumNumber + 1);

        }

        private double GetPrecursorMonoisotopicMZ(int spectrumNumber)
        {
            spectrumNumber--;
            foreach (Generated.CVParamType cv in _mzMLConnection.run.spectrumList.spectrum[spectrumNumber].precursorList.precursor[0].selectedIonList.selectedIon[0].cvParam)
            {
                if (cv.accession.Equals(_precursorMass))
                {
                    return double.Parse(cv.value);
                }
            }
            throw new ArgumentNullException("Could not determine precursor monoisotopic mass for " + spectrumNumber + 1);

        }

        private double GetTotalIonCurrent(int spectrumNumber)
        {
            spectrumNumber--;
            foreach (Generated.CVParamType cv in _mzMLConnection.run.spectrumList.spectrum[spectrumNumber].cvParam)
            {
                if (cv.accession.Equals(_totalIonCurrent))
                {
                    return double.Parse(cv.value);
                }
            }
            throw new ArgumentNullException("Could not determine total ion current " + spectrumNumber + 1);

        }

        private int GetPrecursorScanNumber(int v)
        {
            return Convert.ToInt32(Regex.Match(GetPrecursorID(v), @"\d+$").Value);
        }
    }

    public static class MzmlMethods
    {

        public static void CreateAndWriteMyIndexedMZmlwithCalibratedSpectra(IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile, string outputFile)
        {
            Generated.indexedmzML _indexedmzMLConnection = new Generated.indexedmzML();
            _indexedmzMLConnection.mzML = new Generated.mzMLType();
            _indexedmzMLConnection.mzML.version = "1";

            _indexedmzMLConnection.mzML.cvList = new Generated.CVListType();
            _indexedmzMLConnection.mzML.cvList.count = "1";
            _indexedmzMLConnection.mzML.cvList.cv = new Generated.CVType[1];
            _indexedmzMLConnection.mzML.cvList.cv[0] = new Generated.CVType();
            _indexedmzMLConnection.mzML.cvList.cv[0].URI = @"https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo";
            _indexedmzMLConnection.mzML.cvList.cv[0].fullName = "Proteomics Standards Initiative Mass Spectrometry Ontology";
            _indexedmzMLConnection.mzML.cvList.cv[0].id = "MS";

            _indexedmzMLConnection.mzML.fileDescription = new Generated.FileDescriptionType();
            _indexedmzMLConnection.mzML.fileDescription.fileContent = new Generated.ParamGroupType();
            _indexedmzMLConnection.mzML.fileDescription.fileContent.cvParam = new Generated.CVParamType[2];
            _indexedmzMLConnection.mzML.fileDescription.fileContent.cvParam[0] = new Generated.CVParamType();
            _indexedmzMLConnection.mzML.fileDescription.fileContent.cvParam[0].accession = "MS:1000579"; // MS1 Data
            _indexedmzMLConnection.mzML.fileDescription.fileContent.cvParam[1] = new Generated.CVParamType();
            _indexedmzMLConnection.mzML.fileDescription.fileContent.cvParam[1].accession = "MS:1000580"; // MSn Data

            _indexedmzMLConnection.mzML.softwareList = new Generated.SoftwareListType();
            _indexedmzMLConnection.mzML.softwareList.count = "1";

            _indexedmzMLConnection.mzML.softwareList.software = new Generated.SoftwareType[1];
            // For a RAW file!!!
            // ToDo: read softwareList from mzML file
            //_indexedmzMLConnection.mzML.softwareList.software[1] = new SoftwareType();
            //_indexedmzMLConnection.mzML.softwareList.software[1].id = "ThermoSoftware";
            //_indexedmzMLConnection.mzML.softwareList.software[1].version = rawFile.GetSofwareVersion();
            //_indexedmzMLConnection.mzML.softwareList.software[1].cvParam = new CVParamType[1];
            //_indexedmzMLConnection.mzML.softwareList.software[1].cvParam[0] = new CVParamType();
            //_indexedmzMLConnection.mzML.softwareList.software[1].cvParam[0].accession = "MS:1000693";

            _indexedmzMLConnection.mzML.softwareList.software[0] = new Generated.SoftwareType();
            _indexedmzMLConnection.mzML.softwareList.software[0].id = "StefanSoftware";
            _indexedmzMLConnection.mzML.softwareList.software[0].version = "1";
            _indexedmzMLConnection.mzML.softwareList.software[0].cvParam = new Generated.CVParamType[1];
            _indexedmzMLConnection.mzML.softwareList.software[0].cvParam[0] = new Generated.CVParamType();
            _indexedmzMLConnection.mzML.softwareList.software[0].cvParam[0].accession = "MS:1000799";
            _indexedmzMLConnection.mzML.softwareList.software[0].cvParam[0].value = "StefanSoftware";


            // Leaving empty. Can't figure out the configurations. 
            // ToDo: read instrumentConfigurationList from mzML file
            _indexedmzMLConnection.mzML.instrumentConfigurationList = new Generated.InstrumentConfigurationListType();

            _indexedmzMLConnection.mzML.dataProcessingList = new Generated.DataProcessingListType();
            // Only writing mine! Might have had some other data processing (but not if it is a raw file)
            // ToDo: read dataProcessingList from mzML file
            _indexedmzMLConnection.mzML.dataProcessingList.count = "1";
            _indexedmzMLConnection.mzML.dataProcessingList.dataProcessing = new Generated.DataProcessingType[1];
            _indexedmzMLConnection.mzML.dataProcessingList.dataProcessing[0] = new Generated.DataProcessingType();
            _indexedmzMLConnection.mzML.dataProcessingList.dataProcessing[0].id = "StefanDataProcessing";


            _indexedmzMLConnection.mzML.run = new Generated.RunType();

            // ToDo: Finish the chromatogram writing!
            _indexedmzMLConnection.mzML.run.chromatogramList = new Generated.ChromatogramListType();
            _indexedmzMLConnection.mzML.run.chromatogramList.count = "1";
            _indexedmzMLConnection.mzML.run.chromatogramList.chromatogram = new Generated.ChromatogramType[1];
            _indexedmzMLConnection.mzML.run.chromatogramList.chromatogram[0] = new Generated.ChromatogramType();

            _indexedmzMLConnection.mzML.run.spectrumList = new Generated.SpectrumListType();
            _indexedmzMLConnection.mzML.run.spectrumList.count = (myMsDataFile.LastSpectrumNumber - myMsDataFile.FirstSpectrumNumber + 1).ToString();
            _indexedmzMLConnection.mzML.run.spectrumList.defaultDataProcessingRef = "StefanDataProcessing";
            _indexedmzMLConnection.mzML.run.spectrumList.spectrum = new Generated.SpectrumType[myMsDataFile.LastSpectrumNumber - myMsDataFile.FirstSpectrumNumber + 1];

            // Loop over all spectra
            int totalNumSpectra = myMsDataFile.LastSpectrumNumber - myMsDataFile.FirstSpectrumNumber + 1;
            for (int i = 0; i < totalNumSpectra; i++)
            {
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i] = new Generated.SpectrumType();

                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].defaultArrayLength = myMsDataFile.GetScan(i + myMsDataFile.FirstSpectrumNumber).MassSpectrum.Count;

                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].index = i.ToString();
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].id = myMsDataFile.GetScan(i + myMsDataFile.FirstSpectrumNumber).id;

                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam = new Generated.CVParamType[8];

                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[0] = new Generated.CVParamType();

                if (myMsDataFile.GetScan(i + myMsDataFile.FirstSpectrumNumber).MsnOrder == 1)
                {
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[0].accession = "MS:1000579";
                }
                else if (myMsDataFile.GetScan(i + myMsDataFile.FirstSpectrumNumber).MsnOrder == 2)
                {
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[0].accession = "MS:1000580";

                    // So needs a precursor!
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList = new Generated.PrecursorListType();
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.count = 1.ToString();
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor = new Generated.PrecursorType[1];
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0] = new Generated.PrecursorType();
                    string precursorID;
                    myMsDataFile.GetScan(i + myMsDataFile.FirstSpectrumNumber).TryGetPrecursorID(out precursorID);
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].spectrumRef = precursorID;
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].selectedIonList = new Generated.SelectedIonListType();
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].selectedIonList.count = 1.ToString();
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].selectedIonList.selectedIon = new Generated.ParamGroupType[1];
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].selectedIonList.selectedIon[0] = new Generated.ParamGroupType();
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].selectedIonList.selectedIon[0].cvParam = new Generated.CVParamType[3];
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].selectedIonList.selectedIon[0].cvParam[0] = new Generated.CVParamType();
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].selectedIonList.selectedIon[0].cvParam[0].name = "selected ion m/z";

                    double selectedIonGuesssMonoisotopicMZ;
                    myMsDataFile.GetScan(i + myMsDataFile.FirstSpectrumNumber).TryGetSelectedIonGuessMonoisotopicMZ(out selectedIonGuesssMonoisotopicMZ);
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].selectedIonList.selectedIon[0].cvParam[0].value = selectedIonGuesssMonoisotopicMZ.ToString();

                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].selectedIonList.selectedIon[0].cvParam[0].accession = "MS:1000744";
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].selectedIonList.selectedIon[0].cvParam[1] = new Generated.CVParamType();
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].selectedIonList.selectedIon[0].cvParam[1].name = "charge state";
                    int selectedIonGuessChargeStateGuess;
                    myMsDataFile.GetScan(i + myMsDataFile.FirstSpectrumNumber).TryGetSelectedIonGuessChargeStateGuess(out selectedIonGuessChargeStateGuess);
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].selectedIonList.selectedIon[0].cvParam[1].value = selectedIonGuessChargeStateGuess.ToString();
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].selectedIonList.selectedIon[0].cvParam[1].accession = "MS:1000041";
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].selectedIonList.selectedIon[0].cvParam[2] = new Generated.CVParamType();
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].selectedIonList.selectedIon[0].cvParam[2].name = "peak intensity";
                    double selectedIonGuesssMonoisotopicIntensity;
                    myMsDataFile.GetScan(i + myMsDataFile.FirstSpectrumNumber).TryGetSelectedIonGuessMonoisotopicIntensity(out selectedIonGuesssMonoisotopicIntensity);
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].selectedIonList.selectedIon[0].cvParam[2].value = selectedIonGuesssMonoisotopicIntensity.ToString();
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].selectedIonList.selectedIon[0].cvParam[2].accession = "MS:1000042";


                    MzRange isolationRange;
                    myMsDataFile.GetScan(i + myMsDataFile.FirstSpectrumNumber).TryGetIsolationRange(out isolationRange);
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].isolationWindow = new Generated.ParamGroupType();
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].isolationWindow.cvParam = new Generated.CVParamType[3];
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].isolationWindow.cvParam[0] = new Generated.CVParamType();
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].isolationWindow.cvParam[0].accession = "MS:1000827";
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].isolationWindow.cvParam[0].name = "isolation window target m/z";
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].isolationWindow.cvParam[0].value = isolationRange.Mean.ToString();
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].isolationWindow.cvParam[1] = new Generated.CVParamType();
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].isolationWindow.cvParam[1].accession = "MS:1000828";
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].isolationWindow.cvParam[1].name = "isolation window lower offset";
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].isolationWindow.cvParam[1].value = (isolationRange.Width / 2).ToString();
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].isolationWindow.cvParam[2] = new Generated.CVParamType();
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].isolationWindow.cvParam[2].accession = "MS:1000829";
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].isolationWindow.cvParam[2].name = "isolation window upper offset";
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].isolationWindow.cvParam[2].value = (isolationRange.Width / 2).ToString();


                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].activation = new Generated.ParamGroupType();
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].activation.cvParam = new Generated.CVParamType[1];
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].activation.cvParam[0] = new Generated.CVParamType();

                    DissociationType dissociationType;
                    myMsDataFile.GetScan(i + myMsDataFile.FirstSpectrumNumber).TryGetDissociationType(out dissociationType);
                    switch (dissociationType)
                    {
                        case DissociationType.HCD:
                            _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].activation.cvParam[0].accession = "MS:1000422";
                            _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].activation.cvParam[0].name = "beam-type collision-induced dissociation";
                            break;
                        case DissociationType.CID:
                            _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].activation.cvParam[0].accession = "MS:1000133";
                            _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].activation.cvParam[0].name = "collision-induced dissociation";
                            break;
                        case DissociationType.Unknown:
                            _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].activation.cvParam[0].accession = "MS:1000044";
                            _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].precursorList.precursor[0].activation.cvParam[0].name = "dissociation method";
                            break;
                    }



                }

                // OPTIONAL, but need for CSMSL reader. ms level
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[1] = new Generated.CVParamType();
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[1].name = "ms level";
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[1].accession = "MS:1000511";
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[1].value = myMsDataFile.GetScan(i + myMsDataFile.FirstSpectrumNumber).MsnOrder.ToString();

                // Centroid?
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[2] = new Generated.CVParamType();
                if (myMsDataFile.GetScan(i + myMsDataFile.FirstSpectrumNumber).isCentroid)
                {
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[2].name = "centroid spectrum";
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[2].accession = "MS:1000127";
                }
                else
                {
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[2].name = "profile spectrum";
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[2].accession = "MS:1000128";
                }

                // Polarity
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[3] = new Generated.CVParamType();
                if (myMsDataFile.GetScan(i + myMsDataFile.FirstSpectrumNumber).Polarity == Polarity.Negative)
                {
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[3].name = "negative scan";
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[3].accession = "MS:1000129";
                }
                else if (myMsDataFile.GetScan(i + myMsDataFile.FirstSpectrumNumber).Polarity == Polarity.Positive)
                {
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[3].name = "positive scan";
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[3].accession = "MS:1000130";
                }

                // Spectrum title
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[4] = new Generated.CVParamType();
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[4].name = "spectrum title";
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[4].accession = "MS:1000796";
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[4].value = myMsDataFile.GetScan(i + myMsDataFile.FirstSpectrumNumber).id;

                if ((myMsDataFile.GetScan(i + myMsDataFile.FirstSpectrumNumber).MassSpectrum.Count) > 0)
                {
                    // Lowest observed mz
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[5] = new Generated.CVParamType();
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[5].name = "lowest observed m/z";
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[5].accession = "MS:1000528";
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[5].value = myMsDataFile.GetScan(i + myMsDataFile.FirstSpectrumNumber).MassSpectrum.FirstX.ToString();

                    // Highest observed mz
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[6] = new Generated.CVParamType();
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[6].name = "highest observed m/z";
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[6].accession = "MS:1000527";
                    _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[6].value = myMsDataFile.GetScan(i + myMsDataFile.FirstSpectrumNumber).MassSpectrum.LastX.ToString();
                }


                // Total ion current
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[7] = new Generated.CVParamType();
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[7].name = "total ion current";
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[7].accession = "MS:1000285";
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].cvParam[7].value = myMsDataFile.GetScan(i + myMsDataFile.FirstSpectrumNumber).TotalIonCurrent.ToString();

                // Retention time
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList = new Generated.ScanListType();
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.count = "1";
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan = new Generated.ScanType[1];
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan[0] = new Generated.ScanType();
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan[0].cvParam = new Generated.CVParamType[2];
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan[0].cvParam[0] = new Generated.CVParamType();
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan[0].cvParam[0].name = "scan start time";
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan[0].cvParam[0].accession = "MS:1000016";
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan[0].cvParam[0].value = myMsDataFile.GetScan(i + myMsDataFile.FirstSpectrumNumber).RetentionTime.ToString();
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan[0].cvParam[0].unitCvRef = "UO";
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan[0].cvParam[0].unitAccession = "UO:0000031";
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan[0].cvParam[0].unitName = "minute";
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan[0].cvParam[1] = new Generated.CVParamType();
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan[0].cvParam[1].name = "filter string";
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan[0].cvParam[1].accession = "MS:1000512";
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan[0].cvParam[1].value = myMsDataFile.GetScan(i + myMsDataFile.FirstSpectrumNumber).ScanFilter;

                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan[0].scanWindowList = new Generated.ScanWindowListType();
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan[0].scanWindowList.count = 1;
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan[0].scanWindowList.scanWindow = new Generated.ParamGroupType[1];
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan[0].scanWindowList.scanWindow[0] = new Generated.ParamGroupType();
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan[0].scanWindowList.scanWindow[0].cvParam = new Generated.CVParamType[2];
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan[0].scanWindowList.scanWindow[0].cvParam[0] = new Generated.CVParamType();
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan[0].scanWindowList.scanWindow[0].cvParam[0].name = "scan window lower limit";
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan[0].scanWindowList.scanWindow[0].cvParam[0].accession = "MS:1000501";
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan[0].scanWindowList.scanWindow[0].cvParam[0].value = myMsDataFile.GetScan(i + myMsDataFile.FirstSpectrumNumber).ScanWindowRange.Minimum.ToString();
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan[0].scanWindowList.scanWindow[0].cvParam[1] = new Generated.CVParamType();
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan[0].scanWindowList.scanWindow[0].cvParam[1].name = "scan window upper limit";
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan[0].scanWindowList.scanWindow[0].cvParam[1].accession = "MS:1000500";
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].scanList.scan[0].scanWindowList.scanWindow[0].cvParam[1].value = myMsDataFile.GetScan(i + myMsDataFile.FirstSpectrumNumber).ScanWindowRange.Maximum.ToString();


                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList = new Generated.BinaryDataArrayListType();

                // ONLY WRITING M/Z AND INTENSITY DATA, NOT THE CHARGE! (but can add charge info later)
                // CHARGE (and other stuff) CAN BE IMPORTANT IN ML APPLICATIONS!!!!!
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.count = 2.ToString();
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray = new Generated.BinaryDataArrayType[2];

                // M/Z Data
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[0] = new Generated.BinaryDataArrayType();
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[0].binary = Mzml.ConvertDoublestoBase64(myMsDataFile.GetScan(i + myMsDataFile.FirstSpectrumNumber).MassSpectrum.xArray, false);
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[0].encodedLength = (4 * Math.Ceiling(((double)_indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[0].binary.Length / 3))).ToString();
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[0].cvParam = new Generated.CVParamType[2];
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[0].cvParam[0] = new Generated.CVParamType();
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[0].cvParam[0].accession = "MS:1000514";
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[0].cvParam[0].name = "m/z array";
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[0].cvParam[1] = new Generated.CVParamType();
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[0].cvParam[1].accession = "MS:1000523";
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[0].cvParam[1].name = "64-bit float";
                //_indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[0].cvParam[0] = new CVParamType();
                //_indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[0].cvParam[0].accession = "MS:1000574";
                //_indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[0].cvParam[0].name = "zlib compression";

                // Intensity Data
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[1] = new Generated.BinaryDataArrayType();
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[1].binary = Mzml.ConvertDoublestoBase64(myMsDataFile.GetScan(i + myMsDataFile.FirstSpectrumNumber).MassSpectrum.yArray, false);
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[1].encodedLength = (4 * Math.Ceiling(((double)_indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[1].binary.Length / 3))).ToString();
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[1].cvParam = new Generated.CVParamType[2];
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[1].cvParam[0] = new Generated.CVParamType();
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[1].cvParam[0].accession = "MS:1000515";
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[1].cvParam[0].name = "intensity array";
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[1].cvParam[1] = new Generated.CVParamType();
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[1].cvParam[1].accession = "MS:1000523";
                _indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[1].cvParam[1].name = "64-bit float";
                //_indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[1].cvParam[0] = new CVParamType();
                //_indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[1].cvParam[0].accession = "MS:1000574";
                //_indexedmzMLConnection.mzML.run.spectrumList.spectrum[i].binaryDataArrayList.binaryDataArray[1].cvParam[0].name = "zlib compression";
            }

            Mzml.Write(outputFile, _indexedmzMLConnection);

        }

    }

}

