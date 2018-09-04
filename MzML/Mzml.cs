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
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Security.Cryptography;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace IO.MzML
{
    public class Mzml : MsDataFile
    {
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

        private static readonly Dictionary<string, Polarity> polarityDictionary = new Dictionary<string, Polarity>
        {
            {"MS:1000129",Polarity.Negative},
            {"MS:1000130",Polarity.Positive}
        };

        private static readonly Dictionary<string, MZAnalyzerType> analyzerDictionary = new Dictionary<string, MZAnalyzerType>
        {
            { "MS:1000443", MZAnalyzerType.Unknown},
            { "MS:1000081",MZAnalyzerType.Quadrupole},
            { "MS:1000291",MZAnalyzerType.IonTrap2D},
            { "MS:1000082",MZAnalyzerType.IonTrap3D},
            { "MS:1000484",MZAnalyzerType.Orbitrap},
            { "MS:1000084",MZAnalyzerType.TOF},
            { "MS:1000079",MZAnalyzerType.FTICR},
            { "MS:1000080",MZAnalyzerType.Sector}
        };

        private static readonly Dictionary<string, DissociationType> dissociationDictionary = new Dictionary<string, DissociationType>
        {
            { "MS:1000133",DissociationType.CID},
            { "MS:1001880",DissociationType.ISCID},
            { "MS:1000422",DissociationType.HCD},
            { "MS:1000598",DissociationType.ETD},
            { "MS:1000435",DissociationType.IRMPD},
            { "MS:1000599",DissociationType.PQD},
            { "MS:1000044",DissociationType.Unknown}
        };

        private Mzml(MsDataScan[] scans, SourceFile sourceFile) : base(scans, sourceFile)
        {
        }

        public static Mzml LoadAllStaticData(string filePath, FilteringParams filterParams = null, int maxThreads = -1)
        {
            if (!File.Exists(filePath))
            {
                throw new FileNotFoundException();
            }

            Generated.mzMLType _mzMLConnection;

            try
            {
                using (FileStream fs = new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.Read))
                {
                    var _indexedmzMLConnection = (Generated.indexedmzML)MzmlMethods.indexedSerializer.Deserialize(fs);
                    _mzMLConnection = _indexedmzMLConnection.mzML;
                }
            }
            catch
            {
                using (FileStream fs = new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.Read))
                    _mzMLConnection = (Generated.mzMLType)MzmlMethods.mzmlSerializer.Deserialize(fs);
            }

            SourceFile sourceFile;
            if (_mzMLConnection.fileDescription.sourceFileList != null && _mzMLConnection.fileDescription.sourceFileList.sourceFile != null && _mzMLConnection.fileDescription.sourceFileList.sourceFile[0] != null && _mzMLConnection.fileDescription.sourceFileList.sourceFile[0].cvParam != null)
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
                using (FileStream stream = File.OpenRead(filePath))
                {
                    using (SHA1Managed sha = new SHA1Managed())
                    {
                        byte[] checksum = sha.ComputeHash(stream);
                        sendCheckSum = BitConverter.ToString(checksum)
                            .Replace("-", string.Empty);
                    }
                }
                sourceFile = new SourceFile(
                    @"no nativeID format",
                    @"mzML format",
                    sendCheckSum,
                    @"SHA-1",
                    Path.GetFullPath(filePath),
                    Path.GetFileNameWithoutExtension(filePath));
            }

            var numSpecta = _mzMLConnection.run.spectrumList.spectrum.Length;
            MsDataScan[] scans = new MsDataScan[numSpecta];

            Parallel.ForEach(Partitioner.Create(0, numSpecta), new ParallelOptions { MaxDegreeOfParallelism = maxThreads }, fff =>
            {
                for (int i = fff.Item1; i < fff.Item2; i++)
                {
                    scans[i] = GetMsDataOneBasedScanFromConnection(_mzMLConnection, i + 1, filterParams);
                }
            });

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
                    throw new MzLibException("Scan number " + scan.OneBasedScanNumber.ToString() + " appeared multiple times in " + filePath);
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
                            scans[i].setOneBasedPrecursorScanNumber(scans[j].OneBasedScanNumber);
                            break;
                        }
                    }
                }
            }

            return new Mzml(scans, sourceFile);
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
                else if (analyzerDictionary.TryGetValue(configs[0].componentList.analyzer[0].cvParam[0].accession, out MZAnalyzerType returnVal))
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
                        analyzerDictionary.TryGetValue(configs[i].componentList.analyzer[0].cvParam[0].accession, out MZAnalyzerType returnVal);
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
                    tic = double.Parse(cv.value);
                }
                if (polarity.Equals(Polarity.Unknown))
                {
                    polarityDictionary.TryGetValue(cv.accession, out polarity);
                }
            }

            if (!msOrder.HasValue || !isCentroid.HasValue)
                throw new MzLibException("!msOrder.HasValue || !isCentroid.HasValue");

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

            if (filterParams != null && intensities.Length > 0 && (filterParams.MinimumAllowedIntensityRatioToBasePeakM.HasValue || filterParams.NumberOfPeaksToKeepPerWindow.HasValue)
                && ((filterParams.ApplyTrimmingToMs1 && msOrder.Value == 1) || (filterParams.ApplyTrimmingToMsMs && msOrder.Value > 1)))
            {
                if (filterParams.NumberOfWindows == null)
                {
                    int numPeaks = TopNpeakHelper(ref intensities, ref masses, filterParams);
                    Array.Resize(ref intensities, numPeaks);
                    Array.Resize(ref masses, numPeaks);
                }
                else
                {
                    WindowModeHelper(ref intensities, ref masses, filterParams);
                }
            }
            Array.Sort(masses, intensities);
            var mzmlMzSpectrum = new MzSpectrum(masses, intensities, false);

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
                        rtInMinutes = double.Parse(cv.value);
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
                        injectionTime = double.Parse(cv.value);
                    }
                    if (cv.accession.Equals(_oneBasedScanNumber)) //get the real one based spectrum number (if available), the other assumes they are in order. This is present in .mgf->.mzml conversions from MSConvert
                    {
                        oneBasedScanNumber = int.Parse(cv.value);
                    }
                }
            }

            double high = double.NaN;
            double low = double.NaN;

            if (_mzMLConnection.run.spectrumList.spectrum[oneBasedIndex - 1].scanList.scan[0].scanWindowList != null)
            {
                foreach (Generated.CVParamType cv in _mzMLConnection.run.spectrumList.spectrum[oneBasedIndex - 1].scanList.scan[0].scanWindowList.scanWindow[0].cvParam)
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
            }

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
                    selectedIonMz = double.Parse(cv.value);
                }
                if (cv.accession.Equals(_precursorCharge))
                {
                    selectedIonCharge = int.Parse(cv.value);
                }
                if (cv.accession.Equals(_peakIntensity))
                {
                    selectedIonIntensity = double.Parse(cv.value);
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
                        isolationMz = double.Parse(cv.value);
                    }
                    if (cv.accession.Equals(_isolationWindowLowerOffset))
                    {
                        lowIsolation = double.Parse(cv.value);
                    }
                    if (cv.accession.Equals(_isolationWindowUpperOffset))
                    {
                        highIsolation = double.Parse(cv.value);
                    }
                }
            }

            DissociationType dissociationType = DissociationType.Unknown;
            if (_mzMLConnection.run.spectrumList.spectrum[oneBasedIndex - 1].precursorList.precursor[0].activation.cvParam != null)
            {
                foreach (Generated.CVParamType cv in _mzMLConnection.run.spectrumList.spectrum[oneBasedIndex - 1].precursorList.precursor[0].activation.cvParam)
                {
                    dissociationDictionary.TryGetValue(cv.accession, out dissociationType);
                }
            }
            double? monoisotopicMz = null;
            if (_mzMLConnection.run.spectrumList.spectrum[oneBasedIndex - 1].scanList.scan[0].userParam != null)
            {
                foreach (var userParam in _mzMLConnection.run.spectrumList.spectrum[oneBasedIndex - 1].scanList.scan[0].userParam)
                {
                    if (userParam.name.EndsWith("Monoisotopic M/Z:"))
                    {
                        monoisotopicMz = double.Parse(userParam.value);
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
        private static double[] ConvertBase64ToDoubles(byte[] bytes, bool zlibCompressed = false, bool is32bit = true)
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

        private static int GetOneBasedPrecursorScanNumber(Generated.mzMLType _mzMLConnection, int oneBasedSpectrumNumber)
        {
            string precursorID = _mzMLConnection.run.spectrumList.spectrum[oneBasedSpectrumNumber - 1].precursorList.precursor[0].spectrumRef;
            do
            {
                oneBasedSpectrumNumber--;
            } while (!precursorID.Equals(_mzMLConnection.run.spectrumList.spectrum[oneBasedSpectrumNumber - 1].id));
            return oneBasedSpectrumNumber;
        }
    }
}