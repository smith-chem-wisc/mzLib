using Easy.Common.Extensions;
using MassSpectrometry;
using MzLibUtil;
using System.Collections.Concurrent;
using System.Globalization;
using System.Security.Cryptography;
using ThermoFisher.CommonCore.Data.Business;
using ThermoFisher.CommonCore.Data.FilterEnums;
using ThermoFisher.CommonCore.Data.Interfaces;
using ThermoFisher.CommonCore.RawFileReader;

// old namespace to ensure backwards compatibility
namespace IO.ThermoRawFileReader
{
    public class ThermoRawFileReader : Readers.ThermoRawFileReader
    {
        public ThermoRawFileReader(string path) : base(path) { }
    }

    public class ThermoRawFileReaderLicence : Readers.ThermoRawFileReaderLicence
    {

    }
}

// This .cs file uses:
// RawFileReader reading tool. Copyright © 2016 by Thermo Fisher Scientific, Inc. All rights reserved.
// See the full Software Licence Agreement for detailed requirements for use.
namespace Readers
{
    // I think that ThermoRawFileReader should be used to store the data from the files, 
    // but the actual implementation details should be completely hidden. 
    public class ThermoRawFileReader : MsDataFile
    {
        private IRawDataPlus? dynamicConnection;
        private int[] MsOrdersByScan;
        public ThermoRawFileReader(string path) : base(path) { }

        public override MsDataFile LoadAllStaticData(FilteringParams filteringParams = null, int maxThreads = 1)
        {
            if (!File.Exists(FilePath))
            {
                throw new FileNotFoundException();
            }

            // I don't know why this line needs to be here, but it does...
            var temp = RawFileReaderAdapter.FileFactory(FilePath);

            using (var threadManager = RawFileReaderFactory.CreateThreadManager(FilePath))
            {
                var rawFileAccessor = threadManager.CreateThreadAccessor();

                if (rawFileAccessor.IsError)
                {
                    throw new MzLibException("Error opening RAW file!");
                }

                if (!rawFileAccessor.IsOpen)
                {
                    throw new MzLibException("Unable to access RAW file!");
                }

                if (rawFileAccessor.InAcquisition)
                {
                    throw new MzLibException("RAW file still being acquired!");
                }

                rawFileAccessor.SelectInstrument(Device.MS, 1);
                var msDataScans = new MsDataScan[rawFileAccessor.RunHeaderEx.LastSpectrum];

                Parallel.ForEach(Partitioner.Create(0, msDataScans.Length),
                    new ParallelOptions { MaxDegreeOfParallelism = maxThreads }, (fff, loopState) =>
                    {
                        using (var myThreadDataReader = threadManager.CreateThreadAccessor())
                        {
                            myThreadDataReader.SelectInstrument(Device.MS, 1);

                            for (int s = fff.Item1; s < fff.Item2; s++)
                            {
                                var scan = GetOneBasedScan(myThreadDataReader, filteringParams, s + 1);
                                msDataScans[s] = scan;
                            }
                        }
                    });


                rawFileAccessor.Dispose();
                Scans = msDataScans;
                SourceFile = GetSourceFile();
            }
            
            temp.Dispose();

            return this;
        }

        public override SourceFile GetSourceFile()
        {
            string sendCheckSum;
            using (FileStream stream = File.OpenRead(FilePath))
            {
                SHA1 sha = SHA1.Create();
                byte[] checksum = sha.ComputeHash(stream);
                sendCheckSum = BitConverter.ToString(checksum)
                    .Replace("-", string.Empty);
            }

            SourceFile sourceFile = new SourceFile(
                @"Thermo nativeID format",
                @"Thermo RAW format",
                sendCheckSum,
                @"SHA-1",
                FilePath,
                Path.GetFileNameWithoutExtension(FilePath));

            return sourceFile;
        }

        /// <summary>
        /// Initiates a dynamic connection with a Thermo .raw file. Data can be "streamed" instead of loaded all at once. Use 
        /// GetOneBasedScanFromDynamicConnection to get data from a particular scan. Use CloseDynamicConnection to close the 
        /// dynamic connection after all desired data has been retrieved from the dynamic connection.
        /// </summary>
        public override void InitiateDynamicConnection()
        {
            if (!File.Exists(FilePath))
            {
                throw new FileNotFoundException();
            }

            if (dynamicConnection != null)
            {
                dynamicConnection.Dispose();
            }

            dynamicConnection = RawFileReaderAdapter.FileFactory(FilePath);

            if (!dynamicConnection.IsOpen)
            {
                throw new MzLibException("Unable to access RAW file!");
            }

            if (dynamicConnection.IsError)
            {
                throw new MzLibException("Error opening RAW file!");
            }

            if (dynamicConnection.InAcquisition)
            {
                throw new MzLibException("RAW file still being acquired!");
            }

            dynamicConnection.SelectInstrument(Device.MS, 1);

            GetMsOrdersByScanInDynamicConnection();
        }

        /// <summary>
        /// Allows access to a .raw file one scan at a time via an open dynamic connection. Returns null if the raw file does not contain the 
        /// scan number specified. Use InitiateDynamicConnection to open a dynamic connection before using this method.
        /// </summary>
        public override MsDataScan GetOneBasedScanFromDynamicConnection(int oneBasedScanNumber, IFilteringParams filterParams = null)
        {
            var dymConnection = RawFileReaderAdapter.FileFactory(FilePath);
            dymConnection.SelectInstrument(Device.MS, 1);

            if (dymConnection == null)
            {
                throw new MzLibException("The dynamic connection has not been created yet!");
            }

            if (oneBasedScanNumber > dymConnection.RunHeaderEx.LastSpectrum ||
                oneBasedScanNumber < dymConnection.RunHeaderEx.FirstSpectrum)
            {
                return null;
            }

            return ThermoRawFileReader.GetOneBasedScan(dymConnection, filterParams, oneBasedScanNumber);
        }

        /// <summary>
        /// Disposes of the dynamic connection, if one is open.
        /// </summary>
        public override void CloseDynamicConnection()
        {
            if (dynamicConnection != null)
            {
                dynamicConnection.Dispose();
            }
        }

        public override int[] GetMsOrderByScanInDynamicConnection()
        {

            if (dynamicConnection != null)
            {

                int lastSpectrum = dynamicConnection.RunHeaderEx.LastSpectrum;
                var scanEvents = dynamicConnection.GetScanEvents(1, lastSpectrum);

                int[] msorders = scanEvents.Select(p => (int)p.MSOrder).ToArray();

                return msorders;
            }

            return null;
        }

        /// <summary>
        /// This method ensures backwards compatibility with previous mzLib implementations
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="filteringParams"></param>
        /// <param name="maxThreads"></param>
        /// <returns></returns>
        public static MsDataFile LoadAllStaticData(string filePath, FilteringParams filteringParams = null,
            int maxThreads = 1)
        {
            return MsDataFileReader.GetDataFile(filePath).LoadAllStaticData(filteringParams, maxThreads);
        }

        private static MsDataScan GetOneBasedScan(IRawDataPlus rawFile, IFilteringParams filteringParams,
            int scanNumber)
        {
            var filter = rawFile.GetFilterForScanNumber(scanNumber);

            string scanFilterString = filter.ToString();

            int msOrder = (int)filter.MSOrder;
            if (msOrder < 1 || msOrder > 10)
            {
                throw new MzLibException("Unknown MS Order (" + msOrder + ") for scan number " + scanNumber);
            }

            string nativeId = "controllerType=0 controllerNumber=1 scan=" + scanNumber;
            MzSpectrum spectrum = GetSpectrum(rawFile, filteringParams, scanNumber, scanFilterString, msOrder);

            var scanStats = rawFile.GetScanStatsForScanNumber(scanNumber);
            double scanRangeHigh = scanStats.HighMass;
            double scanRangeLow = scanStats.LowMass;
            MzRange scanWindowRange = new(scanRangeLow, scanRangeHigh);

            double? ionInjectionTime = null;
            double? precursorSelectedMonoisotopicIonMz = null;
            int? selectedIonChargeState = null;
            double? ms2IsolationWidth = null;
            int? precursorScanNumber = null;
            double? isolationMz = null;
            string HcdEnergy = null;
            string scanDescript = null;
            ActivationType activationType = ActivationType.Any; // thermo enum
            DissociationType dissociationType = DissociationType.Unknown; // mzLib enum

            var trailer = rawFile.GetTrailerExtraInformation(scanNumber);
            string[] labels = trailer.Labels;
            string[] values = trailer.Values;

            for (int i = 0; i < trailer.Labels.Length; i++)
            {
                if (labels[i].StartsWith("Ion Injection Time (ms)", StringComparison.Ordinal))
                {
                    ionInjectionTime = double.Parse(values[i], CultureInfo.InvariantCulture) == 0 ?
                        (double?)null :
                        double.Parse(values[i], CultureInfo.InvariantCulture);
                }

                if (msOrder < 2)
                {
                    continue;
                }

                if (labels[i].StartsWith("MS" + msOrder + " Isolation Width", StringComparison.Ordinal))
                {
                    ms2IsolationWidth = double.Parse(values[i], CultureInfo.InvariantCulture) == 0 ?
                        (double?)null :
                        double.Parse(values[i], CultureInfo.InvariantCulture);
                }

                if (labels[i].StartsWith("Monoisotopic M/Z", StringComparison.Ordinal))
                {
                    precursorSelectedMonoisotopicIonMz = double.Parse(values[i], CultureInfo.InvariantCulture) == 0 ?
                        (double?)null :
                        double.Parse(values[i], CultureInfo.InvariantCulture);
                }

                if (labels[i].StartsWith("Charge State", StringComparison.Ordinal))
                {
                    selectedIonChargeState = int.Parse(values[i], CultureInfo.InvariantCulture) == 0 ?
                        (int?)null :
                        int.Parse(values[i], CultureInfo.InvariantCulture);
                }

                if (labels[i].StartsWith("Master Scan Number", StringComparison.Ordinal)
                    || labels[i].StartsWith("Master Index", StringComparison.Ordinal))
                {
                    precursorScanNumber = int.Parse(values[i], CultureInfo.InvariantCulture) <= 1 ?
                        (int?)null :
                        int.Parse(values[i], CultureInfo.InvariantCulture);
                }

                if (labels[i].StartsWith("HCD Energy:", StringComparison.Ordinal))
                {
                    HcdEnergy = values[i];
                }

                if (labels[i].StartsWith("Scan Description", StringComparison.Ordinal))
                {
                    scanDescript = values[i].TrimEnd();
                }
            }

            if (msOrder > 1)
            {
                var scanEvent = rawFile.GetScanEventForScanNumber(scanNumber);
                var reaction = scanEvent.GetReaction(0);
                isolationMz = reaction.PrecursorMass;
                activationType = reaction.ActivationType;

                dissociationType = GetDissociationType(activationType);

                // thermo does not have an enum value for ETHcD, so this needs to be detected from the scan filter
                if (scanFilterString.Contains("@etd", StringComparison.OrdinalIgnoreCase)
                    && scanFilterString.Contains("@hcd", StringComparison.OrdinalIgnoreCase))
                {
                    dissociationType = DissociationType.EThcD;
                }

                if (ms2IsolationWidth == null)
                {
                    ms2IsolationWidth = reaction.IsolationWidth;
                }

                if (precursorScanNumber == null)
                {
                    // we weren't able to get the precursor scan number, so we'll have to guess;
                    // loop back to find precursor scan
                    // (assumed to be the first scan before this scan with an MS order of this scan's MS order - 1)
                    // e.g., if this is an MS2 scan, find the first MS1 scan before this and assume that's the precursor scan
                    for (int i = scanNumber; i >= 1; i--)
                    {
                        var possiblePrecursorScanFilter = rawFile.GetFilterForScanNumber(i);
                        int order = (int)possiblePrecursorScanFilter.MSOrder;
                        if (order == msOrder - 1)
                        {
                            precursorScanNumber = i;
                            break;
                        }
                    }

                    if (precursorScanNumber == null)
                    {
                        throw new MzLibException("Could not get precursor for scan #" + scanNumber);
                    }
                }
            }

            // at this point, we have the m/z value of the species that got fragmented, from the scan header
            // this section of the code finds that peak in the spectrum (it's actual intensity and centroided m/z values)
            // the intention is to remove any rounding issues caused by what is in the scan header and what is observable in the spectrum
            double? selectedIonIntensity = null;

            if (isolationMz.HasValue)
            {
                if (spectrum.Size != 0)
                {
                    int closest = spectrum.GetClosestPeakIndex(isolationMz.Value);

                    double mz = spectrum.XArray[closest];
                    double intensity = spectrum.YArray[closest];

                    if (Math.Abs(mz - isolationMz.Value) < 0.1)
                    {
                        selectedIonIntensity = intensity;
                        isolationMz = mz;
                    }
                }
            }

            return new MsDataScan(
                massSpectrum: spectrum,
                oneBasedScanNumber: scanNumber,
                msnOrder: msOrder,
                isCentroid: true,
                polarity: GetPolarity(filter.Polarity),
                retentionTime: rawFile.RetentionTimeFromScanNumber(scanNumber),
                scanWindowRange: scanWindowRange,
                scanFilter: scanFilterString,
                mzAnalyzer: GetMassAnalyzerType(filter.MassAnalyzer),
                totalIonCurrent: spectrum.SumOfAllY,
                injectionTime: ionInjectionTime,
                noiseData: null, //TODO: implement reading noise data. it's unused right now, so it's just left as null
                nativeId: nativeId,
                selectedIonMz: isolationMz,
                selectedIonChargeStateGuess: selectedIonChargeState,
                selectedIonIntensity: selectedIonIntensity,
                isolationMZ: isolationMz,
                isolationWidth: ms2IsolationWidth,
                dissociationType: dissociationType,
                oneBasedPrecursorScanNumber: precursorScanNumber,
                selectedIonMonoisotopicGuessMz: precursorSelectedMonoisotopicIonMz,
                hcdEnergy: HcdEnergy,
                scanDescription: scanDescript);
        }

        private static MzSpectrum GetSpectrum(IRawDataPlus rawFile, IFilteringParams filterParams,
            int scanNumber, string scanFilter, int scanOrder)
        {
            MzSpectrum spectrum;
            double[] xArray;
            double[] yArray;

            if (string.IsNullOrEmpty(scanFilter))
            {
                return new MzSpectrum(new double[0], new double[0], false);
            }

            var centroidStream = rawFile.GetCentroidStream(scanNumber, false);

            // PreferredMasses should be used if centroidStream data is null; it's probably ITMS data
            if (centroidStream.Masses == null || centroidStream.Intensities == null)
            {
                var scan = Scan.FromFile(rawFile, scanNumber);

                var mzs = scan.PreferredMasses;
                xArray = scan.PreferredMasses;
                yArray = scan.PreferredIntensities;

                if (xArray == null || yArray == null)
                {
                    throw new MzLibException("Could not centroid data from scan " + scanNumber);
                }
            }
            else
            {
                xArray = centroidStream.Masses;
                yArray = centroidStream.Intensities;
            }

            //Remove Zero Intensity Peaks
            double zeroEquivalentIntensity = 0.01;
            int zeroIntensityCount = yArray.Count(i => i < zeroEquivalentIntensity);
            int intensityValueCount = yArray.Count();
            if (zeroIntensityCount > 0 && zeroIntensityCount < intensityValueCount)
            {
                Array.Sort(yArray, xArray);
                double[] nonZeroIntensities = new double[intensityValueCount - zeroIntensityCount];
                double[] nonZeroMzs = new double[intensityValueCount - zeroIntensityCount];
                yArray = yArray.SubArray(zeroIntensityCount, intensityValueCount - zeroIntensityCount);
                xArray = xArray.SubArray(zeroIntensityCount, intensityValueCount - zeroIntensityCount);
                Array.Sort(xArray, yArray);
            }

            if (filterParams != null
                && xArray.Length > 0
                && (filterParams.MinimumAllowedIntensityRatioToBasePeakM.HasValue || filterParams.NumberOfPeaksToKeepPerWindow.HasValue)
                && ((filterParams.ApplyTrimmingToMs1 && scanOrder == 1) || (filterParams.ApplyTrimmingToMsMs && scanOrder > 1)))
            {
                var count = xArray.Length;

                var mzArray = new double[count];
                var intensityArray = new double[count];
                Array.Copy(xArray, mzArray, count);
                Array.Copy(yArray, intensityArray, count);

                var scanStats = rawFile.GetScanStatsForScanNumber(scanNumber);
                double scanRangeHigh = scanStats.HighMass;
                double scanRangeLow = scanStats.LowMass;

                WindowModeHelper.Run(ref intensityArray, ref mzArray, filterParams,
                     scanRangeLow, scanRangeHigh);

                Array.Sort(mzArray, intensityArray);
                spectrum = new MzSpectrum(mzArray, intensityArray, false);
            }
            else
            {
                spectrum = new MzSpectrum(xArray, yArray, false);
            }

            return spectrum;
        }

        private static MZAnalyzerType GetMassAnalyzerType(MassAnalyzerType massAnalyzerType)
        {
            switch (massAnalyzerType)
            {
                case MassAnalyzerType.MassAnalyzerFTMS: return MZAnalyzerType.Orbitrap;
                case MassAnalyzerType.MassAnalyzerITMS: return MZAnalyzerType.IonTrap2D;
                case MassAnalyzerType.MassAnalyzerSector: return MZAnalyzerType.Sector;
                case MassAnalyzerType.MassAnalyzerTOFMS: return MZAnalyzerType.TOF;

                default: return MZAnalyzerType.Unknown;
            }
        }

        private static MassSpectrometry.Polarity GetPolarity(PolarityType polarity)
        {
            switch (polarity)
            {
                case PolarityType.Positive: return MassSpectrometry.Polarity.Positive;
                case PolarityType.Negative: return MassSpectrometry.Polarity.Negative;

                default: throw new MzLibException("Cannot interpret polarity type: " + polarity);
            }
        }

        private static DissociationType GetDissociationType(ActivationType activationType)
        {
            switch (activationType)
            {
                case ActivationType.CollisionInducedDissociation: return DissociationType.CID;
                case ActivationType.ElectronTransferDissociation: return DissociationType.ETD;
                case ActivationType.HigherEnergyCollisionalDissociation: return DissociationType.HCD;
                case ActivationType.ElectronCaptureDissociation: return DissociationType.ECD;
                case ActivationType.MultiPhotonDissociation: return DissociationType.MPD;
                case ActivationType.PQD: return DissociationType.PQD;
                case ActivationType.UltraVioletPhotoDissociation: return DissociationType.UVPD;
                case ActivationType.NegativeElectronTransferDissociation: return DissociationType.NETD;

                default: return DissociationType.Unknown;
            }
        }

        /// <summary>
        /// Gets all the MS orders of all scans in a dynamic connection. This is useful if you want to open all MS1 scans
        /// without loading all of the other MSn scans.
        /// </summary>
        private int[] GetMsOrdersByScanInDynamicConnection()
        {
            if (MsOrdersByScan.IsNotNullOrEmpty())
            {
                return MsOrdersByScan;
            }

            if (dynamicConnection != null)
            {
                int lastSpectrum = dynamicConnection.RunHeaderEx.LastSpectrum;
                var scanEvents = dynamicConnection.GetScanEvents(1, lastSpectrum);

                int[] msorders = scanEvents.Select(p => (int)p.MSOrder).ToArray();

                MsOrdersByScan = msorders;

                return msorders;
            }

            return null;
        }
    }
}
