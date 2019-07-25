﻿using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Concurrent;
using System.Globalization;
using System.IO;
using System.Security.Cryptography;
using System.Threading.Tasks;
using ThermoFisher.CommonCore.Data.Business;
using ThermoFisher.CommonCore.Data.FilterEnums;
using ThermoFisher.CommonCore.Data.Interfaces;
using ThermoFisher.CommonCore.RawFileReader;
using UsefulProteomicsDatabases;

// This .cs file uses:
// RawFileReader reading tool. Copyright © 2016 by Thermo Fisher Scientific, Inc. All rights reserved.
// See the full Software Licence Agreement for detailed requirements for use.

namespace ThermoRawFileReader
{
    public class ThermoRawFileReaderData : MsDataFile
    {
        protected ThermoRawFileReaderData(MsDataScan[] scans, SourceFile sourceFile) : base(scans, sourceFile)
        {
        }

        /// <summary>
        /// Loads all scan data from a Thermo .raw file.
        /// </summary>
        public static ThermoRawFileReaderData LoadAllStaticData(string filePath, IFilteringParams filterParams = null, int maxThreads = -1)
        {
            if (!File.Exists(filePath))
            {
                throw new FileNotFoundException();
            }

            Loaders.LoadElements();

            // I don't know why this line needs to be here, but it does...
            var temp = RawFileReaderAdapter.FileFactory(filePath);

            var threadManager = RawFileReaderFactory.CreateThreadManager(filePath);
            var rawFileAccessor = threadManager.CreateThreadAccessor();

            if (!rawFileAccessor.IsOpen)
            {
                throw new MzLibException("Unable to access RAW file!");
            }

            if (rawFileAccessor.IsError)
            {
                throw new MzLibException("Error opening RAW file!");
            }

            if (rawFileAccessor.InAcquisition)
            {
                throw new MzLibException("RAW file still being acquired!");
            }

            rawFileAccessor.SelectInstrument(Device.MS, 1);
            var msDataScans = new MsDataScan[rawFileAccessor.RunHeaderEx.LastSpectrum];

            Parallel.ForEach(Partitioner.Create(0, msDataScans.Length), new ParallelOptions { MaxDegreeOfParallelism = maxThreads }, (fff, loopState) =>
            {
                IRawDataPlus myThreadDataReader = threadManager.CreateThreadAccessor();
                myThreadDataReader.SelectInstrument(Device.MS, 1);

                for (int s = fff.Item1; s < fff.Item2; s++)
                {
                    try
                    {
                        var scan = GetOneBasedScan(myThreadDataReader, filterParams, s + 1);
                        msDataScans[s] = scan;
                    }
                    catch (Exception ex)
                    {
                        throw new MzLibException("Error reading scan " + (s + 1) + ": " + ex.Message);
                    }
                }
            });

            rawFileAccessor.Dispose();

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

            SourceFile sourceFile = new SourceFile(
                @"Thermo nativeID format",
                @"Thermo RAW format",
                sendCheckSum,
                @"SHA-1",
                filePath,
                Path.GetFileNameWithoutExtension(filePath));

            return new ThermoRawFileReaderData(msDataScans, sourceFile);
        }

        public static MsDataScan GetOneBasedScan(IRawDataPlus rawFile, IFilteringParams filteringParams, int scanNumber)
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
            MzRange scanWindowRange = new MzRange(scanRangeLow, scanRangeHigh);

            double? ionInjectionTime = null;
            double? precursorSelectedMonoisotopicIonMz = null;
            int? selectedIonChargeState = null;
            double? ms2IsolationWidth = null;
            int? precursorScanNumber = null;
            double? isolationMz = null;
            ActivationType activationType = ActivationType.Any;

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
            }

            if (msOrder > 1)
            {
                var scanEvent = rawFile.GetScanEventForScanNumber(scanNumber);
                var reaction = scanEvent.GetReaction(0);
                isolationMz = reaction.PrecursorMass;
                activationType = reaction.ActivationType;

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
                int? closest = spectrum.GetClosestPeakIndex(isolationMz.Value);

                if (closest.HasValue)
                {
                    double mz = spectrum.XArray[closest.Value];
                    double intensity = spectrum.YArray[closest.Value];

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
                dissociationType: GetDissociationType(activationType),
                oneBasedPrecursorScanNumber: precursorScanNumber,
                selectedIonMonoisotopicGuessMz: precursorSelectedMonoisotopicIonMz);
        }

        private static MzSpectrum GetSpectrum(IRawDataPlus rawFile, IFilteringParams filterParams, int scanNumber, string scanFilter, int scanOrder)
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

                WindowModeHelper(ref intensityArray, ref mzArray, filterParams, scanRangeLow, scanRangeHigh);

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

                default: return DissociationType.Unknown;
            }
        }
    }
}
