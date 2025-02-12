using MzLibUtil;
using Readers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ThermoFisher.CommonCore.Data;

namespace MassSpectrometry
{
    public class TimsDataScan : MsDataScan
    {
        public int ScanNumberStart { get; }
        public int ScanNumberEnd { get; }
        public double OneOverK0 { get; }
        public int? PrecursorId { get; }
        public long FrameId { get; }
        /// <summary>
        /// For PASEF Aggregate scans, contains the list of Frames where the same precursor was samples
        /// This is a list of succesive PASEF scans capturing data on the same ion-mobility scan range and quadrupole isolation window
        /// </summary>
        public List<long> FrameIds { get; }
        internal int ComponentSpectraTotalPeaks { get; private set; }

        // Need to incorporate scan range somehow
        public TimsDataScan(MzSpectrum massSpectrum,
            int oneBasedScanNumber,
            int msnOrder,
            bool isCentroid,
            Polarity polarity,
            double retentionTime,
            MzRange scanWindowRange,
            string scanFilter,
            MZAnalyzerType mzAnalyzer,
            double totalIonCurrent,
            double? injectionTime,
            double[,] noiseData,
            string nativeId,
            long frameId,
            int scanNumberStart,
            int scanNumberEnd,
            double medianOneOverK0,
            int? precursorId = null,
            double? selectedIonMz = null,
            int? selectedIonChargeStateGuess = null,
            double? selectedIonIntensity = null,
            double? isolationMZ = null,
            double? isolationWidth = null,
            DissociationType? dissociationType = null,
            int? oneBasedPrecursorScanNumber = null,
            double? selectedIonMonoisotopicGuessMz = null,
            string hcdEnergy = null,
            List<long> frames = null) : 
            base(massSpectrum, oneBasedScanNumber, msnOrder, isCentroid, polarity,
                retentionTime, scanWindowRange, scanFilter, mzAnalyzer, totalIonCurrent,
                injectionTime, noiseData, nativeId, selectedIonMz, selectedIonChargeStateGuess,
                selectedIonIntensity, isolationMZ, isolationWidth, dissociationType,
                oneBasedPrecursorScanNumber, selectedIonMonoisotopicGuessMz, hcdEnergy)
        {
            FrameId = frameId;
            FrameIds = frames;
            ScanNumberStart = scanNumberStart;
            ScanNumberEnd = scanNumberEnd;
            OneOverK0 = medianOneOverK0;
            PrecursorId = precursorId;
            ComponentSpectraTotalPeaks = 0;
        }

        internal void AverageComponentSpectra(FrameProxyFactory proxyFactory, FilteringParams filteringParams = null)
        {
            MassSpectrum = TofSpectraMerger.MergeArraysToMs2Spectrum(MzArrays, IntensityArrays, filteringParams);
            TotalIonCurrent = MassSpectrum.SumOfAllY;
            MzArrays.Clear();
            IntensityArrays.Clear();
        }

        internal List<double[]> MzArrays;
        internal List<int[]> IntensityArrays;

        internal void AddComponentArrays(double[] mzs, int[] intensities)
        {
            if (MzArrays == null)
            {
                MzArrays = new();
                IntensityArrays = new();
            }
            MzArrays.Add(mzs);
            IntensityArrays.Add(intensities);
        }

        public MzSpectrum?[] Ms1SpectraIndexedByZeroBasedScanNumber { get; private set; }

        public void AddMs1Spectrum(MzSpectrum spectrum, int scanNumber)
        {
            if (Ms1SpectraIndexedByZeroBasedScanNumber.IsNullOrEmpty())
            {
                Ms1SpectraIndexedByZeroBasedScanNumber = new MzSpectrum[ScanNumberEnd - ScanNumberStart + 1];
            }
            Ms1SpectraIndexedByZeroBasedScanNumber[scanNumber - ScanNumberStart] = spectrum;
        }
    }
}
