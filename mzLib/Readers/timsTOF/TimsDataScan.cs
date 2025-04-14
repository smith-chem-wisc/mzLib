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

        public List<TimsSpectrum> ComponentSpectra { get; private set; }

        /// <summary>
        /// Average component MS2 spectra to create a single MS2 spectrum
        /// </summary>
        /// <param name="proxyFactory"></param>
        /// <param name="filteringParams"></param>
        internal void AverageComponentSpectra(FrameProxyFactory proxyFactory, TofSpectraMerger spectraMerger, FilteringParams filteringParams = null)
        {
            MassSpectrum = spectraMerger.CreateMzSpectrum(ComponentSpectra, proxyFactory, msnLevel: 2, filteringParams);
            TotalIonCurrent = MassSpectrum.SumOfAllY;
            ComponentSpectraTotalPeaks = ComponentSpectra.Sum(s => s.Size);
            ComponentSpectra = null;
        }

        internal void AddComponentSpectrum(TimsSpectrum spectrum)
        {
            if(spectrum==null) return;
            ComponentSpectra ??= new();
            ComponentSpectra.Add(spectrum);
        }
    }

    /// <summary>
    /// This is similar to an mz spectrum, but much more lightweight
    /// It stores intensities as ints and tof indices instead of mz values
    /// </summary>
    public class TimsSpectrum
    {
        public uint[] XArray { get; init; }
        public int[] YArray { get; init; }

        public int Size => XArray.Length;
        public int ZeroIndexedIonMobilityScanNumber { get; init; }

        public TimsSpectrum(uint[] tofIndices, int[] intensities)
        {
            XArray = tofIndices;
            YArray = intensities;
            ZeroIndexedIonMobilityScanNumber =  -1;
        }

        public TimsSpectrum(uint[] tofIndices, int[] intensities, int zeroIndexedTimsScanIndex)
        {
            XArray = tofIndices;
            YArray = intensities;
            ZeroIndexedIonMobilityScanNumber = zeroIndexedTimsScanIndex;
        }
    }
}
