using MzLibUtil;
using Readers.Bruker.TimsTofReader;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public class TimsDataScan : MsDataScan
    {
        public int ScanNumberStart { get; }
        public int ScanNumberEnd { get; }
        public double OneOverK0 { get; }
        public int PrecursorId { get; }
        public long FrameId { get; }
        /// <summary>
        /// For PASEF Aggregate scans, contains the list of Frames where the same precursor was samples
        /// This is a list of succesive PASEF scans capturing data on the same ion-mobility scan range and quadrupole isolation window
        /// </summary>
        public List<long> FrameIds { get; }
        internal List<MzSpectrum> ComponentSpectra { get; private set; }
        internal List<ListNode<TofPeak>> ComponentSpectraListNodes { get; private set; }
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
            int precursorId,
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
            if(msnOrder > 1)
            {
                ComponentSpectra = new();
            }
        }

        public void AverageComponentSpectra()
        {
            // TODO: Probably need to add, like, checks and stuff. But oh well.
            MassSpectrum = TofSpectraMerger.MergeSpectra(ComponentSpectra);
            ComponentSpectra.Clear();
            //MassSpectrum = TofSpectraMerger.MergeSpectra(ComponentSpectraListNodes, ComponentSpectraTotalPeaks);
            //ComponentSpectraListNodes.Clear();
        }

        public void AddComponentSpectrum(MzSpectrum spectrum)
        {
            if (ComponentSpectra == null) return;
            ComponentSpectra.Add(spectrum);
        }

        //internal void AddComponentSpectrum(ListNode<TofPeak> spectrumHead, int spectrumLength)
        //{
        //    if (ComponentSpectraListNodes == null)
        //    {
        //        ComponentSpectraListNodes = new();
        //    }
        //    ComponentSpectraListNodes.Add(spectrumHead);
        //    ComponentSpectraTotalPeaks += spectrumLength;
        //}

    }
}
