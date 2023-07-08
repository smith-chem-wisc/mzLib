using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ.PeakPicking
{
    /// <summary>
    /// This class handles situations where a cluster of isobaric species elute close together.
    /// </summary>
    public class IsobarCluster
    {
        public readonly double PeakFindingMass;
        public Xic ReferenceXic { get; private set; }
        public List<Xic> OtherXics { get; private set; }
        /// <summary>
        /// The file with the most isobaric identifications is defined as the reference file
        /// </summary>
        public SpectraFileInfo ReferenceFile { get; private set; }
        public List<SpectraFileInfo> OtherFiles { get; }
        public List<Identification> Identifications { get; private set; }
        public DoubleRange RetentionTimeRange { get; private set; }
        public IsobarCluster(List<Identification> identifications, PeakIndexingEngine peakIndexingEngine)
        {
            Identifications = identifications;
            PeakFindingMass = identifications.GroupBy(id => id.PrecursorChargeState)
                .MaxBy(group => group.Count())
                .Average(group => group.PeakfindingMass);
            // Group IDs by file
            //  Define reference file
            //  Find time range for each file
            //  Set Cluster time range (all files) as the min ID time to maxID time, +- Max(15s, time range*0.25)
            ReferenceFile = Identifications.GroupBy(id => id.FileInfo)
                .MaxBy(group => group.Count())
                .Key;
            OtherFiles = Identifications.Select(id => id.FileInfo)
                .Distinct()
                .Where(file => file != ReferenceFile)
                .ToList();
            SetTimeRange();
            // Create XICs
            //  Pull XIC for each file in time range
            //  Set XIC for file w/ most IDs as reference
            //  Align XICs to reference
            PullXics(peakIndexingEngine);
            // Reconcile Extrema
            // Define Regions
            // Assign IDs to regions 
            // Check for discrepancies
        }

        private void SetTimeRange()
        {
            var firstMs2Rt = Identifications.Select(i => i.Ms2RetentionTimeInMinutes).Min();
            var lastMs2Rt = Identifications.Select(i => i.Ms2RetentionTimeInMinutes).Max();
            double ms2RtSpan = lastMs2Rt - firstMs2Rt;
            double timeBuffer = ms2RtSpan < 1
                ? 0.25
                : (ms2RtSpan) / 4.0;
            RetentionTimeRange = new DoubleRange(Math.Max(0.0, firstMs2Rt - timeBuffer), lastMs2Rt + timeBuffer); 
        }

        private void PullXics(PeakIndexingEngine peakIndexingEngine)
        {
            List<IndexedMassSpectralPeak> refPeaks = peakIndexingEngine.ExtractPeaks(PeakFindingMass, ReferenceFile);
            ReferenceXic = new Xic(refPeaks, PeakFindingMass, ReferenceFile, 0, referenceXic: true);
            OtherXics = new();
            foreach(SpectraFileInfo file in OtherFiles)
            {
                List<IndexedMassSpectralPeak> peaks = peakIndexingEngine.ExtractPeaks(PeakFindingMass, file);
                double rtAdjustment = XicProcessing.AlignPeaks(refPeaks, peaks);
                OtherXics.Add(new Xic(peaks, PeakFindingMass, file, rtAdjustment, referenceXic: false));
            }
        }

        public IsobarCluster(List<Identification> identifications, List<Xic> xics)
        {
            Identifications = identifications;
            OtherXics = xics;

            ReferenceFile = Identifications.GroupBy(id => id.FileInfo)
                .ToDictionary(keySelector: group => group.Key, elementSelector: group => group.Count())
                .MaxBy(kvp => kvp.Value)
                .Key;

            // build reference XIC
            ReferenceXic = OtherXics.First(x => x.SpectraFile == ReferenceFile);

            // Foreach loop to build non-reference XICs
        }
    }
}
