using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Easy.Common.Extensions;

namespace FlashLFQ.PeakPicking
{
    /// <summary>
    /// This class handles situations where a cluster of isobaric species elute close together.
    /// </summary>
    public class IsobarCluster
    {
        public readonly double PeakFindingMass;

        public List<ChromatographicPeak> Peaks;
        public List<Identification> Identifications { get; }
        public DoubleRange RetentionTimeRange { get; private set; }
        /// <summary>
        /// The file with the most unambiguous peaks. The xics for all other
        /// files are aligned against the reference file
        /// </summary>
        public SpectraFileInfo ReferenceFile { get; }
        public Dictionary<SpectraFileInfo, Xic> Xics { get; private set; }

        // Unclear if we need this or the pruned extrema can be stored in the xic
        public Dictionary<SpectraFileInfo, List<Extremum>> ConsensusExtrema { get; private set; }
        
        public Dictionary<int, List<Identification>> RegionIdDictionary;

        public IsobarCluster(List<ChromatographicPeak> peaks, PeakIndexingEngine peakIndexingEngine)
        {
            Peaks = peaks;
            Identifications = peaks.SelectMany(peak => peak.Identifications).ToList();

            PeakFindingMass = Identifications
                .Average(group => group.PeakfindingMass);

            // The file with the most unambiguous peaks is chosen as the reference file that all other files are aligned to
            ReferenceFile = Peaks.Where(peak => peak.Identifications.IsNotNullOrEmpty() && peak.Identifications.Count ==1)
                .GroupBy(peak => peak.SpectraFileInfo)
                .MaxBy(group => group.Count())
                .Key;

            SetTimeRange();

            PullXics(peakIndexingEngine);

            FindConsensusExtrema();

            DefinePeakRegions();

            AssignIDs();

        }

        private void SetTimeRange()
        {
            double firstPeakStart = Peaks.Min(peak => peak.IsotopicEnvelopes.Min(e => e.IndexedPeak.RetentionTime));
            double lastPeakEnd = Peaks.Max(peak => peak.IsotopicEnvelopes.Max(e => e.IndexedPeak.RetentionTime));
            RetentionTimeRange = new DoubleRange(firstPeakStart, lastPeakEnd);
        }

        /// <summary>
        /// Create XICs for each spectra file. XICs will be 1.5x longer than the RetentionTimeRange
        /// as to increase accuracy of peak alignment
        /// </summary>
        /// <param name="peakIndexingEngine"></param>
        private void PullXics(PeakIndexingEngine peakIndexingEngine)
        {
            
            double timeBuffer = RetentionTimeRange.Width < 1
                ? 0.25
                : (RetentionTimeRange.Width) / 4.0;

            
            Xics = new Dictionary<SpectraFileInfo, Xic>();

            //TODO: add time span param to Extract peaks and pull rtRange.min - buffer, rtRange.max + buffer
            List<IndexedMassSpectralPeak> refPeaks = peakIndexingEngine.ExtractPeaks(PeakFindingMass, ReferenceFile);
            Xics.Add(ReferenceFile, new Xic(refPeaks, PeakFindingMass, ReferenceFile, 0, referenceXic: true));
            
            List<SpectraFileInfo> otherFiles = Peaks.Select(peak => peak.SpectraFileInfo)
                .Distinct()
                .Where(file => !file.Equals(ReferenceFile))
                .ToList();

            foreach(SpectraFileInfo file in otherFiles)
            {
                List<IndexedMassSpectralPeak> peaks = peakIndexingEngine.ExtractPeaks(PeakFindingMass, file);
                double rtAdjustment = XicProcessing.AlignPeaks(refPeaks, peaks);
                Xics.Add(file, new Xic(peaks, PeakFindingMass, file, rtAdjustment, referenceXic: false));
            }
        }

        /// <summary>
        /// Finds matching extrema between each run/XIC, then writes them to the
        /// ConsensusExtrema dictionary
        /// </summary>
        private void FindConsensusExtrema()
        {
            List<Xic> otherXics = Xics
                .Where(kvp => !kvp.Key.Equals(ReferenceFile))
                .Select(kvp => kvp.Value)
                .ToList();
            Extremum[,] extrema2dArray = XicProcessing
                .ReconcileExtrema(Xics[ReferenceFile].Extrema, otherXics.Select(xic => xic.Extrema).ToList());

            ConsensusExtrema = new();
            List<int> colIndices = Enumerable.Range(0, extrema2dArray.GetLength(1)).ToList();
            List <Extremum> firstRow = colIndices
                .Select(col => extrema2dArray[0, col])
                .Where(e => e != null)
                .ToList();
            ConsensusExtrema.Add(ReferenceFile, firstRow);

            foreach (int row in Enumerable.Range(1, extrema2dArray.GetLength(0)-1) )
            {
                List<Extremum> extrema = colIndices
                    .Select(col => extrema2dArray[row, col])
                    .Where(e => e != null)
                    .ToList();
                ConsensusExtrema.Add(otherXics[row-1].SpectraFile, extrema);
            }
        }
        /// <summary>
        /// Every XIC has a dictionary mapping an integer to a double range called PeakRegions.
        /// A peak region must start with either a minimum or the beginning of an XIC, contain at least
        /// one maximum, and end with either a minimum of the end of the XIC. The time span will be 
        /// slightly different for each XIC, but the integer key will map to the same peak region
        /// for each XIC in the Isobar cluster. Identifications will then be sorted into peak regions, and 
        /// discrepancies will be resolved.
        /// </summary>
        internal void DefinePeakRegions()
        {
            // write a dictionary for each file that maps the region number (key, same for each file)
            // to a double range containing the start and end times of the region (slightly different for each file)
            foreach (var kvp in ConsensusExtrema)
            {
                Dictionary<int, DoubleRange> peakRegions = new();
                int regionNumber = 1;
                double startTime = RetentionTimeRange.Minimum;
                bool maxFound = kvp.Value.First().ExtremumType == ExtremumType.Maximum;

                // Iterating through each set of extrema to define regions for each file
                foreach (Extremum extrema in kvp.Value)
                {
                    if (extrema.ExtremumType == ExtremumType.Minimum)
                    {
                        if (!maxFound)
                            startTime = extrema.RetentionTime;
                        else
                        {
                            peakRegions.Add(regionNumber++, new DoubleRange(startTime, extrema.RetentionTime));
                            maxFound = false;
                            startTime = extrema.RetentionTime;
                        }
                    }
                    else if (extrema.ExtremumType == ExtremumType.Maximum)
                    {
                        maxFound = true;
                    }
                }

                Xics[kvp.Key].PeakRegions = peakRegions;
            }
        }

        internal void AssignIDs()
        {
            RegionIdDictionary = new Dictionary<int, List<Identification>>();

            foreach (Identification id in Identifications)
            {
                // Find the region the id belongs to
                // region boundaries are slightly different for each file, which is why it's 
                // done this way
                Xic idXic = Xics[id.FileInfo];
                foreach (var region in idXic.PeakRegions)
                {
                    if (region.Value.Contains(id.Ms2RetentionTimeInMinutes))
                    {
                        if(!RegionIdDictionary.ContainsKey(region.Key))
                            RegionIdDictionary.Add(region.Key, new List<Identification>());

                        RegionIdDictionary[region.Key].Add(id);
                        break;
                    }
                }
            }

            // This is going to be less than ideal, but here we go
            // take every kvp in the region id dictionary
            // transform into Dict<FullSeq, List<int> regions>
            //  any full seq with multiple regions - handle that case
            //  any regions with multiple IDs - handle that case

            #region OnePeptideMultiRegion

            Dictionary<string, List<int>> fullSeqToRegion = new();

            foreach (var kvp in RegionIdDictionary)
            {
                foreach (var id in kvp.Value)
                {
                    string[] fullSeqs = id.ModifiedSequence.Split('|');
                    foreach (string seq in fullSeqs)
                    {
                        if (!fullSeqToRegion.TryAdd(seq, new List<int> { kvp.Key }))
                        {
                            fullSeqToRegion[seq].Add(kvp.Key);
                        }
                    }
                }
            }

            foreach (var kvp in fullSeqToRegion)
            {
                if (kvp.Value.Distinct().Count() > 1)
                {
                    // Handle case where one peptide just has multiple peaks (e.g., mod on aromatic)
                    // Handle case where peptide is assigned differently in different files
                    // do something
                }
            }

            #endregion

            #region OneRegionMultiPeptide

            foreach (var kvp in RegionIdDictionary)
            {
                if (kvp.Value.Select(id => id.ModifiedSequence).Distinct().Count() > 1)
                {
                    // Does this happen in the same file? e.g. coeluting isobars both observed in MS2
                    // Is there disagreement between files? 
                }
            }

            #endregion

        }

        internal void AssignPeaks()
        {

        }

        public IsobarCluster(List<Identification> identifications, List<Xic> xics)
        {
            Identifications = identifications;
            //OtherXics = xics;

            ReferenceFile = Identifications.GroupBy(id => id.FileInfo)
                .ToDictionary(keySelector: group => group.Key, elementSelector: group => group.Count())
                .MaxBy(kvp => kvp.Value)
                .Key;

            // build reference XIC
            //ReferenceXic = OtherXics.First(x => x.SpectraFile == ReferenceFile);

            // Foreach loop to build non-reference XICs
        }
    }
}
