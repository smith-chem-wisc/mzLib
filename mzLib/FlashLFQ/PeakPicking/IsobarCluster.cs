using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using Chemistry;
using Easy.Common.Extensions;

namespace FlashLFQ.PeakPicking
{
    /// <summary>
    /// This class handles situations where a cluster of isobaric species elute close together.
    /// </summary>
    public class IsobarCluster
    {
        public readonly double PeakFindingMz;
        public readonly int ChargeState;

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

        /// <summary>
        /// Dictionary linking each peak region to a list of chromatographic peaks with
        /// apexes inside that region
        /// </summary>
        private Dictionary<int, List<ChromatographicPeak>> _regionPeakDictionary;

        public Dictionary<int, List<ChromatographicPeak>> RegionPeakDictionary
        {
            get { return _regionPeakDictionary; }
        }

        /// <summary>
        /// This dictionary is used to determine what peaks should be used for peptide quantification
        /// </summary>
        public Dictionary<string, List<int>> SequenceRegionDictionary { get; private set; }

        private FlashLfqEngine _flashLfqEngine;
        private FlashLfqResults _results;

        //TODO: Need to link IDs/Peptides to peaks within the isobar cluster class

        public IsobarCluster(List<ChromatographicPeak> peaks, FlashLfqEngine flashLfqEngine, 
            FlashLfqResults results = null, DoubleRange rtRange = null)
        {
            // I have no idea what causes a peak to be apex-less.
            // TODO: Figure out why some peaks have no apex
            if (peaks.Any(peak => peak.Apex == null))
            {
                throw new ArgumentException("All peaks must have a defined apex");
            }

            Peaks = peaks;
            Identifications = peaks.SelectMany(peak => peak.Identifications).ToList();
            _flashLfqEngine = flashLfqEngine;
            _results = results;

            ChargeState = Peaks
                .Where(p => p?.Apex != null)
                .Select(p => p.Apex.ChargeState)
                .GroupBy(z => z)
                .OrderByDescending(group => group.Count())
                .First().Key;
            PeakFindingMz = Identifications
                .Average(group => group.PeakfindingMass.ToMz(ChargeState));

            // Sometimes there are no peaks with only one ID. In that case, we look for the file with the most peaks with
            // the fewest number of associated IDs
            int minPeakIds = Peaks
                .Where(peak => peak.Identifications.IsNotNullOrEmpty())
                .Min(peak => peak.Identifications.Count);
            // The file with the most unambiguous peaks is chosen as the reference file that all other files are aligned to
            ReferenceFile = Peaks.Where(peak => peak.Identifications.IsNotNullOrEmpty() && peak.Identifications.Count == minPeakIds)
                .GroupBy(peak => peak.SpectraFileInfo)
                .MaxBy(group => group.Count())
                .Key;

            SetTimeRange(rtRange);

            PullXics();

            FindConsensusExtrema();

            DefinePeakRegions();
        }

        private void SetTimeRange(DoubleRange rtRange)
        {
            double firstPeakStart = Peaks.Min(peak => peak.IsotopicEnvelopes.Min(e => e.IndexedPeak.RetentionTime));
            double lastPeakEnd = Peaks.Max(peak => peak.IsotopicEnvelopes.Max(e => e.IndexedPeak.RetentionTime));
            RetentionTimeRange = rtRange ?? new DoubleRange(firstPeakStart, lastPeakEnd);
        }

        /// <summary>
        /// Create XICs for each spectra file. XICs will be 1.5x longer than the RetentionTimeRange
        /// as to increase accuracy of peak alignment
        /// </summary>
        /// <param name="peakIndexingEngine"></param>
        private void PullXics()
        {
            
            double timeBuffer = RetentionTimeRange.Width < 1
                ? .25
                : (RetentionTimeRange.Width) / 4.0;

            Xics = new Dictionary<SpectraFileInfo, Xic>();

            //TODO: This is deserializing each file multiple times, which comes with a performance hit
            List<IndexedMassSpectralPeak> refPeaks = _flashLfqEngine.IndexDict[ReferenceFile]
                .ExtractPeaks(PeakFindingMz,
                    startTime: RetentionTimeRange.Minimum - timeBuffer,
                    endTime: RetentionTimeRange.Maximum + timeBuffer);
            if (refPeaks.Count <= 5)
            {
                throw new Exception("XIC too short");
            }
            Xics.Add(ReferenceFile, new Xic(refPeaks, PeakFindingMz, ReferenceFile, 0, referenceXic: true));
            
            List<SpectraFileInfo> otherFiles = Peaks.Select(peak => peak.SpectraFileInfo)
                .Distinct()
                .Where(file => !file.Equals(ReferenceFile))
                .ToList();

            foreach(SpectraFileInfo file in otherFiles)
            {
                List<IndexedMassSpectralPeak> peaks = _flashLfqEngine.IndexDict[file]
                    .ExtractPeaks(PeakFindingMz,
                        startTime: RetentionTimeRange.Minimum - timeBuffer,
                        endTime: RetentionTimeRange.Maximum + timeBuffer);
                double rtAdjustment = XicProcessing.AlignPeaks(refPeaks, peaks);
                if (peaks.Count <= 5)
                {
                    throw new Exception("XIC too short");
                }
                Xics.Add(file, new Xic(peaks, PeakFindingMz, file, rtAdjustment, referenceXic: false));
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

                //Xics[kvp.Key].PeakRegions = peakRegions;
                Xics[kvp.Key].SetPeakRegions(peakRegions);
            }
        }


        public void ReassignPeakIDs()
        {
            _regionPeakDictionary = new Dictionary<int, List<ChromatographicPeak>>();
            SequenceRegionDictionary = new Dictionary<string, List<int>>();
            foreach(int region in Xics[ReferenceFile].PeakRegions.Keys)
                _regionPeakDictionary.Add(region, new List<ChromatographicPeak>());

            // Link each region to a list of IDs
            foreach (ChromatographicPeak peak in Peaks)
            {
                //TODO: group by spectra file info so we're not pulling the same xic multiple times
                Xic xicContainingPeak = Xics[peak.SpectraFileInfo];
                foreach (var region in xicContainingPeak.PeakRegions)
                {
                    if (region.Value.Contains(peak.ApexRetentionTime))
                    {
                        if (!_regionPeakDictionary.ContainsKey(region.Key))
                            _regionPeakDictionary.Add(region.Key, new List<ChromatographicPeak>());
                        _regionPeakDictionary[region.Key].Add(peak);
                        peak.Region = region.Key;

                        break;
                    }
                }
            }

            // Something went wrong
            // TODO: Make this throw an exception
            if (_regionPeakDictionary.Count == 0)
                return;

            #region OnePeptideMultiRegion

            Dictionary<string, List<int>> msmsIdFullSeqToRegion = new();
            foreach (var kvp in _regionPeakDictionary)
            {
                foreach (var id in kvp.Value
                             .Where(peak => !peak.IsMbrPeak)
                             .SelectMany(peak => peak.Identifications))
                {
                    string[] fullSeqs = id.ModifiedSequence.Split('|');
                    foreach (string seq in fullSeqs)
                    {
                        if (!msmsIdFullSeqToRegion.TryAdd(seq, new List<int> { kvp.Key }))
                        {
                            msmsIdFullSeqToRegion[seq].Add(kvp.Key);
                        }
                    }
                }
            }

            // This isn't gonna work for MBR. There, we could have one MSMS id and still need to run this
            foreach (var seqToRegionsKvp in msmsIdFullSeqToRegion)
            {
                if (seqToRegionsKvp.Value.Distinct().Count() > 1)
                {
                    // This algorithm needs to assign peaks such that ambiguity is minimized

                    // Handle case where consensus is peak ambiguity, e.g. most files have two peptides
                    // assigned to one peak.
                    List<int> regions = seqToRegionsKvp.Value;

                    // For each region (key), contains the number of IDs made in that region
                    var regionCountDictionary = regions
                        .GroupBy(i => i)
                        .ToDictionary(
                            keySelector: group => group.Key,
                            elementSelector: group => group.Count());

                    if (regionCountDictionary.Values.Max() == regionCountDictionary.Values.Min())
                    {
                        //TODO: Come up with tie-breaking mechanism
                        continue;
                    }

                    // TODO: Come up with tie-breaking mechanism here. Like, if two regions each have 3 IDs
                    int bestRegion = regionCountDictionary.MaxBy(kvp => kvp.Value).Key;
                    foreach (string seq in _regionPeakDictionary[bestRegion]
                                .SelectMany(peak => peak.Identifications)
                                .SelectMany(id => id.ModifiedSequence.Split('|')))
                    {
                        //TODO: Need to weight unambiguous identifications higher 
                        if (!SequenceRegionDictionary.ContainsKey(seq))
                            SequenceRegionDictionary.Add(seq, new List<int>());
                        SequenceRegionDictionary[seq].Add(bestRegion);
                    }

                    foreach (int region in seqToRegionsKvp.Value)
                    {
                        if (region == bestRegion) 
                        {
                            continue; // These were accurately assigned already
                        }
                        // Peaks are added to the dictionary within the foreach loop, so it's important to stash the 
                        // peaks before entering the loop
                        foreach (var peak in _regionPeakDictionary[region])
                        {
                            // Remove and reassign the id
                            if (peak.Identifications
                                .SelectMany(id => id.ModifiedSequence.Split('|'))
                                .Contains(seqToRegionsKvp.Key))
                            {
                                int idIndex = peak.Identifications
                                    .FindIndex(id => id.ModifiedSequence == seqToRegionsKvp.Key);
                                Identification id = null;
                                // If an identification with a modified sequence exactly matching
                                // the sequence in question (i.e., not ambiguous and separated by '|'),
                                // Then we remove the id and add it to a different peak
                                if (idIndex >= 0)
                                {
                                    id = peak.Identifications[idIndex];
                                    peak.RemoveIdentification(id);
                                }
                                // Otherwise, we borrow an identification, just like in MBR
                                else
                                {
                                    // If there is an id with an exact sequence match (i.e., non-ambiguous), we use that
                                    id = _regionPeakDictionary[bestRegion]
                                            .SelectMany(peak => peak.Identifications)
                                            .FirstOrDefault(id => id.ModifiedSequence.Equals(seqToRegionsKvp.Key))
                                         // Otherwise, we use an ambiguous id 
                                         ?? _regionPeakDictionary[bestRegion].First().Identifications.First();
                                    
                                }
                                
                                var bestPeak = _regionPeakDictionary[bestRegion]
                                    .FirstOrDefault(p => p.SpectraFileInfo.Equals(peak.SpectraFileInfo));

                                if (bestPeak == null)
                                {
                                    var xic = Xics[peak.SpectraFileInfo];
                                    double rtApex = xic.PeakRegions[bestRegion].Mean;

                                    // Need to find a better way to avoid serializing/deserializing
                                    // XIC contains all the indexes ms peaks you would need (for a given charge state)
                                    // Probably just pass that in

                                    // This is faster, but it feels extremely fragile. Side-effects almost guaranteed.
                                    // There's a refactor a-coming, but not today
                                    _flashLfqEngine.SwapPeakIndexingEngine(peak.SpectraFileInfo);

                                    bestPeak = _flashLfqEngine.GetChromatographicPeak(
                                        id, 
                                        peak.SpectraFileInfo, 
                                        isMbrPeak: peak.SpectraFileInfo != id.FileInfo, // TODO: Add an alternate classication system for ambgiuous peaks. MBR isn't totally accurate 
                                        rtApex, 
                                        xic.PeakRegions[bestRegion]);
                                    _regionPeakDictionary[bestRegion].Add(bestPeak);

                                    //TODO: Clean this up
                                    if(_results != null)
                                        _results.Peaks[peak.SpectraFileInfo].Add(bestPeak);
                                }
                                bestPeak.AddIdentification(id);
                                bestPeak.Region = bestRegion;
                                foreach(string seq in id.ModifiedSequence.Split('|'))
                                {
                                    if (!SequenceRegionDictionary.ContainsKey(seq))
                                        SequenceRegionDictionary.Add(seq, new List<int>());
                                    SequenceRegionDictionary[seq].Add(bestRegion);
                                }
                            }

                            //TODO: Do additional MBR stuff here. If there's a file without a peak in the region, create one
                        }
                    }
                    // Assign every id to a new peak in the best region. 

                    // Differential assignment
                    // One identification is assigned to different peaks in different files 
                    // Handle case where one peptide just has multiple peaks (e.g., mod on aromatic)

                    // TODO: MBR for ambiguous peaks, improved mbr
                }
            }

            #endregion

            #region OneRegionMultiPeptide

            foreach (var kvp in _regionPeakDictionary)
            {
                if (kvp.Value
                        .SelectMany(peak => peak.Identifications)
                        .Select(id => id.ModifiedSequence)
                        .Distinct().Count() > 1)
                {
                    // Does this happen in the same file? e.g. coeluting isobars both observed in MS2
                    // Is there disagreement between files? 
                    // Is this even necessary? Like, this is just co-eluting species
                    // If none of the peptides are assigned to different regions, can't we just leave well enough alone
                }
            }

            #endregion

        }

        // Clusters things with identical peak finding masses
        // Pull chromPeaks from all files, the recursively groups them together
        public static List<IsobarCluster> FindIsobarClusters(List<Identification> identifications, 
            FlashLfqResults results, PeakIndexingEngine indexingEngine, FlashLfqEngine flashLfqEngine, 
            List<Exception> exceptions = null)
        {
            var isobarGroups = identifications
                .GroupBy(id => Math.Round(id.PeakfindingMass, 2));
            List<IsobarCluster> isobarClusters = new();

            foreach (var group in isobarGroups)
            {
                List<Identification> isobaricIds = group.ToList();
                List<ChromatographicPeak> isobaricPeaks = new();
                foreach (var peaksKvp in results.Peaks)
                {
                    isobaricPeaks.AddRange(peaksKvp.Value
                            .Where(peak => peak.Identifications.Intersect(isobaricIds).Any())
                            .ToList());
                }
                // In MBR, one ID can lead to many peaks. However, if there's only one peak, we don't need to be concerned
                if (isobaricPeaks.Count() <= 1) continue;
                isobaricPeaks.Sort();

                List<List<ChromatographicPeak>> clusteredPeaks = ClusterPeaks(isobaricPeaks);
                foreach (var peakCluster in 
                         clusteredPeaks.Where(peakList => peakList.IsNotNullOrEmpty()))
                {
                    if (peakCluster.All(peak => peak.Apex != null))
                    {
                        IsobarCluster isoCluster = null;
                        try
                        {
                            isoCluster = new IsobarCluster(peakCluster, flashLfqEngine, results);
                            isobarClusters.Add(isoCluster);
                        }
                        catch (Exception ex)
                        {
                            if (exceptions != null)
                                exceptions.Add(ex);
                        }
                    }
                        
                }

            }

            return isobarClusters;
        }

        //TODO: SUmmary comment
        /// <summary>
        /// Add summary comment here
        /// </summary>
        /// <param name="orderedPeaks"></param>
        /// <returns></returns>
        internal static List<List<ChromatographicPeak>> ClusterPeaks(List<ChromatographicPeak> orderedPeaks)
        {
            List<List<ChromatographicPeak>> peakClusters = new List<List<ChromatographicPeak>>();
            ClusterPeaks(orderedPeaks, ref peakClusters, minIndex:0, maxIndex: orderedPeaks.Count - 1);
            return peakClusters;
        }

        private static void ClusterPeaks(
            List<ChromatographicPeak> orderedPeaks, 
            ref List<List<ChromatographicPeak>> peakClusters, 
            int minIndex, 
            int maxIndex)
        {
            // Maximum allowable time difference for two peaks to be grouped into the same cluster
            double deltaMax = 1;

            // Define the base case
            if (minIndex == maxIndex || maxIndex < 0)
                return;
            int leftIndex = minIndex;
            int rightIndex = maxIndex;

            int centerIndex = (maxIndex - minIndex + 1) / 2;

            List<ChromatographicPeak> peakCluster =
                new List<ChromatographicPeak>() { orderedPeaks[centerIndex] };
            int[] directions = { -1, 1 };

            foreach (int direction in directions)
            {
                double currentRt = orderedPeaks[centerIndex].ApexRetentionTime;
                for (int i = centerIndex + direction; minIndex <= i && i <= maxIndex; i += direction)
                {
                    double rtDelta = Math.Abs(orderedPeaks[i].ApexRetentionTime - currentRt);
                    if (rtDelta < deltaMax)
                    {
                        peakCluster.Add(orderedPeaks[i]);
                    }
                    else
                    {
                        if (direction < 0) leftIndex = i;
                        else rightIndex = i;

                        break;
                    }
                }
            }

            var peakIds = peakCluster
                .SelectMany(peak => peak.Identifications.Select(id => id.ModifiedSequence))
                .ToList();

            // Check to see if the cluster contains modified aromatic residues and/or multiple species
            if (peakIds.Any(seq => ModifiedAromaticResidueRegex.IsMatch(seq))
                || peakIds.Distinct().Count() > 1)
            {
                peakClusters.Add(peakCluster);
            }
            
            // Recurse left
            ClusterPeaks(orderedPeaks, ref peakClusters, minIndex, leftIndex);
            // Recurse right
            ClusterPeaks(orderedPeaks, ref peakClusters, rightIndex, maxIndex);
        }

        public static Regex ModifiedAromaticResidueRegex
        {
            get
            {
                if (_modifiedAromaticResidueRegex == null)
                    _modifiedAromaticResidueRegex = new Regex(@"[WYF]\[");
                return _modifiedAromaticResidueRegex;
            }
        }

        private static Regex _modifiedAromaticResidueRegex;
    }
}
