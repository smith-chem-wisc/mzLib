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
        public Dictionary<string, int> SequenceRegionDictionary { get; private set; }

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

            var test = Identifications.Select(id => id.ModifiedSequence).ToList();
            if (test.Any(seq => seq.Equals("F[CF3:CF3 on F]KDLGEEHFK")))
            {
                int placeholder = 0;
            }

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
                if (peaks.Count <= 5)
                {
                    throw new Exception("XIC too short");
                }
                double rtAdjustment = XicProcessing.AlignPeaks(refPeaks, peaks);
                // If the rtAdjustment is some huge number, it probably went wrong
                // so, set to 0. Proceed from there
                rtAdjustment = Math.Abs(rtAdjustment) > (RetentionTimeRange.Mean / 10) + 0.1 
                    ? 0
                    : rtAdjustment;
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
                double startTime = Math.Min(RetentionTimeRange.Minimum, kvp.Value.First().RetentionTime);

                if (kvp.Value.Count <= 1)
                {
                    peakRegions.Add(regionNumber, RetentionTimeRange);
                    Xics[kvp.Key].SetPeakRegions(peakRegions);
                    continue;
                }

                bool maxFound = false;
                double lastMaxRt = startTime;
                if (kvp.Value.First().ExtremumType == ExtremumType.Maximum)
                {
                    maxFound = true;
                }

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
                        if (!maxFound || extrema == kvp.Value.First()) // Intended reference comparison. If there are two maxima at the start, the extrema check avoids a 0-width region
                        {
                            maxFound = true;
                            lastMaxRt = extrema.RetentionTime;
                        }
                        else
                        {
                            double pseudoMinimum = (lastMaxRt + extrema.RetentionTime) / 2.0;
                            peakRegions.Add(regionNumber++, new DoubleRange(startTime, pseudoMinimum));
                            startTime = pseudoMinimum;
                            lastMaxRt = extrema.RetentionTime;
                        }
                    }
                }
                // If loop completes on a a maximum, add a final region from last start time - end
                if (maxFound && RetentionTimeRange.Maximum - startTime > 0)
                {
                    peakRegions.Add(regionNumber, new DoubleRange(startTime, RetentionTimeRange.Maximum));
                }

                if(peakRegions.Count == 0) peakRegions.Add(1, RetentionTimeRange);
                Xics[kvp.Key].SetPeakRegions(peakRegions);
            }
        }


        public void ReassignPeakIDs()
        {
            SequenceRegionDictionary = new Dictionary<string, int>();

            _regionPeakDictionary = new Dictionary<int, List<ChromatographicPeak>>();
            foreach(int region in Xics[ReferenceFile].PeakRegions.Keys)
                _regionPeakDictionary.Add(region, new List<ChromatographicPeak>());

            // Populate the _regionPeakDictionary. Each consensus region is linked to peaks
            // with their retention time apex in that region.
            // Also sets the region property for each peak
            foreach (ChromatographicPeak peak in Peaks)
            {
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

            // TODO: Make this throw an exception
            if (_regionPeakDictionary.Count == 0)
                return;


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
                    // This algorithm assigns peaks such that ambiguity is minimized
                    List<int> regions = seqToRegionsKvp.Value;

                    // For each region (key), contains the number of IDs made in that region
                    var regionCountDictionary = regions
                        .GroupBy(i => i)
                        .ToDictionary(
                            keySelector: group => group.Key,
                            elementSelector: group => group.Count());

                    int maxNumberOfPeaks = regionCountDictionary.Max(kvp => kvp.Value);
                    var bestRegionCandidates = regionCountDictionary
                        .Where(kvp => kvp.Value == maxNumberOfPeaks)
                        .ToList();
                    int bestRegion = bestRegionCandidates.First().Key;

                    // if multiple regions have the same number of associated peaks, we select the 
                    // region with the greatest mean peak intensity
                    if (bestRegionCandidates.Count() > 1)
                    {
                        double highestMeanRegionIntensity = 0;
                        foreach (var kvp in regionCountDictionary
                                     .Where(kvp => kvp.Value == maxNumberOfPeaks))
                        {
                            double meanRegionIntensity =
                                RegionPeakDictionary[kvp.Key].Average(peak => peak.Intensity);

                            if (meanRegionIntensity > highestMeanRegionIntensity)
                            {
                                bestRegion = kvp.Key;
                                highestMeanRegionIntensity = meanRegionIntensity;
                            }
                        }
                    }

                    foreach (string seq in _regionPeakDictionary[bestRegion]
                                 .SelectMany(peak => peak.Identifications)
                                 .SelectMany(id => id.ModifiedSequence.Split('|'))
                                 .Distinct())
                    {
                        //TODO: Need to weight unambiguous identifications higher
                        SequenceRegionDictionary.TryAdd(seq, bestRegion);
                    }

                    foreach (int region in seqToRegionsKvp.Value)
                    {
                        if (region == bestRegion || Math.Abs(region - bestRegion) > 1) 
                        {
                            // If region == bestRegion, these were accurately assigned already
                            // If region is more than 2 regions away from the best region, then we don't want to merge them
                            // That check should maybe be based on RT differences instead of region number
                            continue; 
                        }

                        // Peaks are added to the dictionary within the foreach loop, so it's important to stash the 
                        // peaks before entering the loop
                        foreach (var peak in _regionPeakDictionary[region])
                        {
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

                                    _flashLfqEngine.SwapPeakIndexingEngine(peak.SpectraFileInfo);

                                    bestPeak = _flashLfqEngine.GetChromatographicPeak(
                                        id, 
                                        peak.SpectraFileInfo, 
                                        isMbrPeak: peak.SpectraFileInfo != id.FileInfo, // TODO: Add an alternate classification system for ambiguous peaks. MBR isn't totally accurate 
                                        rtApex, 
                                        xic.PeakRegions[bestRegion]);

                                    if (bestPeak == null) continue; //TODO: Figure out why this is happening!
                                    _regionPeakDictionary[bestRegion].Add(bestPeak);

                                    //TODO: Clean this up
                                    if(_results != null)
                                        _results.Peaks[peak.SpectraFileInfo].Add(bestPeak);
                                }
                                bestPeak.AddIdentification(id);
                                bestPeak.Region = bestRegion;

                                // I think this is unnecessary and leads to more ambiguous assignments than is accurate

                                //foreach(string seq in id.ModifiedSequence.Split('|'))
                                //{
                                //    SequenceRegionDictionary.TryAdd(seq, bestRegion);
                                //}
                            }

                        }
                    }
                }
            }
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

            // Maximum allowable time difference for two peaks to be grouped into the same cluster
            // This may need to be adjusted to be larger
            //double meanRt = orderedPeaks.Average(peak => peak.ApexRetentionTime);
            //double deltaMax = 1 + meanRt / 10.0;
            double deltaMax = 1.0;
            ClusterPeaks(orderedPeaks, ref peakClusters, deltaMax, minIndex:0, maxIndex: orderedPeaks.Count - 1);
            MergePeakClusters(ref peakClusters);
            return peakClusters;
        }

        private static void ClusterPeaks(
            List<ChromatographicPeak> orderedPeaks, 
            ref List<List<ChromatographicPeak>> peakClusters, 
            double deltaMax,
            int minIndex, 
            int maxIndex)
        {
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
            ClusterPeaks(orderedPeaks, ref peakClusters, deltaMax, minIndex, leftIndex);
            // Recurse right
            ClusterPeaks(orderedPeaks, ref peakClusters, deltaMax, rightIndex, maxIndex);
        }

        internal static void MergePeakClusters(ref List<List<ChromatographicPeak>> peakClusters)
        {
            peakClusters = peakClusters
                .OrderBy(cluster => cluster.Average(peak => peak.ApexRetentionTime))
                .ToList();
            List<HashSet<string>> seqHashSets = peakClusters
                .Select(cluster => cluster
                    .SelectMany(peak => peak.Identifications)
                    .SelectMany(id => id.ModifiedSequence.Split('|'))
                    .ToHashSet())
                .ToList();

            // This algorithm is O(n^2) but I don't think it can be improved
            // Go through every peak cluster and associated set of sequences
            // If there is an overlap (i.e., same peptide id'd in separate clusters)
            // Then we merge the overlapping clusters + any clusters between them
            int i = 0;
            while (i < seqHashSets.Count - 1)
            {
                int j = i + 1;
                while (j < seqHashSets.Count)
                {
                    if (seqHashSets[i].Overlaps(seqHashSets[j]))
                    {
                        HashSet<string> newSeqSet = new();
                        List<ChromatographicPeak> newPeakCluster = new();
                        for (int k = i; k <= j; k++)
                        {
                            newSeqSet.UnionWith(seqHashSets[k]);
                            newPeakCluster.AddRange(peakClusters[k]);
                        }
                        // Remove everything between i and j, inclusive
                        // Then add new element to start of list
                        // Removal and insertion isn't super efficient for regular lists,
                        // TODO: Optimize to use stack or LinkedList
                        seqHashSets.RemoveRange(i, j-i+1);
                        seqHashSets.Insert(i, newSeqSet);
                        peakClusters.RemoveRange(i, j-i+1);
                        peakClusters.Insert(i, newPeakCluster);

                        i--; // decrement one so that when it is incremented in the outer loop, we end up at the same index we started at
                        break;
                    }

                    j++;
                }
                i++; // i is only incremented in the outer loop. Incrementing when joining elements prevents recursive joining
            }
        }

        public static Regex ModifiedAromaticResidueRegex
        {
            get
            {
                if (_modifiedAromaticResidueRegex == null)
                    _modifiedAromaticResidueRegex = new Regex(@"[WYFH]\[");
                return _modifiedAromaticResidueRegex;
            }
        }

        private static Regex _modifiedAromaticResidueRegex;
    }
}
