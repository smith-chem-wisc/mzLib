﻿using System.Collections.Generic;
using System.Linq;

namespace FlashLFQ
{
    public class FlashLFQResults
    {
        #region Public Fields

        public static List<RawFileInfo> rawFiles;
        public string log;
        public Dictionary<string, Peptide> peptideBaseSequences;
        public Dictionary<string, Peptide> peptideModifiedSequences;
        public Dictionary<string, ProteinGroup> proteinGroups;
        public Dictionary<RawFileInfo, List<ChromatographicPeak>> peaks;

        #endregion Public Fields

        #region Public Constructors

        public FlashLFQResults()
        {
            peptideBaseSequences = new Dictionary<string, Peptide>();
            peptideModifiedSequences = new Dictionary<string, Peptide>();
            proteinGroups = new Dictionary<string, ProteinGroup>();
            peaks = new Dictionary<RawFileInfo, List<ChromatographicPeak>>();
        }

        #endregion Public Constructors

        #region Public Methods

        public void CalculatePeptideResults(bool sumByBaseSequenceNotModifiedSequence)
        {
            foreach (var file in peaks)
            {
                // match peaks to sequences
                var sequenceToPeaksMatch = new Dictionary<string, HashSet<ChromatographicPeak>>();

                foreach (var peak in file.Value)
                {
                    foreach (var id in peak.identifyingScans)
                    {
                        string seq = id.BaseSequence;
                        if (!sumByBaseSequenceNotModifiedSequence)
                            seq = id.ModifiedSequence;

                        if (sequenceToPeaksMatch.TryGetValue(seq, out HashSet<ChromatographicPeak> featuresForThisBaseSeq))
                            featuresForThisBaseSeq.Add(peak);
                        else
                            sequenceToPeaksMatch.Add(seq, new HashSet<ChromatographicPeak>() { peak });
                    }
                }

                // calculate and store peptide intensities
                foreach (var sequence in sequenceToPeaksMatch)
                {
                    // calculate intensity and detection type for this peptide and file
                    double intensity = 0;
                    DetectionType detectionType = DetectionType.NotDetected;
                    var pgs = new HashSet<ProteinGroup>(sequence.Value.SelectMany(p => p.identifyingScans).SelectMany(v => v.proteinGroups));

                    if (sequence.Value.First().isMbrFeature)
                    {
                        intensity = sequence.Value.Max(p => p.intensity);
                        detectionType = DetectionType.MBR;
                    }
                    else
                    {
                        intensity = sequence.Value.Sum(p => p.intensity);
                        detectionType = DetectionType.MSMS;

                        if (intensity == 0)
                            detectionType = DetectionType.MSMSIdentifiedButNotQuantified;

                        if (sequence.Value.Max(p => p.NumIdentificationsByBaseSeq) > 1)
                        {
                            double ambigPeakIntensity = sequence.Value.Where(p => p.NumIdentificationsByBaseSeq > 1).Sum(v => v.intensity);

                            if ((ambigPeakIntensity / intensity) < 0.3)
                                intensity = sequence.Value.Select(p => (p.intensity / p.NumIdentificationsByBaseSeq)).Sum();
                            else
                            {
                                detectionType = DetectionType.MSMSAmbiguousPeakfinding;
                                intensity = 0;
                            }
                        }
                    }

                    // store intensity and detection type for this peptide and file
                    if (sumByBaseSequenceNotModifiedSequence)
                    {
                        peptideBaseSequences[sequence.Key].intensities[file.Key] = intensity;
                        peptideBaseSequences[sequence.Key].detectionTypes[file.Key] = detectionType;
                        peptideBaseSequences[sequence.Key].proteinGroups = pgs;
                    }
                    else
                    {
                        peptideModifiedSequences[sequence.Key].intensities[file.Key] = intensity;
                        peptideModifiedSequences[sequence.Key].detectionTypes[file.Key] = detectionType;
                        peptideModifiedSequences[sequence.Key].proteinGroups = pgs;
                    }
                }
            }
        }

        public void CalculateProteinResults()
        {
            // get all peptides without ambiguous peaks
            var allFeatures = peaks.Values.SelectMany(p => p).ToList();
            var allAmbiguousFeatures = allFeatures.Where(p => p.NumIdentificationsByBaseSeq > 1).ToList();
            var ambiguousFeatureSeqs = new HashSet<string>(allAmbiguousFeatures.SelectMany(p => p.identifyingScans.Select(v => v.BaseSequence)));

            foreach (var feature in allFeatures)
                if (ambiguousFeatureSeqs.Contains(feature.identifyingScans.First().BaseSequence))
                    allAmbiguousFeatures.Add(feature);

            var allUnambiguousFeatures = allFeatures.Except(allAmbiguousFeatures).ToList();

            // match these peaks to proteins
            var proteinsWithFeatures = new Dictionary<ProteinGroup, List<ChromatographicPeak>>();
            foreach (var feature in allUnambiguousFeatures)
            {
                foreach (var proteinGroup in feature.identifyingScans.First().proteinGroups)
                {
                    if (proteinsWithFeatures.TryGetValue(proteinGroup, out List<ChromatographicPeak> featuresForThisProtein))
                        featuresForThisProtein.Add(feature);
                    else
                        proteinsWithFeatures.Add(proteinGroup, new List<ChromatographicPeak> { feature });
                }
            }

            // calculate protein's intensity for this file
            foreach (var proteinGroupsFeatures in proteinsWithFeatures)
            {
                foreach (var feature in proteinGroupsFeatures.Value)
                {
                    int numProteinGroupsClaimingThisFeature = feature.identifyingScans.SelectMany(p => p.proteinGroups).Distinct().Count();
                    proteinGroupsFeatures.Key.intensities[feature.rawFileInfo] += (feature.intensity / numProteinGroupsClaimingThisFeature);
                }
            }
        }

        #endregion Public Methods
    }
}