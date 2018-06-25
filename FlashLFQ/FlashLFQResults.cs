using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace FlashLFQ
{
    public class FlashLFQResults
    {
        public readonly List<SpectraFileInfo> spectraFiles;
        public readonly Dictionary<string, Peptide> peptideBaseSequences;
        public readonly Dictionary<string, Peptide> peptideModifiedSequences;
        public readonly Dictionary<string, ProteinGroup> proteinGroups;
        public readonly Dictionary<SpectraFileInfo, List<ChromatographicPeak>> peaks;

        public FlashLFQResults(List<SpectraFileInfo> rawFiles)
        {
            this.spectraFiles = rawFiles;
            peptideBaseSequences = new Dictionary<string, Peptide>();
            peptideModifiedSequences = new Dictionary<string, Peptide>();
            proteinGroups = new Dictionary<string, ProteinGroup>();
            peaks = new Dictionary<SpectraFileInfo, List<ChromatographicPeak>>();
        }

        public void CalculatePeptideResults(bool sumByBaseSequenceNotModifiedSequence)
        {
            foreach (var file in peaks)
            {
                // match peaks to sequences
                var sequenceToPeaksMatch = new Dictionary<string, HashSet<ChromatographicPeak>>();

                foreach (var peak in file.Value)
                {
                    foreach (var id in peak.identifications)
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
                    var pgs = new HashSet<ProteinGroup>(sequence.Value.SelectMany(p => p.identifications).SelectMany(v => v.proteinGroups));

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
                        if (!peptideBaseSequences.ContainsKey(sequence.Key))
                            peptideBaseSequences.Add(sequence.Key, new Peptide(sequence.Key));

                        peptideBaseSequences[sequence.Key].SetIntensity(file.Key, intensity);
                        peptideBaseSequences[sequence.Key].SetDetectionType(file.Key, detectionType);
                        peptideBaseSequences[sequence.Key].proteinGroups = pgs;
                    }
                    else
                    {
                        if (!peptideModifiedSequences.ContainsKey(sequence.Key))
                            peptideModifiedSequences.Add(sequence.Key, new Peptide(sequence.Key));

                        peptideModifiedSequences[sequence.Key].SetIntensity(file.Key, intensity);
                        peptideModifiedSequences[sequence.Key].SetDetectionType(file.Key, detectionType);
                        peptideModifiedSequences[sequence.Key].proteinGroups = pgs;
                    }
                }
            }
        }

        public void CalculateProteinResults()
        {
            // get all peptides without ambiguous peaks
            var allFeatures = peaks.Values.SelectMany(p => p).ToList();

            foreach (var peak in allFeatures)
                foreach (var id in peak.identifications)
                    foreach (var proteinGroup in id.proteinGroups)
                        if (!proteinGroups.ContainsKey(proteinGroup.ProteinGroupName))
                            proteinGroups.Add(proteinGroup.ProteinGroupName, proteinGroup);

            var allAmbiguousFeatures = allFeatures.Where(p => p.NumIdentificationsByBaseSeq > 1).ToList();
            var ambiguousFeatureSeqs = new HashSet<string>(allAmbiguousFeatures.SelectMany(p => p.identifications.Select(v => v.BaseSequence)));

            foreach (var feature in allFeatures)
                if (ambiguousFeatureSeqs.Contains(feature.identifications.First().BaseSequence))
                    allAmbiguousFeatures.Add(feature);

            var allUnambiguousFeatures = allFeatures.Except(allAmbiguousFeatures).ToList();

            // match these peaks to proteins
            var proteinsWithFeatures = new Dictionary<ProteinGroup, List<ChromatographicPeak>>();
            foreach (var feature in allUnambiguousFeatures)
            {
                foreach (var proteinGroup in feature.identifications.First().proteinGroups)
                {
                    if (proteinsWithFeatures.TryGetValue(proteinGroup, out List<ChromatographicPeak> featuresForThisProtein))
                        featuresForThisProtein.Add(feature);
                    else
                        proteinsWithFeatures.Add(proteinGroup, new List<ChromatographicPeak> { feature });
                }
            }

            // calculate protein's intensity for this file
            foreach (var proteinGroup in proteinsWithFeatures)
            {
                int topNFeatures = 3;
                Dictionary<SpectraFileInfo, List<double>> fileToPepIntensities = new Dictionary<SpectraFileInfo, List<double>>();

                foreach (var feature in proteinGroup.Value)
                {
                    int numProteinGroupsClaimingThisFeature = feature.identifications.SelectMany(p => p.proteinGroups).Distinct().Count();

                    if (fileToPepIntensities.TryGetValue(feature.rawFileInfo, out var featureIntensitiesForThisProtein))
                    {
                        fileToPepIntensities[feature.rawFileInfo].Add(feature.intensity / numProteinGroupsClaimingThisFeature);
                    }
                    else
                    {
                        fileToPepIntensities.Add(feature.rawFileInfo, new List<double> { feature.intensity / numProteinGroupsClaimingThisFeature });
                    }
                }

                foreach (var file in fileToPepIntensities)
                {
                    // need to observe at least one MS2-identified peptide for a protein in a file. if they're all MBR-identified, the protein
                    // intensity is zero. this is to prevent false-positives but will reduce the number of quantified proteins
                    if (proteinGroup.Value.Any(p => !p.isMbrFeature))
                    {
                        proteinGroup.Key.intensities[file.Key] = file.Value.OrderByDescending(p => p).Take(topNFeatures).Sum();
                    }
                    else
                    {
                        proteinGroup.Key.intensities[file.Key] = 0;
                    }
                }
            }
        }

        public void WriteResults(string peaksOutputPath, string modPeptideOutputPath, string baseSeqPeptideOutputPath, string proteinOutputPath)
        {
            if (peaksOutputPath != null)
            {
                using (StreamWriter output = new StreamWriter(peaksOutputPath))
                {
                    output.WriteLine(ChromatographicPeak.TabSeparatedHeader);

                    foreach (var peak in peaks.SelectMany(p => p.Value))
                    {
                        output.WriteLine(peak.ToString());
                    }
                }
            }

            if (modPeptideOutputPath != null)
            {
                using (StreamWriter output = new StreamWriter(modPeptideOutputPath))
                {
                    output.WriteLine(Peptide.TabSeparatedHeader(spectraFiles));

                    foreach (var peptide in peptideModifiedSequences.OrderBy(p => p.Key))
                    {
                        output.WriteLine(peptide.Value.ToString(spectraFiles));
                    }
                }
            }

            if (baseSeqPeptideOutputPath != null)
            {
                using (StreamWriter output = new StreamWriter(baseSeqPeptideOutputPath))
                {
                    output.WriteLine(Peptide.TabSeparatedHeader(spectraFiles));

                    foreach (var peptide in peptideBaseSequences.OrderBy(p => p.Key))
                    {
                        output.WriteLine(peptide.Value.ToString(spectraFiles));
                    }
                }
            }

            if (proteinOutputPath != null)
            {
                using (StreamWriter output = new StreamWriter(proteinOutputPath))
                {
                    output.WriteLine(ProteinGroup.TabSeparatedHeader(spectraFiles));

                    foreach (var protein in proteinGroups.OrderBy(p => p.Key))
                    {
                        output.WriteLine(protein.Value.ToString(spectraFiles));
                    }
                }
            }
        }
    }
}