using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace FlashLFQ
{
    public class FlashLfqResults
    {
        public readonly List<SpectraFileInfo> SpectraFiles;
        public readonly Dictionary<string, Peptide> PeptideBaseSequences;
        public readonly Dictionary<string, Peptide> PeptideModifiedSequences;
        public readonly Dictionary<string, ProteinGroup> ProteinGroups;
        public readonly Dictionary<SpectraFileInfo, List<ChromatographicPeak>> Peaks;

        public FlashLfqResults(List<SpectraFileInfo> rawFiles)
        {
            SpectraFiles = rawFiles;
            PeptideBaseSequences = new Dictionary<string, Peptide>();
            PeptideModifiedSequences = new Dictionary<string, Peptide>();
            ProteinGroups = new Dictionary<string, ProteinGroup>();
            Peaks = new Dictionary<SpectraFileInfo, List<ChromatographicPeak>>();
        }

        public void MergeResultsWith(FlashLfqResults mergeFrom)
        {
            this.SpectraFiles.AddRange(mergeFrom.SpectraFiles);

            foreach (var pep in mergeFrom.PeptideBaseSequences)
            {
                if (this.PeptideBaseSequences.TryGetValue(pep.Key, out var peptide))
                {
                    Peptide mergeFromPep = pep.Value;
                    Peptide mergeToPep = peptide;

                    foreach (SpectraFileInfo file in mergeFrom.SpectraFiles)
                    {
                        mergeToPep.SetIntensity(file, mergeFromPep.GetIntensity(file));
                        mergeToPep.SetDetectionType(file, mergeFromPep.GetDetectionType(file));
                    }
                }
                else
                {
                    this.PeptideBaseSequences.Add(pep.Key, pep.Value);
                }
            }
            foreach (var pep in mergeFrom.PeptideModifiedSequences)
            {
                if (this.PeptideModifiedSequences.TryGetValue(pep.Key, out var peptide))
                {
                    Peptide mergeFromPep = pep.Value;
                    Peptide mergeToPep = peptide;

                    foreach (SpectraFileInfo file in mergeFrom.SpectraFiles)
                    {
                        mergeToPep.SetIntensity(file, mergeFromPep.GetIntensity(file));
                        mergeToPep.SetDetectionType(file, mergeFromPep.GetDetectionType(file));
                    }
                }
                else
                {
                    this.PeptideModifiedSequences.Add(pep.Key, pep.Value);
                }
            }
            foreach (var pg in mergeFrom.ProteinGroups)
            {
                if (this.ProteinGroups.TryGetValue(pg.Key, out var proteinGroup))
                {
                    ProteinGroup mergeFromPg = pg.Value;
                    ProteinGroup mergeToPg = proteinGroup;

                    foreach (SpectraFileInfo file in mergeFrom.SpectraFiles)
                    {
                        mergeToPg.SetIntensity(file, mergeFromPg.GetIntensity(file));
                    }
                }
                else
                {
                    this.ProteinGroups.Add(pg.Key, pg.Value);
                }
            }
            foreach (var fromPeaks in mergeFrom.Peaks)
            {
                if (this.Peaks.TryGetValue(fromPeaks.Key, out var toPeaks))
                {
                    toPeaks.AddRange(fromPeaks.Value);
                }
                else
                {
                    this.Peaks.Add(fromPeaks.Key, fromPeaks.Value);
                }
            }
        }

        public void CalculatePeptideResults(bool sumByBaseSequenceNotModifiedSequence)
        {
            foreach (var file in Peaks)
            {
                // match peaks to sequences
                var sequenceToPeaksMatch = new Dictionary<string, HashSet<ChromatographicPeak>>();

                foreach (var peak in file.Value)
                {
                    foreach (var id in peak.Identifications)
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
                    var pgs = new HashSet<ProteinGroup>(sequence.Value.SelectMany(p => p.Identifications).SelectMany(v => v.proteinGroups));

                    if (sequence.Value.First().IsMbrFeature)
                    {
                        intensity = sequence.Value.Max(p => p.Intensity);
                        detectionType = DetectionType.MBR;
                    }
                    else
                    {
                        intensity = sequence.Value.Sum(p => p.Intensity);
                        detectionType = DetectionType.MSMS;

                        if (intensity == 0)
                            detectionType = DetectionType.MSMSIdentifiedButNotQuantified;

                        if (sequence.Value.Max(p => p.NumIdentificationsByBaseSeq) > 1)
                        {
                            double ambigPeakIntensity = sequence.Value.Where(p => p.NumIdentificationsByBaseSeq > 1).Sum(v => v.Intensity);

                            if ((ambigPeakIntensity / intensity) < 0.3)
                                intensity = sequence.Value.Select(p => (p.Intensity / p.NumIdentificationsByBaseSeq)).Sum();
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
                        if (!PeptideBaseSequences.ContainsKey(sequence.Key))
                            PeptideBaseSequences.Add(sequence.Key, new Peptide(sequence.Key));

                        PeptideBaseSequences[sequence.Key].SetIntensity(file.Key, intensity);
                        PeptideBaseSequences[sequence.Key].SetDetectionType(file.Key, detectionType);
                        PeptideBaseSequences[sequence.Key].proteinGroups = pgs;
                    }
                    else
                    {
                        if (!PeptideModifiedSequences.ContainsKey(sequence.Key))
                            PeptideModifiedSequences.Add(sequence.Key, new Peptide(sequence.Key));

                        PeptideModifiedSequences[sequence.Key].SetIntensity(file.Key, intensity);
                        PeptideModifiedSequences[sequence.Key].SetDetectionType(file.Key, detectionType);
                        PeptideModifiedSequences[sequence.Key].proteinGroups = pgs;
                    }
                }
            }
        }

        public void CalculateProteinResultsTop3()
        {
            // get all peptides without ambiguous peaks
            var allFeatures = Peaks.Values.SelectMany(p => p).ToList();

            foreach (var peak in allFeatures)
            {
                foreach (var id in peak.Identifications)
                {
                    foreach (var proteinGroup in id.proteinGroups)
                    {
                        if (!ProteinGroups.ContainsKey(proteinGroup.ProteinGroupName))
                        {
                            ProteinGroups.Add(proteinGroup.ProteinGroupName, proteinGroup);
                        }
                    }
                }
            }

            var allAmbiguousFeatures = allFeatures.Where(p => p.NumIdentificationsByBaseSeq > 1).ToList();
            var ambiguousFeatureSeqs = new HashSet<string>(allAmbiguousFeatures.SelectMany(p => p.Identifications.Select(v => v.BaseSequence)));

            foreach (var feature in allFeatures)
            {
                if (ambiguousFeatureSeqs.Contains(feature.Identifications.First().BaseSequence))
                {
                    allAmbiguousFeatures.Add(feature);
                }
            }

            var allUnambiguousFeatures = allFeatures.Except(allAmbiguousFeatures).ToList();

            // match these peaks to proteins
            var proteinsWithFeatures = new Dictionary<ProteinGroup, List<ChromatographicPeak>>();
            foreach (var feature in allUnambiguousFeatures)
            {
                // only use unmodified peptides for protein quant
                // ONLY WORKS WITH METAMORPHEUS OUTPUT
                //if (feature.Identifications.First().BaseSequence != feature.Identifications.First().ModifiedSequence)
                //{
                //    continue;
                //}

                foreach (var proteinGroup in feature.Identifications.First().proteinGroups)
                {
                    if (proteinsWithFeatures.TryGetValue(proteinGroup, out List<ChromatographicPeak> featuresForThisProtein))
                    {
                        featuresForThisProtein.Add(feature);
                    }
                    else
                    {
                        proteinsWithFeatures.Add(proteinGroup, new List<ChromatographicPeak> { feature });
                    }
                }
            }

            // calculate protein's intensity for this file
            foreach (var proteinGroup in proteinsWithFeatures)
            {
                int topNFeatures = 3;
                Dictionary<SpectraFileInfo, List<double>> fileToPepIntensities = new Dictionary<SpectraFileInfo, List<double>>();

                foreach (var feature in proteinGroup.Value)
                {
                    int numProteinGroupsClaimingThisFeature = feature.Identifications.SelectMany(p => p.proteinGroups).Distinct().Count();

                    if (fileToPepIntensities.TryGetValue(feature.RawFileInfo, out var featureIntensitiesForThisProtein))
                    {
                        featureIntensitiesForThisProtein.Add(feature.Intensity / numProteinGroupsClaimingThisFeature);
                    }
                    else
                    {
                        fileToPepIntensities.Add(feature.RawFileInfo, new List<double> { feature.Intensity / numProteinGroupsClaimingThisFeature });
                    }
                }

                foreach (var file in fileToPepIntensities)
                {
                    proteinGroup.Key.SetIntensity(file.Key, file.Value.OrderByDescending(p => p).Take(topNFeatures).Sum());
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

                    foreach (var peak in Peaks.SelectMany(p => p.Value))
                    {
                        output.WriteLine(peak.ToString());
                    }
                }
            }

            if (modPeptideOutputPath != null)
            {
                using (StreamWriter output = new StreamWriter(modPeptideOutputPath))
                {
                    output.WriteLine(Peptide.TabSeparatedHeader(SpectraFiles));

                    foreach (var peptide in PeptideModifiedSequences.OrderBy(p => p.Key))
                    {
                        output.WriteLine(peptide.Value.ToString(SpectraFiles));
                    }
                }
            }

            if (baseSeqPeptideOutputPath != null)
            {
                using (StreamWriter output = new StreamWriter(baseSeqPeptideOutputPath))
                {
                    output.WriteLine(Peptide.TabSeparatedHeader(SpectraFiles));

                    foreach (var peptide in PeptideBaseSequences.OrderBy(p => p.Key))
                    {
                        output.WriteLine(peptide.Value.ToString(SpectraFiles));
                    }
                }
            }

            if (proteinOutputPath != null)
            {
                using (StreamWriter output = new StreamWriter(proteinOutputPath))
                {
                    output.WriteLine(ProteinGroup.TabSeparatedHeader(SpectraFiles));

                    foreach (var protein in ProteinGroups.OrderBy(p => p.Key))
                    {
                        output.WriteLine(protein.Value.ToString(SpectraFiles));
                    }
                }
            }
        }
    }
}