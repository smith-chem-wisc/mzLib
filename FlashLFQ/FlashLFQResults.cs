using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace FlashLFQ
{
    public class FlashLfqResults
    {
        public readonly List<SpectraFileInfo> SpectraFiles;
        public readonly Dictionary<string, Peptide> PeptideModifiedSequences;
        public readonly Dictionary<string, ProteinGroup> ProteinGroups;
        public readonly Dictionary<SpectraFileInfo, List<ChromatographicPeak>> Peaks;

        public FlashLfqResults(List<SpectraFileInfo> spectraFiles, List<Identification> identifications)
        {
            SpectraFiles = spectraFiles;
            PeptideModifiedSequences = new Dictionary<string, Peptide>();
            ProteinGroups = new Dictionary<string, ProteinGroup>();
            Peaks = new Dictionary<SpectraFileInfo, List<ChromatographicPeak>>();

            foreach (SpectraFileInfo file in spectraFiles)
            {
                Peaks.Add(file, new List<ChromatographicPeak>());
            }

            foreach (Identification id in identifications)
            {
                if (!PeptideModifiedSequences.TryGetValue(id.ModifiedSequence, out Peptide peptide))
                {
                    PeptideModifiedSequences.Add(id.ModifiedSequence, new Peptide(id.ModifiedSequence, id.BaseSequence, id.UseForProteinQuant, id.ProteinGroups));
                }
                else
                {
                    foreach (ProteinGroup pg in id.ProteinGroups)
                    {
                        peptide.ProteinGroups.Add(pg);
                    }
                }

                foreach (ProteinGroup proteinGroup in id.ProteinGroups)
                {
                    if (!ProteinGroups.ContainsKey(proteinGroup.ProteinGroupName))
                    {
                        ProteinGroups.Add(proteinGroup.ProteinGroupName, proteinGroup);
                    }
                }
            }
        }

        public void MergeResultsWith(FlashLfqResults mergeFrom)
        {
            this.SpectraFiles.AddRange(mergeFrom.SpectraFiles);

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

        public void CalculatePeptideResults()
        {
            foreach (var sequence in PeptideModifiedSequences)
            {
                foreach (SpectraFileInfo file in SpectraFiles)
                {
                    sequence.Value.SetDetectionType(file, DetectionType.NotDetected);
                    sequence.Value.SetIntensity(file, 0);
                }
            }

            foreach (var filePeaks in Peaks)
            {
                var groupedPeaks = filePeaks.Value.Where(p => p.NumIdentificationsByFullSeq == 1)
                    .GroupBy(p => p.Identifications.First().ModifiedSequence).ToList();

                foreach (var sequenceWithPeaks in groupedPeaks)
                {
                    string sequence = sequenceWithPeaks.Key;
                    double intensity = sequenceWithPeaks.Max(p => p.Intensity);
                    ChromatographicPeak bestPeak = sequenceWithPeaks.First(p => p.Intensity == intensity);
                    DetectionType detectionType;

                    if (bestPeak.IsMbrPeak && intensity > 0)
                    {
                        detectionType = DetectionType.MBR;
                    }
                    else if (!bestPeak.IsMbrPeak && intensity > 0)
                    {
                        detectionType = DetectionType.MSMS;
                    }
                    else if (!bestPeak.IsMbrPeak && intensity == 0)
                    {
                        detectionType = DetectionType.MSMSIdentifiedButNotQuantified;
                    }
                    else
                    {
                        detectionType = DetectionType.NotDetected;
                    }

                    PeptideModifiedSequences[sequence].SetIntensity(filePeaks.Key, intensity);
                    PeptideModifiedSequences[sequence].SetDetectionType(filePeaks.Key, detectionType);
                }

                // report ambiguous quantification
                var ambiguousPeaks = filePeaks.Value.Where(p => p.NumIdentificationsByFullSeq > 1).ToList();
                foreach (ChromatographicPeak ambiguousPeak in ambiguousPeaks)
                {
                    foreach (Identification id in ambiguousPeak.Identifications)
                    {
                        string sequence = id.ModifiedSequence;

                        double alreadyRecordedIntensity = PeptideModifiedSequences[sequence].GetIntensity(filePeaks.Key);
                        double fractionAmbiguous = ambiguousPeak.Intensity / (alreadyRecordedIntensity + ambiguousPeak.Intensity);

                        if (fractionAmbiguous > 0.3)
                        {
                            PeptideModifiedSequences[sequence].SetIntensity(filePeaks.Key, 0);
                            PeptideModifiedSequences[sequence].SetDetectionType(filePeaks.Key, DetectionType.MSMSAmbiguousPeakfinding);
                        }
                    }
                }
            }

            HandleAmbiguityInFractions();
        }

        private void HandleAmbiguityInFractions()
        {
            // handle ambiguous peaks in fractionated data
            // if the largest fraction intensity is ambiguous, zero out the other fractions for the sample
            var sampleGroups = SpectraFiles.GroupBy(p => p.Condition);
            foreach (var sampleGroup in sampleGroups)
            {
                var samples = sampleGroup.Select(p => p).GroupBy(p => p.BiologicalReplicate);

                foreach (var sample in samples)
                {
                    // skip unfractionated samples
                    if (sample.Select(p => p.Fraction).Distinct().Count() == 1)
                    {
                        continue;
                    }

                    var peaksForEachSequence = new Dictionary<(SpectraFileInfo, string), List<ChromatographicPeak>>();

                    foreach (SpectraFileInfo file in sample)
                    {
                        foreach (ChromatographicPeak peak in Peaks[file])
                        {
                            foreach (Identification id in peak.Identifications)
                            {
                                if (peaksForEachSequence.TryGetValue((file, id.ModifiedSequence), out var peaks))
                                {
                                    peaks.Add(peak);
                                }
                                else
                                {
                                    peaksForEachSequence.Add((file, id.ModifiedSequence), new List<ChromatographicPeak> { peak });
                                }
                            }
                        }
                    }

                    var peptides = PeptideModifiedSequences.Values.ToList();
                    List<(double, DetectionType)> fractionIntensitiesWithDetectionTypes = new List<(double, DetectionType)>();
                    foreach (var peptide in peptides)
                    {
                        fractionIntensitiesWithDetectionTypes.Clear();
                        bool ambiguityObservedInSample = false;

                        foreach (SpectraFileInfo file in sample)
                        {
                            double fractionIntensity = peptide.GetIntensity(file);
                            DetectionType detectionType = peptide.GetDetectionType(file);

                            if (detectionType == DetectionType.MSMSAmbiguousPeakfinding)
                            {
                                ambiguityObservedInSample = true;

                                fractionIntensity = peaksForEachSequence[(file, peptide.Sequence)].Sum(p => p.Intensity);
                            }

                            fractionIntensitiesWithDetectionTypes.Add((fractionIntensity, detectionType));
                        }

                        if (ambiguityObservedInSample)
                        {
                            (double, DetectionType) highestIntensity = fractionIntensitiesWithDetectionTypes.OrderByDescending(p => p.Item1).First();

                            DetectionType highestFractionIntensityDetectionType = highestIntensity.Item2;
                            if (highestFractionIntensityDetectionType == DetectionType.MSMSAmbiguousPeakfinding)
                            {
                                // highest fraction intensity is ambiguous - zero out the other fractions
                                foreach (SpectraFileInfo file in sample)
                                {
                                    peptide.SetIntensity(file, 0);
                                }
                            }
                        }
                    }
                }
            }
        }

        public void CalculateProteinResultsTop3(bool useSharedPeptides)
        {
            foreach (var proteinGroup in ProteinGroups)
            {
                foreach (SpectraFileInfo file in SpectraFiles)
                {
                    proteinGroup.Value.SetIntensity(file, 0);
                }
            }

            int topNPeaks = 3;

            List<Peptide> peptides = PeptideModifiedSequences.Values.Where(p => p.UnambiguousPeptideQuant()).ToList();
            Dictionary<ProteinGroup, List<Peptide>> proteinGroupToPeptides = new Dictionary<ProteinGroup, List<Peptide>>();

            foreach (Peptide peptide in peptides)
            {
                if (!peptide.UseForProteinQuant)
                {
                    continue;
                }

                foreach (ProteinGroup pg in peptide.ProteinGroups)
                {
                    if (proteinGroupToPeptides.TryGetValue(pg, out var peptidesForThisProtein))
                    {
                        peptidesForThisProtein.Add(peptide);
                    }
                    else
                    {
                        proteinGroupToPeptides.Add(pg, new List<Peptide> { peptide });
                    }
                }
            }

            foreach (ProteinGroup pg in ProteinGroups.Values)
            {
                if (proteinGroupToPeptides.TryGetValue(pg, out var peptidesForThisProtein))
                {
                    foreach (SpectraFileInfo file in SpectraFiles)
                    {
                        // top N peptides in the file
                        double proteinIntensity = peptidesForThisProtein.Where(p => p.ProteinGroups.Count == 1 || useSharedPeptides)
                            .Select(p => p.GetIntensity(file)).OrderByDescending(p => p).Take(topNPeaks).Sum();

                        pg.SetIntensity(file, proteinIntensity);
                    }
                }
            }
        }

        public void WriteResults(string peaksOutputPath, string modPeptideOutputPath, string proteinOutputPath, string bayesianProteinQuantOutput, bool silent)
        {
            if (!silent)
            {
                Console.WriteLine("Writing output...");
            }

            if (peaksOutputPath != null)
            {
                using (StreamWriter output = new StreamWriter(peaksOutputPath))
                {
                    output.WriteLine(ChromatographicPeak.TabSeparatedHeader);

                    foreach (var peak in Peaks.SelectMany(p => p.Value)
                        .OrderBy(p => p.SpectraFileInfo.FilenameWithoutExtension)
                        .ThenByDescending(p => p.Intensity))
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

            if (bayesianProteinQuantOutput != null)
            {
                StringBuilder header = new StringBuilder();
                StringBuilder[] proteinStringBuilders = new StringBuilder[ProteinGroups.Count];

                for (int i = 0; i < proteinStringBuilders.Length; i++)
                {
                    proteinStringBuilders[i] = new StringBuilder();
                }

                using (StreamWriter output = new StreamWriter(bayesianProteinQuantOutput))
                {
                    if (!ProteinGroups.Any())
                    {
                        return;
                    }

                    var firstProteinQuantResults = ProteinGroups.First().Value.ConditionToQuantificationResults;

                    if (!firstProteinQuantResults.Any())
                    {
                        return;
                    }

                    string tabSepHeader = null;

                    if (firstProteinQuantResults.First().Value is PairedProteinQuantResult)
                    {
                        tabSepHeader = PairedProteinQuantResult.TabSeparatedHeader();
                    }
                    else
                    {
                        tabSepHeader = UnpairedProteinQuantResult.TabSeparatedHeader();
                    }

                    foreach (var condition in firstProteinQuantResults.Keys)
                    {
                        header.Append(tabSepHeader);

                        int p = 0;

                        // sort by protein false discovery rate, then by number of measurements
                        foreach (var protein in ProteinGroups
                            .OrderByDescending(v => v.Value.ConditionToQuantificationResults[condition].IsStatisticallyValid)
                            .ThenByDescending(v => v.Value.ConditionToQuantificationResults[condition].BayesFactor)
                            .ThenByDescending(v => v.Value.ConditionToQuantificationResults[condition].Peptides.Count))
                        {
                            proteinStringBuilders[p].Append(
                                protein.Value.ConditionToQuantificationResults[condition].ToString());

                            p++;
                        }
                    }

                    output.WriteLine(header);

                    foreach (var proteinStringBuilder in proteinStringBuilders)
                    {
                        output.WriteLine(proteinStringBuilder);
                    }
                }
            }

            if (!silent)
            {
                Console.WriteLine("Finished writing output");
            }
        }
    }
}