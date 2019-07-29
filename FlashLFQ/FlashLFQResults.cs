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
                if (!PeptideModifiedSequences.ContainsKey(id.ModifiedSequence))
                {
                    PeptideModifiedSequences.Add(id.ModifiedSequence, new Peptide(id.ModifiedSequence, id.UseForProteinQuant, id.proteinGroups));
                }

                foreach (ProteinGroup proteinGroup in id.proteinGroups)
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
                var groupedPeaks = filePeaks.Value.Where(p => p.NumIdentificationsByFullSeq == 1).GroupBy(p => p.Identifications.First().ModifiedSequence).ToList();

                foreach (var sequenceWithPeaks in groupedPeaks)
                {
                    string sequence = sequenceWithPeaks.Key;
                    double intensity = sequenceWithPeaks.Sum(p => p.Intensity);
                    DetectionType detectionType;

                    if (sequenceWithPeaks.First().IsMbrPeak && intensity > 0)
                    {
                        detectionType = DetectionType.MBR;
                    }
                    else if (!sequenceWithPeaks.First().IsMbrPeak && intensity > 0)
                    {
                        detectionType = DetectionType.MSMS;
                    }
                    else if (!sequenceWithPeaks.First().IsMbrPeak && intensity == 0)
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
                        double fractionAmbiguous = (ambiguousPeak.Intensity + alreadyRecordedIntensity) / alreadyRecordedIntensity;

                        if (fractionAmbiguous > 0.3)
                        {
                            PeptideModifiedSequences[sequence].SetIntensity(filePeaks.Key, 0);
                            PeptideModifiedSequences[sequence].SetDetectionType(filePeaks.Key, DetectionType.MSMSAmbiguousPeakfinding);
                        }
                    }
                }
            }
        }

        public void CalculateProteinResultsTop3()
        {
            foreach (var proteinGroup in ProteinGroups)
            {
                foreach (SpectraFileInfo file in SpectraFiles)
                {
                    proteinGroup.Value.SetIntensity(file, 0);
                }
            }

            int topNPeaks = 3;

            List<ProteinGroup> allProteinGroups = Peaks.Values.SelectMany(p => p)
                .SelectMany(p => p.Identifications)
                .SelectMany(p => p.proteinGroups).Distinct().ToList();
            
            List<ChromatographicPeak> unambiguousPeaks = Peaks.Values.SelectMany(p => p).Where(p => p.NumIdentificationsByFullSeq == 1 && p.Intensity > 0).ToList();
            Dictionary<ProteinGroup, List<ChromatographicPeak>> proteinGroupToPeaks = new Dictionary<ProteinGroup, List<ChromatographicPeak>>();

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
                        // top N peaks, prioritizing protein-uniqueness and then intensity
                        double proteinIntensity =
                            peptidesForThisProtein.OrderBy(p => p.ProteinGroups.Distinct().Count())
                            .ThenByDescending(p => p.GetIntensity(file)).Take(topNPeaks).Sum(p => p.GetIntensity(file));

                        pg.SetIntensity(file, proteinIntensity);
                    }
                }
            }
        }

        public void WriteResults(string peaksOutputPath, string modPeptideOutputPath, string proteinOutputPath, string bayesianProteinQuantOutput)
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

                    foreach (var condition in firstProteinQuantResults.Keys)
                    {
                        header.Append(ProteinQuantificationEngineResult.TabSeparatedHeader());

                        int p = 0;

                        // sort by protein false discovery rate, then by number of fold-change measurements
                        foreach (var protein in ProteinGroups
                            .OrderBy(v => v.Value.conditionToQuantificationResults[condition].FalseDiscoveryRate)
                            .ThenByDescending(v => v.Value.conditionToQuantificationResults[condition].peptideFoldChangeMeasurements.SelectMany(b => b.foldChanges).Count()))
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
        }
    }
}