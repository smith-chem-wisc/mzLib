using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace FlashLFQ
{
    public class FlashLfqResults
    {
        public readonly List<SpectraFileInfo> SpectraFiles;
        public readonly Dictionary<string, Peptide> PeptideModifiedSequences;
        public readonly Dictionary<string, ProteinGroup> ProteinGroups;
        public readonly Dictionary<SpectraFileInfo, List<ChromatographicPeak>> Peaks;

        public FlashLfqResults(List<SpectraFileInfo> spectraFiles)
        {
            SpectraFiles = spectraFiles;
            PeptideModifiedSequences = new Dictionary<string, Peptide>();
            ProteinGroups = new Dictionary<string, ProteinGroup>();
            Peaks = new Dictionary<SpectraFileInfo, List<ChromatographicPeak>>();
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
            foreach (var file in Peaks)
            {
                var groupedPeaks = file.Value.Where(p => p.NumIdentificationsByFullSeq == 1).GroupBy(p => p.Identifications.First().ModifiedSequence).ToList();

                foreach (var sequenceWithPeaks in groupedPeaks)
                {
                    string sequence = sequenceWithPeaks.Key;
                    double intensity = sequenceWithPeaks.Sum(p => p.Intensity);
                    DetectionType detectionType;
                    var pgs = new HashSet<ProteinGroup>(sequenceWithPeaks.SelectMany(p => p.Identifications).SelectMany(v => v.ProteinGroups));

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

                    if (!PeptideModifiedSequences.ContainsKey(sequence))
                    {
                        bool useForProteinQuant = sequenceWithPeaks.First().Identifications.First().UseForProteinQuant;
                        PeptideModifiedSequences.Add(sequence, new Peptide(sequence, useForProteinQuant));
                    }

                    PeptideModifiedSequences[sequence].SetIntensity(file.Key, intensity);
                    PeptideModifiedSequences[sequence].SetDetectionType(file.Key, detectionType);
                    PeptideModifiedSequences[sequence].ProteinGroups = pgs;
                }

                // report ambiguous quantification
                var ambiguousPeaks = file.Value.Where(p => p.NumIdentificationsByFullSeq > 1).ToList();
                foreach (ChromatographicPeak ambiguousPeak in ambiguousPeaks)
                {
                    foreach (Identification id in ambiguousPeak.Identifications)
                    {
                        string sequence = id.ModifiedSequence;
                        
                        if (!PeptideModifiedSequences.ContainsKey(sequence))
                        {
                            bool useForProteinQuant = id.UseForProteinQuant;
                            PeptideModifiedSequences.Add(sequence, new Peptide(sequence, useForProteinQuant));
                        }

                        double alreadyRecordedIntensity = PeptideModifiedSequences[sequence].GetIntensity(file.Key);
                        double fractionAmbiguous = (ambiguousPeak.Intensity + alreadyRecordedIntensity) / alreadyRecordedIntensity;

                        if (fractionAmbiguous > 0.3)
                        {
                            PeptideModifiedSequences[sequence].SetIntensity(file.Key, 0);
                            PeptideModifiedSequences[sequence].SetDetectionType(file.Key, DetectionType.MSMSAmbiguousPeakfinding);
                            PeptideModifiedSequences[sequence].ProteinGroups = id.ProteinGroups;
                        }
                    }
                }
            }
        }

        public void CalculateProteinResultsTop3()
        {
            int topNPeaks = 3;

            List<ProteinGroup> allProteinGroups = PeptideModifiedSequences.Values.SelectMany(p => p.ProteinGroups).Distinct().ToList();
            foreach (ProteinGroup pg in allProteinGroups)
            {
                ProteinGroups.Add(pg.ProteinGroupName, pg);
            }

            List<Peptide> peptides = PeptideModifiedSequences.Values.Where(p => p.UnambiguousPeptideQuant()).ToList();
            Dictionary<ProteinGroup, List<Peptide>> proteinGroupToPeptides = new Dictionary<ProteinGroup, List<Peptide>>();

            foreach(Peptide peptide in peptides)
            {
                if(!peptide.UseForProteinQuant)
                {
                    continue;
                }

                foreach(ProteinGroup pg in peptide.ProteinGroups)
                {
                    if(proteinGroupToPeptides.TryGetValue(pg,out var peptidesForThisProtein))
                    {
                        peptidesForThisProtein.Add(peptide);
                    }
                    else
                    {
                        proteinGroupToPeptides.Add(pg, new List<Peptide> { peptide });
                    }
                }
            }

            foreach(ProteinGroup pg in ProteinGroups.Values)
            {
                if(proteinGroupToPeptides.TryGetValue(pg, out var peptidesForThisProtein))
                {
                    foreach(SpectraFileInfo file in SpectraFiles)
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

        public void WriteResults(string peaksOutputPath, string modPeptideOutputPath, string proteinOutputPath)
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
        }
    }
}