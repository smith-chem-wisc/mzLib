using FlashLFQ.IsoTracker;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace FlashLFQ
{
    public class FlashLfqResults
    {
        public bool IsoTracker = false;
        public readonly List<SpectraFileInfo> SpectraFiles;
        public readonly Dictionary<string, Peptide> PeptideModifiedSequences;
        public readonly Dictionary<string, ProteinGroup> ProteinGroups;
        public readonly Dictionary<SpectraFileInfo, List<ChromatographicPeak>> Peaks;
        private readonly HashSet<string> _peptideModifiedSequencesToQuantify;
        public IDictionary<string, Dictionary<PeakRegion, List<ChromatographicPeak>>> IsobaricPeptideDict = null;
        public string PepResultString { get; set; }
        public double MbrQValueThreshold { get; set; }

        public FlashLfqResults(List<SpectraFileInfo> spectraFiles, List<Identification> identifications, double mbrQValueThreshold = 0.05,
            HashSet<string> peptideModifiedSequencesToQuantify = null, bool isIsoTracker = false)
        {
            SpectraFiles = spectraFiles;
            PeptideModifiedSequences = new Dictionary<string, Peptide>();
            ProteinGroups = new Dictionary<string, ProteinGroup>();
            Peaks = new Dictionary<SpectraFileInfo, List<ChromatographicPeak>>();
            MbrQValueThreshold = mbrQValueThreshold;
            _peptideModifiedSequencesToQuantify = peptideModifiedSequencesToQuantify ?? identifications.Where(id => !id.IsDecoy).Select(id => id.ModifiedSequence).ToHashSet();
            IsoTracker = isIsoTracker;

            foreach (SpectraFileInfo file in spectraFiles)
            {
                Peaks.Add(file, new List<ChromatographicPeak>());
            }

            // Only quantify peptides within the set of valid peptide modified (full) sequences. This is done to enable pepitde-level FDR control of reported results
            foreach (Identification id in identifications.Where(id => !id.IsDecoy & _peptideModifiedSequencesToQuantify.Contains(id.ModifiedSequence)))
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
                    if (IsoTracker)
                    {
                        output.WriteLine(Peptide.TabSeparatedHeader_IsoTracker(SpectraFiles));
                        // we want to output with same iso group index followed by peak order.
                        foreach (var peptide in PeptideModifiedSequences
                                     .OrderBy(p => p.Value.IsoGroupIndex ?? int.MaxValue)
                                     .ThenBy(p => p.Value.PeakOrder ?? int.MinValue))
                        {
                            output.WriteLine(peptide.Value.ToString(SpectraFiles, IsoTracker));
                        }
                    }
                    else
                    {
                        output.WriteLine(Peptide.TabSeparatedHeader(SpectraFiles));
                        foreach (var peptide in PeptideModifiedSequences.OrderBy(p => p.Key))
                        {
                            output.WriteLine(peptide.Value.ToString(SpectraFiles, IsoTracker));
                        }
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

        public void ReNormalizeResults(bool integrate = false, int maxThreads = 10, bool useSharedPeptides = false)
        {
            foreach(var peak in Peaks.SelectMany(p => p.Value))
            {
                peak.CalculateIntensityForThisFeature(integrate);
            }
            new IntensityNormalizationEngine(this, integrate, silent: true, maxThreads).NormalizeResults();
            CalculatePeptideResults(quantifyAmbiguousPeptides: false);
            CalculateProteinResultsMedianPolish(useSharedPeptides: useSharedPeptides);
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

        public void CalculatePeptideResults(bool quantifyAmbiguousPeptides)
        {
            foreach (var sequence in PeptideModifiedSequences)
            {
                foreach (SpectraFileInfo file in SpectraFiles)
                {
                    sequence.Value.SetDetectionType(file, DetectionType.NotDetected);
                    sequence.Value.SetIntensity(file, 0);
                    sequence.Value.SetRetentionTime(file,0);
                }
            }


            foreach (var filePeaks in Peaks)
            {
                var groupedPeaks = filePeaks.Value
                    .Where(p => p.NumIdentificationsByFullSeq == 1)
                    .Where(p => !p.Identifications.First().IsDecoy)
                    .Where(p => p.DetectionType != DetectionType.MBR || (p is MbrChromatographicPeak m && m.MbrQValue < MbrQValueThreshold && !m.RandomRt))
                    .GroupBy(p => p.Identifications.First().ModifiedSequence)
                    .Where(group => _peptideModifiedSequencesToQuantify.Contains(group.Key))
                    .ToDictionary(p => p.Key, p => p.ToList());

                foreach (var sequenceWithPeaks in groupedPeaks)
                {
                    string sequence = sequenceWithPeaks.Key;
                    double intensity = sequenceWithPeaks.Value.Max(p => p.Intensity);
                    ChromatographicPeak bestPeak = sequenceWithPeaks.Value.First(p => p.Intensity == intensity);
                    DetectionType detectionType;

                    if (bestPeak.DetectionType == DetectionType.MBR && intensity > 0)
                    {
                        detectionType = DetectionType.MBR;
                    }
                    else if (bestPeak.DetectionType != DetectionType.MBR && intensity > 0)
                    {
                        detectionType = DetectionType.MSMS;
                    }
                    else if (bestPeak.DetectionType != DetectionType.MBR && intensity == 0)
                    {
                        detectionType = DetectionType.MSMSIdentifiedButNotQuantified;
                    }
                    else
                    {
                        detectionType = DetectionType.NotDetected;
                    }

                    PeptideModifiedSequences[sequence].SetIntensity(filePeaks.Key, intensity);
                    PeptideModifiedSequences[sequence].SetRetentionTime(filePeaks.Key, bestPeak.ApexRetentionTime);
                    PeptideModifiedSequences[sequence].SetDetectionType(filePeaks.Key, detectionType);
                }

                // report ambiguous quantification
                var ambiguousPeaks = filePeaks.Value
                    .Where(p => p.NumIdentificationsByFullSeq > 1)
                    .Where(p => !p.Identifications.First().IsDecoy)
                    .Where(p => p.DetectionType != DetectionType.MBR || (p is MbrChromatographicPeak m && m.MbrQValue < MbrQValueThreshold && !m.RandomRt))
                    .ToList();
                foreach (ChromatographicPeak ambiguousPeak in ambiguousPeaks)
                {
                    foreach (Identification id in ambiguousPeak.Identifications.Where(id => !id.IsDecoy))
                    {
                        if (!_peptideModifiedSequencesToQuantify.Contains(id.ModifiedSequence)) continue; // Ignore the ids/sequences we don't want to quantify

                        string sequence = id.ModifiedSequence;

                        double alreadyRecordedIntensity = PeptideModifiedSequences[sequence].GetIntensity(filePeaks.Key);
                        double fractionAmbiguous = ambiguousPeak.Intensity / (alreadyRecordedIntensity + ambiguousPeak.Intensity);

                        if (quantifyAmbiguousPeptides)
                        {
                            // If the peptide intensity hasn't been recorded, the intensity is set equal to the intensity of the ambiguous peak
                            if (Math.Abs(alreadyRecordedIntensity) < 0.01)
                            {
                                PeptideModifiedSequences[sequence].SetDetectionType(filePeaks.Key, DetectionType.MSMSAmbiguousPeakfinding);
                                PeptideModifiedSequences[sequence].SetRetentionTime(filePeaks.Key, ambiguousPeak.ApexRetentionTime);
                                PeptideModifiedSequences[sequence].SetIntensity(filePeaks.Key, ambiguousPeak.Intensity);
                            }
                            // If the peptide intensity has already been recorded, that value is retained. 
                            else if (fractionAmbiguous > 0.3)
                            {
                                PeptideModifiedSequences[sequence].SetDetectionType(filePeaks.Key, DetectionType.MSMSAmbiguousPeakfinding);
                            }
                        }
                        else if (fractionAmbiguous > 0.3)
                        {
                            PeptideModifiedSequences[sequence].SetDetectionType(filePeaks.Key, DetectionType.MSMSAmbiguousPeakfinding);
                            PeptideModifiedSequences[sequence].SetIntensity(filePeaks.Key, 0);
                            PeptideModifiedSequences[sequence].SetRetentionTime(filePeaks.Key, ambiguousPeak.ApexRetentionTime);
                        }
                    }
                }
                
            }

            if (IsoTracker && IsobaricPeptideDict != null)
            {
                // We view each Isobaric peak as an individual peptide, so we need to add them to the peptide list
                RevisedModifiedPeptides();
            }

            if (!quantifyAmbiguousPeptides)
            {
                HandleAmbiguityInFractions();
            }
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
                        foreach (ChromatographicPeak peak in Peaks[file].Where(p => p.DetectionType != DetectionType.MBR 
                            || (p is MbrChromatographicPeak m && m.MbrQValue < MbrQValueThreshold)))
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
            var strategy = new Top3ProteinQuantificationStrategy<Peptide, ProteinGroup, SpectraFileInfo>();
            var engine = new GenericProteinQuantificationEngine<Peptide, ProteinGroup, SpectraFileInfo>(
                strategy,
                PeptideModifiedSequences,
                ProteinGroups,
                SpectraFiles);

            engine.Run(useSharedPeptides);
        }

        /// <summary>
        /// This method uses the median polish algorithm to calculate protein quantities in each biological replicate.
        /// See https://mgimond.github.io/ES218/Week11a.html for an example of the median polish algorithm.
        /// </summary>
        public void CalculateProteinResultsMedianPolish(bool useSharedPeptides)
        {
            var strategy = new MedianPolishProteinQuantificationStrategy<Peptide, ProteinGroup, SpectraFileInfo>();
            var engine = new GenericProteinQuantificationEngine<Peptide, ProteinGroup, SpectraFileInfo>(
                strategy,
                PeptideModifiedSequences,
                ProteinGroups,
                SpectraFiles);

            engine.Run(useSharedPeptides);
        }

        /// <summary>
        /// This method is used to re-edit the peptide List by adding the isobaric peptides and remove the former peptide.
        /// </summary>
        internal void RevisedModifiedPeptides()
        {
            int isoGroupIndex = 1;
            //If the isobaric peptide dictionary is not empty, then we need to revise the peptide list.
            foreach (var isoPeptides in IsobaricPeptideDict.Where(p=>p.Value.Count != 0)) 
            {
                string peptideSequence = isoPeptides.Key;
                Peptide originalPeptide = PeptideModifiedSequences[peptideSequence];

                // Remove the formal peptide from the peptide list
                var allIDs = isoPeptides.Value.Values
                    .SelectMany(p => p)
                    .Where(p => p != null)
                    .SelectMany(p=>p.Identifications)
                    .Where(p=>p.BaseSequence == originalPeptide.BaseSequence) // Avoid to remove any peptide with different base sequence
                    .DistinctBy(p=>p.ModifiedSequence)
                    .Select(p=>p.ModifiedSequence)
                    .ToList();

                foreach (var modSeq in allIDs)
                {
                    if (PeptideModifiedSequences.ContainsKey(modSeq))
                    {
                        PeptideModifiedSequences.Remove(modSeq);
                    }
                }

                // Add the isobaric peptides to the peptide list

                //If there is only one peak for the isobaric peptides, then we don't view them as isobaric peptides.
                if (isoPeptides.Value.Values.Count == 1)
                {
                    var isoPeptidePeaks = isoPeptides.Value.Values.First();
                    var allSeq = isoPeptidePeaks
                        .Where(p => p != null)
                        .SelectMany(p => p.Identifications)
                        .Where(p=>p.BaseSequence == originalPeptide.BaseSequence) // do not output the peptide with different base sequence in the peptide result
                        .Select(p => p.ModifiedSequence)
                        .Distinct()
                        .ToList();
                    Peptide peptide = new Peptide(string.Join(" | ", allSeq), originalPeptide.BaseSequence, originalPeptide.UseForProteinQuant, originalPeptide.ProteinGroups);
                    peptide.SetIsobaricPeptide(isoPeptidePeaks); //When we set the peptide as IsobaricPeptide, then the retention time, intensity and detectionType will be set from the chromPeak automatically.
                    PeptideModifiedSequences[peptide.Sequence] = peptide;
                }
                //If there are multiple peaks for the isobaric peptides, then we view them as isobaric peptides.
                else
                {
                    int peakIndex = 1;
                    foreach (var isoPeptidePeaks in isoPeptides.Value.Values.ToList())
                    {
                        var allSeq = isoPeptidePeaks
                            .Where(p => p != null)
                            .SelectMany(p => p.Identifications)
                            .Where(p=>p.BaseSequence == originalPeptide.BaseSequence)// do not output the peptide with different base sequence that was merged in RunErrorCheck
                            .Select(p => p.ModifiedSequence)
                            .Distinct()
                            .ToList();
                        Peptide peptide = new Peptide(string.Join(" | ", allSeq) + " Isopeptide_peak" + peakIndex, originalPeptide.BaseSequence, originalPeptide.UseForProteinQuant, originalPeptide.ProteinGroups, isoGroupIndex, peakIndex);
                        peptide.SetIsobaricPeptide(isoPeptidePeaks); //When we set the peptide as IsobaricPeptide, then the retention time, intensity and detectionType will be set from the chromPeak automatically.
                        PeptideModifiedSequences[peptide.Sequence] = peptide;
                        peakIndex++;
                    }
                    isoGroupIndex++;
                }
            }
        }
    }
}