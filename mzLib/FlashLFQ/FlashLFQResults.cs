﻿using Easy.Common.Extensions;
using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using FlashLFQ.IsoTracker;

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

        /// <summary>
        /// This method uses the median polish algorithm to calculate protein quantities in each biological replicate.
        /// See https://mgimond.github.io/ES218/Week11a.html for an example of the median polish algorithm.
        /// </summary>
        public void CalculateProteinResultsMedianPolish(bool useSharedPeptides)
        {
            // reset protein intensities to 0
            foreach (var proteinGroup in ProteinGroups)
            {
                foreach (SpectraFileInfo file in SpectraFiles)
                {
                    proteinGroup.Value.SetIntensity(file, 0);
                }
            }

            // associate peptide w/ proteins in a dictionary for easy lookup
            List<Peptide> peptides = PeptideModifiedSequences.Values.Where(p => p.UnambiguousPeptideQuant()).ToList();
            Dictionary<ProteinGroup, List<Peptide>> proteinGroupToPeptides = new Dictionary<ProteinGroup, List<Peptide>>();

            foreach (Peptide peptide in peptides)
            {
                if (!peptide.UseForProteinQuant || (peptide.ProteinGroups.Count > 1 && !useSharedPeptides))
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

            var filesGroupedByCondition = SpectraFiles.GroupBy(p => p.Condition).ToList();

            // quantify each protein
            foreach (ProteinGroup proteinGroup in ProteinGroups.Values)
            {
                if (proteinGroupToPeptides.TryGetValue(proteinGroup, out var peptidesForThisProtein))
                {
                    // set up peptide intensity table
                    // top row is the column effects, left column is the row effects
                    // the other cells are peptide intensity measurements
                    int numSamples = SpectraFiles.Select(p => p.Condition + p.BiologicalReplicate).Distinct().Count();
                    double[][] peptideIntensityMatrix = new double[peptidesForThisProtein.Count + 1][];
                    for (int i = 0; i < peptideIntensityMatrix.Length; i++)
                    {
                        peptideIntensityMatrix[i] = new double[numSamples + 1];
                    }

                    // populate matrix w/ log2-transformed peptide intensities
                    // if a value is missing, it will be filled with NaN
                    int sampleN = 0;
                    foreach (var group in SpectraFiles.GroupBy(p => p.Condition).OrderBy(p => p.Key))
                    {
                        foreach (var sample in group.GroupBy(p => p.BiologicalReplicate).OrderBy(p => p.Key))
                        {
                            foreach (Peptide peptide in peptidesForThisProtein)
                            {
                                double sampleIntensity = 0;
                                double highestFractionIntensity = 0;

                                // the fraction w/ the highest intensity is used as the sample intensity for this peptide.
                                // if there is more than one replicate of the fraction, then the replicate intensities are averaged
                                foreach (var fraction in sample.GroupBy(p => p.Fraction))
                                {
                                    double fractionIntensity = 0;
                                    int replicatesWithValidValues = 0;

                                    foreach (SpectraFileInfo replicate in fraction.OrderBy(p => p.TechnicalReplicate))
                                    {
                                        double replicateIntensity = peptide.GetIntensity(replicate);

                                        if (replicateIntensity > 0)
                                        {
                                            fractionIntensity += replicateIntensity;
                                            replicatesWithValidValues++;
                                        }
                                    }

                                    if (replicatesWithValidValues > 0)
                                    {
                                        fractionIntensity /= replicatesWithValidValues;
                                    }

                                    if (fractionIntensity > highestFractionIntensity)
                                    {
                                        highestFractionIntensity = fractionIntensity;
                                        sampleIntensity = highestFractionIntensity;
                                    }
                                }

                                int sampleNumber = sample.Key;

                                if (sampleIntensity == 0)
                                {
                                    sampleIntensity = double.NaN;
                                }
                                else
                                {
                                    sampleIntensity = Math.Log(sampleIntensity, 2);
                                }

                                peptideIntensityMatrix[peptidesForThisProtein.IndexOf(peptide) + 1][sampleN + 1] = sampleIntensity;
                            }

                            sampleN++;
                        }
                    }

                    // if there are any peptides that have only one measurement, mark them as NaN
                    // unless we have ONLY peptides with one measurement
                    var peptidesWithMoreThanOneMmt = peptideIntensityMatrix.Skip(1).Count(row => row.Skip(1).Count(cell => !double.IsNaN(cell)) > 1);
                    if (peptidesWithMoreThanOneMmt > 0)
                    {
                        for (int i = 1; i < peptideIntensityMatrix.Length; i++)
                        {
                            int validValueCount = peptideIntensityMatrix[i].Count(p => !double.IsNaN(p) && p != 0);

                            if (validValueCount < 2 && numSamples >= 2)
                            {
                                for (int j = 1; j < peptideIntensityMatrix[0].Length; j++)
                                {
                                    peptideIntensityMatrix[i][j] = double.NaN;
                                }
                            }
                        }
                    }

                    // do median polish protein quantification
                    // row effects in a protein can be considered ~ relative ionization efficiency
                    // column effects are differences between conditions
                    MedianPolish(peptideIntensityMatrix);

                    double overallEffect = peptideIntensityMatrix[0][0];
                    double[] columnEffects = peptideIntensityMatrix[0].Skip(1).ToArray();
                    double referenceProteinIntensity = Math.Pow(2, overallEffect) * peptidesForThisProtein.Count;

                    // check for unquantifiable proteins; these are proteins w/ quantified peptides, but
                    // the protein is still not quantifiable because there are not peptides to compare across runs
                    List<string> possibleUnquantifiableSample = new List<string>();
                    sampleN = 0;
                    foreach (var group in SpectraFiles.GroupBy(p => p.Condition).OrderBy(p => p.Key))
                    {
                        foreach (var sample in group.GroupBy(p => p.BiologicalReplicate).OrderBy(p => p.Key))
                        {
                            bool isMissingValue = true;

                            foreach (SpectraFileInfo spectraFile in sample)
                            {
                                if (peptidesForThisProtein.Any(p => p.GetIntensity(spectraFile) != 0))
                                {
                                    isMissingValue = false;
                                    break;
                                }
                            }

                            if (!isMissingValue && columnEffects[sampleN] == 0)
                            {
                                possibleUnquantifiableSample.Add(group.Key + "_" + sample.Key);
                            }

                            sampleN++;
                        }
                    }

                    // set the sample protein intensities
                    sampleN = 0;
                    foreach (var group in SpectraFiles.GroupBy(p => p.Condition).OrderBy(p => p.Key))
                    {
                        foreach (var sample in group.GroupBy(p => p.BiologicalReplicate).OrderBy(p => p.Key))
                        {
                            // this step un-logs the protein "intensity". in reality this value is more like a fold-change 
                            // than an intensity, but unlike a fold-change it's not relative to a particular sample.
                            // by multiplying this value by the reference protein intensity calculated earlier, then we get 
                            // a protein intensity value
                            double columnEffect = columnEffects[sampleN];
                            double sampleProteinIntensity = Math.Pow(2, columnEffect) * referenceProteinIntensity;

                            // the column effect can be 0 in some cases. sometimes it's a valid value and sometimes it's not.
                            // so we need to check to see if it is actually a valid value
                            bool isMissingValue = true;

                            foreach (SpectraFileInfo spectraFile in sample)
                            {
                                if (peptidesForThisProtein.Any(p => p.GetIntensity(spectraFile) != 0))
                                {
                                    isMissingValue = false;
                                    break;
                                }
                            }

                            if (!isMissingValue)
                            {
                                if (possibleUnquantifiableSample.Count > 1 && possibleUnquantifiableSample.Contains(group.Key + "_" + sample.Key))
                                {
                                    proteinGroup.SetIntensity(sample.First(), double.NaN);
                                }
                                else
                                {
                                    proteinGroup.SetIntensity(sample.First(), sampleProteinIntensity);
                                }
                            }

                            sampleN++;
                        }
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

        public static void MedianPolish(double[][] table, int maxIterations = 10, double improvementCutoff = 0.0001)
        {
            // technically, this is weighted mean polish and not median polish.
            // but it should give similar results while being more robust to issues
            // arising from missing values.
            // the weights are inverse square difference to median.

            // subtract overall effect
            List<double> allValues = table.SelectMany(p => p.Where(p => !double.IsNaN(p) && p != 0)).ToList();

            if (allValues.Any())
            {
                double overallEffect = allValues.Median();
                table[0][0] += overallEffect;

                for (int r = 1; r < table.Length; r++)
                {
                    for (int c = 1; c < table[0].Length; c++)
                    {
                        table[r][c] -= overallEffect;
                    }
                }
            }

            double sumAbsoluteResiduals = double.MaxValue;

            for (int i = 0; i < maxIterations; i++)
            {
                // subtract row effects
                for (int r = 0; r < table.Length; r++)
                {
                    List<double> rowValues = table[r].Skip(1).Where(p => !double.IsNaN(p)).ToList();

                    if (rowValues.Any())
                    {
                        double rowMedian = rowValues.Median();
                        double[] weights = rowValues.Select(p => 1.0 / Math.Max(0.0001, Math.Pow(p - rowMedian, 2))).ToArray();
                        double rowEffect = rowValues.Sum(p => p * weights[rowValues.IndexOf(p)]) / weights.Sum();
                        table[r][0] += rowEffect;

                        for (int c = 1; c < table[0].Length; c++)
                        {
                            table[r][c] -= rowEffect;
                        }
                    }
                }

                // subtract column effects
                for (int c = 0; c < table[0].Length; c++)
                {
                    List<double> colValues = table.Skip(1).Select(p => p[c]).Where(p => !double.IsNaN(p)).ToList();

                    if (colValues.Any())
                    {
                        double colMedian = colValues.Median();
                        double[] weights = colValues.Select(p => 1.0 / Math.Max(0.0001, Math.Pow(p - colMedian, 2))).ToArray();
                        double colEffect = colValues.Sum(p => p * weights[colValues.IndexOf(p)]) / weights.Sum();
                        table[0][c] += colEffect;

                        for (int r = 1; r < table.Length; r++)
                        {
                            table[r][c] -= colEffect;
                        }
                    }
                }

                // calculate sum of absolute residuals and end the algorithm if it is not improving
                double iterationSumAbsoluteResiduals = table.Skip(1).SelectMany(p => p.Skip(1)).Where(p => !double.IsNaN(p)).Sum(p => Math.Abs(p));

                if (Math.Abs((iterationSumAbsoluteResiduals - sumAbsoluteResiduals) / sumAbsoluteResiduals) < improvementCutoff)
                {
                    break;
                }

                sumAbsoluteResiduals = iterationSumAbsoluteResiduals;
            }
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