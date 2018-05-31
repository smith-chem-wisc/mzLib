using FlashLFQ.BoundedNelderMeadOptimizer;
using MathNet.Numerics.Statistics;
using SharpLearning.Optimization;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace FlashLFQ
{
    public class IntensityNormalizationEngine
    {
        private const int numPeptidesDesiredFromEachFraction = 500;
        private const int numPeptidesDesiredInMatrix = 5000;
        private readonly FlashLFQResults results;
        private readonly bool integrate;
        private readonly bool silent;

        #region Public Constructors

        public IntensityNormalizationEngine(FlashLFQResults results, bool integrate, bool silent)
        {
            this.results = results;
            this.integrate = integrate;
            this.silent = silent;
        }

        #endregion Public Constructors

        #region Public Methods

        public void NormalizeResults()
        {
            results.CalculatePeptideResults(false);
            
            // run normalization functions, recalculating intensity between each function
            if (!silent)
            {
                Console.WriteLine("Normalizing techreps");
            }
            NormalizeTechreps();
            results.CalculatePeptideResults(false);

            if (!silent)
            {
                Console.WriteLine("Normalizing bioreps and fractions");
            }
            NormalizeBiorepsAndFractions();
            results.CalculatePeptideResults(false);

            if (!silent)
            {
                Console.WriteLine("Normalizing conditions");
            }
            NormalizeConditions();
            results.CalculatePeptideResults(false);
        }

        #endregion Public Methods

        #region Private Methods

        private void NormalizeTechreps()
        {
            List<double> nonZeroIntensities = new List<double>();

            var peptides = results.peptideModifiedSequences.Select(v => v.Value).ToList();
            var conditions = results.spectraFiles.GroupBy(p => p.condition);

            foreach (var condition in conditions)
            {
                var bioreps = condition.GroupBy(p => p.biologicalReplicate);

                foreach (var biorep in bioreps)
                {
                    var fractions = biorep.GroupBy(p => p.fraction);

                    foreach (var fraction in fractions)
                    {
                        int numReplicates = fraction.Max(p => p.technicalReplicate);
                        var seqsToPeaksPerFile = new Dictionary<string, List<ChromatographicPeak>>[numReplicates + 1];

                        for (int i = 0; i < seqsToPeaksPerFile.Length; i++)
                        {
                            var rep = fraction.Where(p => p.technicalReplicate == i).FirstOrDefault();

                            if (rep != null)
                            {
                                seqsToPeaksPerFile[i] = new Dictionary<string, List<ChromatographicPeak>>();

                                foreach (var peak in results.peaks[rep])
                                {
                                    if (seqsToPeaksPerFile[i].TryGetValue(peak.identifications.First().ModifiedSequence, out var list))
                                    {
                                        list.Add(peak);
                                    }
                                    else
                                    {
                                        seqsToPeaksPerFile[i].Add(peak.identifications.First().ModifiedSequence, new List<ChromatographicPeak> { peak });
                                    }
                                }
                            }
                        }

                        for (int p = 0; p < peptides.Count; p++)
                        {
                            double avgIntensity = 0;
                            nonZeroIntensities.Clear();

                            foreach (var technicalRep in fraction)
                            {
                                if (peptides[p].GetIntensity(technicalRep) > 0)
                                {
                                    nonZeroIntensities.Add(peptides[p].GetIntensity(technicalRep));
                                }
                            }
                            if (nonZeroIntensities.Any())
                            {
                                avgIntensity = nonZeroIntensities.Average();

                                foreach (var technicalRep in fraction)
                                {
                                    double normFactorForThisTechRepAndPeptide = avgIntensity / peptides[p].GetIntensity(technicalRep);

                                    if (seqsToPeaksPerFile[technicalRep.technicalReplicate].TryGetValue(peptides[p].Sequence, out var peaksForPepAndTechrep))
                                    {
                                        //TODO: set intensity to 0 if very different between technical reps

                                        foreach (var peak in peaksForPepAndTechrep)
                                        {
                                            foreach (var isotopeEnvelope in peak.isotopicEnvelopes)
                                            {
                                                isotopeEnvelope.Normalize(normFactorForThisTechRepAndPeptide);
                                            }

                                            // recalculate intensity after normalization
                                            peak.CalculateIntensityForThisFeature(integrate);
                                        }
                                    }
                                    else
                                    {
                                        peptides[p].SetIntensity(technicalRep, avgIntensity);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        private void NormalizeBiorepsAndFractions()
        {
            var peptides = results.peptideModifiedSequences.Select(v => v.Value).ToList();
            var conditions = results.spectraFiles.Select(p => p.condition).Distinct().OrderBy(p => p).ToList();
            var filesForCond1Biorep1 = results.spectraFiles.Where(p => p.condition == conditions[0] && p.biologicalReplicate == 0).ToList();

            foreach (var condition in conditions)
            {
                var filesForThisCondition = results.spectraFiles.Where(p => p.condition.Equals(condition)).ToList();

                int numB = filesForThisCondition.Select(p => p.biologicalReplicate).Distinct().Count();

                for (int b = 0; b < numB; b++)
                {
                    // condition 1 biorep 1 is the reference, don't normalize it
                    if (b == 0 && conditions.IndexOf(condition) == 0)
                        continue;

                    // run the normalization function
                    if (!silent)
                    {
                        Console.WriteLine("Normalizing condition \"" + condition + "\" biorep " + (b + 1));
                    }

                    var filesForThisBiorep = filesForThisCondition.Where(p => p.biologicalReplicate == b);

                    int numF = Math.Max(filesForCond1Biorep1.Max(p => p.fraction), filesForThisBiorep.Max(p => p.fraction)) + 1;

                    // only normalize on peptides seen in both bioreps
                    List<Peptide> seenInBothBioreps = new List<Peptide>();
                    for (int p = 0; p < peptides.Count; p++)
                    {
                        bool seenInBiorep1 = false;
                        bool seenInBiorep2 = false;

                        foreach (var file in filesForCond1Biorep1)
                        {
                            if (peptides[p].GetIntensity(file) > 0)
                            {
                                seenInBiorep1 = true;
                            }
                        }

                        foreach (var file in filesForThisBiorep)
                        {
                            if (peptides[p].GetIntensity(file) > 0)
                            {
                                seenInBiorep2 = true;
                            }
                        }

                        if (seenInBiorep1 && seenInBiorep2)
                        {
                            seenInBothBioreps.Add(peptides[p]);
                        }
                    }

                    // we're only normalizing on a subset of data here because it might take too long to do the whole set
                    seenInBothBioreps = SubsetData(seenInBothBioreps, filesForThisBiorep.Concat(filesForCond1Biorep1).ToList());

                    // add the data to the array to set up for the normalization function
                    int numP = seenInBothBioreps.Count;
                    double[,,] myIntensityArray = new double[numP, 2, numF];

                    for (int p = 0; p < numP; p++)
                    {
                        var peptide = seenInBothBioreps[p];

                        foreach (var file in filesForCond1Biorep1)
                        {
                            myIntensityArray[p, 0, file.fraction] = peptide.GetIntensity(file);
                        }

                        foreach (var file in filesForThisBiorep)
                        {
                            myIntensityArray[p, 1, file.fraction] = peptide.GetIntensity(file);
                        }
                    }

                    // solve for normalization factors
                    var normFactors = GetNormalizationFactors(myIntensityArray, numP, 2, numF, b, condition);
                    if (normFactors.All(p => p == 1.0) && !silent)
                    {
                        Console.WriteLine("Warning: Could not solve for optimal normalization factors for condition \"" + condition + "\" biorep " + (b + 1));
                    }

                    // multiply each peak's isotope envelopes by its file's normalization factor
                    foreach (var spectraFile in filesForThisBiorep)
                    {
                        foreach (var peak in results.peaks[spectraFile])
                        {
                            foreach (var isotopeEnvelope in peak.isotopicEnvelopes)
                            {
                                isotopeEnvelope.Normalize(normFactors[spectraFile.fraction]);
                            }

                            // recalculate intensity after normalization
                            peak.CalculateIntensityForThisFeature(integrate);
                        }
                    }
                }
            }
        }

        private void NormalizeConditions()
        {
            var peptides = results.peptideModifiedSequences.Select(v => v.Value).ToList();
            var conditions = results.spectraFiles.Select(p => p.condition).Distinct().OrderBy(p => p).ToList();

            // calculate intensities for the first condition (all other conditions are normalized to this condition)
            var filesForConditionOne = results.spectraFiles.Where(p => p.condition.Equals(conditions[0])).ToList();

            double[] peptideIntensitiesForConditionOne = new double[peptides.Count];
            for (int p = 0; p < peptides.Count; p++)
            {
                List<double> biorepIntensities = new List<double>();

                foreach (var file in filesForConditionOne)
                {
                    if (peptides[p].GetIntensity(file) > 0)
                    {
                        biorepIntensities.Add(peptides[p].GetIntensity(file));
                    }
                }

                if (biorepIntensities.Any())
                {
                    peptideIntensitiesForConditionOne[p] = biorepIntensities.Average();
                }
            }

            // calculate intensities for other conditions and normalize to the first condition
            for (int c = 1; c < conditions.Count; c++)
            {
                // calculate intensities for this condition
                var filesForThisCondition = results.spectraFiles.Where(p => p.condition.Equals(conditions[c])).ToList();
                double[] peptideIntensitiesForThisCondition = new double[peptides.Count];

                for (int p = 0; p < peptides.Count; p++)
                {
                    List<double> biorepIntensities = new List<double>();

                    foreach (var file in filesForThisCondition)
                    {
                        if (peptides[p].GetIntensity(file) > 0)
                            biorepIntensities.Add(peptides[p].GetIntensity(file));
                    }

                    if (biorepIntensities.Any())
                        peptideIntensitiesForThisCondition[p] = biorepIntensities.Average();
                }

                // calculate fold-changes between this condition and the first condition
                List<double> logFoldChanges = new List<double>();
                for (int p = 0; p < peptides.Count; p++)
                {
                    if (peptideIntensitiesForConditionOne[p] > 0 && peptideIntensitiesForThisCondition[p] > 0)
                    {
                        logFoldChanges.Add(Math.Log(peptideIntensitiesForThisCondition[p], 2) - Math.Log(peptideIntensitiesForConditionOne[p], 2));
                    }
                }

                double medianLogFoldChangeBetweenConditions = Statistics.Median(logFoldChanges);
                double unloggedMedianLogFoldChange = Math.Pow(2, medianLogFoldChangeBetweenConditions);
                double conditionNormalizationFactor = 1.0 / unloggedMedianLogFoldChange;

                // normalize to median fold-change
                foreach (var file in filesForThisCondition)
                {
                    foreach (var peak in results.peaks[file])
                    {
                        foreach (var isotopeEnvelope in peak.isotopicEnvelopes)
                        {
                            isotopeEnvelope.Normalize(conditionNormalizationFactor);
                        }

                        // recalculate intensity after normalization
                        peak.CalculateIntensityForThisFeature(integrate);
                    }
                }
            }
        }

        private List<Peptide> SubsetData(List<Peptide> initialList, List<SpectraFileInfo> spectraFiles)
        {
            List<SpectraFileInfo>[] bothBioreps = new List<SpectraFileInfo>[2];
            var temp1 = spectraFiles.GroupBy(p => p.condition).ToList();
            if (temp1.Count() == 1)
            {
                // normalizing bioreps within a condition
                var temp2 = spectraFiles.GroupBy(p => p.biologicalReplicate).ToList();
                bothBioreps[0] = temp2[0].ToList();
                bothBioreps[1] = temp2[1].ToList();
            }
            else
            {
                // normalizing bioreps between conditions
                bothBioreps[0] = temp1[0].ToList();
                bothBioreps[1] = temp1[1].ToList();
            }

            HashSet<Peptide> subsetList = new HashSet<Peptide>();
            int maxFractionIndex = bothBioreps.SelectMany(p => p).Max(v => v.fraction);

            foreach (var biorep in bothBioreps)
            {
                List<int> fractions = biorep.Select(p => p.fraction).Distinct().ToList();

                int numToAddPerFraction = numPeptidesDesiredInMatrix / fractions.Count;
                if (numToAddPerFraction < numPeptidesDesiredFromEachFraction)
                {
                    numToAddPerFraction = numPeptidesDesiredFromEachFraction;
                }

                int[] peptidesAddedPerFraction = new int[fractions.Count];
                Queue<Peptide>[] peptidesObservedInEachFraction = new Queue<Peptide>[fractions.Count];

                foreach (var file in biorep)
                {
                    if (peptidesObservedInEachFraction[file.fraction] == null)
                    {
                        peptidesObservedInEachFraction[file.fraction] = new Queue<Peptide>(initialList.Where(p => p.GetIntensity(file) > 0)
                            .OrderByDescending(p => p.GetIntensity(file)));
                    }
                }

                foreach (var fraction in fractions)
                {
                    while (peptidesAddedPerFraction[fraction] < numToAddPerFraction && peptidesObservedInEachFraction[fraction].Any())
                    {
                        var peptide = peptidesObservedInEachFraction[fraction].Dequeue();

                        // don't need to check if the return list already contains the peptide because it's a HashSet (no duplicates are allowed)
                        subsetList.Add(peptide);

                        // this peptide is in the return list regardless of whether or not it was actually just added;
                        // we just want to guarantee this fraction has 500 peptides in the return list to normalize with
                        peptidesAddedPerFraction[fraction]++;
                    }
                }
            }

            return subsetList.ToList();
        }

        private static double[] GetNormalizationFactors(double[,,] proteinIntensities, int numP, int numB, int numF, int bb, string cond)
        {
            double step = 0.01;
            object locker = new object();

            // initialize normalization factors to 1.0
            // normalization factor optimization must improve on these to be valid
            double bestError = 0;
            double[] bestNormFactors = new double[numF];
            for (int i = 0; i < bestNormFactors.Length; i++)
            {
                bestNormFactors[i] = 1.0;
            }

            // constraint (normalization factors must be >0.3 and <3
            var parameterArray = new ParameterBounds[numF];
            for (int f = 0; f < numF; f++)
            {
                parameterArray[f] = new ParameterBounds(0.3, 3, Transform.Linear);
            }

            //TODO: put this into a helper method to avoid code repetition...
            // calculate the error between bioreps if all norm factors are 1 (initial error)
            {
                double[,] originalBiorepIntensities = new double[numP, numB];
                double[] temp = new double[2];

                for (int p = 0; p < numP; p++)
                {
                    for (int b = 0; b < numB; b++)
                    {
                        for (int f = 0; f < numF; f++)
                        {
                            originalBiorepIntensities[p, b] += proteinIntensities[p, b, f];
                        }
                    }
                }

                for (int p = 0; p < numP; p++)
                {
                    for (int b2 = 1; b2 < numB; b2++)
                    {
                        temp[0] = originalBiorepIntensities[p, 0];
                        temp[1] = originalBiorepIntensities[p, b2];

                        // calculate initial error (error if all norm factors are 1)
                        // error metric is sum square error of coefficient of variation of each peptide
                        double coefficientOfVariation = Statistics.StandardDeviation(temp) / temp.Average();
                        bestError += Math.Pow(coefficientOfVariation, 2);
                    }
                }
            }

            int startN = (int)(parameterArray[0].Min / step);
            int endN = (int)(parameterArray[0].Max / step);

            Parallel.For(startN, endN, n =>
            {
                double startPosition = n * step;

                double[,] biorepIntensities = new double[numP, numB];
                double[] temp = new double[2];

                // define minimization metric
                Func<double[], OptimizerResult> minimize = v =>
                {
                    // sum the intensities with the current normalization factors
                    Array.Clear(biorepIntensities, 0, biorepIntensities.Length);

                    for (int p = 0; p < numP; p++)
                    {
                        for (int b = 0; b < numB; b++)
                        {
                            for (int f = 0; f < numF; f++)
                            {
                                if (b == 0)
                                {
                                    biorepIntensities[p, b] += proteinIntensities[p, b, f];
                                }
                                else
                                {
                                    biorepIntensities[p, b] += proteinIntensities[p, b, f] * v[f];
                                }
                            }
                        }
                    }

                    // calculate the error for these normalization factors
                    double candidateError = 0;
                    for (int p = 0; p < numP; p++)
                    {
                        for (int b2 = 1; b2 < numB; b2++)
                        {
                            temp[0] = biorepIntensities[p, 0];
                            temp[1] = biorepIntensities[p, b2];

                            // error metric is sum square error of coefficient of variation of each peptide
                            double coefficientOfVariation = Statistics.StandardDeviation(temp) / temp.Average();
                            candidateError += Math.Pow(coefficientOfVariation, 2);
                        }
                    }

                    return new OptimizerResult(v, candidateError);
                };

                // create optimizer
                OptimizerResult result = new NelderMeadWithStartPoints(
                    parameters: parameterArray,
                    startingValue: startPosition,
                    maxIterationsWithoutImprovement: 10
                ).OptimizeBest(minimize);

                double error = result.Error;
                double[] normFactors = result.ParameterSet;

                if (error < bestError)
                {
                    lock (locker)
                    {
                        if (error < bestError)
                        {
                            bestError = error;
                            bestNormFactors = normFactors;
                        }
                    }
                }
            });
            
            return bestNormFactors;
        }

        #endregion Private Methods
    }
}
