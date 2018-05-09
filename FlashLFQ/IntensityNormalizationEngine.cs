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
        private FlashLFQResults results;
        private bool integrate;
        private bool silent;

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
            // run normalization functions, recalculating intensity between each function
            if (!silent)
            {
                Console.WriteLine("Normalizing techreps");
            }
            results.CalculatePeptideResults(false);
            NormalizeTechreps();

            if (!silent)
            {
                Console.WriteLine("Normalizing bioreps and fractions");
            }
            results.CalculatePeptideResults(false);
            NormalizeBiorepsAndFractions();

            if (!silent)
            {
                Console.WriteLine("Normalizing conditions");
            }
            results.CalculatePeptideResults(false);
            NormalizeConditions();
        }

        #endregion Public Methods

        #region Private Methods

        private void NormalizeTechreps()
        {
            List<double> nonZeroIntensities = new List<double>();

            var peptides = results.peptideModifiedSequences.Select(v => v.Value).ToList();
            var conditions = results.rawFiles.Select(p => p.condition).Distinct().OrderBy(p => p).ToList();

            foreach (string condition in conditions)
            {
                var filesForThisCondition = results.rawFiles.Where(p => p.condition == condition);
                var biorepsForThisCondition = filesForThisCondition.Select(p => p.biologicalReplicate).Distinct();

                foreach (var biorep in biorepsForThisCondition)
                {
                    var filesForThisBiorep = filesForThisCondition.Where(p => p.biologicalReplicate == biorep);
                    var fractionsForThisBiorep = filesForThisCondition.Select(p => p.fraction).Distinct();

                    foreach (var fraction in fractionsForThisBiorep)
                    {
                        var replicatesForThisFraction = filesForThisBiorep.Where(p => p.fraction == fraction);
                        int numReplicates = replicatesForThisFraction.Max(p => p.technicalReplicate);

                        var seqsToPeaksPerFile = new Dictionary<string, List<ChromatographicPeak>>[numReplicates + 1];
                        for(int i = 0; i < seqsToPeaksPerFile.Length; i++)
                        {
                            var rep = replicatesForThisFraction.Where(p => p.technicalReplicate == i).FirstOrDefault();

                            if (rep != null)
                            {
                                seqsToPeaksPerFile[i] = new Dictionary<string, List<ChromatographicPeak>>();

                                foreach(var peak in results.peaks[rep])
                                {
                                    if (seqsToPeaksPerFile[i].TryGetValue(peak.identifyingScans.First().ModifiedSequence, out var list))
                                        list.Add(peak);
                                    else
                                        seqsToPeaksPerFile[i].Add(peak.identifyingScans.First().ModifiedSequence, new List<ChromatographicPeak> { peak });
                                }
                            }
                        }

                        for (int p = 0; p < peptides.Count; p++)
                        {
                            double avgIntensity = 0;
                            nonZeroIntensities.Clear();

                            foreach (var technicalRep in replicatesForThisFraction)
                            {
                                if (peptides[p].intensities[technicalRep] > 0)
                                {
                                    nonZeroIntensities.Add(peptides[p].intensities[technicalRep]);
                                }
                            }
                            if (nonZeroIntensities.Any())
                            {
                                avgIntensity = nonZeroIntensities.Average();

                                foreach (var technicalRep in replicatesForThisFraction)
                                {
                                    double normFactorForThisTechRepAndPeptide = avgIntensity / peptides[p].intensities[technicalRep];

                                    if (seqsToPeaksPerFile[technicalRep.technicalReplicate].TryGetValue(peptides[p].Sequence, out var peaksForPepAndTechrep))
                                    {
                                        //TODO: set intensity to 0 if very different between technical reps
                                        if (normFactorForThisTechRepAndPeptide > 10.0 || normFactorForThisTechRepAndPeptide < 0.1)
                                        {

                                        }

                                        foreach (var peak in peaksForPepAndTechrep)
                                        {
                                            foreach (var isotopeEnvelope in peak.isotopeClusters)
                                            {
                                                isotopeEnvelope.isotopeClusterIntensity *= normFactorForThisTechRepAndPeptide;
                                            }

                                            // recalculate intensity after normalization
                                            peak.CalculateIntensityForThisFeature(integrate);
                                        }
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
            var conditions = results.rawFiles.Select(p => p.condition).Distinct().OrderBy(p => p).ToList();

            Parallel.ForEach(conditions, condition =>
            {
                var rawFilesForThisCondition = results.rawFiles.Where(p => p.condition.Equals(condition)).ToList();
                int numB = rawFilesForThisCondition.Select(p => p.biologicalReplicate).Distinct().Count();

                if (numB > 1)
                {
                    int numF = rawFilesForThisCondition.Select(p => p.fraction).Distinct().Count();
                    //int numT = rawFilesForThisCondition.Select(p => p.techrep).Distinct().Count();

                    // only normalize on peptides seen in every biorep
                    List<Peptide> seenInAllBiorepsForThisCondition = new List<Peptide>();
                    HashSet<int> bioreps = new HashSet<int>(rawFilesForThisCondition.Select(p => p.biologicalReplicate));
                    for (int p = 0; p < peptides.Count; p++)
                    {
                        int numBiorepsQuantifiedIn = 0;
                        foreach (var file in rawFilesForThisCondition)
                        {
                            if (peptides[p].intensities[file] > 0)
                            {
                                numBiorepsQuantifiedIn++;
                            }
                        }
                        if (numBiorepsQuantifiedIn == bioreps.Count)
                            seenInAllBiorepsForThisCondition.Add(peptides[p]);
                    }


                    // we're only normalizing on a subset of data here because it might take too long to do the whole set
                    seenInAllBiorepsForThisCondition = SubsetData(seenInAllBiorepsForThisCondition, condition);

                    // add the data to the array to set up for the normalization function
                    int numP = seenInAllBiorepsForThisCondition.Count;
                    double[,,] myIntensityArray = new double[numP, numB, numF];

                    for (int p = 0; p < numP; p++)
                    {
                        var peptide = seenInAllBiorepsForThisCondition[p];

                        foreach (var file in rawFilesForThisCondition)
                        {
                            myIntensityArray[p, file.biologicalReplicate, file.fraction] = peptide.intensities[file];
                        }
                    }

                    // run the normalization function (one biorep at a time, referenced to biorep 0)
                    Parallel.For(1, numB, b =>
                    {
                        var spectraFilesForThisBiorep = rawFilesForThisCondition.Where(p => p.biologicalReplicate == b);

                        double[,,] intensityArrayForTheseTwoBioreps = new double[numP, 2, numF];
                        for (int p = 0; p < numP; p++)
                        {
                            for (int f = 0; f < numF; f++)
                            {
                                intensityArrayForTheseTwoBioreps[p, 0, f] = myIntensityArray[p, 0, f];
                                intensityArrayForTheseTwoBioreps[p, 1, f] = myIntensityArray[p, b, f];
                            }
                        }

                        // solve the array to generate normalization factors
                        var normFactors = OptimizeNormalizationFactorsWithSharplearning(intensityArrayForTheseTwoBioreps, numP, 2, numF);
                        if (normFactors.All(p => p == 1.0) && !silent)
                        {
                            Console.WriteLine("Warning: Could not solve for optimal normalization factors for condition \"" + condition + "\" biorep " + b);
                        }

                        // multiply each peak's isotope envelopes by its file's normalization factor
                        foreach (var biorepFile in spectraFilesForThisBiorep)
                        {
                            if (biorepFile.biologicalReplicate == 0)
                            {
                                // since we're only normalizing one biorep at a time to the reference biorep,
                                // don't need to normalize for biorep 0, this is the reference biorep (what we're normalizing to)
                                continue;
                            }

                            double normFactorForThisFraction = normFactors[0 * numF + biorepFile.fraction];

                            foreach (var peak in results.peaks[biorepFile])
                            {
                                foreach (var isotopeEnvelope in peak.isotopeClusters)
                                {
                                    isotopeEnvelope.isotopeClusterIntensity *= normFactorForThisFraction;
                                }

                                // recalculate intensity after normalization
                                peak.CalculateIntensityForThisFeature(integrate);
                            }
                        }
                    });
                }
            });
        }

        private void NormalizeConditions()
        {
            var peptides = results.peptideModifiedSequences.Select(v => v.Value).ToList();
            var conditions = results.rawFiles.Select(p => p.condition).Distinct().OrderBy(p => p).ToList();

            // calculate intensities for the first condition (all other conditions are normalized to this condition)
            var filesForConditionOne = results.rawFiles.Where(p => p.condition.Equals(conditions[0])).ToList();

            double[] peptideIntensitiesForConditionOne = new double[peptides.Count];
            for (int p = 0; p < peptides.Count; p++)
            {
                List<double> biorepIntensities = new List<double>();

                foreach (var file in filesForConditionOne)
                {
                    if (peptides[p].intensities[file] > 0)
                        biorepIntensities.Add(peptides[p].intensities[file]);
                }

                if (biorepIntensities.Any())
                    peptideIntensitiesForConditionOne[p] = Statistics.Mean(biorepIntensities);
            }

            // calculate intensities for other conditions and normalize to the first condition
            for (int c = 1; c < conditions.Count; c++)
            {
                // calculate intensities for this condition
                var filesForThisCondition = results.rawFiles.Where(p => p.condition.Equals(conditions[c])).ToList();
                double[] peptideIntensitiesForThisCondition = new double[peptides.Count];

                for (int p = 0; p < peptides.Count; p++)
                {
                    List<double> biorepIntensities = new List<double>();

                    foreach (var file in filesForThisCondition)
                    {
                        if (peptides[p].intensities[file] > 0)
                            biorepIntensities.Add(peptides[p].intensities[file]);
                    }

                    if (biorepIntensities.Any())
                        peptideIntensitiesForThisCondition[p] = Statistics.Mean(biorepIntensities);
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
                        foreach (var isotopeEnvelope in peak.isotopeClusters)
                        {
                            isotopeEnvelope.isotopeClusterIntensity *= conditionNormalizationFactor;
                        }

                        // recalculate intensity after normalization
                        peak.CalculateIntensityForThisFeature(integrate);
                    }
                }
            }
        }

        private List<Peptide> SubsetData(List<Peptide> initialList, string condition)
        {
            int initialListCount = initialList.Count;
            HashSet<Peptide> subsetList = new HashSet<Peptide>();
            List<int> bioreps = results.rawFiles.Where(p => p.condition.Equals(condition)).Select(p => p.biologicalReplicate).Distinct().ToList();
            int maxFractionIndex = results.rawFiles.Where(p => p.condition.Equals(condition)).Max(p => p.fraction);

            RawFileInfo[,] files = new RawFileInfo[bioreps.Count, maxFractionIndex + 1];
            foreach (var file in results.rawFiles)
            {
                if (file.condition.Equals(condition))
                {
                    files[file.biologicalReplicate, file.fraction] = file;
                }
            }

            foreach (var biorep in bioreps)
            {
                List<int> fractions = results.rawFiles.Where(p => p.biologicalReplicate.Equals(biorep)).Select(p => p.fraction).Distinct().ToList();

                int numToAddPerFraction = numPeptidesDesiredInMatrix / fractions.Count;
                if (numToAddPerFraction < numPeptidesDesiredFromEachFraction)
                {
                    numToAddPerFraction = numPeptidesDesiredFromEachFraction;
                }

                int[] peptidesAddedPerFraction = new int[fractions.Count];
                Queue<Peptide>[] peptidesObservedInEachFraction = new Queue<Peptide>[fractions.Count];

                foreach (var fraction in fractions)
                {
                    peptidesObservedInEachFraction[fraction] = new Queue<Peptide>(initialList.Where(p => p.intensities[files[biorep, fraction]] > 0)
                        .OrderByDescending(p => p.intensities[files[biorep, fraction]]));
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

        private static double[] OptimizeNormalizationFactorsWithSharplearning(double[,,] proteinIntensities, int numP, int numB, int numF)
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
                    double b1Intensity = originalBiorepIntensities[p, 0];

                    for (int b2 = 1; b2 < numB; b2++)
                    {
                        double b2Intensity = originalBiorepIntensities[p, b2];
                        temp[0] = b1Intensity;
                        temp[1] = b2Intensity;

                        // calculate initial error (error if all norm factors are 1)
                        // add square error for this peptide/biorep
                        bestError += Math.Pow(Statistics.StandardDeviation(temp) / temp.Average(), 2);
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
                                    biorepIntensities[p, b] += proteinIntensities[p, b, f];
                                else
                                    biorepIntensities[p, b] += proteinIntensities[p, b, f] * v[f];
                            }
                        }
                    }

                    // calculate the error for these normalization factors
                    double candidateError = 0;
                    for (int p = 0; p < numP; p++)
                    {
                        double b1Intensity = biorepIntensities[p, 0];

                        for (int b2 = 1; b2 < numB; b2++)
                        {
                            double b2Intensity = biorepIntensities[p, b2];
                            Array.Clear(temp, 0, temp.Length);
                            temp[0] = b1Intensity;
                            temp[1] = b2Intensity;
                            double e = Math.Pow(Statistics.StandardDeviation(temp) / temp.Average(), 2);
                            candidateError += e;
                        }
                    }

                    return new OptimizerResult(v, candidateError);
                };

                // create optimizer
                var result = new NelderMeadWithStartPoints(parameterArray, 8, 0.0001, 5, 0, 0, 1, 2, -0.5, 0.5, startPosition).OptimizeBest(minimize);
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
