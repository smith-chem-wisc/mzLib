using FlashLFQ.BoundedNelderMeadOptimizer;
using MathNet.Numerics.Statistics;
using SharpLearning.Optimization;
using System;
using System.Collections.Generic;
using System.Linq;

namespace FlashLFQ
{
    public class IntensityNormalizationEngine
    {
        private const int numPeptidesDesiredFromEachFraction = 500;
        private const int numPeptidesDesiredInMatrix = 5000;
        private readonly FlashLfqResults results;
        private readonly bool integrate;
        private readonly bool silent;
        private readonly bool quantifyAmbiguousPeptides;
        private readonly int maxThreads;

        public IntensityNormalizationEngine(FlashLfqResults results, bool integrate, bool silent, int maxThreads, bool quantifyAmbiguousPeptides = false)
        {
            this.results = results;
            this.integrate = integrate;
            this.silent = silent;
            this.maxThreads = maxThreads;
        }

        /// <summary>
        /// Runs the normalization functions.
        /// </summary>
        public void NormalizeResults()
        {
            results.CalculatePeptideResults(quantifyAmbiguousPeptides);

            // run normalization functions, recalculating intensity between each function
            if (!silent)
            {
                Console.WriteLine("Normalizing fractions");
            }
            NormalizeFractions();
            results.CalculatePeptideResults(quantifyAmbiguousPeptides);

            if (!silent)
            {
                Console.WriteLine("Normalizing bioreps and conditions");
            }
            NormalizeBioreps();
            results.CalculatePeptideResults(quantifyAmbiguousPeptides);

            if (!silent)
            {
                Console.WriteLine("Normalizing techreps");
            }
            NormalizeTechreps();
            results.CalculatePeptideResults(quantifyAmbiguousPeptides);
        }

        /// <summary>
        /// This method normalizes peptide intensities so that the median fold-change between technical replicates
        /// is zero. The median is used instead of the average because it is more robust to outliers (i.e., if there are
        /// many changing peptides, the median will perform better than using the average).
        /// </summary>
        private void NormalizeTechreps()
        {
            var peptides = results.PeptideModifiedSequences.Select(v => v.Value).ToList();
            var conditions = results.SpectraFiles.GroupBy(v => v.Condition);

            foreach (var condition in conditions)
            {
                var bioreps = condition.GroupBy(v => v.BiologicalReplicate);

                foreach (var biorep in bioreps)
                {
                    var fractions = biorep.GroupBy(v => v.Fraction);

                    foreach (var fraction in fractions)
                    {
                        var techreps = fraction.ToList();

                        for (int t = 1; t < techreps.Count; t++)
                        {
                            List<double> foldChanges = new List<double>();

                            for (int p = 0; p < peptides.Count; p++)
                            {
                                double techrep1Intensity = peptides[p].GetIntensity(techreps[0]);
                                double techrepTIntensity = peptides[p].GetIntensity(techreps[t]);

                                if (techrep1Intensity > 0 && techrepTIntensity > 0)
                                {
                                    foldChanges.Add(techrepTIntensity / techrep1Intensity);
                                }
                            }

                            if (!foldChanges.Any())
                            {
                                // TODO: throw an exception?
                                return;
                            }

                            double medianFoldChange = foldChanges.Median();
                            double normalizationFactor = 1.0 / medianFoldChange;

                            // normalize to median fold-change
                            foreach (var peak in results.Peaks[techreps[t]])
                            {
                                foreach (var isotopeEnvelope in peak.IsotopicEnvelopes)
                                {
                                    isotopeEnvelope.Normalize(normalizationFactor);
                                }

                                // recalculate intensity after normalization
                                peak.CalculateIntensityForThisFeature(integrate);
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        /// This method uses a bounded Nelder-Mead optimization algorithm to find linear normalization factors
        /// so that coefficient of variation of peptide intensity between two biological replicates will be minimized.
        /// Calls "GetNormalizationFactors" to calculate the normalization factors.
        /// </summary>
        private void NormalizeFractions()
        {
            if (results.SpectraFiles.Max(p => p.Fraction) == 0)
            {
                return;
            }

            var peptides = results.PeptideModifiedSequences.Select(v => v.Value).ToList();
            var conditions = results.SpectraFiles.Select(p => p.Condition).Distinct().OrderBy(p => p).ToList();
            var filesForCond1Biorep1 = results.SpectraFiles.Where(p => p.Condition == conditions[0] && p.BiologicalReplicate == 0 && p.TechnicalReplicate == 0).ToList();

            foreach (var condition in conditions)
            {
                var filesForThisCondition = results.SpectraFiles.Where(p => p.Condition.Equals(condition)).ToList();

                int numB = filesForThisCondition.Select(p => p.BiologicalReplicate).Distinct().Count();

                for (int b = 0; b < numB; b++)
                {
                    // condition 1 biorep 1 is the reference, don't normalize it
                    if (b == 0 && conditions.IndexOf(condition) == 0)
                    {
                        continue;
                    }

                    // run the normalization function
                    if (!silent)
                    {
                        Console.WriteLine("Normalizing condition \"" + condition + "\" biorep " + (b + 1));
                    }

                    var filesForThisBiorep = filesForThisCondition.Where(p => p.BiologicalReplicate == b && p.TechnicalReplicate == 0);

                    int numF = Math.Max(filesForCond1Biorep1.Max(p => p.Fraction), filesForThisBiorep.Max(p => p.Fraction)) + 1;

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
                    //seenInBothBioreps = SubsetData(seenInBothBioreps, filesForThisBiorep.Concat(filesForCond1Biorep1).ToList());

                    // add the data to the array to set up for the normalization function
                    int numP = seenInBothBioreps.Count;
                    double[,,] myIntensityArray = new double[numP, 2, numF];

                    for (int p = 0; p < numP; p++)
                    {
                        var peptide = seenInBothBioreps[p];

                        foreach (var file in filesForCond1Biorep1)
                        {
                            myIntensityArray[p, 0, file.Fraction] = peptide.GetIntensity(file);
                        }

                        foreach (var file in filesForThisBiorep)
                        {
                            myIntensityArray[p, 1, file.Fraction] = peptide.GetIntensity(file);
                        }
                    }

                    // solve for normalization factors
                    var normFactors = GetNormalizationFactors(myIntensityArray, numP, numF);
                    if (normFactors.All(p => p == 1.0) && !silent)
                    {
                        Console.WriteLine("Warning: Could not solve for optimal normalization factors for condition \"" + condition + "\" biorep " + (b + 1));
                    }

                    // multiply each peak's isotope envelopes by its file's normalization factor
                    foreach (var spectraFile in filesForThisBiorep)
                    {
                        foreach (var peak in results.Peaks[spectraFile])
                        {
                            foreach (var isotopeEnvelope in peak.IsotopicEnvelopes)
                            {
                                isotopeEnvelope.Normalize(normFactors[spectraFile.Fraction]);
                            }

                            // recalculate intensity after normalization
                            peak.CalculateIntensityForThisFeature(integrate);
                        }
                    }
                }
            }
        }

        /// <summary>
        /// This method normalizes peptide intensities so that the median fold-change between any two biological replicates
        /// (regardless of condition) is ~zero. The median is used instead of the average because it is more robust to outliers.
        /// The assumption in this method is that the median fold-change between bioreps of different conditions
        /// is zero (i.e., that most peptides do not change in abundance between conditions).
        /// </summary>
        private void NormalizeBioreps()
        {
            var peptides = results.PeptideModifiedSequences.Select(v => v.Value).ToList();
            var conditions = results.SpectraFiles.GroupBy(v => v.Condition).ToList();

            double[,] biorepIntensityPair = new double[peptides.Count, 2];

            var firstConditionFirstBiorep = conditions.First().Where(v => v.BiologicalReplicate == 0 && v.TechnicalReplicate == 0);

            foreach (var file in firstConditionFirstBiorep)
            {
                for (int p = 0; p < peptides.Count; p++)
                {
                    biorepIntensityPair[p, 0] += peptides[p].GetIntensity(file);
                }
            }

            foreach (var condition in conditions)
            {
                var bioreps = condition.GroupBy(v => v.BiologicalReplicate);

                foreach (var biorep in bioreps)
                {
                    for (int p = 0; p < peptides.Count; p++)
                    {
                        biorepIntensityPair[p, 1] = 0;
                    }

                    var fractions = biorep.GroupBy(v => v.Fraction);

                    foreach (var fraction in fractions)
                    {
                        var firstTechrep = fraction.Where(v => v.TechnicalReplicate == 0).First();

                        for (int p = 0; p < peptides.Count; p++)
                        {
                            biorepIntensityPair[p, 1] += peptides[p].GetIntensity(firstTechrep);
                        }
                    }

                    List<double> foldChanges = new List<double>();

                    for (int p = 0; p < peptides.Count; p++)
                    {
                        if (biorepIntensityPair[p, 0] > 0 && biorepIntensityPair[p, 1] > 0)
                        {
                            foldChanges.Add(biorepIntensityPair[p, 1] / biorepIntensityPair[p, 0]);
                        }
                    }

                    if (!foldChanges.Any())
                    {
                        // TODO: throw an exception?
                        return;
                    }

                    double medianFoldChange = foldChanges.Median();
                    double normalizationFactor = 1.0 / medianFoldChange;

                    // normalize to median fold-change
                    foreach (var file in biorep)
                    {
                        foreach (var peak in results.Peaks[file])
                        {
                            foreach (var isotopeEnvelope in peak.IsotopicEnvelopes)
                            {
                                isotopeEnvelope.Normalize(normalizationFactor);
                            }

                            // recalculate intensity after normalization
                            peak.CalculateIntensityForThisFeature(integrate);
                        }
                    }
                }
            }
        }

        /// <summary>
        /// This method takes a list of peptides and creates a subset list of peptides to normalize with, to avoid
        /// excessive computation time in normalization functions.
        /// </summary>
        //private List<Peptide> SubsetData(List<Peptide> initialList, List<SpectraFileInfo> spectraFiles)
        //{
        //    List<SpectraFileInfo>[] bothBioreps = new List<SpectraFileInfo>[2];
        //    var temp1 = spectraFiles.GroupBy(p => p.Condition).ToList();
        //    if (temp1.Count() == 1)
        //    {
        //        // normalizing bioreps within a condition
        //        var temp2 = spectraFiles.GroupBy(p => p.BiologicalReplicate).ToList();
        //        bothBioreps[0] = temp2[0].ToList();
        //        bothBioreps[1] = temp2[1].ToList();
        //    }
        //    else
        //    {
        //        // normalizing bioreps between conditions
        //        bothBioreps[0] = temp1[0].ToList();
        //        bothBioreps[1] = temp1[1].ToList();
        //    }

        //    HashSet<Peptide> subsetList = new HashSet<Peptide>();
        //    int maxFractionIndex = bothBioreps.SelectMany(p => p).Max(v => v.Fraction);

        //    foreach (var biorep in bothBioreps)
        //    {
        //        List<int> fractions = biorep.Select(p => p.Fraction).Distinct().ToList();

        //        int numToAddPerFraction = numPeptidesDesiredInMatrix / fractions.Count;
        //        if (numToAddPerFraction < numPeptidesDesiredFromEachFraction)
        //        {
        //            numToAddPerFraction = numPeptidesDesiredFromEachFraction;
        //        }

        //        int[] peptidesAddedPerFraction = new int[fractions.Count];
        //        Queue<Peptide>[] peptidesObservedInEachFraction = new Queue<Peptide>[fractions.Count];

        //        foreach (var file in biorep)
        //        {
        //            if (peptidesObservedInEachFraction[file.Fraction] == null)
        //            {
        //                peptidesObservedInEachFraction[file.Fraction] = new Queue<Peptide>(initialList.Where(p => p.GetIntensity(file) > 0)
        //                    .OrderByDescending(p => p.GetIntensity(file)));
        //            }
        //        }

        //        foreach (var fraction in fractions)
        //        {
        //            while (peptidesAddedPerFraction[fraction] < numToAddPerFraction && peptidesObservedInEachFraction[fraction].Any())
        //            {
        //                var peptide = peptidesObservedInEachFraction[fraction].Dequeue();

        //                // don't need to check if the return list already contains the peptide because it's a HashSet (no duplicates are allowed)
        //                subsetList.Add(peptide);

        //                // this peptide is in the return list regardless of whether or not it was actually just added;
        //                // we just want to guarantee this fraction has 500 peptides in the return list to normalize with
        //                peptidesAddedPerFraction[fraction]++;
        //            }
        //        }
        //    }

        //    return subsetList.ToList();
        //}

        /// <summary>
        /// Calculates normalization factors for fractionated data using a bounded Nelder-Mead optimization algorithm.
        /// Called by NormalizeFractions().
        /// </summary>
        private static double[] GetNormalizationFactors(double[,,] peptideIntensities, int numP, int numF)
        {
            object locker = new object();

            double[] referenceSample = new double[numP];
            double[,] sampleToNormalize = new double[numP, numF];

            // populate the peptide sample quantity array for normalization calculations
            for (int p = 0; p < numP; p++)
            {
                for (int f = 0; f < numF; f++)
                {
                    referenceSample[p] += peptideIntensities[p, 0, f];
                    sampleToNormalize[p, f] = peptideIntensities[p, 1, f];
                }
            }

            // initialize normalization factors to 1.0
            // normalization factor optimization must improve on these to be valid
            double[] bestNormFactors = new double[numF];
            for (int i = 0; i < bestNormFactors.Length; i++)
            {
                bestNormFactors[i] = 1.0;
            }

            // calculate the error between bioreps if all normalization factors are 1 (initial error)
            double[] initialErrors = new double[numP];
            double bestError = CalculateNormalizationFactorError(ref referenceSample, ref sampleToNormalize, ref bestNormFactors);

            // constraint (normalization factors must be >0.3 and <3
            var parameterArray = new ParameterBounds[numF];
            for (int f = 0; f < numF; f++)
            {
                parameterArray[f] = new ParameterBounds(0.3, 3, Transform.Linear);
            }

            // find approximate best starting area for each fraction normalization factor
            for (int f = 0; f < numF; f++)
            {
                double bestFractionError = double.PositiveInfinity;
                double start = parameterArray[0].Min;
                double end = parameterArray[0].Max;
                double[] factors = new double[numF];
                for (int i = 0; i < factors.Length; i++)
                {
                    factors[i] = 1.0;
                }

                for (double n = start; n <= end; n += 0.01)
                {
                    factors[f] = Math.Round(n, 2);

                    double error = CalculateNormalizationFactorError(ref referenceSample, ref sampleToNormalize, ref factors);

                    if (error < bestFractionError)
                    {
                        bestFractionError = error;
                        bestNormFactors[f] = factors[f];
                    }
                }
            }

            // find the best normalization factors (minimize error)
            double[] errors = new double[numP];

            // define minimization metric
            Func<double[], OptimizerResult> minimize = v =>
            {
                // calculate error with these normalization factors
                double candidateError = CalculateNormalizationFactorError(ref referenceSample, ref sampleToNormalize, ref v);

                return new OptimizerResult(v, candidateError);
            };

            // create optimizer
            OptimizerResult result = new NelderMeadWithStartPoints(
                parameters: parameterArray,
                startingValue: bestNormFactors,
                maxIterationsWithoutImprovement: 10
            ).OptimizeBest(minimize);

            double sampleError = result.Error;
            double[] normalizationFactors = result.ParameterSet;

            if (sampleError < bestError)
            {
                lock (locker)
                {
                    if (sampleError < bestError)
                    {
                        bestError = sampleError;
                        bestNormFactors = normalizationFactors;
                    }
                }
            }

            return bestNormFactors;
        }

        private static double CalculateNormalizationFactorError(ref double[] reference, ref double[,] sampleToNormalize, ref double[] normalizationFactors)
        {
            int numP = reference.Length;
            int numF = normalizationFactors.Length;

            double totalError = 0;
            
            for (int p = 0; p < numP; p++)
            {
                // sum the intensities with the current normalization factors
                double normalizedReplicateIntensity = 0;
                for (int f = 0; f < numF; f++)
                {
                    normalizedReplicateIntensity += sampleToNormalize[p, f] * normalizationFactors[f];
                }
                
                double peptideError = Math.Log(normalizedReplicateIntensity) - Math.Log(reference[p]);

                totalError += peptideError;
            }

            return totalError;
        }
    }
}