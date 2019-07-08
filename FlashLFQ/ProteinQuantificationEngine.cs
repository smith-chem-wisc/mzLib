using BayesianEstimation;
using MathNet.Numerics.Random;
using MathNet.Numerics.Statistics;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace FlashLFQ
{
    /// <summary>
    /// This is the Bayesian protein quantification engine used by FlashLFQ.
    /// </summary>
    public class ProteinQuantificationEngine
    {
        private readonly FlashLfqResults results;
        private readonly int maxThreads;
        private bool imputeMissing;
        private bool useSkepticalPrior;
        private bool UseSharedPeptides;
        private bool PairedSamplesOnly;
        private MersenneTwister random;
        private bool conditionsAreDefined;
        private string baseConditionString;
        private double? foldChangeCutoff;

        /// <summary>
        /// Constructs the protein quantification engine
        /// </summary>
        public ProteinQuantificationEngine(FlashLfqResults results, int maxThreads, string baseCondition, double? foldChangeCutoff = null, int? randomSeed = null)
        {
            this.maxThreads = maxThreads;
            this.results = results;
            imputeMissing = false;
            UseSharedPeptides = false;
            useSkepticalPrior = true;
            conditionsAreDefined = results.SpectraFiles.Any(p => !string.IsNullOrWhiteSpace(p.Condition));
            this.baseConditionString = baseCondition;
            this.foldChangeCutoff = foldChangeCutoff;

            if (randomSeed == null)
            {
                random = new MersenneTwister();
            }
            else
            {
                random = new MersenneTwister(randomSeed.Value);
            }

            if (!conditionsAreDefined)
            {
                SetTemporaryConditions();
            }
        }

        /// <summary>
        /// Runs the protein quantification engine
        /// </summary>
        public void Run()
        {
            HashSet<Peptide> sharedPeptides = new HashSet<Peptide>();

            // determine shared peptides
            if (!UseSharedPeptides)
            {
                sharedPeptides = new HashSet<Peptide>(results.PeptideModifiedSequences.Where(p => p.Value.proteinGroups.Count > 1).Select(p => p.Value));
            }

            // match proteins to peptides
            Dictionary<ProteinGroup, List<Peptide>> proteinsToPeptides = new Dictionary<ProteinGroup, List<Peptide>>();
            foreach (var peptide in results.PeptideModifiedSequences)
            {
                if (!peptide.Value.UseForProteinQuant || (!UseSharedPeptides && sharedPeptides.Contains(peptide.Value)))
                {
                    continue;
                }

                foreach (ProteinGroup protein in peptide.Value.proteinGroups)
                {
                    if (proteinsToPeptides.TryGetValue(protein, out var peptides))
                    {
                        peptides.Add(peptide.Value);
                    }
                    else
                    {
                        proteinsToPeptides.Add(protein, new List<Peptide> { peptide.Value });
                    }
                }
            }

            var proteinList = proteinsToPeptides.ToList();

            // generate a random seed for each protein for the MCMC sampler
            Dictionary<ProteinGroup, int> randomSeedsForEachProtein = new Dictionary<ProteinGroup, int>();
            foreach (var protein in proteinList)
            {
                randomSeedsForEachProtein.Add(protein.Key, random.NextFullRangeInt32());
            }

            // figure out which conditions to compare
            var conditions = results.SpectraFiles.Select(p => p.Condition).Distinct().ToList();
            conditions.Remove(baseConditionString);

            foreach (var condition in conditions)
            {
                var proteinToFoldChangeMeasurements = GetProteinFoldChangeMeasurements(proteinsToPeptides, condition);

                // calculate fold-change for each protein, for each condition
                Parallel.ForEach(Partitioner.Create(0, proteinList.Count), new ParallelOptions { MaxDegreeOfParallelism = maxThreads }, partitionRange =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        ProteinGroup protein = proteinList[i].Key;
                        List<Peptide> peptides = proteinList[i].Value;

                        // get the fold-change measurements for the peptides assigned to this protein
                        List<double> foldChanges = proteinToFoldChangeMeasurements[protein];

                        double mean = foldChanges.Mean();

                        // the Math.Max is here because in some edge cases the SD can be 0, which causes a crash
                        // also, the prior distribution for SD would be very small if SD < 0.001
                        double sd = Math.Max(0.001, foldChanges.StandardDeviation());

                        double nullHypothesis = foldChangeCutoff == null ? DetermineExperimentalNull() : foldChangeCutoff.Value;
                        double priorMuSd = useSkepticalPrior ? DetermineExperimentalNull() : sd * 100000;
                        double priorMuMean = useSkepticalPrior ? 0 : mean;

                        AdaptiveMetropolisWithinGibbs sampler = new AdaptiveMetropolisWithinGibbs(
                            foldChanges.ToArray(),
                            new StudentTDistributionModel(
                                priorMuMean: priorMuMean, priorMuSd: priorMuSd, muInitialGuess: mean,
                                priorSdStart: sd / 1000, priorSdEnd: sd * 1000, sdInitialGuess: sd,
                                priorNuExponent: 1.0 / 29.0, nuInitialGuess: 5),
                            seed: randomSeedsForEachProtein[protein]);

                        // burn in and then sample the MCMC chain
                        sampler.Run(2000, 2000);

                        var chain = sampler.MarkovChain;

                        // store the result
                        var result = new BayesianProteinQuantificationResult(
                            baseConditionString,
                            condition,
                            chain.Select(p => p[0]).ToArray(),
                            chain.Select(p => p[1]).ToArray(),
                            chain.Select(p => p[2]).ToArray(),
                            nullHypothesis,
                            foldChanges
                        );

                        protein.conditionToQuantificationResult.Add(condition, result);
                    }
                });
            }

            if (!conditionsAreDefined)
            {
                RemoveTemporaryConditions();
            }
        }

        private Dictionary<ProteinGroup, List<double>> GetProteinFoldChangeMeasurements(Dictionary<ProteinGroup, List<Peptide>> proteinsToPeptides, string peturbedCondition)
        {
            var proteinsWithFoldChanges = new Dictionary<ProteinGroup, List<double>>();

            var keyValuePairs = new Dictionary<(ProteinGroup protein, Peptide peptide, string condition), double[]>();

            foreach (var condition in results.SpectraFiles.GroupBy(p => p.Condition))
            {
                if (condition.Key != peturbedCondition && condition.Key != baseConditionString)
                {
                    continue;
                }

                var bioreps = condition.GroupBy(p => p.BiologicalReplicate);
                int numB = bioreps.Count();

                foreach (var biorep in bioreps)
                {
                    foreach (var fraction in biorep.GroupBy(p => p.Fraction))
                    {
                        foreach (var protein in proteinsToPeptides)
                        {
                            foreach (var peptide in protein.Value)
                            {
                                //todo: techreps
                                var fractionIntensity = peptide.GetIntensity(fraction.First(p => p.TechnicalReplicate == 0));
                                var tuple = (protein.Key, peptide, condition.Key);

                                if (keyValuePairs.ContainsKey(tuple))
                                {
                                    keyValuePairs[tuple][biorep.Key] += fractionIntensity;
                                }
                                else
                                {
                                    double[] biorepIntensities = new double[numB];
                                    biorepIntensities[biorep.Key] = fractionIntensity;

                                    keyValuePairs.Add(tuple, biorepIntensities);
                                }
                            }
                        }
                    }
                }
            }
            
            foreach (var protein in proteinsToPeptides)
            {
                proteinsWithFoldChanges.Add(protein.Key, new List<double>());
                var res = keyValuePairs.Where(p => p.Key.protein == protein.Key);

                var cond1 = res.Where(p => p.Key.condition == baseConditionString);
                var cond2 = res.Where(p => p.Key.condition == peturbedCondition);

                foreach (var peptide in protein.Value)
                {
                    List<double> peptideFoldChangeMeasurements = new List<double>();
                    var cond1Intensities = cond1.Where(p => p.Key.peptide == peptide).Select(p => p.Value).First();
                    var cond2Intensities = cond2.Where(p => p.Key.peptide == peptide).Select(p => p.Value).First();

                    int n = Math.Min(cond1Intensities.Length, cond2Intensities.Length);

                    int valid = 0;
                    for (int i = 0; i < n; i++)
                    {
                        double basemmt = Math.Log(cond1Intensities[i], 2);
                        double mmt = Math.Log(cond2Intensities[i], 2);

                        double foldChange = mmt - basemmt;

                        if (!double.IsNaN(foldChange) && !double.IsInfinity(foldChange))
                        {
                            peptideFoldChangeMeasurements.Add(foldChange);
                            valid++;
                        }
                    }

                    if (valid < n && !PairedSamplesOnly)
                    {
                        // there is at least one missing value... try to maximize the number of valid fold-change measurements
                        List<int> condition1Used = new List<int>();
                        List<int> condition2Used = new List<int>();
                        peptideFoldChangeMeasurements.Clear();

                        for (int j = 0; j < cond1Intensities.Length; j++)
                        {
                            for (int k = 0; k < cond2Intensities.Length; k++)
                            {
                                double condition1SampleIntensity = Math.Log(cond1Intensities[j], 2);
                                double condition2SampleIntensity = Math.Log(cond2Intensities[k], 2);

                                double foldChange = condition2SampleIntensity - condition1SampleIntensity;

                                if (!double.IsNaN(foldChange) && !double.IsInfinity(foldChange) && !condition1Used.Contains(j) && !condition2Used.Contains(k))
                                {
                                    peptideFoldChangeMeasurements.Add(foldChange);
                                    condition1Used.Add(j);
                                    condition2Used.Add(k);
                                }
                            }
                        }
                    }

                    proteinsWithFoldChanges[protein.Key].AddRange(peptideFoldChangeMeasurements);
                }
            }

            return proteinsWithFoldChanges;
        }

        private double DetermineExperimentalNull()
        {
            // todo: determine this from bioreps, or from conditions if bioreps are not defined
            return 0.6;
        }

        private void SetTemporaryConditions()
        {
            foreach (var file in results.SpectraFiles)
            {
                file.Condition = file.FilenameWithoutExtension;
            }
        }

        private void RemoveTemporaryConditions()
        {
            foreach (var file in results.SpectraFiles)
            {
                file.Condition = "";
            }
        }
    }
}
