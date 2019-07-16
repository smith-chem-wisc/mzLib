using BayesianEstimation;
using MathNet.Numerics.Random;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace FlashLFQ
{
    /// <summary>
    /// This is the Bayesian protein quantification engine used by FlashLFQ. It uses a Markov Chain Monte Carlo (MCMC)
    /// algorithm to estimate the fold-change of a protein between conditions, given its constituent peptides fold-change measurements.
    /// This estimation is performed by fitting Student's t-distribution to the peptide fold-change data.
    /// 
    /// The posterior error probability (PEP) of a protein is the probability that the protein's fold-change is below a certain value ("cutoff").
    /// This "cutoff" can be defined by the user or can be estimated from the data using an "experimental null" strategy.
    /// The PEP is estimated using a "skeptical prior probability", which means that the initial assumption (given no data) is that the 
    /// protein is not changing in abundance between conditions. The number of data points (peptide fold-changes for a protein), the 
    /// magnitude of the fold-change relative to the cutoff, and the precision of the data points all factor in to this PEP calculation.
    /// 
    /// In other words, the default assumption is that the protein is not changing; it must be "proven" to be changing by examining the evidence.
    /// 
    /// See the BayesianEstimation project for details of the MCMC algorithm.
    /// </summary>
    public class ProteinQuantificationEngine
    {
        private readonly FlashLfqResults results;
        private readonly int MaxThreads;
        private readonly bool UseSharedPeptides;
        private readonly bool PairedSamples;
        private readonly string BaseCondition;
        private readonly double? FoldChangeCutoff;
        private readonly int BurnInSteps;
        private readonly int NSteps;
        private readonly MersenneTwister Rng;
        private bool conditionsAreDefined;
        private Dictionary<ProteinGroup, List<Peptide>> ProteinsWithConstituentPeptides;
        private Dictionary<(Peptide, string, int), double> PeptideToSampleQuantity;
        private Dictionary<ProteinGroup, double> ProteinToBaseConditionIntensity;

        /// <summary>
        /// Constructs the protein quantification engine. 
        /// 
        /// "baseCondition" refers to the condition that the protein fold-changes could be calculated relative to.
        /// 
        /// If "foldChangeCutoff" is null, the fold-change cutoff will be estimated from the data ("experimental null").
        /// 
        /// </summary>
        public ProteinQuantificationEngine(FlashLfqResults results, int maxThreads, string baseCondition, bool useSharedPeptides = false,
            double? foldChangeCutoff = null, int? randomSeed = null, int mcmcBurninSteps = 1000, int mcmcSteps = 3000, bool pairedSamples = false)
        {
            if (string.IsNullOrWhiteSpace(baseCondition))
            {
                throw new MzLibException("The protein quantification base condition must be defined");
            }

            this.MaxThreads = maxThreads;
            this.results = results;
            UseSharedPeptides = useSharedPeptides;
            this.BaseCondition = baseCondition;
            this.FoldChangeCutoff = foldChangeCutoff;
            this.BurnInSteps = mcmcBurninSteps;
            this.NSteps = mcmcSteps;
            this.PairedSamples = pairedSamples;

            conditionsAreDefined = results.SpectraFiles.All(p => !string.IsNullOrWhiteSpace(p.Condition));

            if (randomSeed == null)
            {
                Rng = new MersenneTwister();
            }
            else
            {
                Rng = new MersenneTwister(randomSeed.Value);
            }
        }

        /// <summary>
        /// Runs the protein quantification engine.
        /// </summary>
        public void Run()
        {
            Setup();

            // generate a random seed for each protein for the MCMC sampler
            var randomSeedsForEachProtein = new Dictionary<ProteinGroup, List<int>>();
            foreach (var protein in ProteinsWithConstituentPeptides)
            {
                randomSeedsForEachProtein.Add(protein.Key, new List<int> { Rng.NextFullRangeInt32(), Rng.NextFullRangeInt32() });
            }

            var proteinList = ProteinsWithConstituentPeptides.ToList();

            // figure out which conditions to compare
            var conditions = results.SpectraFiles.Select(p => p.Condition).Distinct().ToList();
            conditions.Remove(BaseCondition);

            double baseConditionExperimentalNull = 0;
            if (FoldChangeCutoff == null)
            {
                baseConditionExperimentalNull = DetermineExperimentalNull(BaseCondition);
            }

            // calculate fold-change for each protein
            foreach (string condition in conditions)
            {
                double nullHyp = FoldChangeCutoff == null ?
                    Math.Max(DetermineExperimentalNull(condition), baseConditionExperimentalNull) : FoldChangeCutoff.Value;

                Parallel.ForEach(Partitioner.Create(0, proteinList.Count), new ParallelOptions { MaxDegreeOfParallelism = MaxThreads }, partitionRange =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        ProteinGroup protein = proteinList[i].Key;

                        // get the fold-change measurements for the peptides assigned to this protein
                        List<double> peptideFoldChanges = GetPeptideFoldChangeMeasurements(protein, BaseCondition, condition, null, null);

                        ProteinQuantificationEngineResult result;

                        if (peptideFoldChanges.Count > 1)
                        {
                            // run the Bayesian analysis
                            result = RunBayesianProteinQuant(peptideFoldChanges, protein, condition,
                                randomSeedsForEachProtein[protein][0], null, false, BurnInSteps, NSteps, out var mus);

                            var skepticalResult = RunBayesianProteinQuant(peptideFoldChanges, protein, condition,
                                randomSeedsForEachProtein[protein][1], nullHyp, true, BurnInSteps, NSteps, out var skepticalMus);

                            result.cutoff = nullHyp;
                            result.CalculatePosteriorErrorProbability(skepticalMus);
                        }
                        else if (peptideFoldChanges.Count == 1)
                        {
                            result = new ProteinQuantificationEngineResult(protein, BaseCondition, condition,
                                new double[] { peptideFoldChanges.First() }, new double[] { double.NaN },
                                new double[] { double.NaN }, peptideFoldChanges, ProteinToBaseConditionIntensity[protein]);

                            result.cutoff = nullHyp;
                            result.CalculatePosteriorErrorProbability(new double[] { nullHyp });
                        }
                        else
                        {
                            result = new ProteinQuantificationEngineResult(protein, BaseCondition, condition,
                                new double[] { 0.0 }, new double[] { double.NaN }, new double[] { double.NaN },
                                peptideFoldChanges, ProteinToBaseConditionIntensity[protein]);

                            result.cutoff = nullHyp;
                            result.CalculatePosteriorErrorProbability(new double[] { nullHyp });
                        }

                        protein.conditionToQuantificationResults.Add(condition, result);
                    }
                });
            }

            // calculate FDR for the condition
            foreach (var condition in conditions)
            {
                var res = proteinList.Select(p => p.Key.conditionToQuantificationResults[condition]).ToList();
                CalculateFalseDiscoveryRates(res);
            }

            // remove temporary conditions
            if (!conditionsAreDefined)
            {
                foreach (var file in results.SpectraFiles)
                {
                    file.Condition = "";
                }
            }
        }

        /// <summary>
        /// Pairs up proteins with constituent peptides and calculates peptide biorep intensities for easy lookup.
        /// </summary>
        private void Setup()
        {
            // temporarily set conditions if they aren't defined
            if (!conditionsAreDefined)
            {
                foreach (var file in results.SpectraFiles)
                {
                    file.Condition = file.FullFilePathWithExtension;
                }
            }

            HashSet<Peptide> sharedPeptides = new HashSet<Peptide>();

            // determine shared peptides
            if (!UseSharedPeptides)
            {
                sharedPeptides = new HashSet<Peptide>(results.PeptideModifiedSequences.Where(p => p.Value.proteinGroups.Count > 1).Select(p => p.Value));
            }

            // match proteins to peptides
            ProteinsWithConstituentPeptides = new Dictionary<ProteinGroup, List<Peptide>>();
            foreach (var peptide in results.PeptideModifiedSequences)
            {
                if (!peptide.Value.UseForProteinQuant || (!UseSharedPeptides && sharedPeptides.Contains(peptide.Value)))
                {
                    continue;
                }

                foreach (ProteinGroup protein in peptide.Value.proteinGroups)
                {
                    if (ProteinsWithConstituentPeptides.TryGetValue(protein, out var peptides))
                    {
                        peptides.Add(peptide.Value);
                    }
                    else
                    {
                        ProteinsWithConstituentPeptides.Add(protein, new List<Peptide> { peptide.Value });
                    }
                }
            }

            // calculate peptide biorep intensities
            var allpeptides = ProteinsWithConstituentPeptides.Values.SelectMany(p => p).Distinct().ToList();
            PeptideToSampleQuantity = new Dictionary<(Peptide, string, int), double>();
            foreach (var condition in results.SpectraFiles.GroupBy(p => p.Condition))
            {
                foreach (var biorep in condition.GroupBy(p => p.BiologicalReplicate))
                {
                    foreach (Peptide peptide in allpeptides)
                    {
                        double biorepIntensity = 0;

                        foreach (var fraction in biorep.GroupBy(p => p.Fraction))
                        {
                            double fractionIntensity = 0;
                            int nonZeroTechrepCount = 0;

                            foreach (var techrep in fraction.GroupBy(p => p.TechnicalReplicate))
                            {
                                double techrepIntensity = peptide.GetIntensity(techrep.First());

                                if (techrepIntensity > 0)
                                {
                                    fractionIntensity += techrepIntensity;
                                    nonZeroTechrepCount++;
                                }
                            }

                            if (nonZeroTechrepCount > 0)
                            {
                                fractionIntensity /= nonZeroTechrepCount;
                            }

                            biorepIntensity += fractionIntensity;
                        }

                        PeptideToSampleQuantity.Add((peptide, condition.Key, biorep.Key), biorepIntensity);
                    }
                }
            }

            // calculate reference intensities
            ProteinToBaseConditionIntensity = new Dictionary<ProteinGroup, double>();
            var numBRef = results.SpectraFiles.Where(p => p.Condition == BaseCondition).Select(p => p.BiologicalReplicate).Distinct().Count();
            foreach (var protein in ProteinsWithConstituentPeptides)
            {
                List<double> intensities = new List<double>();

                foreach (Peptide peptide in protein.Value)
                {
                    double avgIntensity = 0;
                    int n = 0;

                    for (int b = 0; b < numBRef; b++)
                    {
                        double intensity = PeptideToSampleQuantity[(peptide, BaseCondition, b)];
                        avgIntensity += intensity;

                        if (intensity > 0)
                        {
                            n++;
                        }
                    }

                    if (n == 0)
                    {
                        continue;
                    }

                    avgIntensity /= n;
                    intensities.Add(avgIntensity);
                }

                double referenceProteinIntensity = 0;

                if (intensities.Any())
                {
                    double top3 = intensities.OrderByDescending(p => p).Take(3).Sum();
                    referenceProteinIntensity = top3;
                }

                ProteinToBaseConditionIntensity.Add(protein.Key, referenceProteinIntensity);
            }
        }

        /// <summary>
        /// Computes a list of fold-change measurements between the constituent peptides of a protein.
        /// 
        /// One of the following must be true:
        /// ---> condition1 and condition2 are different, and biorep1 and biorep2 are null.
        /// (This returns a list of peptide fold-changes between the conditions.)
        /// 
        /// ---> biorep1 and biorep2 are not null.
        /// (This returns a list of peptide fold-changes between the defined bioreps of (a) condition(s).)
        /// </summary>
        private List<double> GetPeptideFoldChangeMeasurements(ProteinGroup protein, string condition1, string condition2, int? biorep1, int? biorep2)
        {
            List<double> allPeptideFoldChanges = new List<double>();

            List<Peptide> peptides = ProteinsWithConstituentPeptides[protein];

            int numB1 = results.SpectraFiles.Where(p => p.Condition == condition1).Max(p => p.BiologicalReplicate) + 1;
            int numB2 = results.SpectraFiles.Where(p => p.Condition == condition2).Max(p => p.BiologicalReplicate) + 1;

            List<double> peptideFoldChanges = new List<double>();
            foreach (var peptide in peptides)
            {
                peptideFoldChanges.Clear();

                // fold-changes between different conditions
                if (biorep1 == null && biorep2 == null && condition1 != condition2)
                {
                    double[] cond1Intensities = new double[numB1];
                    double[] cond2Intensities = new double[numB2];

                    for (int b = 0; b < cond1Intensities.Length; b++)
                    {
                        cond1Intensities[b] = PeptideToSampleQuantity[(peptide, condition1, b)];
                    }

                    for (int b = 0; b < cond2Intensities.Length; b++)
                    {
                        cond2Intensities[b] = PeptideToSampleQuantity[(peptide, condition2, b)];
                    }

                    int n = Math.Min(cond1Intensities.Length, cond2Intensities.Length);

                    int valid = 0;
                    for (int i = 0; i < n; i++)
                    {
                        double? foldChange = GetLogFoldChange(cond1Intensities[i], cond2Intensities[i]);

                        if (foldChange != null)
                        {
                            peptideFoldChanges.Add(foldChange.Value);
                            valid++;
                        }
                    }
                    if (valid < n && !PairedSamples)
                    {
                        // there is at least one missing value... try to maximize the number of valid fold-change measurements
                        List<int> condition1Used = new List<int>();
                        List<int> condition2Used = new List<int>();
                        peptideFoldChanges.Clear();

                        for (int j = 0; j < cond1Intensities.Length; j++)
                        {
                            for (int k = 0; k < cond2Intensities.Length; k++)
                            {
                                double? foldChange = GetLogFoldChange(cond1Intensities[j], cond2Intensities[k]);

                                if (foldChange != null && !condition1Used.Contains(j) && !condition2Used.Contains(k))
                                {
                                    peptideFoldChanges.Add(foldChange.Value);
                                    condition1Used.Add(j);
                                    condition2Used.Add(k);
                                }
                            }
                        }
                    }
                }
                // fold-changes between two bioreps of the same condition
                else if (biorep1 != null && biorep2 != null && condition1 == condition2)
                {
                    double intensity1 = PeptideToSampleQuantity[(peptide, condition1, biorep1.Value)];
                    double intensity2 = PeptideToSampleQuantity[(peptide, condition2, biorep2.Value)];

                    double? lfc = GetLogFoldChange(intensity1, intensity2);

                    if (lfc != null)
                    {
                        peptideFoldChanges.Add(lfc.Value);
                    }
                }
                else
                {
                    throw new MzLibException("Problem getting peptide fold-change measurements for protein " + protein.ProteinGroupName);
                }

                allPeptideFoldChanges.AddRange(peptideFoldChanges);
            }

            return allPeptideFoldChanges;
        }

        /// <summary>
        /// Computes the log-fold change between two intensities. If there is an error (e.g., one of the intensities is zero), 
        /// null is returned.
        /// </summary>
        private double? GetLogFoldChange(double intensity1, double intensity2)
        {
            double logFoldChange = Math.Log(intensity2, 2) - Math.Log(intensity1, 2);

            if (!double.IsNaN(logFoldChange) && !double.IsInfinity(logFoldChange))
            {
                return logFoldChange;
            }

            return null;
        }

        /// <summary>
        /// Estimates the fold-change between conditions for a protein, given its constituent peptides' fold change measurements. 
        /// If "skepticalPrior" is set to true, then the result's estimate of the fold-change between conditions is drawn towards 
        /// zero. This is useful for estimating the posterior error probability (PEP) that a fold-change is below a certain cutoff.
        /// </summary>
        private ProteinQuantificationEngineResult RunBayesianProteinQuant(List<double> foldChanges, ProteinGroup protein, string condition, int randomSeed,
            double? nullHypothesisCutoff, bool skepticalPrior, int burnin, int n, out double[] mus)
        {
            if (skepticalPrior && nullHypothesisCutoff == null)
            {
                throw new MzLibException("Invalid protein quantification parameters given for protein " + protein.ProteinGroupName);
            }

            // the Math.Max is here because in some edge cases the SD can be 0, which causes a crash
            // also, the prior distribution for SD would be very narrow if SD < 0.001
            double sd = Math.Max(0.001, foldChanges.StandardDeviation());
            double mean = foldChanges.Mean();

            // the Math.max is here for cases where the cutoff is very small (e.g., 0)
            double priorMuSd = Math.Max(0.05, skepticalPrior ? nullHypothesisCutoff.Value : sd * 100);
            double priorMuMean = skepticalPrior ? 0 : mean;

            AdaptiveMetropolisWithinGibbs sampler = new AdaptiveMetropolisWithinGibbs(
                foldChanges.ToArray(),
                new ProteinFoldChangeEstimationModel(
                    priorMuMean: priorMuMean,
                    priorMuSd: priorMuSd,
                    muInitialGuess: mean,
                    priorSdStart: sd / 1000,
                    priorSdEnd: sd * 1000,
                    sdInitialGuess: sd,
                    priorNuExponent: 1.0 / 29.0,
                    nuInitialGuess: 5),
                    seed: randomSeed
                );

            // burn in and then sample the MCMC chain
            sampler.Run(burnin, n);

            mus = sampler.MarkovChain.Select(p => p[0]).ToArray();
            double[] sds = sampler.MarkovChain.Select(p => p[1]).ToArray();
            double[] nus = sampler.MarkovChain.Select(p => p[2]).ToArray();

            // store the result
            var result = new ProteinQuantificationEngineResult(
                protein,
                BaseCondition,
                condition,
                mus,
                sds,
                nus,
                foldChanges,
                ProteinToBaseConditionIntensity[protein]
            );

            return result;
        }

        /// <summary>
        /// This function determines a fold-change cutoff from the data. The cutoff is the estimate of the biological variance. 
        /// If bioreps are defined, the standard deviation of the protein fold-change between bioreps within a condition is used 
        /// as the cutoff. If the condition only has 1 biorep, then the cutoff is the standard deviation of the protein fold-change 
        /// between conditions.
        /// </summary>
        private double DetermineExperimentalNull(string treatmentCondition)
        {
            double experimentalNull;

            int biorepCount = results.SpectraFiles.Where(p => p.Condition == treatmentCondition).Select(p => p.BiologicalReplicate).Distinct().Count();

            if ((PairedSamples || biorepCount == 1) && treatmentCondition == BaseCondition)
            {
                return 0;
            }

            var proteinList = ProteinsWithConstituentPeptides.ToList();

            // generate a random seed for each protein for the MCMC sampler
            var randomSeedsForEachProtein = new Dictionary<ProteinGroup, List<int>>();
            foreach (var protein in ProteinsWithConstituentPeptides)
            {
                randomSeedsForEachProtein.Add(protein.Key, new List<int> { Rng.NextFullRangeInt32() });
            }

            var experimentalNullResults = new Dictionary<string, List<ProteinQuantificationEngineResult>>();

            Parallel.ForEach(Partitioner.Create(0, proteinList.Count), new ParallelOptions { MaxDegreeOfParallelism = MaxThreads }, partitionRange =>
            {
                for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                {
                    ProteinGroup protein = proteinList[i].Key;

                    List<double> peptideFoldChanges = new List<double>();

                    if (biorepCount > 1)
                    {
                        if (!PairedSamples)
                        {
                            // get the fold-change measurements for the peptides assigned to this protein between bioreps of this condition
                            for (int b = 0; b < biorepCount - 1; b++)
                            {
                                peptideFoldChanges.AddRange(GetPeptideFoldChangeMeasurements(protein, treatmentCondition, treatmentCondition, b, b + 1));
                            }
                        }
                        else if (PairedSamples)
                        {
                            peptideFoldChanges.AddRange(GetPeptideFoldChangeMeasurements(protein, BaseCondition, treatmentCondition, null, null));
                        }
                    }
                    else
                    {
                        peptideFoldChanges.AddRange(GetPeptideFoldChangeMeasurements(protein, BaseCondition, treatmentCondition, null, null));
                    }

                    ProteinQuantificationEngineResult result;

                    if (peptideFoldChanges.Count > 1)
                    {
                        // run the Bayesian analysis
                        result = RunBayesianProteinQuant(peptideFoldChanges, protein, treatmentCondition,
                            randomSeedsForEachProtein[protein][0], null, false, 500, 500, out var mus);

                        lock (experimentalNullResults)
                        {
                            experimentalNullResults.Add(protein.ProteinGroupName, new List<ProteinQuantificationEngineResult> { result });
                        }
                    }
                    else
                    {
                        continue;
                    }
                }
            });

            double stdDev = 0;

            // this is a weird case where there are very few proteins in the data...
            if (experimentalNullResults.Count < 10)
            {
                if (experimentalNullResults.Any())
                {
                    stdDev = experimentalNullResults.SelectMany(p => p.Value).Select(p => p.StandardDeviationPointEstimate).Average();
                }
                else
                {
                    // TODO: not sure what to do here.. there were no protein quant results for this condition
                    // throw an exception... or just make the cutoff 1.0?
                    stdDev = 0.3333;

                    //throw new MzLibException("Could not determine experimental null for condition " + treatmentCondition);
                }
            }
            else
            {
                stdDev = experimentalNullResults.SelectMany(p => p.Value).Select(p => p.FoldChangePointEstimate).InterquartileRange() / 1.34896;
            }

            experimentalNull = 3.0 * stdDev;

            return experimentalNull;
        }

        /// <summary>
        /// Calculates the false discovery rate of each protein from the Bayesian-estimated PEP values.
        /// </summary>
        private void CalculateFalseDiscoveryRates(List<ProteinQuantificationEngineResult> bayesianProteinQuantificationResults)
        {
            var bayesianQuantResults = bayesianProteinQuantificationResults
                .OrderByDescending(p => p.DegreeOfEvidence)
                .ThenByDescending(p => p.foldChangeMeasurements.Count)
                .ToList();

            double runningPEP = 0;
            double qValue = 0;

            for (int p = 0; p < bayesianQuantResults.Count; p++)
            {
                runningPEP += bayesianQuantResults[p].PosteriorErrorProbability;
                qValue = runningPEP / (p + 1);
                bayesianQuantResults[p].FalseDiscoveryRate = qValue;
            }
        }
    }
}
