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
        private readonly FlashLfqResults Results;
        private readonly int MaxThreads;
        private readonly bool UseSharedPeptides;
        private readonly bool PairedSamples;
        private readonly string ControlCondition;
        private readonly double? NullHypothesisInterval;
        private readonly int BurnInSteps;
        private readonly int McmcSteps;
        private readonly MersenneTwister Rng;
        private Dictionary<ProteinGroup, List<Peptide>> ProteinsWithConstituentPeptides;
        private Dictionary<(Peptide, string, int), double> PeptideToSampleQuantity;
        public readonly int RandomSeed;

        /// <summary>
        /// Constructs the protein quantification engine. 
        /// 
        /// "controlCondition" refers to the condition that the protein fold-changes could be calculated relative to.
        /// 
        /// If "nullHypothesisInterval" is null, the null hypothesis interval will be estimated from the data ("experimental null").
        /// 
        /// </summary>
        public ProteinQuantificationEngine(FlashLfqResults results, int maxThreads, string controlCondition, bool useSharedPeptides = false,
            double? nullHypothesisInterval = null, int? randomSeed = null, int mcmcBurninSteps = 1000, int mcmcSteps = 3000, bool pairedSamples = false)
        {
            if (string.IsNullOrWhiteSpace(controlCondition))
            {
                throw new MzLibException("A control condition must be defined to run the Bayesian protein quantification");
            }

            bool conditionsAreDefined = results.SpectraFiles.All(p => !string.IsNullOrWhiteSpace(p.Condition));
            if (!conditionsAreDefined)
            {
                throw new MzLibException("Conditions must be defined to run the Bayesian protein quantification");
            }

            this.MaxThreads = maxThreads;
            this.Results = results;
            UseSharedPeptides = useSharedPeptides;
            this.ControlCondition = controlCondition;
            this.NullHypothesisInterval = nullHypothesisInterval;
            this.BurnInSteps = mcmcBurninSteps;
            this.McmcSteps = mcmcSteps;
            this.PairedSamples = pairedSamples;

            if (randomSeed == null)
            {
                Rng = new MersenneTwister();
                RandomSeed = Rng.NextFullRangeInt32();
            }
            else
            {
                RandomSeed = randomSeed.Value;
            }

            Rng = new MersenneTwister(RandomSeed);
        }

        /// <summary>
        /// Runs the protein quantification engine.
        /// </summary>
        public void Run()
        {
            Setup();
            CalculatePeptideIonizationEfficiencies();

            // generate a random seed for each protein for the MCMC sampler
            var randomSeedsForEachProtein = new Dictionary<ProteinGroup, List<int>>();
            foreach (var protein in ProteinsWithConstituentPeptides)
            {
                randomSeedsForEachProtein.Add(protein.Key, new List<int> { Rng.NextFullRangeInt32(), Rng.NextFullRangeInt32() });
            }

            var proteinList = ProteinsWithConstituentPeptides.ToList();

            // figure out which conditions to compare
            var conditions = Results.SpectraFiles.Select(p => p.Condition).Distinct().ToList();
            conditions.Remove(ControlCondition);

            // calculate fold-change for each protein (diffuse prior)
            foreach (string condition in conditions)
            {
                Parallel.ForEach(Partitioner.Create(0, proteinList.Count), new ParallelOptions { MaxDegreeOfParallelism = MaxThreads }, partitionRange =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        ProteinGroup protein = proteinList[i].Key;

                        // get the fold-change measurements for the peptides assigned to this protein
                        List<(Peptide, List<double>)> peptideFoldChanges = GetPeptideFoldChangeMeasurements(protein, ControlCondition, condition, null, null);
                        int measurementCount = peptideFoldChanges.SelectMany(p => p.Item2).Count();

                        ProteinQuantificationEngineResult result;

                        if (measurementCount > 1)
                        {
                            if (PairedSamples)
                            {
                                result = RunPairedBayesianProteinQuant(peptideFoldChanges, protein, condition,
                                    randomSeedsForEachProtein[protein][0], null, false, BurnInSteps, McmcSteps, out var mus);
                            }
                            else
                            {
                                result = RunUnpairedBayesianProteinQuant(protein, condition,
                                    randomSeedsForEachProtein[protein][0], null, false, BurnInSteps, McmcSteps, out var mus);
                            }
                        }
                        else if (measurementCount == 1)
                        {
                            result = new ProteinQuantificationEngineResult(protein, ControlCondition, condition,
                                new double[] { peptideFoldChanges.SelectMany(p => p.Item2).First() }, new double[] { double.NaN },
                                new double[] { double.NaN });
                        }
                        else
                        {
                            result = new ProteinQuantificationEngineResult(protein, ControlCondition, condition,
                                new double[] { 0.0 }, new double[] { double.NaN }, new double[] { double.NaN });
                        }

                        protein.ConditionToQuantificationResults.Add(condition, result);
                    }
                });
            }

            // calculate PEPs for each protein (concentrated/skeptical prior)
            foreach (string condition in conditions)
            {
                double nullHyp = NullHypothesisInterval == null ? DetermineExperimentalNull(condition) : NullHypothesisInterval.Value;

                Parallel.ForEach(Partitioner.Create(0, proteinList.Count), new ParallelOptions { MaxDegreeOfParallelism = MaxThreads }, partitionRange =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        ProteinGroup protein = proteinList[i].Key;

                        // get the fold-change measurements for the peptides assigned to this protein
                        List<(Peptide, List<double>)> peptideFoldChanges = GetPeptideFoldChangeMeasurements(protein, ControlCondition, condition, null, null);
                        int measurementCount = peptideFoldChanges.SelectMany(p => p.Item2).Count();

                        ProteinQuantificationEngineResult result = protein.ConditionToQuantificationResults[condition];

                        if (measurementCount > 1)
                        {
                            double[] skepticalMus;

                            if (PairedSamples)
                            {
                                RunPairedBayesianProteinQuant(peptideFoldChanges, protein, condition,
                                    randomSeedsForEachProtein[protein][1], nullHyp, true, BurnInSteps, McmcSteps, out skepticalMus);
                            }
                            else
                            {
                                RunUnpairedBayesianProteinQuant(protein, condition,
                                    randomSeedsForEachProtein[protein][1], nullHyp, true, BurnInSteps, McmcSteps, out skepticalMus);
                            }

                            result.NullHypothesisCutoff = nullHyp;
                            result.CalculatePosteriorErrorProbability(skepticalMus);
                        }
                        else if (measurementCount == 1)
                        {
                            result.NullHypothesisCutoff = nullHyp;
                            result.CalculatePosteriorErrorProbability(new double[] { nullHyp });
                        }
                        else
                        {
                            result.NullHypothesisCutoff = nullHyp;
                            result.CalculatePosteriorErrorProbability(new double[] { nullHyp });
                        }
                    }
                });
            }

            // calculate FDR for the condition
            foreach (var condition in conditions)
            {
                var res = proteinList.Select(p => p.Key.ConditionToQuantificationResults[condition]).ToList();
                CalculateFalseDiscoveryRates(res);
            }
        }

        /// <summary>
        /// Pairs up proteins with constituent peptides and calculates peptide biorep intensities for easy lookup.
        /// </summary>
        private void Setup()
        {
            HashSet<Peptide> sharedPeptides = new HashSet<Peptide>();

            // determine shared peptides
            if (!UseSharedPeptides)
            {
                sharedPeptides = new HashSet<Peptide>(Results.PeptideModifiedSequences.Where(p => p.Value.ProteinGroups.Count > 1).Select(p => p.Value));
            }

            // match proteins to peptides
            ProteinsWithConstituentPeptides = Results.PeptideModifiedSequences.Values
                .SelectMany(p => p.ProteinGroups)
                .Distinct()
                .ToDictionary(p => p, p => new List<Peptide>());

            foreach (var peptide in Results.PeptideModifiedSequences)
            {
                if (!peptide.Value.UseForProteinQuant || (!UseSharedPeptides && sharedPeptides.Contains(peptide.Value)))
                {
                    continue;
                }

                foreach (ProteinGroup protein in peptide.Value.ProteinGroups)
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
            foreach (var condition in Results.SpectraFiles.GroupBy(p => p.Condition))
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
        private List<(Peptide, List<double>)> GetPeptideFoldChangeMeasurements(ProteinGroup protein, string condition1, string condition2, int? biorep1, int? biorep2)
        {
            List<(Peptide, List<double>)> allPeptideFoldChanges = new List<(Peptide, List<double>)>();

            List<Peptide> peptides = ProteinsWithConstituentPeptides[protein];

            int numB1 = Results.SpectraFiles.Where(p => p.Condition == condition1).Max(p => p.BiologicalReplicate) + 1;
            int numB2 = Results.SpectraFiles.Where(p => p.Condition == condition2).Max(p => p.BiologicalReplicate) + 1;

            foreach (var peptide in peptides)
            {
                List<double> peptideFoldChanges = new List<double>();

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
                    //if (valid < n && !PairedSamples)
                    //{
                    //    // there is at least one missing value... try to maximize the number of valid fold-change measurements
                    //    List<int> condition1Used = new List<int>();
                    //    List<int> condition2Used = new List<int>();
                    //    peptideFoldChanges.Clear();

                    //    for (int j = 0; j < cond1Intensities.Length; j++)
                    //    {
                    //        for (int k = 0; k < cond2Intensities.Length; k++)
                    //        {
                    //            double? foldChange = GetLogFoldChange(cond1Intensities[j], cond2Intensities[k]);

                    //            if (foldChange != null && !condition1Used.Contains(j) && !condition2Used.Contains(k))
                    //            {
                    //                peptideFoldChanges.Add(foldChange.Value);
                    //                condition1Used.Add(j);
                    //                condition2Used.Add(k);
                    //            }
                    //        }
                    //    }
                    //}
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

                allPeptideFoldChanges.Add((peptide, peptideFoldChanges));
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
        /// Estimates a relative ionization efficiency for each peptide.
        /// </summary>
        private void CalculatePeptideIonizationEfficiencies()
        {
            var peptides = PeptideToSampleQuantity.Where(v => v.Key.Item2 == ControlCondition).GroupBy(p => p.Key.Item1).ToList();

            Dictionary<int, int> randomSeeds = new Dictionary<int, int>();
            for (int i = 0; i < peptides.Count; i++)
            {
                randomSeeds.Add(i, Rng.NextFullRangeInt32());
            }

            Parallel.ForEach(Partitioner.Create(0, peptides.Count), new ParallelOptions { MaxDegreeOfParallelism = MaxThreads }, partitionRange =>
            {
                List<double> mmts = new List<double>();

                for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                {
                    var peptide = peptides[i];
                    mmts.Clear();

                    foreach (var mmt in peptide.Where(p => p.Value > 0))
                    {
                        mmts.Add(mmt.Value);
                    }

                    double ionizationEfficiencyEstimate = 0;

                    if (mmts.Count > 1)
                    {
                        double median = mmts.Median();
                        double sd = mmts.StandardDeviation();

                        AdaptiveMetropolisWithinGibbs sampler = new AdaptiveMetropolisWithinGibbs(
                            mmts.ToArray(),
                            new ProteinFoldChangeEstimationModel(
                                priorMuMean: median,
                                priorMuSd: sd,
                                muInitialGuess: median,
                                priorSdStart: sd / 1000,
                                priorSdEnd: sd * 1000,
                                sdInitialGuess: sd,
                                priorNuExponent: 1.0 / 29.0,
                                nuInitialGuess: 5),
                                seed: randomSeeds[i]
                            );

                        // burn in and then sample the MCMC chain
                        sampler.Run(1000, 500);

                        ionizationEfficiencyEstimate = sampler.MarkovChain.Select(p => p[0]).Median();
                    }
                    else if (mmts.Count == 1)
                    {
                        ionizationEfficiencyEstimate = mmts.First();
                    }

                    peptide.Key.IonizationEfficiency = ionizationEfficiencyEstimate;
                }
            });
        }

        private ProteinQuantificationEngineResult RunUnpairedBayesianProteinQuant(ProteinGroup protein, string condition, int randomSeed,
            double? nullHypothesisInterval, bool skepticalPrior, int burnin, int n, out double[] mus)
        {
            var peptides = ProteinsWithConstituentPeptides[protein];

            List<double> sample1Intensities = new List<double>();
            List<double> sample2Intensities = new List<double>();

            int samplesInControlCondition = Results.SpectraFiles.Where(p => p.Condition == ControlCondition).Max(p => p.BiologicalReplicate) + 1;
            int samplesInTreatmentCondition = Results.SpectraFiles.Where(p => p.Condition == condition).Max(p => p.BiologicalReplicate) + 1;

            foreach (var peptide in peptides)
            {
                for (int s = 0; s < samplesInControlCondition; s++)
                {
                    if (PeptideToSampleQuantity.TryGetValue((peptide, ControlCondition, s), out double intensity))
                    {
                        intensity /= peptide.IonizationEfficiency;
                        sample1Intensities.Add(intensity);
                    }
                }

                for (int s = 0; s < samplesInTreatmentCondition; s++)
                {
                    if (PeptideToSampleQuantity.TryGetValue((peptide, condition, s), out double intensity))
                    {
                        intensity /= peptide.IonizationEfficiency;
                        sample2Intensities.Add(intensity);
                    }
                }
            }

            // sample1
            double sd1 = Math.Max(0.001, sample1Intensities.StandardDeviation());
            double mean1 = sample1Intensities.Mean();

            double priorMuSd1 = Math.Max(0.05, skepticalPrior ? nullHypothesisInterval.Value : sd1 * 100);
            double priorMuMean1 = skepticalPrior ? 1.0 : mean1;

            AdaptiveMetropolisWithinGibbs sampler = new AdaptiveMetropolisWithinGibbs(
                sample1Intensities.ToArray(),
                new ProteinFoldChangeEstimationModel(
                    priorMuMean: priorMuMean1,
                    priorMuSd: priorMuSd1,
                    muInitialGuess: mean1,
                    priorSdStart: sd1 / 1000,
                    priorSdEnd: sd1 * 1000,
                    sdInitialGuess: sd1,
                    priorNuExponent: 1.0 / 29.0,
                    nuInitialGuess: 5),
                    seed: randomSeed
                );

            // burn in and then sample the MCMC chain
            sampler.Run(burnin, n);

            double[] mus1 = sampler.MarkovChain.Select(p => p[0]).ToArray();
            double[] sds1 = sampler.MarkovChain.Select(p => p[1]).ToArray();
            double[] nus1 = sampler.MarkovChain.Select(p => p[2]).ToArray();

            //sample2
            double sd2 = Math.Max(0.001, sample2Intensities.StandardDeviation());
            double mean2 = sample2Intensities.Mean();

            double priorMuSd2 = Math.Max(0.05, skepticalPrior ? nullHypothesisInterval.Value : sd2 * 100);
            double priorMuMean2 = skepticalPrior ? 1.0 : mean2;

            sampler = new AdaptiveMetropolisWithinGibbs(
                sample2Intensities.ToArray(),
                new ProteinFoldChangeEstimationModel(
                    priorMuMean: priorMuMean2,
                    priorMuSd: priorMuSd2,
                    muInitialGuess: mean2,
                    priorSdStart: sd2 / 1000,
                    priorSdEnd: sd2 * 1000,
                    sdInitialGuess: sd2,
                    priorNuExponent: 1.0 / 29.0,
                    nuInitialGuess: 5),
                    seed: randomSeed
                );

            // burn in and then sample the MCMC chain
            sampler.Run(burnin, n);

            double[] mus2 = sampler.MarkovChain.Select(p => p[0]).ToArray();
            double[] sds2 = sampler.MarkovChain.Select(p => p[1]).ToArray();
            double[] nus2 = sampler.MarkovChain.Select(p => p[2]).ToArray();

            // calculate fold-change
            double[] muDiffs = new double[mus1.Length];
            for (int i = 0; i < mus1.Length; i++)
            {
                muDiffs[i] = Math.Log(mus2[i] / mus1[i], 2);
            }

            // TODO: log-sd and not linear SD
            double[] jointSd = new double[mus1.Length];
            for (int i = 0; i < sds1.Length; i++)
            {
                double a = mus2[i];
                double b = mus1[i];
                double a_1 = sds2[i];
                double b_1 = sds1[i];

                jointSd[i] = Math.Sqrt((Math.Pow(a, 2) * Math.Pow(b_1, 2) + Math.Pow(a_1, 2) * Math.Pow(b, 2)) / Math.Pow(b, 4));
            }

            var result = new ProteinQuantificationEngineResult(
                protein,
                ControlCondition,
                condition,
                muDiffs,
                jointSd,
                nus1);

            mus = muDiffs;

            return result;
        }

        /// <summary>
        /// Estimates the fold-change between conditions for a protein, given its constituent peptides' fold change measurements. 
        /// If "skepticalPrior" is set to true, then the result's estimate of the fold-change between conditions is drawn towards 
        /// zero. This is useful for estimating the posterior error probability (PEP) that a fold-change is below a certain cutoff.
        /// </summary>
        private ProteinQuantificationEngineResult RunPairedBayesianProteinQuant(List<(Peptide, List<double>)> peptideFoldChanges, ProteinGroup protein, string condition, int randomSeed,
            double? nullHypothesisCutoff, bool skepticalPrior, int burnin, int n, out double[] mus)
        {
            if (skepticalPrior && nullHypothesisCutoff == null)
            {
                throw new MzLibException("Invalid protein quantification parameters given for protein " + protein.ProteinGroupName);
            }

            List<double> foldChanges = peptideFoldChanges.SelectMany(p => p.Item2).ToList();

            // the Math.Max is here because in some edge cases the SD can be 0, which causes a crash
            // also, the prior distribution for SD would be very narrow if SD < 0.001
            double sd = Math.Max(0.001, foldChanges.StandardDeviation());
            double mean = foldChanges.Mean();

            // the Math.max is here for cases where the cutoff is very small (e.g., 0)

            // "nullHypothesisCutoff.Value" means that the standard deviation of the Student's t prior probability distribution is 
            // the null hypothesis (e.g., the fold-change cutoff). so if the fold-change cutoff is 0.3, then the standard deviation 
            // of the prior probability distribution for mu is 0.3. the "degrees of freedom" (nu) of the prior will always be 1.0.
            // the reason for this is explained in the ProteinFoldChangeEstimationModel.cs file.

            // The standard deviation of the prior for mu is equal to the null hypothesis because the cumulative probability density of a 
            // Student's t distribution with nu=1 is 0.5, between -1 * SD and 1 * SD. In other words, 1 standard deviation from the mean is 50%
            // of the probability. This was chosen because the probability of the null hypothesis and the alternative hypothesis, given no data,
            // should each be 50%. So in summary, the prior here was chosen because given no data, the initial assumption is that the 
            // null hypothesis and the alternative hypothesis are equally likely. The datapoints are then used to "convince" the algorithm 
            // one way or another (towards the null or the alternative), with some estimated probability that either the null or alternative hypothesis
            // is true.
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
                ControlCondition,
                condition,
                mus,
                sds,
                nus
            );

            return result;
        }

        /// <summary>
        /// This function determines a fold-change cutoff (a null hypothesis) from the data. Each protein has an estimate of the
        /// standard deviation of its peptide fold-changes. The protein's standard deviations are pooled and the null hypothesis
        /// is estimated as the most common standard deviation. This is difficult to do because the distribution is usually skewed,
        /// so median and average are pretty poor estimates of the mode. Instead, the smallest interval with 10% of the distribution's 
        /// density is found, and the median of this range is an estimate of the mode. This estimate of the mode standard deviation 
        /// is used as the null hypothesis.
        /// </summary>
        private double DetermineExperimentalNull(string treatmentCondition)
        {
            double experimentalNull = 1.0;
            //var res = Results.ProteinGroups.Values.Select(p => p.ConditionToQuantificationResults[treatmentCondition]).ToList();

            //List<double> stdDevs = res.Where(p => p.PeptideFoldChangeMeasurements
            //    .SelectMany(v => v.foldChanges).Count() > 1).Select(p => p.StandardDeviationPointEstimate).ToList();

            //if (stdDevs.Any())
            //{
            //    var tenPercentHdi = Util.GetHighestDensityInterval(stdDevs.ToArray(), 0.1);
            //    var sdsWithinHdi = stdDevs.Where(p => p <= tenPercentHdi.hdi_end && p >= tenPercentHdi.hdi_start).ToList();

            //    var modeStdDevPointEstimate = sdsWithinHdi.Median();
            //    experimentalNull = modeStdDevPointEstimate;

            //    if (double.IsNaN(experimentalNull))
            //    {
            //        experimentalNull = stdDevs.Median();
            //    }
            //}

            return experimentalNull;
        }

        /// <summary>
        /// Calculates the false discovery rate of each protein from the Bayesian-estimated PEP values.
        /// </summary>
        private void CalculateFalseDiscoveryRates(List<ProteinQuantificationEngineResult> bayesianProteinQuantificationResults)
        {
            var bayesianQuantResults = bayesianProteinQuantificationResults
                .OrderByDescending(p => 1.0 - p.PosteriorErrorProbability)
                //.ThenByDescending(p => p.PeptideFoldChangeMeasurements.SelectMany(v => v.foldChanges).Count())
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
