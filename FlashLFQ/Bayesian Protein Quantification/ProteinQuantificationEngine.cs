using BayesianEstimation;
using MathNet.Numerics.Distributions;
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
        private readonly double NullHypothesisInterval;
        private readonly int BurnInSteps;
        private readonly int McmcSteps;
        private readonly MersenneTwister Rng;
        private Dictionary<ProteinGroup, List<Peptide>> ProteinsWithConstituentPeptides;
        private Dictionary<(Peptide, string, int), double> PeptideToSampleQuantity;
        private Dictionary<string, List<(double, double)>> PeptideIntensityBiasEstimates;
        private List<string> TreatmentConditions;

        public Dictionary<(ProteinGroup, string), Gamma> ProteinAndConditionToSigmaPrior { get; private set; }
        public Dictionary<(ProteinGroup, string), Weibull> ProteinAndConditionToNuPrior { get; private set; }
        public readonly int RandomSeed;

        public ProteinQuantificationEngine(FlashLfqResults results, int maxThreads, string controlCondition, bool useSharedPeptides = false,
            double nullHypothesisInterval = 0.1375, int? randomSeed = null, int mcmcBurninSteps = 1000, int mcmcSteps = 3000, bool pairedSamples = false)
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
            this.UseSharedPeptides = useSharedPeptides;
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

            // figure out which conditions to compare
            TreatmentConditions = Results.SpectraFiles.Select(p => p.Condition).Distinct().ToList();
            TreatmentConditions.Remove(ControlCondition);
        }

        /// <summary>
        /// Runs the protein quantification engine.
        /// </summary>
        public void Run()
        {
            Setup();
            EstimateProteinFoldChanges();
            EstimateBiasAndProteinSigmaPriors();
            CreateProteinNuPriors();
            PerformHypothesisTesting();
            CalculateFalseDiscoveryRates();
            AssignProteinIntensities();
        }

        public static (double[] mus, double[] sds, double[] nus) FitProteinQuantModel(List<double> measurements, bool skepticalPrior, bool paired,
            int? randomSeed, int burnin, int n, double nullHypothesisInterval, IContinuousDistribution sdPrior, IContinuousDistribution nuPrior)
        {
            if (measurements.Count < 1)
            {
                return (new double[] { 0 }, new double[] { double.NaN }, new double[] { double.NaN });
            }
            else if (measurements.Count == 1)
            {
                double mmt = skepticalPrior ? 0 : measurements.First();

                return (new double[] { mmt }, new double[] { double.NaN }, new double[] { double.NaN });
            }

            // the Math.Max is here because in some edge cases the SD can be 0, which causes a crash
            double sd = Math.Max(0.001, measurements.StandardDeviation());

            double mean = measurements.Mean();

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
            double priorMuSd = skepticalPrior ? nullHypothesisInterval : 20.0;
            double priorMuMean = 0;

            AdaptiveMetropolisWithinGibbs sampler = new AdaptiveMetropolisWithinGibbs(
                measurements.ToArray(),
                new ProteinFoldChangeEstimationModel(
                    priorMuMean: priorMuMean,
                    priorMuSd: priorMuSd,
                    muInitialGuess: mean,
                    sdPrior: sdPrior,
                    sdInitialGuess: sd,
                    nuPrior: nuPrior,
                    nuInitialGuess: nuPrior == null ? 5 : nuPrior.Median),
                    seed: randomSeed
                );

            // burn in and then sample the MCMC chain
            sampler.Run(burnin, n);

            double[] mus = sampler.MarkovChain.Select(p => p[0]).ToArray();
            double[] sds = sampler.MarkovChain.Select(p => p[1]).ToArray();
            double[] nus = sampler.MarkovChain.Select(p => p[2]).ToArray();

            return (mus, sds, nus);
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

            DetermineIonizationEfficiencies();
        }

        private void DetermineIonizationEfficiencies()
        {
            var proteinList = ProteinsWithConstituentPeptides.ToList();

            //DEBUG
            //proteinList = proteinList.Where(p => p.Key.ProteinGroupName == "Q9Y6J9").ToList();

            // calculate fold-change for each protein (diffuse prior)
            foreach (string treatmentCondition in TreatmentConditions)
            {
                // generate a random seed for each protein for the MCMC sampler
                var randomSeedsForEachProtein = new Dictionary<ProteinGroup, List<int>>();
                foreach (var protein in proteinList)
                {
                    randomSeedsForEachProtein.Add(protein.Key, new List<int>());

                    foreach (var peptide in protein.Value)
                    {
                        randomSeedsForEachProtein[protein.Key].Add(Rng.NextFullRangeInt32());
                    }
                }

                Parallel.ForEach(Partitioner.Create(0, proteinList.Count), new ParallelOptions { MaxDegreeOfParallelism = MaxThreads }, partitionRange =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        ProteinGroup protein = proteinList[i].Key;

                        ProteinQuantificationEngineResult result;

                        if (PairedSamples)
                        {
                            result = new PairedProteinQuantResult(protein, ProteinsWithConstituentPeptides[protein], ControlCondition, treatmentCondition, UseSharedPeptides,
                                Results, PeptideToSampleQuantity);
                        }
                        else
                        {
                            result = new UnpairedProteinQuantResult(protein, ProteinsWithConstituentPeptides[protein], ControlCondition, treatmentCondition, UseSharedPeptides,
                                Results, PeptideToSampleQuantity, randomSeedsForEachProtein[protein], BurnInSteps);
                        }

                        protein.ConditionToQuantificationResults.Add(treatmentCondition, result);
                    }
                });
            }
        }

        private void EstimateProteinFoldChanges()
        {
            var proteinList = ProteinsWithConstituentPeptides.ToList();

            // calculate fold-change for each protein (diffuse prior)
            foreach (string treatmentCondition in TreatmentConditions)
            {
                // generate a random seed for each protein for the MCMC sampler
                var randomSeedsForEachProtein = new Dictionary<ProteinGroup, List<int>>();
                foreach (var protein in proteinList)
                {
                    randomSeedsForEachProtein.Add(protein.Key, new List<int> { Rng.NextFullRangeInt32(), Rng.NextFullRangeInt32() });
                }

                Parallel.ForEach(Partitioner.Create(0, proteinList.Count), new ParallelOptions { MaxDegreeOfParallelism = MaxThreads }, partitionRange =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        ProteinGroup protein = proteinList[i].Key;

                        ProteinQuantificationEngineResult result = protein.ConditionToQuantificationResults[treatmentCondition];

                        if (PairedSamples)
                        {
                            ((PairedProteinQuantResult)result).EstimateProteinFoldChange(randomSeedsForEachProtein[protein][0],
                                BurnInSteps, McmcSteps);
                        }
                        else
                        {
                            ((UnpairedProteinQuantResult)result).EstimateProteinFoldChange(randomSeedsForEachProtein[protein][0], randomSeedsForEachProtein[protein][1],
                                BurnInSteps, McmcSteps);
                        }
                    }
                });
            }
        }

        private void PerformHypothesisTesting()
        {
            var proteinList = ProteinsWithConstituentPeptides.ToList();

            //DEBUG
            //proteinList = proteinList.Where(p => p.Key.ProteinGroupName == "Q9Y6J9").ToList();

            foreach (string treatmentCondition in TreatmentConditions)
            {
                // generate a random seed for each protein for the MCMC sampler
                var randomSeedsForEachProtein = new Dictionary<ProteinGroup, List<int>>();
                foreach (var protein in proteinList)
                {
                    randomSeedsForEachProtein.Add(protein.Key, new List<int> { Rng.NextFullRangeInt32(), Rng.NextFullRangeInt32() });
                }

                Parallel.ForEach(Partitioner.Create(0, proteinList.Count), new ParallelOptions { MaxDegreeOfParallelism = MaxThreads }, partitionRange =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        ProteinGroup protein = proteinList[i].Key;
                        double proteinSpecificNullHypothesis = DetermineNullHypothesisWidth(protein, treatmentCondition);

                        ProteinQuantificationEngineResult result = protein.ConditionToQuantificationResults[treatmentCondition];

                        var controlSigmaPrior = ProteinAndConditionToSigmaPrior[(protein, ControlCondition)];
                        var treatmentSigmaPrior = ProteinAndConditionToSigmaPrior[(protein, treatmentCondition)];

                        var controlNuPrior = ProteinAndConditionToNuPrior[(protein, ControlCondition)];
                        var treatmentNuPrior = ProteinAndConditionToNuPrior[(protein, treatmentCondition)];

                        if (PairedSamples)
                        {
                            var pairedResult = (PairedProteinQuantResult)result;

                            pairedResult.EstimateProteinFoldChange(randomSeedsForEachProtein[protein][0],
                                BurnInSteps, McmcSteps, proteinSpecificNullHypothesis);
                        }
                        else
                        {
                            double controlNullInterval = Math.Sqrt(Math.Pow(proteinSpecificNullHypothesis, 2) / 2);
                            double treatmentNullInterval = Math.Sqrt(Math.Pow(proteinSpecificNullHypothesis, 2) / 2);

                            var unpairedResult = (UnpairedProteinQuantResult)result;

                            unpairedResult.EstimateProteinFoldChange(randomSeedsForEachProtein[protein][0], randomSeedsForEachProtein[protein][1],
                                BurnInSteps, McmcSteps, controlNullInterval, treatmentNullInterval,
                                controlSigmaPrior, treatmentSigmaPrior, controlNuPrior, treatmentNuPrior);
                        }
                    }
                });
            }
        }

        private double GetUnbiasedSigma(double sigma, int N)
        {
            if (N < 2)
            {
                return double.NaN;
            }

            double cn = 1 + (1.0 / (4 * (N - 1)));

            return sigma * cn;
        }

        private double GetBiasFromPeptide(Peptide peptide, string condition)
        {
            double biasEstimate = double.NaN;

            int numSamples = Results.SpectraFiles.Where(p => p.Condition == condition).Max(p => p.BiologicalReplicate) + 1;

            List<double> peptideIntensities = new List<double>();
            for (int s = 0; s < numSamples; s++)
            {
                double intensity = PeptideToSampleQuantity[(peptide, condition, s)];

                if (intensity > 0)
                {
                    peptideIntensities.Add(intensity);
                }
            }

            if (peptideIntensities.Count < 1)
            {
                return biasEstimate;
            }

            double medianIntensity = peptideIntensities.Median();

            if (double.IsNaN(medianIntensity) || double.IsInfinity(medianIntensity))
            {
                return biasEstimate;
            }

            double logIntensity = Math.Log(medianIntensity, 2);

            var biasEstimatesForCondition = PeptideIntensityBiasEstimates[condition];

            for (int i = 0; i < biasEstimatesForCondition.Count; i++)
            {
                biasEstimate = biasEstimatesForCondition[i].Item2;

                if (i == biasEstimatesForCondition.Count - 1 || logIntensity > biasEstimatesForCondition[i + 1].Item1)
                {
                    if (i < biasEstimatesForCondition.Count - 1)
                    {
                        biasEstimate = biasEstimatesForCondition[i + 1].Item2;
                    }

                    break;
                }
            }

            return biasEstimate;
        }

        private double DetermineNullHypothesisWidth(ProteinGroup protein, string treatmentCondition)
        {
            double proteinSpecificNullHypothesis = NullHypothesisInterval;
            double intensityDependentBiasForThisProtein = 0;
            List<(double bias, int weight)> biasesWithWeights = new List<(double, int)>();

            var peptides = ProteinsWithConstituentPeptides[protein];

            foreach (var peptide in peptides)
            {
                double bias = GetBiasFromPeptide(peptide, treatmentCondition);

                if (!double.IsNaN(bias))
                {
                    // TODO: weight?
                    biasesWithWeights.Add((bias, 1));
                }
            }

            if (biasesWithWeights.Any())
            {
                intensityDependentBiasForThisProtein = biasesWithWeights.Sum(p => p.bias * p.weight) / biasesWithWeights.Sum(p => p.weight);
            }

            proteinSpecificNullHypothesis += Math.Abs(intensityDependentBiasForThisProtein);

            return proteinSpecificNullHypothesis;
        }

        private void EstimateBiasAndProteinSigmaPriors()
        {
            PeptideIntensityBiasEstimates = new Dictionary<string, List<(double, double)>>();
            var conditions = Results.SpectraFiles.Select(p => p.Condition).Distinct().ToList();

            foreach (string condition in conditions)
            {
                PeptideIntensityBiasEstimates.Add(condition, new List<(double, double)>());

                int numSamples = Results.SpectraFiles.Where(p => p.Condition == condition).Max(p => p.BiologicalReplicate) + 1;

                // used for calculating 1. variance of means, and 2. bias
                List<(double, double)> intensityToPeptideDiffToProtein = new List<(double, double)>();

                // used for calculating mean of variances
                List<(double, double)> intensityToStdev = new List<(double, double)>();

                foreach (var proteinWithPeptides in ProteinsWithConstituentPeptides)
                {
                    var protein = proteinWithPeptides.Key;
                    var peptides = proteinWithPeptides.Value;

                    List<(double logIntensity, double logAbundance)> peptideMedianLogAbundances = new List<(double, double)>();
                    List<double> peptideLogAbundances = new List<double>();

                    int nPeptidesMeasured = 0;

                    foreach (var peptide in peptides)
                    {
                        for (int s = 0; s < numSamples; s++)
                        {
                            if (PeptideToSampleQuantity[(peptide, condition, s)] > 0)
                            {
                                nPeptidesMeasured++;
                                break;
                            }
                        }
                    }

                    //if (nPeptidesMeasured < 3)
                    //{
                    //    continue;
                    //}

                    foreach (var peptide in peptides)
                    {
                        if (peptide.IonizationEfficiency == 0)
                        {
                            continue;
                        }

                        List<double> intensities = new List<double>();
                        for (int s = 0; s < numSamples; s++)
                        {
                            double intensity = PeptideToSampleQuantity[(peptide, condition, s)];

                            if (intensity > 0)
                            {
                                intensities.Add(intensity);
                                peptideLogAbundances.Add(Math.Log(intensity / peptide.IonizationEfficiency, 2));
                            }
                        }

                        if (intensities.Count > 0)
                        {
                            double peptideMedianIntensity = intensities.Median();
                            double logMedianIntensity = Math.Log(peptideMedianIntensity, 2);
                            peptideMedianLogAbundances.Add((logMedianIntensity, Math.Log(peptideMedianIntensity / peptide.IonizationEfficiency, 2)));

                            if (intensities.Count > 1)
                            {
                                var stdev = intensities.Select(p => Math.Log(p / peptide.IonizationEfficiency, 2)).StandardDeviation();
                                stdev = GetUnbiasedSigma(stdev, intensities.Count);

                                intensityToStdev.Add((logMedianIntensity, stdev));
                            }
                        }
                    }

                    double proteinLogAbundance = peptideLogAbundances.Median();

                    foreach (var peptideLogAbundance in peptideMedianLogAbundances)
                    {
                        double diff = peptideLogAbundance.logAbundance - proteinLogAbundance;
                        double logIntensity = peptideLogAbundance.logIntensity;

                        if (diff != 0)
                        {
                            intensityToPeptideDiffToProtein.Add((logIntensity, diff));
                        }
                    }
                }

                intensityToPeptideDiffToProtein.Sort((x, y) => y.Item1.CompareTo(x.Item1));
                intensityToStdev.Sort((x, y) => y.Item1.CompareTo(x.Item1));

                // measure variance of means and mean of means
                List<(double, double[], double[])> meansArrays = new List<(double, double[], double[])>();

                var queue = new Queue<(double, double)>(intensityToPeptideDiffToProtein);

                int numPeptidesPerBin = Math.Min(100, intensityToPeptideDiffToProtein.Count / 10);

                var peptideIntensitiesAndBiasesBinned = new List<List<(double, double)>>();
                var bin = new List<(double, double)>();

                while (queue.Any())
                {
                    var element = queue.Dequeue();
                    bin.Add(element);

                    if (bin.Count == numPeptidesPerBin)
                    {
                        peptideIntensitiesAndBiasesBinned.Add(bin);
                        bin = new List<(double, double)>();
                    }
                }

                Dictionary<int, List<int>> randomSeeds = new Dictionary<int, List<int>>();
                for (int i = 0; i < peptideIntensitiesAndBiasesBinned.Count; i++)
                {
                    randomSeeds.Add(i, new List<int> { Rng.NextFullRangeInt32(), Rng.NextFullRangeInt32() });
                }

                Parallel.ForEach(Partitioner.Create(0, peptideIntensitiesAndBiasesBinned.Count),
                    new ParallelOptions { MaxDegreeOfParallelism = MaxThreads }, partitionRange =>
                    {
                        for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                        {
                            var theBin = peptideIntensitiesAndBiasesBinned[i];
                            var biases = theBin.Select(p => p.Item2).ToList();
                            double intensityEstimate = theBin.Select(p => p.Item1).Median();

                            double mean = biases.Median();
                            double sd = Math.Max(0.1, biases.StandardDeviation());

                            // gamma prior for SD
                            double mean2 = biases.StandardDeviation() / 2;
                            double sd2 = Math.Pow(2 * biases.StandardDeviation(), 2);

                            var reparams = ReparameterizeGamma(mean2, Math.Pow(sd2, 2));

                            AdaptiveMetropolisWithinGibbs sampler = new AdaptiveMetropolisWithinGibbs(
                                biases.ToArray(),
                                new ProteinFoldChangeEstimationModel(
                                    priorMuMean: mean,
                                    priorMuSd: sd * 100,
                                    muInitialGuess: mean,
                                    sdPrior: new Gamma(reparams.shape, reparams.rate), //TODO: prior for this?
                                    sdInitialGuess: sd,
                                    nuPrior: null,
                                    nuInitialGuess: 5),
                                    seed: randomSeeds[i][0]
                                );

                            // burn in and then sample the MCMC chain
                            sampler.Run(BurnInSteps, 1000);

                            double[] mus = sampler.MarkovChain.Select(p => p[0]).ToArray();
                            double[] sds = sampler.MarkovChain.Select(p => p[1]).ToArray();

                            // mean of means
                            double biasEstimate = mus.Median();

                            lock (meansArrays)
                            {
                                meansArrays.Add((intensityEstimate, mus, sds));
                            }

                            var biasEstimates = PeptideIntensityBiasEstimates[condition];

                            lock (biasEstimates)
                            {
                                biasEstimates.Add((intensityEstimate, biasEstimate));
                            }
                        }
                    });

                PeptideIntensityBiasEstimates[condition].Sort((x, y) => y.Item1.CompareTo(x.Item1));

                meansArrays.Sort((x, y) => y.Item1.CompareTo(x.Item1));

                // fit gamma to peptide variances
                List<(double, double[])> peptideSigmasArray = new List<(double, double[])>();
                queue = new Queue<(double, double)>(intensityToStdev);

                numPeptidesPerBin = Math.Min(100, queue.Count / 10);

                var peptideIntensitiesAndSigmasBinned = new List<List<(double, double)>>();
                bin = new List<(double, double)>();

                while (queue.Any())
                {
                    var element = queue.Dequeue();
                    bin.Add(element);

                    if (bin.Count == numPeptidesPerBin)
                    {
                        peptideIntensitiesAndSigmasBinned.Add(bin);
                        bin = new List<(double, double)>();
                    }
                }

                for (int i = 0; i < peptideIntensitiesAndSigmasBinned.Count; i++)
                {
                    if (!randomSeeds.ContainsKey(i))
                    {
                        randomSeeds.Add(i, new List<int> { Rng.NextFullRangeInt32(), Rng.NextFullRangeInt32() });
                    }
                }

                Parallel.ForEach(Partitioner.Create(0, peptideIntensitiesAndSigmasBinned.Count),
                    new ParallelOptions { MaxDegreeOfParallelism = MaxThreads }, partitionRange =>
                    {
                        for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                        {
                            //DEBUG
                            var rng2 = new MersenneTwister(seed: randomSeeds[i][1]);

                            var theBin = peptideIntensitiesAndSigmasBinned[i];
                            var sigmas = theBin.Select(p => p.Item2).ToList();
                            double intensityEstimate = theBin.Select(p => p.Item1).Median();

                            double mean = sigmas.StandardDeviation() / 2;
                            double sd = Math.Pow(2 * sigmas.StandardDeviation(), 2);

                            var reparams = ReparameterizeGamma(mean, Math.Pow(sd, 2));

                            AdaptiveMetropolisWithinGibbs sampler = new AdaptiveMetropolisWithinGibbs(
                                sigmas.ToArray(),
                                new GammaDistributionModel(
                                    new Gamma(reparams.shape, reparams.rate), 1, new Gamma(reparams.shape, reparams.rate), 1),
                                    seed: randomSeeds[i][1]
                                );

                            // burn in and then sample the MCMC chain
                            sampler.Run(BurnInSteps, 1000);

                            double[] shapes = sampler.MarkovChain.Select(p => p[0]).ToArray();
                            double[] rates = sampler.MarkovChain.Select(p => p[0]).ToArray();
                            double[] sigmas2 = new double[shapes.Length];

                            for (int j = 0; j < shapes.Length; j++)
                            {
                                var gamma = new Gamma(shapes[j], rates[j], randomSource: rng2);
                                //sigmas2[j] = gamma.Mode;
                                gamma.Samples(sigmas2);
                            }

                            lock (peptideSigmasArray)
                            {
                                peptideSigmasArray.Add((intensityEstimate, sigmas2));
                            }
                        }
                    });

                peptideSigmasArray.Sort((x, y) => y.Item1.CompareTo(x.Item1));
                CreateProteinSigmaPriors(condition, meansArrays, peptideSigmasArray);
            }
        }

        private void CreateProteinSigmaPriors(string condition, List<(double, double[], double[])> meansArrays, List<(double, double[])> peptideSigmasArray)
        {
            if (ProteinAndConditionToSigmaPrior == null)
            {
                ProteinAndConditionToSigmaPrior = new Dictionary<(ProteinGroup, string), Gamma>();
                //ProteinAndConditionToNuPrior = new Dictionary<(ProteinGroup, string), Gamma>();
            }

            int numSamples = Results.SpectraFiles.Where(p => p.Condition == condition).Max(p => p.BiologicalReplicate) + 1;

            Dictionary<int, List<int>> randomSeeds = new Dictionary<int, List<int>>();
            for (int i = 0; i < ProteinsWithConstituentPeptides.Count; i++)
            {
                randomSeeds.Add(i, new List<int> { Rng.NextFullRangeInt32(), Rng.NextFullRangeInt32() });
            }

            var list = ProteinsWithConstituentPeptides.ToList();

            foreach (var protein in list)
            {
                ProteinAndConditionToSigmaPrior.Add((protein.Key, condition), null);
            }

            //DEBUG
            //list = list.Where(p => p.Key.ProteinGroupName == "Q9H488").ToList();

            Parallel.ForEach(Partitioner.Create(0, list.Count),
                    new ParallelOptions { MaxDegreeOfParallelism = MaxThreads }, partitionRange =>
                    {
                        for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                        {
                            var protein = list[i].Key;
                            var peptides = list[i].Value;
                            var rng2 = new MersenneTwister(seed: randomSeeds[i][1]);

                            // should be an array for each peptide in these 2 lists
                            List<double[]> listOfMus = new List<double[]>();
                            List<double[]> listOfSigmas = new List<double[]>();
                            List<double[]> sdOfMeansInBin = new List<double[]>();

                            double[] pooledVarianceArray = new double[meansArrays.First().Item2.Length];

                            foreach (var peptide in peptides)
                            {
                                if (peptide.IonizationEfficiency == 0)
                                {
                                    continue;
                                }

                                List<double> intensities = new List<double>();
                                for (int s = 0; s < numSamples; s++)
                                {
                                    double intensity = PeptideToSampleQuantity[(peptide, condition, s)];

                                    if (intensity > 0)
                                    {
                                        intensities.Add(intensity);
                                    }
                                }

                                if (intensities.Count == 0)
                                {
                                    continue;
                                }

                                double logIntensity = Math.Log(intensities.Median(), 2);

                                var bestMu = meansArrays.FirstOrDefault(p => p.Item1 < logIntensity);
                                var bestSd = peptideSigmasArray.FirstOrDefault(p => p.Item1 < logIntensity);

                                if (bestMu.Item2 == null)
                                {
                                    bestMu = meansArrays.Last();
                                }
                                if (bestSd.Item2 == null)
                                {
                                    bestSd = peptideSigmasArray.Last();
                                }

                                double[] muArray = bestMu.Item2.OrderBy(x => rng2.Next()).ToArray();
                                double[] sdArray = bestSd.Item2;

                                listOfMus.Add(muArray);
                                listOfSigmas.Add(sdArray);
                                sdOfMeansInBin.Add(bestMu.Item3);
                            }

                            if (listOfMus.Count == 0 || listOfSigmas.Count == 0)
                            {
                                continue;
                            }

                            for (int j = 0; j < pooledVarianceArray.Length; j++)
                            {
                                double averageOfVariances = listOfSigmas.Select(sigmas => Math.Pow(sigmas[j], 2)).Average();
                                double varianceOfAverages = 0;

                                if (listOfMus.Count == 1)
                                {
                                    double add = condition == ControlCondition ?
                                    ((UnpairedProteinQuantResult)protein.ConditionToQuantificationResults.First().Value).SigmaControl :
                                        ((UnpairedProteinQuantResult)protein.ConditionToQuantificationResults[condition]).SigmaTreatment;

                                    if (double.IsNaN(add))
                                    {
                                        add = 0;
                                    }

                                    varianceOfAverages = Math.Pow(sdOfMeansInBin[0][j], 2);
                                }
                                else
                                {
                                    varianceOfAverages = listOfMus.Select(mus => mus[j]).Variance();
                                }

                                double pooledVariance = averageOfVariances + varianceOfAverages;
                                pooledVarianceArray[j] = pooledVariance;
                            }

                            List<double> thinnedPooledSigma = new List<double>();
                            for (int j = 0; j < pooledVarianceArray.Length; j += 30)
                            {
                                thinnedPooledSigma.Add(Math.Sqrt(pooledVarianceArray[j]));
                            }

                            double mean = thinnedPooledSigma.StandardDeviation() / 2;
                            double sd = Math.Pow(2 * thinnedPooledSigma.StandardDeviation(), 2);

                            var reparams = ReparameterizeGamma(mean, Math.Pow(sd, 2));
                            AdaptiveMetropolisWithinGibbs sampler = new AdaptiveMetropolisWithinGibbs(
                                            thinnedPooledSigma.ToArray(),
                                            new GammaDistributionModel(
                                                new Gamma(reparams.shape, reparams.rate), 1, new Gamma(reparams.shape, reparams.rate), 1),
                                                seed: randomSeeds[i][0]
                                            );

                            // burn in and then sample the MCMC chain
                            sampler.Run(BurnInSteps, 500);

                            double[] shapes = sampler.MarkovChain.Select(p => p[0]).ToArray();
                            double[] rates = sampler.MarkovChain.Select(p => p[1]).ToArray();

                            double shape = shapes.Median();
                            double rate = rates.Median();

                            var gamma = new Gamma(shape, rate);

                            lock (ProteinAndConditionToSigmaPrior)
                            {
                                ProteinAndConditionToSigmaPrior[(protein, condition)] = gamma;
                            }
                        }
                    });
        }

        private void CreateProteinNuPriors()
        {
            ProteinAndConditionToNuPrior = new Dictionary<(ProteinGroup, string), Weibull>();
            var proteins = ProteinsWithConstituentPeptides.ToList();
            var conditions = Results.SpectraFiles.Select(p => p.Condition).Distinct().ToList();

            foreach (var protein in proteins)
            {
                foreach (var condition in conditions)
                {
                    ProteinAndConditionToNuPrior.Add((protein.Key, condition), null);
                }
            }

            foreach (var treatmentCondition in TreatmentConditions)
            {
                var results = proteins.Select(p => p.Key.ConditionToQuantificationResults[treatmentCondition]);

                var grouped = results.GroupBy(p => p.Peptides.Count).Where(v => v.Key < 5 && v.Count() > 20);

                foreach (var group in grouped)
                {
                    int numPeptides = group.Key;

                    if (!PairedSamples)
                    {
                        var unpairedResults = group.Select(v => (UnpairedProteinQuantResult)v).ToList();
                        var controlNus = unpairedResults.Select(p => p.NuControl).Where(p => !double.IsNaN(p)).ToList();

                        // control
                        AdaptiveMetropolisWithinGibbs sampler = new AdaptiveMetropolisWithinGibbs(
                                            controlNus.ToArray(),
                                            new WeibullDistributionModel(
                                                0, 100, 5, 0, 100, 5),
                                                seed: Rng.NextFullRangeInt32()
                                            );

                        sampler.Run(BurnInSteps, 1000);

                        double[] shapes = sampler.MarkovChain.Select(p => p[0]).ToArray();
                        double[] scales = sampler.MarkovChain.Select(p => p[1]).ToArray();

                        double shape = shapes.Median();
                        double scale = scales.Median();

                        foreach (var protein in group)
                        {
                            if (ProteinAndConditionToNuPrior[(protein.Protein, ControlCondition)] == null)
                            {
                                ProteinAndConditionToNuPrior[(protein.Protein, ControlCondition)] = new Weibull(shape, scale);
                            }
                        }

                        // treatment
                        var treatmentNus = unpairedResults.Select(p => p.NuTreatment).Where(p => !double.IsNaN(p)).ToList();
                        sampler = new AdaptiveMetropolisWithinGibbs(
                                            treatmentNus.ToArray(),
                                            new WeibullDistributionModel(
                                                0, 100, 5, 0, 100, 5),
                                                seed: Rng.NextFullRangeInt32()
                                            );

                        sampler.Run(BurnInSteps, 3000);

                        shapes = sampler.MarkovChain.Select(p => p[0]).ToArray();
                        scales = sampler.MarkovChain.Select(p => p[1]).ToArray();

                        shape = shapes.Median();
                        scale = scales.Median();

                        foreach (var protein in group)
                        {
                            if (ProteinAndConditionToNuPrior[(protein.Protein, treatmentCondition)] == null)
                            {
                                ProteinAndConditionToNuPrior[(protein.Protein, treatmentCondition)] = new Weibull(shape, scale);
                            }
                        }
                    }
                    else
                    {
                        throw new NotImplementedException();
                    }
                }
            }
        }

        private (double shape, double rate) ReparameterizeGamma(double mean, double variance)
        {
            double shape = Math.Pow(mean, 2) / variance;
            double rate = mean / variance;

            return (shape, rate);
        }

        /// <summary>
        /// Calculates the false discovery rate of each protein from the Bayesian-estimated PEP values.
        /// </summary>
        private void CalculateFalseDiscoveryRates()
        {
            var proteinList = ProteinsWithConstituentPeptides.ToList();

            // calculate FDR for the condition
            foreach (string treatmentCondition in TreatmentConditions)
            {
                var bayesianProteinQuantificationResults = proteinList.Select(p => p.Key.ConditionToQuantificationResults[treatmentCondition]).ToList();

                var bayesianQuantResults = bayesianProteinQuantificationResults
                    .OrderByDescending(p => p.IsStatisticallyValid)
                    .ThenByDescending(p => 1.0 - p.PosteriorErrorProbability)
                    .ThenByDescending(p => p.NMeasurements)
                    .ToList();

                double runningPEP = 0;
                double qValue = 0;

                //DEBUG
                double runningRealPEP = 0;
                double FDP = 0;

                for (int p = 0; p < bayesianQuantResults.Count; p++)
                {
                    runningPEP += bayesianQuantResults[p].PosteriorErrorProbability;
                    qValue = runningPEP / (p + 1);
                    bayesianQuantResults[p].FalseDiscoveryRate = qValue;

                    //DEBUG
                    double falsePositive = 0;
                    if (bayesianQuantResults[p].Protein.Organism == "HUMAN")
                    {
                        falsePositive = 1;
                    }

                    runningRealPEP += falsePositive;
                    FDP = runningRealPEP / (p + 1);

                    if (!PairedSamples)
                    {
                        var unpaired = (UnpairedProteinQuantResult)bayesianQuantResults[p];
                        unpaired.FalsePositive = falsePositive;
                        unpaired.FDP = FDP;
                    }
                }
            }
        }

        private void AssignProteinIntensities()
        {
            // calculate reference intensities
            var ProteinToControlConditionIntensity = new Dictionary<ProteinGroup, double>();

            var numBRef = Results.SpectraFiles.Where(p => p.Condition == ControlCondition).Select(p => p.BiologicalReplicate).Distinct().Count();

            foreach (var protein in ProteinsWithConstituentPeptides)
            {
                List<double> intensities = new List<double>();

                foreach (Peptide peptide in protein.Value)
                {
                    double avgIntensity = 0;
                    int n = 0;

                    for (int b = 0; b < numBRef; b++)
                    {
                        double intensity = PeptideToSampleQuantity[(peptide, ControlCondition, b)];
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

                ProteinToControlConditionIntensity.Add(protein.Key, referenceProteinIntensity);

                foreach (var res in protein.Key.ConditionToQuantificationResults)
                {
                    res.Value.ControlConditionIntensity = referenceProteinIntensity;
                    res.Value.TreatmentConditionIntensity = referenceProteinIntensity * Math.Pow(2, res.Value.FoldChangePointEstimate);
                }
            }

            //TODO: assign per-sample protein intensities
            foreach (var protein in Results.ProteinGroups)
            {
                foreach (var file in Results.SpectraFiles)
                {
                    protein.Value.SetIntensity(file, 0);
                }

                if (ProteinsWithConstituentPeptides.TryGetValue(protein.Value, out var peptides))
                {
                    foreach (var file in Results.SpectraFiles)
                    {
                        List<double> abundances = new List<double>();

                        foreach (var peptide in peptides)
                        {
                            if (peptide.IonizationEfficiency == 0)
                            {
                                continue;
                            }

                            double peptideIntensity = peptide.GetIntensity(file);
                            double abundance = peptideIntensity / peptide.IonizationEfficiency;
                            abundances.Add(abundance);
                        }

                        double proteinAbundance = abundances.Any() ? abundances.Median() : 0;
                        double proteinIntensity = proteinAbundance * ProteinToControlConditionIntensity[protein.Value];
                        protein.Value.SetIntensity(file, proteinIntensity);
                    }
                }
            }
        }
    }
}
