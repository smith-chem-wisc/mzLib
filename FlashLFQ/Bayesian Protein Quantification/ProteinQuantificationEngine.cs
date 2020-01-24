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
        private Dictionary<string, List<(double, Gamma)>> PeptideIntensityNoiseEstimates;
        private Dictionary<string, List<(double, double)>> PeptideIntensityBiasEstimates;
        private List<string> TreatmentConditions;

        public Dictionary<(ProteinGroup, string), MultiGammaDistribution> ProteinAndConditionToSigmaPrior { get; private set; }
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
            DeterminePeptideSigmaDistributions();
            DetermineProteinSigmaPriors();
            EstimateProteinFoldChanges();
            DetermineIntensityDependentBias();
            PerformHypothesisTesting();
            CalculateFalseDiscoveryRates();
            AssignProteinIntensities();
        }

        public static (double[] mus, double[] sds, double[] nus) FitProteinQuantModel(List<double> measurements, bool skepticalPrior, bool paired,
            int? randomSeed, int burnin, int n, double nullHypothesisInterval, MultiGammaDistribution sdPrior)
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
            // also, the prior distribution for SD would be very narrow if SD < 0.001
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
            double priorMuSd = skepticalPrior ? nullHypothesisInterval : 3.0;
            double priorMuMean = 0;

            AdaptiveMetropolisWithinGibbs sampler = new AdaptiveMetropolisWithinGibbs(
                measurements.ToArray(),
                new ProteinFoldChangeEstimationModel(
                    priorMuMean: priorMuMean,
                    priorMuSd: priorMuSd,
                    muInitialGuess: mean,
                    sdPrior: sdPrior,
                    sdInitialGuess: sd,
                    priorNuExponent: 1.0 / 29.0,
                    nuInitialGuess: 5),
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
                                BurnInSteps, McmcSteps, controlNullInterval, treatmentNullInterval, controlSigmaPrior, treatmentSigmaPrior);
                        }
                    }
                });
            }
        }

        private void DeterminePeptideSigmaDistributions()
        {
            var peptides = Results.PeptideModifiedSequences.Values.ToList();
            PeptideIntensityNoiseEstimates = new Dictionary<string, List<(double, Gamma)>>();

            var conditions = Results.SpectraFiles.Select(p => p.Condition).Distinct().ToList();

            foreach (var condition in conditions)
            {
                PeptideIntensityNoiseEstimates.Add(condition, new List<(double, Gamma)>());
            }

            if (peptides.Count < 100)
            {
                foreach (var condition in conditions)
                {
                    PeptideIntensityNoiseEstimates[condition] = null;
                }

                return;
            }

            foreach (string condition in conditions)
            {
                List<(double, double)> peptideIntensityToStdev = new List<(double, double)>();
                var samples = Results.SpectraFiles.Where(p => p.Condition == condition).Max(p => p.BiologicalReplicate) + 1;

                List<double> intensities = new List<double>();
                foreach (Peptide peptide in peptides)
                {
                    intensities.Clear();

                    if (peptide.IonizationEfficiency == 0)
                    {
                        continue;
                    }

                    for (int s = 0; s < samples; s++)
                    {
                        if (PeptideToSampleQuantity.TryGetValue((peptide, condition, s), out double intensity) && intensity > 0)
                        {
                            intensities.Add(intensity);
                        }
                    }

                    double median = intensities.Select(p => Math.Log(p, 2)).Median();
                    double stdev = intensities.Select(p => Math.Log(p / peptide.IonizationEfficiency, 2)).StandardDeviation();

                    if (!double.IsNaN(median) && !double.IsNaN(stdev) && intensities.Count > 1 && stdev > 0)
                    {
                        stdev = GetUnbiasedSigma(stdev, intensities.Count);

                        peptideIntensityToStdev.Add((median, stdev));
                    }
                }

                if (peptideIntensityToStdev.Count < 100)
                {
                    PeptideIntensityNoiseEstimates[condition] = null;
                    return;
                }

                peptideIntensityToStdev.Sort((x, y) => y.Item1.CompareTo(x.Item1));
                var queue = new Queue<(double, double)>(peptideIntensityToStdev);

                int numPeptidesPerBin = Math.Min(100, peptideIntensityToStdev.Count / 10);

                List<List<(double, double)>> peptideIntensitiesAndStdevBinned = new List<List<(double, double)>>();
                List<(double, double)> bin = new List<(double, double)>();

                while (queue.Any())
                {
                    var element = queue.Dequeue();
                    bin.Add(element);

                    if (bin.Count == numPeptidesPerBin)
                    {
                        peptideIntensitiesAndStdevBinned.Add(bin);
                        bin = new List<(double, double)>();
                    }
                }

                Dictionary<int, int> randomSeeds = new Dictionary<int, int>();
                for (int i = 0; i < peptideIntensitiesAndStdevBinned.Count; i++)
                {
                    randomSeeds.Add(i, Rng.NextFullRangeInt32());
                }

                Parallel.ForEach(Partitioner.Create(0, peptideIntensitiesAndStdevBinned.Count),
                    new ParallelOptions { MaxDegreeOfParallelism = MaxThreads }, partitionRange =>
                    {
                        for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                        {
                            var theBin = peptideIntensitiesAndStdevBinned[i];
                            int randomSeed = randomSeeds[i];
                            var stdevs = theBin.Select(p => p.Item2).ToList();

                            AdaptiveMetropolisWithinGibbs sampler = new AdaptiveMetropolisWithinGibbs(
                                stdevs.ToArray(),
                                new GammaDistributionModel(
                                    priorShapeStart: 0.01,
                                    priorShapeEnd: 30,
                                    shapeInitialGuess: 2.0,
                                    priorRateStart: 0.01,
                                    priorRateEnd: 10,
                                    rateInitialGuess: 2.0),
                                    seed: randomSeed
                                );

                            // burn in and then sample the MCMC chain
                            sampler.Run(BurnInSteps, McmcSteps);

                            double[] shapes = sampler.MarkovChain.Select(p => p[0]).ToArray();
                            double[] rates = sampler.MarkovChain.Select(p => p[1]).ToArray();

                            double shapeEstimate = shapes.Median();
                            double rateEstimate = rates.Median();

                            var gamma = new Gamma(shapeEstimate, rateEstimate);

                            double intensityEstimate = theBin.Select(p => p.Item1).Median();

                            var estimates = PeptideIntensityNoiseEstimates[condition];

                            lock (estimates)
                            {
                                estimates.Add((intensityEstimate, gamma));
                            }
                        }
                    });

                PeptideIntensityNoiseEstimates[condition].Sort((x, y) => y.Item1.CompareTo(x.Item1));
            }
        }

        private void DetermineProteinSigmaPriors()
        {
            ProteinAndConditionToSigmaPrior = new Dictionary<(ProteinGroup, string), MultiGammaDistribution>();
            var conditions = Results.SpectraFiles.Select(p => p.Condition).Distinct().ToList();
            List<double> intensities = new List<double>();

            foreach (string condition in conditions)
            {
                var samples = Results.SpectraFiles.Where(p => p.Condition == condition).Max(p => p.BiologicalReplicate) + 1;

                foreach (var proteinWithPeptides in ProteinsWithConstituentPeptides)
                {
                    intensities.Clear();
                    ProteinGroup protein = proteinWithPeptides.Key;
                    List<Peptide> peptides = proteinWithPeptides.Value;
                    List<Gamma> gammas = new List<Gamma>();
                    List<int> weights = new List<int>();

                    foreach (var peptide in peptides)
                    {
                        if (peptide.IonizationEfficiency == 0)
                        {
                            continue;
                        }

                        int n = 0;
                        for (int s = 0; s < samples; s++)
                        {
                            if (PeptideToSampleQuantity.TryGetValue((peptide, condition, s), out double intensity) && intensity > 0)
                            {
                                n++;
                                intensities.Add(Math.Log(intensity / peptide.IonizationEfficiency, 2));
                            }
                        }

                        var peptideSigmaDistribution = GetNoiseFromPeptide(peptide, condition);

                        if (peptideSigmaDistribution != null)
                        {
                            gammas.Add(peptideSigmaDistribution);
                            //weights.Add(n);
                            weights.Add(1);
                        }
                    }

                    ProteinAndConditionToSigmaPrior.Add((protein, condition), new MultiGammaDistribution(gammas, weights));
                }
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

        private Gamma GetNoiseFromPeptide(Peptide peptide, string condition)
        {
            Gamma noiseEstimate = null;

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
                return noiseEstimate;
            }

            double medianIntensity = peptideIntensities.Median();

            if (double.IsNaN(medianIntensity) || double.IsInfinity(medianIntensity))
            {
                return noiseEstimate;
            }

            double logIntensity = Math.Log(medianIntensity, 2);

            var noiseEstimatesForCondition = PeptideIntensityNoiseEstimates[condition];

            for (int i = 0; i < noiseEstimatesForCondition.Count; i++)
            {
                noiseEstimate = noiseEstimatesForCondition[i].Item2;

                if (i == noiseEstimatesForCondition.Count - 1 || logIntensity > noiseEstimatesForCondition[i + 1].Item1)
                {
                    if (i < noiseEstimatesForCondition.Count - 1)
                    {
                        noiseEstimate = noiseEstimatesForCondition[i + 1].Item2;
                    }

                    break;
                }
            }

            return noiseEstimate;
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
                    break;
                }
            }

            return biasEstimate;
        }

        private void DetermineIntensityDependentBias()
        {
            PeptideIntensityBiasEstimates = new Dictionary<string, List<(double, double)>>();
            List<double> intensitiesControl = new List<double>();
            List<double> intensitiesTreatment = new List<double>();

            foreach (string treatmentCondition in TreatmentConditions)
            {
                var treatmentIntensityToBias = new List<(double, double, int, int)>();
                PeptideIntensityBiasEstimates.Add(treatmentCondition, new List<(double, double)>());

                foreach (var proteinWithPeptides in ProteinsWithConstituentPeptides)
                {
                    var protein = proteinWithPeptides.Key;
                    var peptides = proteinWithPeptides.Value;

                    int nPeptidesMeasuredControl = 0;
                    int nPeptidesMeasuredTreatment = 0;

                    int numSamplesControl = Results.SpectraFiles.Where(p => p.Condition == ControlCondition).Max(p => p.BiologicalReplicate) + 1;
                    int numSamplesTreatment = Results.SpectraFiles.Where(p => p.Condition == treatmentCondition).Max(p => p.BiologicalReplicate) + 1;

                    foreach (var peptide in peptides)
                    {
                        for (int s = 0; s < numSamplesControl; s++)
                        {
                            if (PeptideToSampleQuantity[(peptide, ControlCondition, s)] > 0)
                            {
                                nPeptidesMeasuredControl++;
                                break;
                            }
                        }

                        for (int s = 0; s < numSamplesTreatment; s++)
                        {
                            if (PeptideToSampleQuantity[(peptide, treatmentCondition, s)] > 0)
                            {
                                nPeptidesMeasuredTreatment++;
                                break;
                            }
                        }
                    }

                    if (nPeptidesMeasuredControl < 3 || nPeptidesMeasuredTreatment < 3)
                    {
                        continue;
                    }

                    var quantResult = protein.ConditionToQuantificationResults[treatmentCondition];

                    double proteinFoldChangeEstimate = quantResult.FoldChangePointEstimate;

                    foreach (var peptide in peptides)
                    {
                        intensitiesControl.Clear();
                        intensitiesTreatment.Clear();

                        for (int s = 0; s < numSamplesControl; s++)
                        {
                            double intensity = PeptideToSampleQuantity[(peptide, ControlCondition, s)];

                            if (intensity > 0)
                            {
                                intensitiesControl.Add(Math.Log(intensity, 2));
                            }
                        }

                        for (int s = 0; s < numSamplesTreatment; s++)
                        {
                            double intensity = PeptideToSampleQuantity[(peptide, treatmentCondition, s)];

                            if (intensity > 0)
                            {
                                intensitiesTreatment.Add(Math.Log(intensity, 2));
                            }
                        }

                        if (intensitiesControl.Count == 0 || intensitiesTreatment.Count == 0)
                        {
                            continue;
                        }

                        double bias = double.NaN;
                        double peptideTreatmentIntensity = intensitiesTreatment.Median();

                        if (PairedSamples)
                        {
                            List<double> biases = new List<double>();

                            //TODO: fix, these aren't paired up correctly
                            for (int i = 0; i < Math.Min(intensitiesTreatment.Count, intensitiesControl.Count); i++)
                            {
                                biases.Add(intensitiesTreatment[i] - intensitiesControl[i]);
                            }

                            bias = biases.Median();
                        }
                        else
                        {
                            double peptideChange = peptideTreatmentIntensity - intensitiesControl.Median();
                            bias = peptideChange - proteinFoldChangeEstimate;
                        }

                        if (!double.IsNaN(bias))
                        {
                            treatmentIntensityToBias.Add((peptideTreatmentIntensity, bias, intensitiesControl.Count, intensitiesTreatment.Count));
                        }
                    }
                }

                treatmentIntensityToBias.Sort((x, y) => y.Item1.CompareTo(x.Item1));

                var queue = new Queue<(double, double, int, int)>(treatmentIntensityToBias);

                int numPeptidesPerBin = Math.Min(100, treatmentIntensityToBias.Count / 10);

                var peptideIntensitiesAndBiasesBinned = new List<List<(double, double, int, int)>>();
                var bin = new List<(double, double, int, int)>();

                while (queue.Any())
                {
                    var element = queue.Dequeue();
                    bin.Add(element);

                    if (bin.Count == numPeptidesPerBin)
                    {
                        peptideIntensitiesAndBiasesBinned.Add(bin);
                        bin = new List<(double, double, int, int)>();
                    }
                }

                Dictionary<int, int> randomSeeds = new Dictionary<int, int>();
                for (int i = 0; i < peptideIntensitiesAndBiasesBinned.Count; i++)
                {
                    randomSeeds.Add(i, Rng.NextFullRangeInt32());
                }

                Parallel.ForEach(Partitioner.Create(0, peptideIntensitiesAndBiasesBinned.Count),
                    new ParallelOptions { MaxDegreeOfParallelism = MaxThreads }, partitionRange =>
                    {
                        for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                        {
                            var theBin = peptideIntensitiesAndBiasesBinned[i];
                            int randomSeed = randomSeeds[i];
                            var biases = theBin.Select(p => p.Item2).ToList();

                            double mean = biases.Mean();
                            double sd = Math.Max(0.1, biases.StandardDeviation());

                            AdaptiveMetropolisWithinGibbs sampler = new AdaptiveMetropolisWithinGibbs(
                                biases.ToArray(),
                                new ProteinFoldChangeEstimationModel(
                                    priorMuMean: mean,
                                    priorMuSd: sd * 100,
                                    muInitialGuess: mean,
                                    sdPrior: null,
                                    sdInitialGuess: sd,
                                    priorNuExponent: 1.0 / 29.0,
                                    nuInitialGuess: 5),
                                    seed: randomSeed
                                );

                            // burn in and then sample the MCMC chain
                            sampler.Run(BurnInSteps, McmcSteps);

                            double[] mus = sampler.MarkovChain.Select(p => p[0]).ToArray();
                            var hdi = Util.GetHighestDensityInterval(mus, 0.5);

                            double biasEstimate = Math.Max(Math.Abs(hdi.hdi_start), Math.Abs(hdi.hdi_end));

                            double intensityEstimate = theBin.Select(p => p.Item1).Median();

                            var estimates = PeptideIntensityBiasEstimates[treatmentCondition];

                            lock (estimates)
                            {
                                estimates.Add((intensityEstimate, biasEstimate));
                            }
                        }
                    });

                PeptideIntensityBiasEstimates[treatmentCondition].Sort((x, y) => y.Item1.CompareTo(x.Item1));

                ////DEBUG
                //List<string> output = new List<string> { "Log Intensity\tBias\tNControl\tNTreatment" };
                //foreach (var biasMeasurement in treatmentIntensityToBias)
                //{
                //    output.Add(biasMeasurement.Item1 + "\t" + biasMeasurement.Item2 + "\t" + biasMeasurement.Item3 + "\t" + biasMeasurement.Item4);
                //}
                //File.WriteAllLines(@"C:\Data\ionstarSample\biasMeasurement_" + treatmentCondition + ".tsv", output);
            }
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

                for (int p = 0; p < bayesianQuantResults.Count; p++)
                {
                    runningPEP += bayesianQuantResults[p].PosteriorErrorProbability;
                    qValue = runningPEP / (p + 1);
                    bayesianQuantResults[p].FalseDiscoveryRate = qValue;
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


            ////TODO: assign per-sample protein intensities
            //foreach (var protein in Results.ProteinGroups)
            //{
            //    foreach (var file in Results.SpectraFiles)
            //    {
            //        protein.Value.SetIntensity(file, 0);
            //    }

            //    if (ProteinsWithConstituentPeptides.TryGetValue(protein.Value, out var peptides))
            //    {
            //        foreach (var file in Results.SpectraFiles)
            //        {
            //            List<double> abundances = new List<double>();

            //            foreach (var peptide in peptides)
            //            {
            //                if (peptide.IonizationEfficiency == 0)
            //                {
            //                    continue;
            //                }

            //                double peptideIntensity = peptide.GetIntensity(file);
            //                double abundance = peptideIntensity / peptide.IonizationEfficiency;
            //                abundances.Add(abundance);
            //            }

            //            double proteinAbundance = abundances.Any() ? abundances.Median() : 0;
            //            double proteinIntensity = proteinAbundance * ProteinToControlConditionIntensity[protein.Value];
            //            protein.Value.SetIntensity(file, proteinIntensity);
            //        }
            //    }
            //}
        }
    }
}
