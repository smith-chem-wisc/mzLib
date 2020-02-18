using BayesianEstimation;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Random;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
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
        private readonly double FoldChangeCutoff;
        private readonly int BurnInSteps;
        private readonly int McmcSteps;
        private readonly MersenneTwister Rng;
        private Dictionary<ProteinGroup, List<Peptide>> ProteinsWithConstituentPeptides;
        private Dictionary<(Peptide, string, int), double> PeptideToSampleQuantity;
        private Dictionary<string, List<(double, double)>> PeptideIntensityBiasEstimates;
        private List<string> TreatmentConditions;
        private Dictionary<(ProteinGroup, string), Gamma> ProteinAndConditionToSigmaPrior;
        public readonly int RandomSeed;

        public ProteinQuantificationEngine(FlashLfqResults results, int maxThreads, string controlCondition, bool useSharedPeptides = false,
            double foldChangeCutoff = 0.585, int? randomSeed = null, int mcmcBurninSteps = 1000, int mcmcSteps = 3000, bool pairedSamples = false)
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
            this.FoldChangeCutoff = Math.Abs(foldChangeCutoff);
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

            //TODO: implement paired samples
            if (PairedSamples)
            {
                throw new NotImplementedException();
            }
        }

        /// <summary>
        /// Runs the protein quantification engine.
        /// </summary>
        public void Run()
        {
            Setup();
            EstimateIntensityDependentPeptideBiasAndUncertainty();
            EstimateProteinFoldChanges();
            PerformHypothesisTesting();
            CalculateFalseDiscoveryRates();
            AssignProteinIntensities();
        }

        public static (double[] mus, double[] sds, double[] nus) FitProteinQuantModel(List<Datum> measurements, bool skepticalPrior, bool paired,
            int? randomSeed, int burnin, int n, double nullHypothesisInterval, IContinuousDistribution sdPrior, IContinuousDistribution nuPrior,
            double minimumNu)
        {
            if (measurements.Count < 1)
            {
                return (new double[] { 0 }, new double[] { double.NaN }, new double[] { double.NaN });
            }
            else if (measurements.Count == 1)
            {
                double mmt = skepticalPrior ? 0 : measurements.First().DataValue;

                return (new double[] { mmt }, new double[] { double.NaN }, new double[] { double.NaN });
            }

            // the Math.Max is here because in some edge cases the SD can be 0, which causes a crash
            double sd = Math.Max(0.001, measurements.Select(p => p.DataValue).StandardDeviation());

            double meanOfData = measurements.Select(p => p.DataValue).Mean();

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

            double priorMuSd = skepticalPrior ? nullHypothesisInterval : 12;
            double priorMuMean = 0;

            AdaptiveMetropolisWithinGibbs sampler = new AdaptiveMetropolisWithinGibbs(
                measurements.ToArray(),
                new ProteinFoldChangeEstimationModel(
                    priorMuMean: priorMuMean,
                    priorMuSd: priorMuSd,
                    muInitialGuess: meanOfData,
                    sdPrior: sdPrior,
                    sdInitialGuess: sd,
                    nuPrior: nuPrior,
                    nuInitialGuess: 5,
                    minimumNu: minimumNu),
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
                                var techrepFile = techrep.First();
                                double techrepIntensity = peptide.GetIntensity(techrepFile);

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

                            //biorepIntensity += fractionIntensity;
                            biorepIntensity = Math.Max(biorepIntensity, fractionIntensity);
                        }

                        PeptideToSampleQuantity.Add((peptide, condition.Key, biorep.Key), biorepIntensity);
                    }
                }
            }

            if (Results.SpectraFiles.Where(p => p.Condition == ControlCondition).Select(p => p.BiologicalReplicate).Distinct().Count() == 1
                && !PairedSamples)
            {
                EstimateIonizationEfficiencies();
            }

            DetermineIonizationEfficiencies();
        }

        private void EstimateIonizationEfficiencies()
        {
            List<string> conditions = new List<string> { ControlCondition };//Results.SpectraFiles.Select(p => p.Condition).Distinct().ToList();
            conditions.AddRange(TreatmentConditions);

            var proteinList = ProteinsWithConstituentPeptides.ToList();

            Parallel.ForEach(Partitioner.Create(0, proteinList.Count), new ParallelOptions { MaxDegreeOfParallelism = MaxThreads }, partitionRange =>
            {
                for (int e = partitionRange.Item1; e < partitionRange.Item2; e++)
                {
                    ProteinGroup protein = proteinList[e].Key;
                    List<Peptide> peptides = proteinList[e].Value;

                    List<List<List<double>>> intensitiesByGroup = new List<List<List<double>>>();

                    foreach (var condition in conditions)
                    {
                        int numSamples = Results.SpectraFiles.Where(p => p.Condition == condition).Max(v => v.BiologicalReplicate) + 1;

                        List<List<double>> conditionIntensities = new List<List<double>>();

                        foreach (Peptide peptide in peptides)
                        {
                            List<double> peptideIntensities = new List<double>();

                            for (int s = 0; s < numSamples; s++)
                            {
                                var intensity = PeptideToSampleQuantity[(peptide, condition, s)];

                                if (intensity > 0)
                                {
                                    intensity = Math.Log(intensity, 2);
                                    peptideIntensities.Add(intensity);
                                }
                            }

                            conditionIntensities.Add(peptideIntensities);
                        }

                        intensitiesByGroup.Add(conditionIntensities);
                    }

                    if (!intensitiesByGroup.SelectMany(p => p.SelectMany(v => v.Select(k => k))).Any())
                    {
                        continue;
                    }

                    double minIntensity = intensitiesByGroup.SelectMany(p => p.SelectMany(v => v.Select(k => k))).Min();
                    double maxIntensity = intensitiesByGroup.SelectMany(p => p.SelectMany(v => v.Select(k => k))).Max();

                    double[] ionizationEfficiencyEstimations = new double[peptides.Count];
                    for (int i = 0; i < ionizationEfficiencyEstimations.Length; i++)
                    {
                        ionizationEfficiencyEstimations[i] = minIntensity - 5;
                    }

                    double oldBestError = double.PositiveInfinity;
                    double bestError = double.PositiveInfinity;
                    double[] bestIonizationEfficiencyEstimations = new double[peptides.Count];

                    for (int steps = 0; steps < 4; steps++)
                    {
                        double step = 0.5;

                        if (steps == 1)
                        {
                            step = 0.2;
                        }
                        if (steps == 2)
                        {
                            step = 0.05;
                        }
                        if (steps == 3)
                        {
                            step = 0.01;
                        }

                        for (int peptide = 0; peptide < peptides.Count; peptide++)
                        {
                            double start = bestIonizationEfficiencyEstimations[peptide] - (step * 10);
                            double end = steps == 0 ? maxIntensity + (step * 10) : bestIonizationEfficiencyEstimations[peptide] + (step * 10);

                            for (double i = start; i < end; i += step)
                            {
                                ionizationEfficiencyEstimations[peptide] = i;
                                double error = CalculateIonizationEfficiencyEstimationError(intensitiesByGroup, ionizationEfficiencyEstimations);

                                if (error < bestError)
                                {
                                    bestError = error;
                                    bestIonizationEfficiencyEstimations = ionizationEfficiencyEstimations.ToArray();
                                }
                            }

                            ionizationEfficiencyEstimations[peptide] = bestIonizationEfficiencyEstimations[peptide];
                        }

                        //if (oldBestError <= bestError)
                        //{
                        //    break;
                        //}

                        oldBestError = bestError;
                    }

                    double globalNormFactor = 0;

                    List<double> group1normintensities = new List<double>();

                    for (int p = 0; p < peptides.Count; p++)
                    {
                        var group1 = intensitiesByGroup[0][p].Select(v => v - bestIonizationEfficiencyEstimations[p]);
                        group1normintensities.AddRange(group1);
                    }

                    globalNormFactor = group1normintensities.Median();
                    if (double.IsNaN(globalNormFactor))
                    {
                        globalNormFactor = 0;
                    }

                    for (int i = 0; i < bestIonizationEfficiencyEstimations.Length; i++)
                    {
                        bestIonizationEfficiencyEstimations[i] += globalNormFactor;
                    }

                    for (int p = 0; p < peptides.Count; p++)
                    {
                        if (Results.SpectraFiles.All(s => peptides[p].GetIntensity(s) == 0))
                        {
                            peptides[p].IonizationEfficiency = 0;
                            continue;
                        }

                        peptides[p].IonizationEfficiency = Math.Pow(2, bestIonizationEfficiencyEstimations[p]);
                    }
                }
            });
        }

        private double CalculateIonizationEfficiencyEstimationError(List<List<List<double>>> intensitiesByGroup,
            double[] ionizationEfficiencyEstimations)
        {
            double totalError = 0;
            
            List<double> temp = new List<double>();

            // error is the sum of squared errors
            foreach (var group in intensitiesByGroup)
            {
                temp.Clear();

                for (int p = 0; p < group.Count; p++)
                {
                    var peptideIntensities = group[p];
                    var ionizationEfficiencyEstimationGuess = ionizationEfficiencyEstimations[p];

                    temp.AddRange(peptideIntensities.Select(v => v - ionizationEfficiencyEstimationGuess));

                    //double variance = peptideIntensities.Count < 10 ? temp.Variance() : Math.Pow(temp.InterquartileRange() * 1.3, 2);

                    //double mean = peptideIntensities.Count < 10 ? temp.Mean() : temp.Median();

                    //if (!double.IsNaN(variance))
                    //{
                    //    variances.Add(variance);
                    //}

                    //if (!double.IsNaN(mean))
                    //{
                    //    means.Add(mean);
                    //}
                }

                double median = temp.Median();

                double errorForGroup = 0;
                foreach (var element in temp)
                {
                    double error = Math.Abs(element - median);
                    errorForGroup += error;
                }

                //double meanOfVariance = variances.Mean();
                //double varianceOfMeans = means.Variance();

                //if (double.IsNaN(meanOfVariance))
                //{
                //    meanOfVariance = 0;
                //}

                //if (double.IsNaN(varianceOfMeans))
                //{
                //    varianceOfMeans = 0;
                //}

                //double pooledVariance = meanOfVariance + varianceOfMeans;
                totalError += errorForGroup;
            }

            return totalError;
        }

        private void DetermineIonizationEfficiencies()
        {
            var proteinList = ProteinsWithConstituentPeptides.ToList();

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
                            result = new PairedProteinQuantResult(protein, ProteinsWithConstituentPeptides[protein], ControlCondition,
                                treatmentCondition, UseSharedPeptides, Results, PeptideToSampleQuantity);
                        }
                        else
                        {
                            result = new UnpairedProteinQuantResult(protein, ProteinsWithConstituentPeptides[protein], ControlCondition,
                                treatmentCondition, UseSharedPeptides, Results, PeptideToSampleQuantity,
                                randomSeedsForEachProtein[protein], BurnInSteps);
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

                        var controlSigmaPrior = ProteinAndConditionToSigmaPrior[(protein, ControlCondition)];
                        var treatmentSigmaPrior = ProteinAndConditionToSigmaPrior[(protein, treatmentCondition)];

                        if (PairedSamples)
                        {
                            ((PairedProteinQuantResult)result).EstimateProteinFoldChange(randomSeedsForEachProtein[protein][0],
                                BurnInSteps, McmcSteps);
                        }
                        else
                        {

                            ((UnpairedProteinQuantResult)result).EstimateProteinFoldChange(randomSeedsForEachProtein[protein][0], randomSeedsForEachProtein[protein][1],
                                BurnInSteps, McmcSteps, controlSigmaPrior: controlSigmaPrior, treatmentSigmaPrior: treatmentSigmaPrior);
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
                        (double, double) proteinSpecificNullHypothesis = DetermineNullHypothesisWidth(protein, treatmentCondition);

                        ProteinQuantificationEngineResult result = protein.ConditionToQuantificationResults[treatmentCondition];

                        var controlSigmaPrior = ProteinAndConditionToSigmaPrior[(protein, ControlCondition)];
                        var treatmentSigmaPrior = ProteinAndConditionToSigmaPrior[(protein, treatmentCondition)];

                        if (PairedSamples)
                        {
                            var pairedResult = (PairedProteinQuantResult)result;

                            pairedResult.EstimateProteinFoldChange(randomSeedsForEachProtein[protein][0],
                                BurnInSteps, McmcSteps, proteinSpecificNullHypothesis.Item1);
                        }
                        else
                        {
                            double controlNullInterval = proteinSpecificNullHypothesis.Item1;
                            double treatmentNullInterval = proteinSpecificNullHypothesis.Item2;

                            var unpairedResult = (UnpairedProteinQuantResult)result;

                            unpairedResult.EstimateProteinFoldChange(randomSeedsForEachProtein[protein][0], randomSeedsForEachProtein[protein][1],
                                BurnInSteps, McmcSteps, controlNullInterval, treatmentNullInterval,
                                controlSigmaPrior, treatmentSigmaPrior);
                        }
                    }
                });
            }
        }

        public static double GetUnbiasedSigma(double sigma, int N)
        {
            //https://www.jstor.org/stable/pdf/2682923.pdf
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

        private (double, double) DetermineNullHypothesisWidth(ProteinGroup protein, string treatmentCondition)
        {
            double proteinSpecificNullHypothesisControl = Math.Sqrt(Math.Pow(FoldChangeCutoff, 2) / 2);
            double proteinSpecificNullHypothesisTreatment = Math.Sqrt(Math.Pow(FoldChangeCutoff, 2) / 2);

            var peptides = ProteinsWithConstituentPeptides[protein];
            var res = protein.ConditionToQuantificationResults[treatmentCondition];

            // add uncertainty to null hypothesis width
            if (PairedSamples)
            {
                //TODO
            }
            else
            {
                var unpairedRes = (UnpairedProteinQuantResult)res;
                var controlUncertaintyWidth = unpairedRes.hdi95Control.hdi_end - unpairedRes.hdi95Control.hdi_start;
                var treatmentUncertaintyWidth = unpairedRes.hdi95Treatment.hdi_end - unpairedRes.hdi95Treatment.hdi_start;

                if (!double.IsNaN(controlUncertaintyWidth))
                {
                    proteinSpecificNullHypothesisControl += controlUncertaintyWidth / 2;
                }
                if (!double.IsNaN(treatmentUncertaintyWidth))
                {
                    proteinSpecificNullHypothesisTreatment += treatmentUncertaintyWidth / 2;
                }
            }

            // add bias to null hypothesis width
            double intensityDependentBiasForThisProtein = 0;
            List<(double bias, int weight)> biasesWithWeights = new List<(double, int)>();

            foreach (var peptide in peptides)
            {
                double bias = GetBiasFromPeptide(peptide, treatmentCondition);

                if (!double.IsNaN(bias))
                {
                    // TODO: weight based on # observations per peptide?
                    biasesWithWeights.Add((bias, 1));
                }
            }

            if (biasesWithWeights.Any())
            {
                intensityDependentBiasForThisProtein = biasesWithWeights.Sum(p => p.bias * p.weight) / biasesWithWeights.Sum(p => p.weight);
            }

            proteinSpecificNullHypothesisControl += Math.Abs(Math.Sqrt(Math.Pow(intensityDependentBiasForThisProtein, 2) / 2));
            proteinSpecificNullHypothesisTreatment += Math.Abs(Math.Sqrt(Math.Pow(intensityDependentBiasForThisProtein, 2) / 2));

            return (proteinSpecificNullHypothesisControl, proteinSpecificNullHypothesisTreatment);
        }

        private void EstimateIntensityDependentPeptideBiasAndUncertainty()
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

                        if (diff != 0 && nPeptidesMeasured > 3)
                        {
                            intensityToPeptideDiffToProtein.Add((logIntensity, diff));
                        }
                    }
                }

                intensityToPeptideDiffToProtein.Sort((x, y) => y.Item1.CompareTo(x.Item1));
                intensityToStdev.Sort((x, y) => y.Item1.CompareTo(x.Item1));

                // measure variance of means and mean of means
                List<(double, double[])> meansAndVarianceOfMeans = new List<(double, double[])>();

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

                List<string> output = new List<string>();

                if (peptideIntensitiesAndBiasesBinned.Any())
                {
                    Parallel.ForEach(Partitioner.Create(0, peptideIntensitiesAndBiasesBinned.Count),
                        new ParallelOptions { MaxDegreeOfParallelism = MaxThreads }, partitionRange =>
                        {
                            for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                            {
                                var theBin = peptideIntensitiesAndBiasesBinned[i];
                                var biases = theBin.Select(p => p.Item2).ToList();
                                double intensityEstimate = theBin.Select(p => p.Item1).Median();

                                double mean = biases.Median();

                                var hdi = Util.GetHighestDensityInterval(biases.ToArray(), 0.682);

                                double sd = Math.Max(0.0001, hdi.hdi_end - hdi.hdi_start);
                                //var temp = biases.Where(p => p < hdi.hdi_end && p > hdi.hdi_end).ToList();

                                //mean = temp.Median();
                                double percentOutliers = biases.Count(p => Math.Abs(p - mean) > 3 * sd) / ((double)biases.Count);
                                //sd = Math.Max(0.0001, biases.StandardDeviation()); //DEBUG

                                // gamma prior for SD
                                double mean2 = sd / 2;
                                double sd2 = Math.Pow(2 * sd, 2);

                                var reparams = ReparameterizeGamma(mean2, Math.Pow(sd2, 2));

                                AdaptiveMetropolisWithinGibbs sampler = new AdaptiveMetropolisWithinGibbs(
                                    biases.ToArray(),
                                    new ProteinFoldChangeEstimationModel(
                                        priorMuMean: mean,
                                        priorMuSd: sd * 100,
                                        muInitialGuess: mean,
                                        sdPrior: new Gamma(reparams.shape, reparams.rate),
                                        sdInitialGuess: sd,
                                        nuPrior: null,
                                        nuInitialGuess: 5),
                                        seed: randomSeeds[i][0]
                                    );

                                // burn in and then sample the MCMC chain
                                sampler.Run(BurnInSteps, McmcSteps);

                                double[] mus = sampler.MarkovChain.Select(p => p[0]).ToArray();
                                double[] sds = sampler.MarkovChain.Select(p => p[1]).ToArray();
                                double[] nus = sampler.MarkovChain.Select(p => p[2]).ToArray();

                                double[] vars = new double[sds.Length];

                                for (int j = 0; j < vars.Length; j++)
                                {
                                    double nu = nus[j];
                                    double multiplier = (nus[j] / (nus[j] - 2));
                                    vars[j] = sds[j] * sds[j] * (nus[j] / (nus[j] / 2));
                                }

                                // mean of means
                                double biasEstimate = mus.Median();

                                lock (output)
                                {
                                    output.Add(intensityEstimate + "\t" + percentOutliers);
                                }

                                lock (meansAndVarianceOfMeans)
                                {
                                    meansAndVarianceOfMeans.Add((intensityEstimate, vars));
                                }

                                var biasEstimates = PeptideIntensityBiasEstimates[condition];

                                lock (biasEstimates)
                                {
                                    biasEstimates.Add((intensityEstimate, biasEstimate));
                                }
                            }
                        });
                }

                File.WriteAllLines(@"C:\Data\Ecoli_Human_Spikein\fractionated_outliers.tsv", output);
                PeptideIntensityBiasEstimates[condition].Sort((x, y) => y.Item1.CompareTo(x.Item1));

                meansAndVarianceOfMeans.Sort((x, y) => y.Item1.CompareTo(x.Item1));

                // fit gamma to peptide variances
                List<(double, double[])> peptideVariancesArray = new List<(double, double[])>();
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

                //DEBUG
                //List<string> output2 = new List<string>();

                if (peptideIntensitiesAndSigmasBinned.Any())
                {
                    Parallel.ForEach(Partitioner.Create(0, peptideIntensitiesAndSigmasBinned.Count),
                        new ParallelOptions { MaxDegreeOfParallelism = MaxThreads }, partitionRange =>
                        {
                            for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                            {
                                var theBin = peptideIntensitiesAndSigmasBinned[i];
                                var sigmas = theBin.Select(p => p.Item2).ToArray();
                                double intensityEstimate = theBin.Select(p => p.Item1).Median();

                                //double mean = sigmas.Mean();
                                //double variance = Math.Max(sigmas.Variance(), 0.01);

                                //var reparams = ReparameterizeGamma(mean, variance);

                                var rng2 = new MersenneTwister(seed: randomSeeds[i][1]);
                                //var gamma = new Gamma(reparams.shape, reparams.rate, rng2);
                                double[] sigmas2 = new double[McmcSteps];

                                for (int j = 0; j < sigmas2.Length; j++)
                                {
                                    int random = rng2.Next(0, sigmas.Length);
                                    sigmas2[j] = sigmas[random];
                                }

                                //gamma.Samples(sigmas2);

                                //DEBUG
                                //var t = gamma.InverseCumulativeDistribution(0.9);
                                //var prop = ((double)sigmas.Count(p => p >= t)) / sigmas.Length;
                                //Array.Sort(sigmas2);


                                lock (peptideVariancesArray)
                                {
                                    peptideVariancesArray.Add((intensityEstimate, sigmas2.Select(p => Math.Pow(p, 2)).ToArray()));
                                }
                            }
                        });
                }

                //File.WriteAllLines(@"C:\Data\Ecoli_Human_Spikein\gammaOutliers.tsv", output2);

                peptideVariancesArray.Sort((x, y) => y.Item1.CompareTo(x.Item1));
                CreateProteinSigmaPriors(condition, meansAndVarianceOfMeans, peptideVariancesArray);
            }
        }

        private void CreateProteinSigmaPriors(string condition, List<(double, double[])> varianceOfPeptideMeansArrays, List<(double, double[])> peptideVariancesArrays)
        {
            // meansAndVarianceOfMeans will be empty if every protein has only 1 peptide (e.g., top-down proteoform quant)
            // peptideVariancesArray will be empty if there is only 1 sample per condition, because within-group variance is unknown

            var list = ProteinsWithConstituentPeptides.ToList();

            if (ProteinAndConditionToSigmaPrior == null)
            {
                ProteinAndConditionToSigmaPrior = new Dictionary<(ProteinGroup, string), Gamma>();
            }

            foreach (var protein in list)
            {
                ProteinAndConditionToSigmaPrior.Add((protein.Key, condition), null);
            }

            int numSamples = Results.SpectraFiles.Where(p => p.Condition == condition).Max(p => p.BiologicalReplicate) + 1;

            Dictionary<int, List<int>> randomSeeds = new Dictionary<int, List<int>>();
            for (int i = 0; i < ProteinsWithConstituentPeptides.Count; i++)
            {
                randomSeeds.Add(i, new List<int> { Rng.NextFullRangeInt32(), Rng.NextFullRangeInt32() });
            }

            Parallel.ForEach(Partitioner.Create(0, list.Count),
                    new ParallelOptions { MaxDegreeOfParallelism = MaxThreads }, partitionRange =>
                     {
                         for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                         {
                             var protein = list[i].Key;
                             var peptides = list[i].Value;
                             var rng2 = new MersenneTwister(seed: randomSeeds[i][1]);

                             // should be an array for each peptide in these 2 lists
                             //List<double[]> means = new List<double[]>();
                             List<double[]> peptideVariances = new List<double[]>();
                             List<double[]> variancesOfPeptideMeans = new List<double[]>();

                             double[] pooledVarianceArray = new double[McmcSteps];

                             List<double> intensities = new List<double>();

                             foreach (var peptide in peptides)
                             {
                                 intensities.Clear();

                                 if (peptide.IonizationEfficiency == 0)
                                 {
                                     continue;
                                 }

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
                                 double[] bestVarianceOfMeans = null;
                                 double[] bestPeptideVariances = null;

                                 for (int j = 0; j < varianceOfPeptideMeansArrays.Count; j++)
                                 {
                                     bestVarianceOfMeans = varianceOfPeptideMeansArrays[j].Item2;

                                     if (j == varianceOfPeptideMeansArrays.Count - 1 || logIntensity > varianceOfPeptideMeansArrays[j + 1].Item1)
                                     {
                                         if (j < varianceOfPeptideMeansArrays.Count - 1)
                                         {
                                             bestVarianceOfMeans = varianceOfPeptideMeansArrays[j + 1].Item2;
                                         }

                                         break;
                                     }
                                 }

                                 for (int j = 0; j < peptideVariancesArrays.Count; j++)
                                 {
                                     bestPeptideVariances = peptideVariancesArrays[j].Item2;

                                     if (j == peptideVariancesArrays.Count - 1 || logIntensity > peptideVariancesArrays[j + 1].Item1)
                                     {
                                         if (j < peptideVariancesArrays.Count - 1)
                                         {
                                             bestPeptideVariances = peptideVariancesArrays[j + 1].Item2;
                                         }

                                         break;
                                     }
                                 }

                                 if (bestVarianceOfMeans != null)
                                 {
                                     variancesOfPeptideMeans.Add(bestVarianceOfMeans);
                                 }

                                 if (bestPeptideVariances != null)
                                 {
                                     peptideVariances.Add(bestPeptideVariances);
                                 }
                             }

                             for (int j = 0; j < pooledVarianceArray.Length; j++)
                             {
                                 double averageOfVariances = peptideVariances.Count > 0 ? peptideVariances.Select(v => v[j]).Average() : 0;
                                 double varianceOfMeans = 0;

                                 if (variancesOfPeptideMeans.Count == 0)
                                 {
                                     varianceOfMeans = 0;
                                 }
                                 else if (variancesOfPeptideMeans.Count == 1)
                                 {
                                     varianceOfMeans = variancesOfPeptideMeans[0][j];
                                 }
                                 else
                                 {
                                     varianceOfMeans = variancesOfPeptideMeans.Select(v => v[j]).Average() + variancesOfPeptideMeans.Select(v => v[j]).Variance();
                                 }

                                 double pooledVariance = averageOfVariances + varianceOfMeans;
                                 pooledVarianceArray[j] = pooledVariance;
                             }

                             // SD has no prior prob distribution for this protein
                             if (pooledVarianceArray.All(p => p == 0))
                             {
                                 continue;
                             }

                             var pooledSigmaArray = pooledVarianceArray.Select(p => Math.Sqrt(p)).ToArray();
                             double mean = pooledSigmaArray.Mean();


                             //double variance = Math.Max(pooledSigmaArray.Variance(), 0.1);

                             //(double shape, double rate) = ReparameterizeGamma(mean, variance);

                             // gamma prior for SD
                             double mean2 = mean / 2;
                             double sd2 = Math.Pow(2 * mean, 2);

                             var reparams = ReparameterizeGamma(mean2, Math.Pow(sd2, 2));

                             var proteinSigmaPrior = new Gamma(reparams.shape, reparams.rate);

                             lock (ProteinAndConditionToSigmaPrior)
                             {
                                 ProteinAndConditionToSigmaPrior[(protein, condition)] = proteinSigmaPrior;
                             }
                         }
                     });
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
                    .ThenByDescending(p => Math.Abs(p.FoldChangePointEstimate))
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

            //TODO: assign per-sample protein intensities
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
