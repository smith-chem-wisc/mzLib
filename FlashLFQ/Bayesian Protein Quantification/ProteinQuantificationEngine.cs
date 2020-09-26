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
        private double FoldChangeCutoff;
        private readonly int BurnInSteps;
        private readonly int McmcSteps;
        private readonly MersenneTwister Rng;
        private Dictionary<ProteinGroup, List<Peptide>> ProteinsWithConstituentPeptides;
        private Dictionary<(Peptide, string, int), (double, DetectionType)> PeptideToSampleQuantity;
        private List<string> TreatmentConditions;
        private Dictionary<(ProteinGroup, string), StudentT> ProteinAndConditionToSigmaPrior;
        public readonly int RandomSeed;

        public ProteinQuantificationEngine(FlashLfqResults results, int maxThreads, string controlCondition, bool useSharedPeptides = false,
            double foldChangeCutoff = 0.1, int? randomSeed = null, int mcmcBurninSteps = 1000, int mcmcSteps = 3000, bool pairedSamples = false)
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
            this.FoldChangeCutoff = foldChangeCutoff;
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

            FoldChangeCutoff = Math.Abs(foldChangeCutoff);
        }

        /// <summary>
        /// Runs the protein quantification engine.
        /// </summary>
        public void Run()
        {
            Setup();
            EstimateIntensityDependentUncertainty();
            EstimateProteinFoldChanges();
            PerformHypothesisTesting();
            CalculateFalseDiscoveryRates();
            AssignProteinIntensities();
        }

        public static (double[] mus, double[] sds, double[] nus) FitProteinQuantModel(List<Datum> measurements, bool skepticalPrior, bool paired,
            int? randomSeed, int burnin, int n, double nullHypothesisInterval, IContinuousDistribution sdPrior)
        {
            int validMeasurements = measurements.Count(p => p.Weight > 0);

            if (validMeasurements < 1)
            {
                return (new double[] { 0 }, new double[] { double.NaN }, new double[] { double.NaN });
            }
            else if (validMeasurements == 1)
            {
                double mmt = skepticalPrior ? 0 : measurements.First(p => p.Weight > 0).X;

                return (new double[] { mmt }, new double[] { double.NaN }, new double[] { double.NaN });
            }

            // the Math.Max is here because in some edge cases the SD can be 0, which causes a crash
            double sd = Math.Max(0.001, measurements.Select(p => p.X).StandardDeviation());

            double meanOfData = measurements.Select(p => p.X).Mean();

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

            double priorMuSd = skepticalPrior ? nullHypothesisInterval : 20;
            double priorMuMean = 0;

            AdaptiveMetropolisWithinGibbs sampler = new AdaptiveMetropolisWithinGibbs(
                measurements.ToArray(),
                new ProteinFoldChangeEstimationModel(
                    priorMuMean: priorMuMean,
                    priorMuSd: priorMuSd,
                    muInitialGuess: meanOfData,
                    sdPrior: sdPrior,
                    sdInitialGuess: sd,
                    nuPrior: new Exponential(1.0 / 29.0),
                    nuInitialGuess: 5,
                    minimumNu: 1,
                    minimumSd: sdPrior == null ? 0 : sdPrior.Mode),
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
            PeptideToSampleQuantity = new Dictionary<(Peptide, string, int), (double, DetectionType)>();

            foreach (var condition in Results.SpectraFiles.GroupBy(p => p.Condition))
            {
                foreach (var sample in condition.GroupBy(p => p.BiologicalReplicate))
                {
                    foreach (Peptide peptide in allpeptides)
                    {
                        double sampleIntensity = 0;
                        double highestFractionIntensity = 0;

                        DetectionType sampleDetectionType = DetectionType.NotDetected;

                        foreach (var fraction in sample.GroupBy(p => p.Fraction))
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

                            sampleIntensity += fractionIntensity;

                            if (fractionIntensity > highestFractionIntensity)
                            {
                                highestFractionIntensity = fractionIntensity;

                                DetectionType fractionDetectionType = peptide.GetDetectionType(fraction.First());
                                sampleDetectionType = fractionDetectionType;
                            }
                        }

                        PeptideToSampleQuantity.Add((peptide, condition.Key, sample.Key), (sampleIntensity, sampleDetectionType));
                    }
                }
            }

            if (Results.SpectraFiles.Where(p => p.Condition == ControlCondition).Select(p => p.BiologicalReplicate).Distinct().Count() == 1
                && !PairedSamples)
            {
                EstimateIonizationEfficienciesIfControlConditionHasOneSample();
            }

            EstimateIonizationEfficiencies();
        }

        /// <summary>
        /// This method is only used if the control condition has 1 sample. In more typical cases, "EstimateIonizationEfficiencies" is called instead.
        /// </summary>
        private void EstimateIonizationEfficienciesIfControlConditionHasOneSample()
        {
            List<string> conditions = new List<string> { ControlCondition };
            conditions.AddRange(TreatmentConditions);

            var proteinList = ProteinsWithConstituentPeptides.ToList();

            Parallel.ForEach(Partitioner.Create(0, proteinList.Count), new ParallelOptions { MaxDegreeOfParallelism = MaxThreads }, partitionRange =>
            {
                for (int e = partitionRange.Item1; e < partitionRange.Item2; e++)
                {
                    ProteinGroup protein = proteinList[e].Key;
                    List<Peptide> peptides = proteinList[e].Value;

                    List<List<List<(double, DetectionType)>>> intensitiesByGroup = new List<List<List<(double, DetectionType)>>>();

                    foreach (var condition in conditions)
                    {
                        int numSamples = Results.SpectraFiles.Where(p => p.Condition == condition).Max(v => v.BiologicalReplicate) + 1;

                        List<List<(double, DetectionType)>> conditionIntensities = new List<List<(double, DetectionType)>>();

                        foreach (Peptide peptide in peptides)
                        {
                            List<(double, DetectionType)> peptideIntensities = new List<(double, DetectionType)>();

                            for (int s = 0; s < numSamples; s++)
                            {
                                var intensity = PeptideToSampleQuantity[(peptide, condition, s)].Item1;

                                if (intensity > 0)
                                {
                                    intensity = Math.Log(intensity, 2);
                                    peptideIntensities.Add((intensity, PeptideToSampleQuantity[(peptide, condition, s)].Item2));
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

                    double minIntensity = intensitiesByGroup.SelectMany(p => p.SelectMany(v => v.Select(k => k.Item1))).Min();
                    double maxIntensity = intensitiesByGroup.SelectMany(p => p.SelectMany(v => v.Select(k => k.Item1))).Max();

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

                        oldBestError = bestError;
                    }

                    double globalNormFactor = 0;

                    List<double> group1normintensities = new List<double>();

                    for (int p = 0; p < peptides.Count; p++)
                    {
                        var group1 = intensitiesByGroup[0][p].Select(v => v.Item1 - bestIonizationEfficiencyEstimations[p]);
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

        /// <summary>
        /// Called by the EstimateIonizationEfficienciesIfControlConditionHasOneSample method.
        /// </summary>
        private double CalculateIonizationEfficiencyEstimationError(List<List<List<(double, DetectionType)>>> intensitiesByGroup,
            double[] ionizationEfficiencyEstimations)
        {
            double totalError = 0;

            List<(double, DetectionType)> normalizedAbundances = new List<(double, DetectionType)>();

            // error is the sum of squared errors
            foreach (var group in intensitiesByGroup)
            {
                normalizedAbundances.Clear();

                for (int p = 0; p < group.Count; p++)
                {
                    var peptideIntensities = group[p];
                    var ionizationEfficiencyEstimationGuess = ionizationEfficiencyEstimations[p];

                    normalizedAbundances.AddRange(peptideIntensities.Select(v => (v.Item1 - ionizationEfficiencyEstimationGuess, v.Item2)));
                }

                double median = normalizedAbundances.Select(p => p.Item1).Median();

                double errorForGroup = 0;
                foreach (var element in normalizedAbundances)
                {
                    double error = Math.Abs(element.Item1 - median);

                    errorForGroup += error;
                }

                totalError += errorForGroup;
            }

            return totalError;
        }

        private void EstimateIonizationEfficiencies()
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
                                treatmentCondition, UseSharedPeptides, Results, PeptideToSampleQuantity);
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

        private static double GetUnbiasedSigma(double sigma, int N)
        {
            //https://www.jstor.org/stable/pdf/2682923.pdf
            if (N < 2)
            {
                return double.NaN;
            }

            double cn = 1 + (1.0 / (4 * (N - 1)));

            return sigma * cn;
        }

        private (double, double) DetermineNullHypothesisWidth(ProteinGroup protein, string treatmentCondition)
        {
            var proteinRes = (UnpairedProteinQuantResult)protein.ConditionToQuantificationResults[treatmentCondition];

            double nullHypothesisControl = Math.Sqrt(Math.Pow(FoldChangeCutoff, 2) / 2);
            double nullHypothesisTreatment = Math.Sqrt(Math.Pow(FoldChangeCutoff, 2) / 2);

            return (nullHypothesisControl, nullHypothesisTreatment);
        }

        private void EstimateIntensityDependentUncertainty()
        {
            var conditions = Results.SpectraFiles.Select(p => p.Condition).Distinct().ToList();

            foreach (string condition in conditions)
            {
                int numSamples = Results.SpectraFiles.Where(p => p.Condition == condition).Max(p => p.BiologicalReplicate) + 1;

                List<(double, double)> intensityToPeptideDiffToProtein = new List<(double, double)>();
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
                            if (PeptideToSampleQuantity[(peptide, condition, s)].Item1 > 0)
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
                            double intensity = PeptideToSampleQuantity[(peptide, condition, s)].Item1;

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
                List<(double, double)> meansAndVarianceOfMeans = new List<(double, double)>();

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

                if (peptideIntensitiesAndBiasesBinned.Any())
                {
                    Parallel.ForEach(Partitioner.Create(0, peptideIntensitiesAndBiasesBinned.Count),
                        new ParallelOptions { MaxDegreeOfParallelism = MaxThreads }, partitionRange =>
                        {
                            for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                            {
                                var theBin = peptideIntensitiesAndBiasesBinned[i];
                                var peptideMeans = theBin.Select(p => p.Item2).ToList();
                                double intensityEstimate = theBin.Min(p => p.Item1);

                                double mean = peptideMeans.Median();

                                double sd = Math.Min(peptideMeans.InterquartileRange() / 1.34896, peptideMeans.StandardDeviation());
                                sd = Math.Max(0.01, sd);

                                AdaptiveMetropolisWithinGibbs sampler = new AdaptiveMetropolisWithinGibbs(
                                    peptideMeans.ToArray(),
                                    new ProteinFoldChangeEstimationModel(
                                        priorMuMean: mean,
                                        priorMuSd: sd * 100,
                                        muInitialGuess: mean,
                                        sdPrior: new StudentT(sd, sd * 100, 1),
                                        sdInitialGuess: sd,
                                        nuPrior: new Exponential(1.0 / 29.0),
                                        nuInitialGuess: 5,
                                        minimumNu: 1,
                                        minimumSd: 0),
                                        seed: randomSeeds[i][0]
                                    );

                                // burn in and then sample the MCMC chain
                                sampler.Run(BurnInSteps, McmcSteps);

                                double[] mus = sampler.MarkovChain.Select(p => p[0]).ToArray();
                                double[] sds = sampler.MarkovChain.Select(p => p[1]).ToArray();
                                double[] nus = sampler.MarkovChain.Select(p => p[2]).ToArray();


                                double scale = sds.Median();

                                lock (meansAndVarianceOfMeans)
                                {
                                    meansAndVarianceOfMeans.Add((intensityEstimate, Math.Pow(scale, 2)));
                                }
                            }
                        });
                }

                meansAndVarianceOfMeans.Sort((x, y) => y.Item1.CompareTo(x.Item1));

                List<(double, double, double)> peptideVariances = new List<(double, double, double)>();
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

                if (peptideIntensitiesAndSigmasBinned.Any())
                {
                    Parallel.ForEach(Partitioner.Create(0, peptideIntensitiesAndSigmasBinned.Count),
                        new ParallelOptions { MaxDegreeOfParallelism = MaxThreads }, partitionRange =>
                        {
                            for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                            {
                                var theBin = peptideIntensitiesAndSigmasBinned[i];
                                var sigmas = theBin.Select(p => p.Item2).ToArray();
                                double intensityEstimate = theBin.Min(p => p.Item1);

                                double iqr15 = sigmas.InterquartileRange() * 1.5;
                                double median = sigmas.Median();
                                double q1 = sigmas.Quantile(0.25);
                                double q3 = sigmas.Quantile(0.75);
                                var trimmedSds = sigmas.Where(p => p > q1 - iqr15 && p < q3 + iqr15).ToList();

                                double medianSd = trimmedSds.Median();
                                double medianVariance = Math.Pow(medianSd, 2);
                                double binVarianceOfSds = trimmedSds.Variance();

                                lock (peptideVariances)
                                {
                                    peptideVariances.Add((intensityEstimate, medianVariance, binVarianceOfSds));
                                }
                            }
                        });
                }

                peptideVariances.Sort((x, y) => y.Item1.CompareTo(x.Item1));
                CreateProteinSigmaPriors(condition, meansAndVarianceOfMeans, peptideVariances);
            }
        }

        private void CreateProteinSigmaPriors(string condition, List<(double, double)> varianceOfPeptideMeans,
            List<(double intensity, double variance, double varianceOfSd)> peptideVariances)
        {
            // meansAndVarianceOfMeans will be empty if every protein has only 1 peptide (e.g., top-down proteoform quant)
            // peptideVariancesArray will be empty if there is only 1 sample per condition, because within-group variance is unknown

            var list = ProteinsWithConstituentPeptides.ToList();

            if (ProteinAndConditionToSigmaPrior == null)
            {
                ProteinAndConditionToSigmaPrior = new Dictionary<(ProteinGroup, string), StudentT>();
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
                             List<double> peptideVariancesForThisProtein = new List<double>();
                             List<double> variancesOfPeptideMeansForThisProtein = new List<double>();
                             List<double> variancesOfPeptideSdForThisProtein = new List<double>();

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
                                     double intensity = PeptideToSampleQuantity[(peptide, condition, s)].Item1;

                                     if (intensity > 0)
                                     {
                                         intensities.Add(intensity);
                                     }
                                 }

                                 if (intensities.Count == 0)
                                 {
                                     continue;
                                 }

                                 foreach (var intensity in intensities)
                                 {
                                     double logIntensity = Math.Log(intensity, 2);
                                     double? bestVarianceOfMeans = null;
                                     double? bestPeptideVariances = null;
                                     double? bestVarianceOfVariances = null;

                                     for (int j = 0; j < varianceOfPeptideMeans.Count; j++)
                                     {
                                         bestVarianceOfMeans = varianceOfPeptideMeans[j].Item2;

                                         if (j == varianceOfPeptideMeans.Count - 1 || logIntensity > varianceOfPeptideMeans[j + 1].Item1)
                                         {
                                             break;
                                         }
                                     }

                                     for (int j = 0; j < peptideVariances.Count; j++)
                                     {
                                         bestPeptideVariances = peptideVariances[j].Item2;
                                         bestVarianceOfVariances = peptideVariances[j].varianceOfSd;

                                         if (j == peptideVariances.Count - 1 || logIntensity > peptideVariances[j + 1].Item1)
                                         {
                                             break;
                                         }
                                     }

                                     if (bestVarianceOfMeans != null)
                                     {
                                         variancesOfPeptideMeansForThisProtein.Add(bestVarianceOfMeans.Value);
                                     }

                                     if (bestPeptideVariances != null)
                                     {
                                         peptideVariancesForThisProtein.Add(bestPeptideVariances.Value);
                                         variancesOfPeptideSdForThisProtein.Add(bestVarianceOfVariances.Value);
                                     }
                                 }
                             }

                             double averageOfVariances = peptideVariancesForThisProtein.Count > 0 ? peptideVariancesForThisProtein.Average() : 0;
                             double varianceOfMeans = 0;

                             if (variancesOfPeptideMeansForThisProtein.Count == 0)
                             {
                                 varianceOfMeans = 0;
                             }
                             else
                             {
                                 varianceOfMeans = variancesOfPeptideMeansForThisProtein.Average();
                             }

                             double pooledVariance = averageOfVariances + varianceOfMeans;

                             if (pooledVariance == 0)
                             {
                                 continue;
                             }

                             double pooledVarianceOfSd = variancesOfPeptideSdForThisProtein.Count > 1
                                ? variancesOfPeptideSdForThisProtein.Average() + variancesOfPeptideSdForThisProtein.Variance()
                                : variancesOfPeptideSdForThisProtein.Count > 0 ? variancesOfPeptideSdForThisProtein.Average() : pooledVariance;

                             double pooledSigma = Math.Sqrt(pooledVariance);
                             double sigmaOfSigma = Math.Max(Math.Sqrt(pooledVarianceOfSd), 0.1);

                             var proteinSigmaPrior = new StudentT(pooledSigma, sigmaOfSigma, 1);

                             lock (ProteinAndConditionToSigmaPrior)
                             {
                                 ProteinAndConditionToSigmaPrior[(protein, condition)] = proteinSigmaPrior;
                             }
                         }
                     });
        }

        /// <summary>
        /// Calculates the false discovery rate of each protein from the Bayes factors.
        /// https://arxiv.org/pdf/1311.3981.pdf
        /// </summary>
        private void CalculateFalseDiscoveryRates()
        {
            var proteinList = ProteinsWithConstituentPeptides.ToList();

            // calculate FDR for the condition
            foreach (string treatmentCondition in TreatmentConditions)
            {
                var bayesianQuantResults = proteinList
                    .Select(p => p.Key.ConditionToQuantificationResults[treatmentCondition])
                    .OrderByDescending(p => p.IsStatisticallyValid)
                    .ThenByDescending(p => p.BayesFactor)
                    .ThenByDescending(p => p.Peptides.Count)
                    .ToList();

                var bayesFactorsAscending = bayesianQuantResults
                    .Where(p => p.IsStatisticallyValid)
                    .Select(p => p.BayesFactor)
                    .OrderBy(p => p)
                    .ToList();

                var validResults = bayesianQuantResults.Where(p => p.IsStatisticallyValid).ToList();

                // pi_0 is the proportion of false-positives in the set of hypothesis tests
                // it can be estimated by finding the largest set of tests where the average bayes factor is less than one (https://arxiv.org/pdf/1311.3981.pdf)
                double pi_0_bayesFactors = 0;

                var runningListOfBayesFactors = new List<double>();

                for (int i = 0; i < bayesFactorsAscending.Count; i++)
                {
                    double numTestsSoFar = i + 1;
                    runningListOfBayesFactors.Add(bayesFactorsAscending[i]);
                    double averageBayesFactor = runningListOfBayesFactors.Average();

                    if (averageBayesFactor < 1)
                    {
                        pi_0_bayesFactors = numTestsSoFar / bayesFactorsAscending.Count;
                    }
                }

                // the above pi_0 can be an underestimate if uncertainty in the data is high
                // pi_0 can be estimated by the proportion of data where the mean estimate is less than the uncertainty in the mean
                double pi_0_uncertainty = validResults.Count(p => Math.Abs(p.FoldChangePointEstimate) < p.UncertaintyInFoldChangeEstimate || Math.Abs(p.FoldChangePointEstimate) < FoldChangeCutoff)
                    / (double)validResults.Count;

                // pi_0 is the maximum of the two above pi_0 estimates
                double pi_0 = Math.Max(pi_0_uncertainty, pi_0_bayesFactors);

                foreach (var bayesianQuantResult in bayesianQuantResults)
                {
                    double bayesFactor = bayesianQuantResult.BayesFactor;

                    // vi is the probability that the alternative hypothesis is true, if accepted
                    double vi = ((1 - pi_0) * bayesFactor) / (pi_0 + (1 - pi_0) * bayesFactor);

                    bayesianQuantResult.PosteriorErrorProbability = 1 - vi;
                }

                List<double> PosteriorErrorProbabilities = new List<double>();

                for (int p = 0; p < bayesianQuantResults.Count; p++)
                {
                    PosteriorErrorProbabilities.Add(bayesianQuantResults[p].PosteriorErrorProbability);
                    bayesianQuantResults[p].FalseDiscoveryRate = PosteriorErrorProbabilities.Average();
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
                        double intensity = PeptideToSampleQuantity[(peptide, ControlCondition, b)].Item1;
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