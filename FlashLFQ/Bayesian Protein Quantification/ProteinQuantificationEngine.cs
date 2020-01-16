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
        private readonly double? NullHypothesisInterval;
        private readonly int BurnInSteps;
        private readonly int McmcSteps;
        private readonly MersenneTwister Rng;
        private Dictionary<ProteinGroup, List<Peptide>> ProteinsWithConstituentPeptides;
        private Dictionary<(Peptide, string, int), double> PeptideToSampleQuantity;
        private Dictionary<string, List<(double, double)>> PeptideIntensityNoiseEstimates;
        public Dictionary<(ProteinGroup, string), (double, double)> ProteinAndConditionToNullHypothesisInterval { get; private set; }
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
        }

        /// <summary>
        /// Runs the protein quantification engine.
        /// </summary>
        public void Run()
        {
            Setup();
            DeterminePeptideNoiseValues();

            // generate a random seed for each protein for the MCMC sampler
            var randomSeedsForEachProtein = new Dictionary<ProteinGroup, List<int>>();
            foreach (var protein in ProteinsWithConstituentPeptides)
            {
                randomSeedsForEachProtein.Add(protein.Key, new List<int> { Rng.NextFullRangeInt32(), Rng.NextFullRangeInt32() });

                foreach (var peptide in protein.Value)
                {
                    randomSeedsForEachProtein[protein.Key].Add(Rng.NextFullRangeInt32());
                }
            }

            var proteinList = ProteinsWithConstituentPeptides.ToList();

            // figure out which conditions to compare
            var treatmentConditions = Results.SpectraFiles.Select(p => p.Condition).Distinct().ToList();
            treatmentConditions.Remove(ControlCondition);

            // calculate fold-change for each protein (diffuse prior)
            foreach (string treatmentCondition in treatmentConditions)
            {
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

                            ((PairedProteinQuantResult)result).EstimateProteinFoldChange(randomSeedsForEachProtein[protein][0], BurnInSteps, McmcSteps);
                        }
                        else
                        {
                            result = new UnpairedProteinQuantResult(protein, ProteinsWithConstituentPeptides[protein], ControlCondition, treatmentCondition, UseSharedPeptides,
                                Results, PeptideToSampleQuantity, randomSeedsForEachProtein[protein], BurnInSteps);
                            ((UnpairedProteinQuantResult)result).EstimateProteinFoldChange(randomSeedsForEachProtein[protein][0], BurnInSteps, McmcSteps);
                        }

                        protein.ConditionToQuantificationResults.Add(treatmentCondition, result);
                    }
                });
            }

            DetermineExperimentalNull();

            foreach (string treatmentCondition in treatmentConditions)
            {
                //double nullHypothesisControl = double.NaN;
                //double nullHypothesisTreatment;

                //if (!PairedSamples)
                //{
                //    if (NullHypothesisInterval == null)
                //    {
                //        var expNull = DetermineExperimentalNull(treatmentCondition);
                //        nullHypothesisControl = expNull.Item1;
                //        nullHypothesisTreatment = expNull.Item2;
                //    }
                //    else
                //    {
                //        nullHypothesisControl = Math.Sqrt(Math.Pow(NullHypothesisInterval.Value, 2) / 2.0);
                //        nullHypothesisTreatment = Math.Sqrt(Math.Pow(NullHypothesisInterval.Value, 2) / 2.0);
                //    }
                //}
                //else
                //{
                //    nullHypothesisTreatment = NullHypothesisInterval == null ? DetermineExperimentalNull(treatmentCondition).Item2 : NullHypothesisInterval.Value;
                //}

                Parallel.ForEach(Partitioner.Create(0, proteinList.Count), new ParallelOptions { MaxDegreeOfParallelism = MaxThreads }, partitionRange =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        ProteinGroup protein = proteinList[i].Key;

                        (double, double) nullHypInterval = ProteinAndConditionToNullHypothesisInterval[(protein, treatmentCondition)];
                        double nullHypControl = nullHypInterval.Item1;
                        double nullHypTreatment = nullHypInterval.Item2;

                        ProteinQuantificationEngineResult result = protein.ConditionToQuantificationResults[treatmentCondition];

                        if (PairedSamples)
                        {
                            var pairedResult = (PairedProteinQuantResult)result;
                            pairedResult.EstimateProteinFoldChange(randomSeedsForEachProtein[protein][1], BurnInSteps, McmcSteps, nullHypTreatment);
                        }
                        else
                        {
                            var unpairedResult = (UnpairedProteinQuantResult)result;
                            unpairedResult.EstimateProteinFoldChange(randomSeedsForEachProtein[protein][1], BurnInSteps, McmcSteps, nullHypControl, nullHypTreatment);
                        }
                    }
                });
            }

            foreach (string treatmentCondition in treatmentConditions)
            {
                Parallel.ForEach(Partitioner.Create(0, proteinList.Count), new ParallelOptions { MaxDegreeOfParallelism = MaxThreads }, partitionRange =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        ProteinGroup protein = proteinList[i].Key;

                        (double, double) nullHypInterval = ProteinAndConditionToNullHypothesisInterval[(protein, treatmentCondition)];
                        double nullHypControl = nullHypInterval.Item1;
                        double nullHypTreatment = nullHypInterval.Item2;

                        ProteinQuantificationEngineResult result = protein.ConditionToQuantificationResults[treatmentCondition];

                        if (PairedSamples)
                        {
                            var pairedResult = (PairedProteinQuantResult)result;
                            pairedResult.EstimateProteinFoldChange(randomSeedsForEachProtein[protein][1], BurnInSteps, McmcSteps, nullHypTreatment);
                        }
                        else
                        {
                            var unpairedResult = (UnpairedProteinQuantResult)result;
                            unpairedResult.EstimateProteinFoldChange(randomSeedsForEachProtein[protein][1], BurnInSteps, McmcSteps, nullHypControl, nullHypTreatment);
                        }
                    }
                });
            }

            // calculate FDR for the condition
            foreach (var treatmentCondition in treatmentConditions)
            {
                var res = proteinList.Select(p => p.Key.ConditionToQuantificationResults[treatmentCondition]).ToList();
                CalculateFalseDiscoveryRates(res);
            }
        }

        public static (double[] mus, double[] sds, double[] nus) FitProteinQuantModel(List<double> measurements, bool skepticalPrior, bool paired,
            int? randomSeed, int burnin, int n, double? nullHypothesisCutoff)
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
            double priorMuSd = Math.Max(0.05, skepticalPrior ? nullHypothesisCutoff.Value : sd * 100);
            double priorMuMean = skepticalPrior ? 0 : mean;

            AdaptiveMetropolisWithinGibbs sampler = new AdaptiveMetropolisWithinGibbs(
                measurements.ToArray(),
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

            double[] mus = sampler.MarkovChain.Select(p => p[0]).ToArray();
            double[] sds = sampler.MarkovChain.Select(p => p[1]).ToArray();
            double[] nus = sampler.MarkovChain.Select(p => p[2]).ToArray();

            //DEBUG
            double muEst = mus.Median();
            double sdEst = sds.Median();
            double nuEst = nus.Median();

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
        }

        /// <summary>
        /// This function determines a fold-change cutoff (a null hypothesis) from the data. If enough peptide-level
        /// data is available, a null hypothesis interval is determined for each protein from intensity-dependent 
        /// peptide noise estimates.
        /// </summary>
        private void DetermineExperimentalNull()
        {
            ProteinAndConditionToNullHypothesisInterval = new Dictionary<(ProteinGroup, string), (double, double)>();
            var conditions = this.Results.SpectraFiles.Select(p => p.Condition).Distinct().ToList();

            // if a user-defined null hypothesis has been provided, use that
            if (NullHypothesisInterval != null)
            {
                (double, double) globalNull;

                if (PairedSamples)
                {
                    globalNull = (NullHypothesisInterval.Value, double.NaN);
                }
                else
                {
                    double nullHyp = Math.Sqrt(Math.Pow(NullHypothesisInterval.Value, 2) / 2.0);
                    globalNull = (nullHyp, nullHyp);
                }

                foreach (string condition in conditions.Where(p => p != ControlCondition))
                {
                    foreach (ProteinGroup protein in ProteinsWithConstituentPeptides.Keys)
                    {
                        ProteinAndConditionToNullHypothesisInterval.Add((protein, condition), globalNull);
                    }
                }

                return;
            }

            // per-protein null interval based on peptide-level noise data (preferred)
            if (PeptideIntensityNoiseEstimates.All(p => p.Value != null))
            {
                foreach (string treatmentCondition in conditions.Where(p => p != ControlCondition))
                {
                    foreach (var proteinWithPeptides in ProteinsWithConstituentPeptides)
                    {
                        ProteinGroup protein = proteinWithPeptides.Key;
                        List<Peptide> peptides = proteinWithPeptides.Value;

                        (double, double) proteinSpecificNull = (0, 0);

                        if (PairedSamples)
                        {
                            //TODO
                        }
                        else
                        {
                            double controlNullInterval = 0;
                            double treatmentNullInterval = 0;
                            var controlPeptides = ((UnpairedProteinQuantResult)protein.ConditionToQuantificationResults[treatmentCondition]).ConditionToPeptidesUsedForQuant[ControlCondition];
                            var treatmentPeptides = ((UnpairedProteinQuantResult)protein.ConditionToQuantificationResults[treatmentCondition]).ConditionToPeptidesUsedForQuant[treatmentCondition];

                            // treatment
                            List<StudentT> normals = new List<StudentT>();

                            foreach (var peptide in treatmentPeptides)
                            {
                                int nMeasurements = 0;
                                int samples = Results.SpectraFiles.Where(p => p.Condition == treatmentCondition).Max(p => p.BiologicalReplicate) + 1;
                                for (int i = 0; i < samples; i++)
                                {
                                    if (PeptideToSampleQuantity.TryGetValue((peptide, treatmentCondition, i), out double intensity) && intensity > 0)
                                    {
                                        nMeasurements++;
                                    }
                                }

                                var log2Noise = GetNoiseFromPeptide(peptide, treatmentCondition);

                                if (log2Noise == null)
                                {
                                    continue;
                                }

                                StudentT n = new StudentT(0, log2Noise.Value, 1.0);

                                for (int i = 0; i < nMeasurements; i++)
                                {
                                    normals.Add(n);
                                }
                            }

                            if (!normals.Any())
                            {
                                treatmentNullInterval = 1.0;
                            }
                            else
                            {
                                for (double n = 0.01; n < 100; n += 0.01)
                                {
                                    double total = 0;

                                    foreach (StudentT normal in normals)
                                    {
                                        double area = normal.CumulativeDistribution(n) - normal.CumulativeDistribution(-n);
                                        total += area;
                                    }

                                    double fraction = total / normals.Count;

                                    if (fraction >= 0.5)
                                    {
                                        treatmentNullInterval = n;
                                        break;
                                    }
                                }
                            }

                            // control
                            normals = new List<StudentT>();

                            foreach (var peptide in controlPeptides)
                            {
                                int nMeasurements = 0;
                                int samples = Results.SpectraFiles.Where(p => p.Condition == ControlCondition).Max(p => p.BiologicalReplicate) + 1;
                                for (int i = 0; i < samples; i++)
                                {
                                    if (PeptideToSampleQuantity.TryGetValue((peptide, ControlCondition, i), out double intensity) && intensity > 0)
                                    {
                                        nMeasurements++;
                                    }
                                }

                                var log2Noise = GetNoiseFromPeptide(peptide, ControlCondition);

                                if (log2Noise == null)
                                {
                                    continue;
                                }

                                StudentT n = new StudentT(0, log2Noise.Value, 1.0);

                                for (int i = 0; i < nMeasurements; i++)
                                {
                                    normals.Add(n);
                                }
                            }

                            if (!normals.Any())
                            {
                                controlNullInterval = 1.0;
                            }
                            else
                            {
                                for (double n = 0.01; n < 100; n += 0.01)
                                {
                                    double total = 0;

                                    foreach (StudentT normal in normals)
                                    {
                                        double area = normal.CumulativeDistribution(n) - normal.CumulativeDistribution(-n);
                                        total += area;
                                    }

                                    double fraction = total / normals.Count;

                                    if (fraction >= 0.5)
                                    {
                                        controlNullInterval = n;
                                        break;
                                    }
                                }
                            }

                            proteinSpecificNull = (controlNullInterval, treatmentNullInterval);
                        }

                        ProteinAndConditionToNullHypothesisInterval.Add((protein, treatmentCondition), proteinSpecificNull);
                    }
                }

                return;
            }

            // if not enough peptide-level data (or not enough samples) to determine per-protein null interval (not preferred)
            foreach (string condition in conditions.Where(p => p != ControlCondition))
            {
                (double, double) globalNull;

                if (PairedSamples)
                {
                    double treatmentExperimentalNull = 1.0;
                    double controlExperimentalNull = double.NaN;

                    var res = Results.ProteinGroups.Values.Select(p => p.ConditionToQuantificationResults[condition]).ToList();

                    List<double> stdDevs = res
                        .Where(p => p.StandardDeviationPointEstimate > 0 && !double.IsNaN(p.StandardDeviationPointEstimate))
                        .Select(p => p.StandardDeviationPointEstimate)
                        .ToList();

                    if (stdDevs.Any())
                    {
                        var tenPercentHdi = Util.GetHighestDensityInterval(stdDevs.ToArray(), 0.1);
                        var sdsWithinHdi = stdDevs.Where(p => p <= tenPercentHdi.hdi_end && p >= tenPercentHdi.hdi_start).ToList();

                        var modeStdDevPointEstimate = sdsWithinHdi.Median();
                        treatmentExperimentalNull = modeStdDevPointEstimate;

                        if (double.IsNaN(treatmentExperimentalNull))
                        {
                            treatmentExperimentalNull = stdDevs.Median();
                        }
                    }

                    globalNull = (controlExperimentalNull, treatmentExperimentalNull);
                }
                else
                {
                    double treatmentExperimentalNull = 1.0;
                    double controlExperimentalNull = 1.0;

                    var res = Results.ProteinGroups.Values.Select(p => (UnpairedProteinQuantResult)p.ConditionToQuantificationResults[condition]).ToList();

                    List<double> stdDevsControl = res
                        .Where(p => p.StandardDeviationPointEstimate > 0 && !double.IsNaN(p.StandardDeviationPointEstimate))
                        .Select(p => p.ConditionToSdEstimates[ControlCondition].Median())
                        .ToList();

                    List<double> stdDevsTreatment = res
                        .Where(p => p.StandardDeviationPointEstimate > 0 && !double.IsNaN(p.StandardDeviationPointEstimate))
                        .Select(p => p.ConditionToSdEstimates[condition].Median())
                        .ToList();

                    if (stdDevsControl.Any())
                    {
                        var tenPercentHdi = Util.GetHighestDensityInterval(stdDevsControl.ToArray(), 0.1);
                        var sdsWithinHdi = stdDevsControl.Where(p => p <= tenPercentHdi.hdi_end && p >= tenPercentHdi.hdi_start).ToList();

                        var modeStdDevPointEstimate = sdsWithinHdi.Median();
                        controlExperimentalNull = modeStdDevPointEstimate;

                        if (double.IsNaN(controlExperimentalNull))
                        {
                            controlExperimentalNull = stdDevsControl.Median();
                        }
                    }

                    if (stdDevsTreatment.Any())
                    {
                        var tenPercentHdi = Util.GetHighestDensityInterval(stdDevsTreatment.ToArray(), 0.1);
                        var sdsWithinHdi = stdDevsTreatment.Where(p => p <= tenPercentHdi.hdi_end && p >= tenPercentHdi.hdi_start).ToList();

                        var modeStdDevPointEstimate = sdsWithinHdi.Median();
                        treatmentExperimentalNull = modeStdDevPointEstimate;

                        if (double.IsNaN(treatmentExperimentalNull))
                        {
                            treatmentExperimentalNull = stdDevsTreatment.Median();
                        }
                    }

                    globalNull = (controlExperimentalNull, treatmentExperimentalNull);
                }

                foreach (ProteinGroup protein in ProteinsWithConstituentPeptides.Keys)
                {
                    ProteinAndConditionToNullHypothesisInterval.Add((protein, condition), globalNull);
                }
            }
        }

        private void DeterminePeptideNoiseValues()
        {
            var peptides = Results.PeptideModifiedSequences.Values.ToList();
            PeptideIntensityNoiseEstimates = new Dictionary<string, List<(double, double)>>();

            var conditions = Results.SpectraFiles.Select(p => p.Condition).Distinct().ToList();

            foreach (var condition in conditions)
            {
                PeptideIntensityNoiseEstimates.Add(condition, new List<(double, double)>());
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

                    for (int s = 0; s < samples; s++)
                    {
                        if (PeptideToSampleQuantity.TryGetValue((peptide, condition, s), out double intensity) && intensity > 0)
                        {
                            intensities.Add(Math.Log(intensity, 2));
                        }
                    }

                    double median = intensities.Median();
                    double stdev = intensities.StandardDeviation();

                    if (!double.IsNaN(median) && !double.IsNaN(stdev))
                    {
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

                            var stdevs = theBin.Select(p => p.Item2).ToList();
                            var res = FitProteinQuantModel(stdevs, false, false, randomSeeds[i], BurnInSteps, 1000, 0);
                            double stdevEstimate = res.mus.Median();
                            double intensityEstimate = theBin.Select(p => p.Item1).Median();

                            var estimates = PeptideIntensityNoiseEstimates[condition];

                            lock (estimates)
                            {
                                estimates.Add((intensityEstimate, stdevEstimate));
                            }
                        }
                    });

                PeptideIntensityNoiseEstimates[condition].Sort((x, y) => y.Item1.CompareTo(x.Item1));
            }
        }

        private double? GetNoiseFromPeptide(Peptide peptide, string condition)
        {
            int numSamples = Results.SpectraFiles.Where(p => p.Condition == condition).Max(p => p.BiologicalReplicate) + 1;

            List<double> peptideIntensities = new List<double>();
            for (int s = 0; s < numSamples; s++)
            {
                peptideIntensities.Add(PeptideToSampleQuantity[(peptide, condition, s)]);
            }

            if (peptideIntensities.Count < 1)
            {
                return null;
            }

            double medianIntensity = peptideIntensities.Median();

            if (double.IsNaN(medianIntensity) || double.IsInfinity(medianIntensity))
            {
                return null;
            }

            double logIntensity = Math.Log(medianIntensity, 2);

            double noiseEstimate = 0;

            var noiseEstimatesForCondition = PeptideIntensityNoiseEstimates[condition];

            for (int i = 0; i < noiseEstimatesForCondition.Count; i++)
            {
                noiseEstimate = noiseEstimatesForCondition[i].Item2;

                if (i == noiseEstimatesForCondition.Count - 1 || logIntensity > noiseEstimatesForCondition[i + 1].Item1)
                {
                    break;
                }
            }

            return noiseEstimate;
        }

        /// <summary>
        /// Calculates the false discovery rate of each protein from the Bayesian-estimated PEP values.
        /// </summary>
        private void CalculateFalseDiscoveryRates(List<ProteinQuantificationEngineResult> bayesianProteinQuantificationResults)
        {
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
}
