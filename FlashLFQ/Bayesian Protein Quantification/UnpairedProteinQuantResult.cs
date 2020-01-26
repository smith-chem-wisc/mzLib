using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace FlashLFQ
{
    public class UnpairedProteinQuantResult : ProteinQuantificationEngineResult
    {
        public readonly Dictionary<string, List<double>> ConditionsWithPeptideSampleQuantities;
        public readonly Dictionary<string, double[]> ConditionToMeanEstimates;
        public readonly Dictionary<string, double[]> ConditionToSdEstimates;
        public readonly Dictionary<string, double[]> ConditionToSkepticalPriorMuEstimates;
        public readonly Dictionary<string, HashSet<Peptide>> ConditionToPeptidesUsedForQuant;
        public double SigmaControl;
        public double SigmaTreatment;
        public double NuControl;
        public double NuTreatment;
        public double FractionMbrControl;
        public double FractionMbrTreatment;
        
        public UnpairedProteinQuantResult(ProteinGroup protein, List<Peptide> peptides, string controlCondition, string treatmentCondition,
            bool useSharedPeptides, FlashLfqResults results, Dictionary<(Peptide, string, int), double> peptideToSampleQuantity, List<int> randomSeeds, int nBurnin)
            : base(protein, peptides, controlCondition, treatmentCondition)
        {
            ConditionsWithPeptideSampleQuantities = new Dictionary<string, List<double>>();
            ConditionToMeanEstimates = new Dictionary<string, double[]>();
            ConditionToSdEstimates = new Dictionary<string, double[]>();
            ConditionToSkepticalPriorMuEstimates = new Dictionary<string, double[]>();
            ConditionToPeptidesUsedForQuant = new Dictionary<string, HashSet<Peptide>>();
            
            GetPeptideSampleQuantities(useSharedPeptides, results, peptideToSampleQuantity, randomSeeds, nBurnin);
        }

        public void EstimateProteinFoldChange(int? randomSeed, int? randomSeed2, int nBurnin, int nSteps,
            double controlNullHypothesisInterval = double.NaN, double treatmentNullHypothesisInterval = double.NaN,
            MultiGammaDistribution controlSigmaPrior = null, MultiGammaDistribution treatmentSigmaPrior = null)
        {
            IsStatisticallyValid = this.ConditionsWithPeptideSampleQuantities[ControlCondition].Count > 1
                && this.ConditionsWithPeptideSampleQuantities[TreatmentCondition].Count > 1;

            bool skepticalPrior = !double.IsNaN(controlNullHypothesisInterval) || !double.IsNaN(treatmentNullHypothesisInterval);

            var conditionToMus = skepticalPrior ? ConditionToSkepticalPriorMuEstimates : ConditionToMeanEstimates;

            double[] controlSds = null;
            double[] treatmentSds = null;

            if (!conditionToMus.TryGetValue(ControlCondition, out double[] controlMus))
            {
                var samples = ConditionsWithPeptideSampleQuantities[ControlCondition];
                
                var (mus, sds, nus) = ProteinQuantificationEngine.FitProteinQuantModel(samples, skepticalPrior, false, randomSeed, nBurnin, nSteps,
                    controlNullHypothesisInterval, controlSigmaPrior);

                controlMus = mus;
                conditionToMus.Add(ControlCondition, controlMus);

                if (conditionToMus == ConditionToMeanEstimates)
                {
                    controlSds = sds;
                    ConditionToSdEstimates.Add(ControlCondition, controlSds);

                    NuControl = nus.Median();
                }
            }

            if (!conditionToMus.TryGetValue(TreatmentCondition, out double[] treatmentMus))
            {
                var samples = ConditionsWithPeptideSampleQuantities[TreatmentCondition];
                
                var (mus, sds, nus) = ProteinQuantificationEngine.FitProteinQuantModel(samples, skepticalPrior, false, randomSeed2, nBurnin, nSteps,
                    treatmentNullHypothesisInterval, treatmentSigmaPrior);

                treatmentMus = mus;
                conditionToMus.Add(TreatmentCondition, treatmentMus);

                if (conditionToMus == ConditionToMeanEstimates)
                {
                    treatmentSds = sds;
                    ConditionToSdEstimates.Add(TreatmentCondition, treatmentSds);

                    NuTreatment = nus.Median();
                }
            }

            double[] diffs = new double[Math.Min(controlMus.Length, treatmentMus.Length)];

            for (int i = 0; i < diffs.Length; i++)
            {
                diffs[i] = treatmentMus[i] - controlMus[i];
            }

            if (!IsStatisticallyValid)
            {
                for (int i = 0; i < diffs.Length; i++)
                {
                    diffs[i] = 0;
                }
            }

            if (!skepticalPrior)
            {
                this.FoldChangePointEstimate = diffs.Median();
                double a_1Estimate = controlSds.Median();
                double b_1Estimate = treatmentSds.Median();

                SigmaControl = a_1Estimate;
                SigmaTreatment = b_1Estimate;

                double propagatedSd = Math.Sqrt(Math.Pow(a_1Estimate, 2) + Math.Pow(b_1Estimate, 2));
                this.StandardDeviationPointEstimate = propagatedSd;
            }
            else
            {
                this.NullHypothesisInterval = Math.Sqrt(Math.Pow(controlNullHypothesisInterval, 2) + Math.Pow(treatmentNullHypothesisInterval, 2));
                this.CalculatePosteriorErrorProbability(diffs);
            }
        }

        public override string ToString()
        {
            int nPeptides = Peptides.Count;
            int nMeasurementsControl = ConditionsWithPeptideSampleQuantities[ControlCondition].Count;
            int nMeasurementsTreatment = ConditionsWithPeptideSampleQuantities[TreatmentCondition].Count;

            var measurementsString = "NOT ENABLED";

            if (measurementsString.Length > 32000)
            {
                measurementsString = "Too long for Excel";
            }

            return
                Protein.ProteinGroupName + "\t" +
                Protein.GeneName + "\t" +
                Protein.Organism + "\t" +
                ControlCondition + "\t" +
                TreatmentCondition + "\t" +
                NullHypothesisInterval.Value + "\t" +
                //PropagatedUncertaintyInDiff + "\t" +
                FoldChangePointEstimate + "\t" +
                StandardDeviationPointEstimate + "\t" +
                NuControl + "\t" +
                NuTreatment + "\t" +
                //PredictedSigmaControl + "\t" +
                SigmaControl + "\t" +
                //PredictedSigmaTreatment + "\t" +
                SigmaTreatment + "\t" +
                FractionMbrControl + "\t" +
                FractionMbrTreatment + "\t" +
                //UncertaintyInMean + "\t" +
                ControlConditionIntensity + "\t" +
                TreatmentConditionIntensity + "\t" +
                nPeptides + "\t" +
                nMeasurementsControl + "\t" +
                nMeasurementsTreatment + "\t" +
                measurementsString + "\t" +
                PosteriorErrorProbability + "\t" +
                //DiffExpectedPep + "\t" +
                FalseDiscoveryRate + "\t\t";
        }

        public static new string TabSeparatedHeader()
        {
            return
                "Protein Group" + "\t" +
                "Gene" + "\t" +
                "Organism" + "\t" +
                "Control Condition" + "\t" +
                "Treatment Condition" + "\t" +
                "Log2 Fold-Change Cutoff" + "\t" +
                //"Propagated Uncertainty" + "\t" +
                "Protein Log2 Fold-Change" + "\t" +
                "Standard Deviation of Peptide Log2 Fold-Changes" + "\t" +
                "Nu Control" + "\t" +
                "Nu Treatment" + "\t" +
                //"Predicted STDEV Control" + "\t" +
                "Sigma Control" + "\t" +
                //"Predicted STDEV Treatment" + "\t" +
                "Sigma Treatment" + "\t" +
                "Fraction MBR Control" + "\t" +
                "Fraction MBR Treatment" + "\t" +
                //"Uncertainty in log2 differences" + "\t" +
                "Protein Intensity in Control Condition" + "\t" +
                "Protein Intensity in Treatment Condition" + "\t" +
                "Number of Peptides" + "\t" +
                "Number of Control Condition Measurements" + "\t" +
                "Number of Treatment Condition Measurements" + "\t" +
                "List of Fold-Change Measurements Grouped by Peptide" + "\t" +
                "Posterior Error Probability" + "\t" +
                //"PEP Error" + "\t" +
                "False Discovery Rate" + "\t\t";
        }

        private void GetPeptideSampleQuantities(bool useSharedPeptides, FlashLfqResults flashLfqResults,
            Dictionary<(Peptide, string, int), double> PeptideToSampleQuantity, List<int> randomSeeds, int nBurnin)
        {
            int numMsms = 0;
            int numMbr = 0;

            if (!ConditionsWithPeptideSampleQuantities.ContainsKey(ControlCondition))
            {
                ConditionToPeptidesUsedForQuant.Add(ControlCondition, new HashSet<Peptide>());

                int numSamples = flashLfqResults.SpectraFiles.Where(p => p.Condition == ControlCondition).Max(p => p.BiologicalReplicate) + 1;
                ConditionsWithPeptideSampleQuantities.Add(ControlCondition, new List<double>());
                List<double> peptideSampleLogIntensities = new List<double>();

                for (int i = 0; i < Peptides.Count; i++)
                {
                    Peptide peptide = Peptides[i];
                    peptideSampleLogIntensities.Clear();
                    int randomSeed = randomSeeds[i];

                    if (!peptide.UseForProteinQuant || (!useSharedPeptides && peptide.ProteinGroups.Count > 1))
                    {
                        continue;
                    }

                    for (int sample = 0; sample < numSamples; sample++)
                    {
                        double intensity = PeptideToSampleQuantity[(peptide, ControlCondition, sample)];

                        if (intensity > 0)
                        {
                            peptideSampleLogIntensities.Add(Math.Log(intensity, 2));

                            //DEBUG
                            var detectionType = flashLfqResults.PeptideModifiedSequences[peptide.Sequence].GetDetectionType(
                                flashLfqResults.SpectraFiles.First(p => p.Condition == ControlCondition && p.BiologicalReplicate == sample));
                            if (detectionType == DetectionType.MSMS)
                            {
                                numMsms++;
                            }
                            else if (detectionType == DetectionType.MBR)
                            {
                                numMbr++;
                            }
                        }
                    }

                    // estimate ionization efficiency
                    if (peptide.IonizationEfficiency == 0)
                    {
                        if (peptideSampleLogIntensities.Count > 1)
                        {
                            var (mus, sds, nus) = ProteinQuantificationEngine.FitProteinQuantModel(peptideSampleLogIntensities, false, false, randomSeed, nBurnin, 1000, 0, null);
                            peptide.IonizationEfficiency = Math.Pow(2, mus.Median());
                        }
                        else if (peptideSampleLogIntensities.Count == 1)
                        {
                            peptide.IonizationEfficiency = Math.Pow(2, peptideSampleLogIntensities.First());
                        }
                    }

                    if (peptide.IonizationEfficiency != 0)
                    {
                        foreach (double logIntensity in peptideSampleLogIntensities)
                        {
                            double linearIntensity = Math.Pow(2, logIntensity);

                            ConditionToPeptidesUsedForQuant[ControlCondition].Add(peptide);
                            ConditionsWithPeptideSampleQuantities[ControlCondition].Add(Math.Log(linearIntensity / peptide.IonizationEfficiency, 2));
                            NMeasurements++;
                        }
                    }
                }
            }

            FractionMbrControl = (double)numMbr / (numMbr + numMsms);
            numMsms = 0;
            numMbr = 0;

            if (!ConditionsWithPeptideSampleQuantities.ContainsKey(TreatmentCondition))
            {
                ConditionToPeptidesUsedForQuant.Add(TreatmentCondition, new HashSet<Peptide>());
                ConditionsWithPeptideSampleQuantities.Add(TreatmentCondition, new List<double>());

                for (int i = 0; i < Peptides.Count; i++)
                {
                    Peptide peptide = Peptides[i];

                    int samples = flashLfqResults.SpectraFiles.Where(p => p.Condition == TreatmentCondition).Max(p => p.BiologicalReplicate) + 1;

                    if (peptide.IonizationEfficiency > 0)
                    {
                        for (int sample = 0; sample < samples; sample++)
                        {
                            if (PeptideToSampleQuantity.TryGetValue((peptide, TreatmentCondition, sample), out double sampleIntensity) && sampleIntensity > 0)
                            {
                                ConditionToPeptidesUsedForQuant[TreatmentCondition].Add(peptide);
                                ConditionsWithPeptideSampleQuantities[TreatmentCondition].Add(Math.Log(sampleIntensity / peptide.IonizationEfficiency, 2));
                                NMeasurements++;

                                //DEBUG
                                var detectionType = flashLfqResults.PeptideModifiedSequences[peptide.Sequence].GetDetectionType(
                                    flashLfqResults.SpectraFiles.First(p => p.Condition == TreatmentCondition && p.BiologicalReplicate == sample));
                                if (detectionType == DetectionType.MSMS)
                                {
                                    numMsms++;
                                }
                                else if (detectionType == DetectionType.MBR)
                                {
                                    numMbr++;
                                }
                            }
                        }
                    }
                }
            }

            FractionMbrTreatment = (double)numMbr / (numMbr + numMsms);
        }
    }
}
