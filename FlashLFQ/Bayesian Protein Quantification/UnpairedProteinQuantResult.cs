using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace FlashLFQ
{
    public class UnpairedProteinQuantResult : ProteinQuantificationEngineResult
    {
        public readonly Dictionary<string, List<double>> ConditionsWithPeptideSampleQuantities;
        public double MuControl;
        public double MuTreatment;
        public double SigmaControl;
        public double SigmaTreatment;
        public double NuControl;
        public double NuTreatment;
        public int NMeasurementsControl;
        public int NMeasurementsTreatment;

        //DEBUG
        public double FalsePositive;
        public double FDP;

        public UnpairedProteinQuantResult(ProteinGroup protein, List<Peptide> peptides, string controlCondition, string treatmentCondition,
            bool useSharedPeptides, FlashLfqResults results, Dictionary<(Peptide, string, int), double> peptideToSampleQuantity, List<int> randomSeeds, int nBurnin)
            : base(protein, peptides, controlCondition, treatmentCondition)
        {
            ConditionsWithPeptideSampleQuantities = new Dictionary<string, List<double>>();
            IsStatisticallyValid = DetermineIfStatisticallyValid(peptideToSampleQuantity, results);

            GetPeptideSampleQuantities(useSharedPeptides, results, peptideToSampleQuantity, randomSeeds, nBurnin);
        }

        public void EstimateProteinFoldChange(int? randomSeed, int? randomSeed2, int nBurnin, int nSteps,
            double controlNullHypothesisInterval = double.NaN, double treatmentNullHypothesisInterval = double.NaN,
            IContinuousDistribution controlSigmaPrior = null, IContinuousDistribution treatmentSigmaPrior = null,
            IContinuousDistribution controlNuPrior = null, IContinuousDistribution treatmentNuPrior = null)
        {
            bool skepticalPrior = !double.IsNaN(controlNullHypothesisInterval) || !double.IsNaN(treatmentNullHypothesisInterval);

            var samples = ConditionsWithPeptideSampleQuantities[ControlCondition];

            var (controlMus, controlSds, controlNus) = ProteinQuantificationEngine.FitProteinQuantModel(samples, skepticalPrior, false, randomSeed, nBurnin, nSteps,
                controlNullHypothesisInterval, controlSigmaPrior, controlNuPrior);

            if (!skepticalPrior)
            {
                MuControl = controlMus.Median();
                SigmaControl = controlSds.Median();
                NuControl = controlNus.Median();
            }

            samples = ConditionsWithPeptideSampleQuantities[TreatmentCondition];

            var (treatmentMus, treatmentSds, treatmentNus) = ProteinQuantificationEngine.FitProteinQuantModel(samples, skepticalPrior, false, randomSeed2, nBurnin, nSteps,
                treatmentNullHypothesisInterval, treatmentSigmaPrior, treatmentNuPrior);

            if (!skepticalPrior)
            {
                MuTreatment = treatmentMus.Median();
                SigmaTreatment = treatmentSds.Median();
                NuTreatment = treatmentNus.Median();
            }

            double[] diffs = new double[Math.Min(controlMus.Length, treatmentMus.Length)];

            if (!IsStatisticallyValid)
            {
                for (int i = 0; i < diffs.Length; i++)
                {
                    diffs[i] = 0;
                }
            }
            else
            {
                for (int i = 0; i < diffs.Length; i++)
                {
                    diffs[i] = treatmentMus[i] - controlMus[i];
                }
            }

            if (!skepticalPrior)
            {
                this.FoldChangePointEstimate = diffs.Median();

                // propagate SD
                double[] propagatedSds = new double[Math.Min(controlSds.Length, treatmentSds.Length)];
                if (!IsStatisticallyValid)
                {
                    for (int i = 0; i < propagatedSds.Length; i++)
                    {
                        propagatedSds[i] = double.NaN;
                    }
                }
                else
                {
                    for (int i = 0; i < propagatedSds.Length; i++)
                    {
                        propagatedSds[i] = Math.Sqrt(Math.Pow(controlSds[i], 2) + Math.Pow(treatmentSds[i], 2));
                    }
                }

                double propagatedSd = propagatedSds.Median();
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
                FoldChangePointEstimate + "\t" +
                StandardDeviationPointEstimate + "\t" +
                NuControl + "\t" +
                NuTreatment + "\t" +
                SigmaControl + "\t" +
                SigmaTreatment + "\t" +
                ControlConditionIntensity + "\t" +
                TreatmentConditionIntensity + "\t" +
                Peptides.Count + "\t" +
                NMeasurementsControl + "\t" +
                NMeasurementsTreatment + "\t" +
                measurementsString + "\t" +
                FalsePositive + "\t" +
                FDP + "\t" +
                PosteriorErrorProbability + "\t" +
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
                "Protein Log2 Fold-Change" + "\t" +
                "Standard Deviation of Peptide Log2 Fold-Changes" + "\t" +
                "Nu Control" + "\t" +
                "Nu Treatment" + "\t" +
                "Sigma Control" + "\t" +
                "Sigma Treatment" + "\t" +
                "Protein Intensity in Control Condition" + "\t" +
                "Protein Intensity in Treatment Condition" + "\t" +
                "Number of Peptides" + "\t" +
                "Number of Control Condition Measurements" + "\t" +
                "Number of Treatment Condition Measurements" + "\t" +
                "List of Fold-Change Measurements Grouped by Peptide" + "\t" +
                "FalsePositive" + "\t" +
                "FDP" + "\t" +
                "Posterior Error Probability" + "\t" +
                "False Discovery Rate" + "\t\t";
        }

        private void GetPeptideSampleQuantities(bool useSharedPeptides, FlashLfqResults flashLfqResults,
            Dictionary<(Peptide, string, int), double> PeptideToSampleQuantity, List<int> randomSeeds, int nBurnin)
        {
            if (!ConditionsWithPeptideSampleQuantities.ContainsKey(ControlCondition))
            {
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
                        }
                    }

                    // estimate ionization efficiency
                    if (peptide.IonizationEfficiency == 0)
                    {
                        if (peptideSampleLogIntensities.Count > 3)
                        {
                            var (mus, sds, nus) = ProteinQuantificationEngine.FitProteinQuantModel(peptideSampleLogIntensities, false, false, randomSeed,
                                nBurnin, 1000, 0, null, null);
                            peptide.IonizationEfficiency = Math.Pow(2, mus.Median());
                        }
                        else if (peptideSampleLogIntensities.Count > 0)
                        {
                            peptide.IonizationEfficiency = Math.Pow(2, peptideSampleLogIntensities.Median());
                        }
                    }

                    if (peptide.IonizationEfficiency != 0)
                    {
                        foreach (double logIntensity in peptideSampleLogIntensities)
                        {
                            double linearIntensity = Math.Pow(2, logIntensity);

                            ConditionsWithPeptideSampleQuantities[ControlCondition].Add(Math.Log(linearIntensity / peptide.IonizationEfficiency, 2));
                        }
                    }
                }
            }

            if (!ConditionsWithPeptideSampleQuantities.ContainsKey(TreatmentCondition))
            {
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
                                ConditionsWithPeptideSampleQuantities[TreatmentCondition].Add(Math.Log(sampleIntensity / peptide.IonizationEfficiency, 2));
                            }
                        }
                    }
                }
            }
        }

        private bool DetermineIfStatisticallyValid(Dictionary<(Peptide, string, int), double> peptideToSampleQuantity, FlashLfqResults results)
        {
            NMeasurementsControl = 0;
            NMeasurementsTreatment = 0;

            int numControlSamples = results.SpectraFiles.Where(p => p.Condition == ControlCondition).Max(p => p.BiologicalReplicate) + 1;
            int numTreatmentSamples = results.SpectraFiles.Where(p => p.Condition == TreatmentCondition).Max(p => p.BiologicalReplicate) + 1;

            foreach (var peptide in Peptides)
            {
                for (int s = 0; s < numControlSamples; s++)
                {
                    if (peptideToSampleQuantity[(peptide, ControlCondition, s)] > 0)
                    {
                        NMeasurementsControl++;
                    }
                }
            }

            foreach (var peptide in Peptides)
            {
                for (int s = 0; s < numTreatmentSamples; s++)
                {
                    if (peptideToSampleQuantity[(peptide, TreatmentCondition, s)] > 0)
                    {
                        NMeasurementsTreatment++;
                    }
                }
            }

            NMeasurements = NMeasurementsControl + NMeasurementsTreatment;

            return NMeasurementsControl > 1 && NMeasurementsTreatment > 1;
        }
    }
}
