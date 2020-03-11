using BayesianEstimation;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace FlashLFQ
{
    public class UnpairedProteinQuantResult : ProteinQuantificationEngineResult
    {
        public readonly Dictionary<string, List<Datum>> ConditionsWithPeptideSampleQuantities;
        public double MuControl;
        public double MuTreatment;
        public double SigmaControl;
        public double SigmaTreatment;
        public double NuControl;
        public double NuTreatment;
        public int NMeasurementsControl;
        public int NMeasurementsTreatment;
        public (double hdi_start, double hdi_end) hdi95Control;
        public (double hdi_start, double hdi_end) hdi95Treatment;

        public UnpairedProteinQuantResult(ProteinGroup protein, List<Peptide> peptides, string controlCondition, string treatmentCondition,
            bool useSharedPeptides, FlashLfqResults results, Dictionary<(Peptide, string, int), (double, DetectionType)> peptideToSampleQuantity,
            List<int> randomSeeds, int nBurnin)
            : base(protein, peptides, controlCondition, treatmentCondition)
        {
            ConditionsWithPeptideSampleQuantities = new Dictionary<string, List<Datum>>();
            GetPeptideSampleQuantities(useSharedPeptides, results, peptideToSampleQuantity, randomSeeds, nBurnin);
            IsStatisticallyValid = DetermineIfStatisticallyValid(peptideToSampleQuantity, results);
        }

        public void EstimateProteinFoldChange(int? randomSeed, int? randomSeed2, int nBurnin, int nSteps,
            double controlNullHypothesisInterval = double.NaN, double treatmentNullHypothesisInterval = double.NaN,
            StudentT controlSigmaPrior = null, StudentT treatmentSigmaPrior = null)
        {
            bool skepticalPrior = !double.IsNaN(controlNullHypothesisInterval) || !double.IsNaN(treatmentNullHypothesisInterval);

            var samples = ConditionsWithPeptideSampleQuantities[ControlCondition];

            var (controlMus, controlSds, controlNus) = ProteinQuantificationEngine.FitProteinQuantModel(samples, skepticalPrior, false, randomSeed, nBurnin, nSteps,
                controlNullHypothesisInterval, controlSigmaPrior);

            if (!skepticalPrior)
            {
                MuControl = controlMus.Median();
                SigmaControl = controlSds.Median();
                NuControl = controlNus.Median();
                hdi95Control = Util.GetHighestDensityInterval(controlMus, 0.95);
            }

            samples = ConditionsWithPeptideSampleQuantities[TreatmentCondition];

            var (treatmentMus, treatmentSds, treatmentNus) = ProteinQuantificationEngine.FitProteinQuantModel(samples, skepticalPrior, false, randomSeed2, nBurnin, nSteps,
                treatmentNullHypothesisInterval, treatmentSigmaPrior);

            if (!skepticalPrior)
            {
                MuTreatment = treatmentMus.Median();
                SigmaTreatment = treatmentSds.Median();
                NuTreatment = treatmentNus.Median();
                hdi95Treatment = Util.GetHighestDensityInterval(treatmentMus, 0.95);
            }

            double[] diffs = new double[Math.Min(controlMus.Length, treatmentMus.Length)];

            for (int i = 0; i < diffs.Length; i++)
            {
                diffs[i] = treatmentMus[i] - controlMus[i];
            }

            if (!skepticalPrior)
            {
                this.FoldChangePointEstimate = diffs.Median();

                double uncertaintyInControl = (hdi95Control.hdi_end - hdi95Control.hdi_start) / 2;
                double uncertaintyInTreatment = (hdi95Treatment.hdi_end - hdi95Treatment.hdi_start) / 2;

                UncertaintyInFoldChangeEstimate = Math.Sqrt(Math.Pow(uncertaintyInControl, 2) + Math.Pow(uncertaintyInTreatment, 2));
                
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
            var controlMmts = this.ConditionsWithPeptideSampleQuantities[ControlCondition];
            var treatmentMmts = this.ConditionsWithPeptideSampleQuantities[TreatmentCondition];

            var controlMmtsString = string.Join(", ", controlMmts.Select(p => Math.Round(p.DataValue, 2)));
            var treatmentMmtsString = string.Join(", ", treatmentMmts.Select(p => Math.Round(p.DataValue, 2)));

            NMeasurementsControl = controlMmts.Count;
            NMeasurementsTreatment = treatmentMmts.Count;

            if (controlMmtsString.Length > 32000)
            {
                controlMmtsString = "Too long for Excel";
            }

            if (treatmentMmtsString.Length > 32000)
            {
                treatmentMmtsString = "Too long for Excel";
            }

            return
                Protein.ProteinGroupName + "\t" +
                Protein.GeneName + "\t" +
                Protein.Organism + "\t" +
                ControlCondition + "\t" +
                TreatmentCondition + "\t" +
                NullHypothesisInterval.Value + "\t" +
                FoldChangePointEstimate + "\t" +
                UncertaintyInFoldChangeEstimate + "\t" +
                StandardDeviationPointEstimate + "\t" +
                ControlConditionIntensity + "\t" +
                TreatmentConditionIntensity + "\t" +
                Peptides.Count + "\t" +
                NMeasurementsControl + "\t" +
                NMeasurementsTreatment + "\t" +
                controlMmtsString + "\t" +
                treatmentMmtsString + "\t" +
                BayesFactor + "\t" +
                ProbabilityOfFalsePositive + "\t" +
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
                "Null Hypothesis Width" + "\t" +
                "Protein Log2 Fold-Change" + "\t" +
                "Uncertainty in Protein Log2 Fold-Change" + "\t" +
                "Standard Deviation of Peptide Log2 Fold-Changes" + "\t" +
                "Protein Intensity in Control Condition" + "\t" +
                "Protein Intensity in Treatment Condition" + "\t" +
                "Number of Peptides" + "\t" +
                "Number of Control Condition Measurements" + "\t" +
                "Number of Treatment Condition Measurements" + "\t" +
                "Control Measurements" + "\t" +
                "Treatment Measurements" + "\t" +
                "Bayes Factor" + "\t" +
                "Probability of False-Positive" + "\t" +
                "False Discovery Rate" + "\t\t";
        }

        private void GetPeptideSampleQuantities(bool useSharedPeptides, FlashLfqResults flashLfqResults,
            Dictionary<(Peptide, string, int), (double, DetectionType)> PeptideToSampleQuantity, List<int> randomSeeds, int nBurnin)
        {
            List<string> conditions = new List<string> { ControlCondition, TreatmentCondition };
            foreach (var condition in conditions)
            {
                ConditionsWithPeptideSampleQuantities.Add(condition, new List<Datum>());
            }

            int peptidesIdentifiedByMsmsInBothConditions = 0;

            for (int i = 0; i < Peptides.Count; i++)
            {
                Peptide peptide = Peptides[i];
                bool identifiedInAllConditions = true;

                foreach (var condition in conditions)
                {
                    bool identifiedInThisCondition = false;
                    int numSamplesInGroup = flashLfqResults.SpectraFiles.Where(p => p.Condition == condition).Max(p => p.BiologicalReplicate) + 1;

                    for (int sample = 0; sample < numSamplesInGroup; sample++)
                    {
                        double intensity = PeptideToSampleQuantity[(peptide, condition, sample)].Item1;
                        DetectionType detectionType = PeptideToSampleQuantity[(peptide, condition, sample)].Item2;

                        if (intensity > 0 && detectionType == DetectionType.MSMS)
                        {
                            identifiedInThisCondition = true;
                            break;
                        }
                    }

                    if (!identifiedInThisCondition)
                    {
                        identifiedInAllConditions = false;
                    }
                }

                if (identifiedInAllConditions)
                {
                    peptidesIdentifiedByMsmsInBothConditions++;
                }
            }

            for (int i = 0; i < Peptides.Count; i++)
            {
                Peptide peptide = Peptides[i];
                int randomSeed = randomSeeds[i];

                if (!peptide.UseForProteinQuant || (!useSharedPeptides && peptide.ProteinGroups.Count > 1))
                {
                    continue;
                }
                
                foreach (string condition in conditions)
                {
                    if (!ConditionsWithPeptideSampleQuantities.ContainsKey(condition))
                    {
                        ConditionsWithPeptideSampleQuantities.Add(condition, new List<Datum>());
                    }

                    int numSamplesInGroup = flashLfqResults.SpectraFiles.Where(p => p.Condition == condition).Max(p => p.BiologicalReplicate) + 1;
                    List<(double, DetectionType)> peptideLogIntensities = new List<(double, DetectionType)>();

                    for (int sample = 0; sample < numSamplesInGroup; sample++)
                    {
                        double intensity = PeptideToSampleQuantity[(peptide, condition, sample)].Item1;
                        DetectionType detectionType = PeptideToSampleQuantity[(peptide, condition, sample)].Item2;

                        if (intensity > 0)
                        {
                            peptideLogIntensities.Add((Math.Log(intensity, 2), detectionType));
                        }
                    }
                    
                    // estimate ionization efficiency
                    if (condition == ControlCondition && peptide.IonizationEfficiency == 0)
                    {
                        if (peptideLogIntensities.Count > 1)
                        {
                            var intensities = peptideLogIntensities.Select(p => p.Item1).ToList();

                            double medianIntensity = intensities.Median();

                            peptide.IonizationEfficiency = Math.Pow(2, medianIntensity);
                        }
                    }

                    if (peptide.IonizationEfficiency != 0)
                    {
                        List<Datum> ionizationEfficiencyNormalizedData = new List<Datum>();

                        foreach (var logIntensity in peptideLogIntensities)
                        {
                            double linearIntensity = Math.Pow(2, logIntensity.Item1);

                            ionizationEfficiencyNormalizedData.Add(
                                new Datum(Math.Log(linearIntensity / peptide.IonizationEfficiency, 2), weight: 1));
                        }

                        ConditionsWithPeptideSampleQuantities[condition].AddRange(ionizationEfficiencyNormalizedData);
                    }
                }
            }
        }

        private bool DetermineIfStatisticallyValid(Dictionary<(Peptide, string, int), (double, DetectionType)> peptideToSampleQuantity, FlashLfqResults results)
        {
            NMeasurementsControl = ConditionsWithPeptideSampleQuantities[ControlCondition].Count(p => p.Weight > 0);
            NMeasurementsTreatment = ConditionsWithPeptideSampleQuantities[TreatmentCondition].Count(p => p.Weight > 0);

            NMeasurements = NMeasurementsControl + NMeasurementsTreatment;

            return NMeasurementsControl > 1 && NMeasurementsTreatment > 1;
        }
    }
}
