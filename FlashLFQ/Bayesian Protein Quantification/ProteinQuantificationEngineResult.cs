using BayesianEstimation;
using System;
using System.Collections.Generic;
using System.Linq;

namespace FlashLFQ
{
    public class ProteinQuantificationEngineResult
    {
        public readonly ProteinGroup protein;
        public readonly string condition1;
        public readonly string condition2;
        public readonly double FoldChangePointEstimate;
        public readonly double ConditionIntensityPointEstimate;
        public readonly double StandardDeviationPointEstimate;
        public readonly double NuPointEstimate;
        public readonly (double, double) hdi_95;
        public readonly List<double> foldChangeMeasurements;
        public double PosteriorErrorProbability { get; private set; }
        public double FalseDiscoveryRate { get; set; }
        public double DegreeOfEvidence { get; private set; }
        public double? cutoff { get; set; }

        public ProteinQuantificationEngineResult(ProteinGroup protein, string condition1, string condition2, double[] mus, double[] sds, double[] nus,
            List<double> fcs, double referenceIntensity)
        {
            this.protein = protein;
            this.condition1 = condition1;
            this.condition2 = condition2;
            this.foldChangeMeasurements = fcs;
            this.FoldChangePointEstimate = mus.Average();
            this.StandardDeviationPointEstimate = sds.Average();
            this.NuPointEstimate = nus.Average();
            this.hdi_95 = Util.GetHighestDensityInterval(mus, 0.95);
            this.ConditionIntensityPointEstimate = Math.Pow(2, FoldChangePointEstimate) * referenceIntensity;
        }

        public void CalculatePosteriorErrorProbability(double[] skepticalMus)
        {
            double pep = 1.0;

            if (cutoff == null || skepticalMus == null)
            {
                pep = double.NaN;
            }

            int numIncreasing = skepticalMus.Count(p => p > cutoff);
            int numDecreasing = skepticalMus.Count(p => p < -cutoff);

            double nullHypothesisCount = 0;
            double alternativeHypothesisCount = 0;

            if (numIncreasing > numDecreasing)
            {
                nullHypothesisCount = skepticalMus.Count(p => p < cutoff);
                alternativeHypothesisCount = numIncreasing;
            }
            else if (numDecreasing > numIncreasing)
            {
                nullHypothesisCount = skepticalMus.Count(p => p > -cutoff);
                alternativeHypothesisCount = numDecreasing;
            }
            else if (skepticalMus.All(p => Math.Abs(p) <= cutoff))
            {
                nullHypothesisCount = skepticalMus.Length;
                alternativeHypothesisCount = 0;
            }
            else if (numIncreasing == numDecreasing)
            {
                nullHypothesisCount = 1.0;
                alternativeHypothesisCount = 1.0;
            }
            else
            {

            }

            pep = nullHypothesisCount / alternativeHypothesisCount;
            DegreeOfEvidence = alternativeHypothesisCount / nullHypothesisCount;

            pep = Math.Min(1.0, pep);

            this.PosteriorErrorProbability = pep;
        }

        public override string ToString()
        {
            return condition1 + "\t" + condition2 + "\t" + cutoff.Value + "\t" + FoldChangePointEstimate + "\t" + ConditionIntensityPointEstimate 
                + "\t" + foldChangeMeasurements.Count() + "\t" + PosteriorErrorProbability + "\t" + FalseDiscoveryRate;
        }
    }
}
