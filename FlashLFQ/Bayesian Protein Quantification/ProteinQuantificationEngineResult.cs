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

        /// <summary>
        /// The posterior error probability (PEP) is the probability that the alternative hypothesis is not true.
        /// 
        /// This is estimated by dividing the number of Markov Chain Monte Carlo iterations that show a fold-change
        /// less than the cutoff by the number of MCMC iterations that show a fold-change greater than the cutoff.
        /// 
        /// This means that if the null and alternative hypotheses are equally likely, the PEP equals 1.
        /// If the null hypothesis is more likely than the alternative, the PEP is also equal to 1.
        /// If the alternative hypothesis is more likely than the null hypothesis, then the PEP is less than 1.
        /// As the alternative hypothesis becomes increasingly likely, the PEP decreases towards zero.
        /// 
        /// If there are no iterations of the MCMC algorithm that show a protein fold-change less than the cutoff,
        /// then the PEP equals zero. 
        /// 
        /// If there are no iterations of the MCMC algorithm that show a protein fold-change greater than the absolute
        /// value of the cutoff, then the PEP equals 1 and the null hypothesis can be accepted, rather than just
        /// rejecting the alternative hypothesis.
        /// 
        /// Additional MCMC iterations can be performed to determine a more precise estimate of the PEP.
        /// 
        /// Numerically:
        /// F = fold change
        /// C = cutoff
        /// 
        /// if F > C, then P(H0|D) = P(F < C) / P (F > C)
        /// if F < -C, then P(H0|D) = P(F > -C) / P (F < -C)
        /// </summary>
        public void CalculatePosteriorErrorProbability(double[] skepticalMus)
        {
            double pep;

            if (cutoff == null || skepticalMus == null)
            {
                pep = double.NaN;
            }

            int numIncreasing = skepticalMus.Count(p => p > cutoff);
            int numDecreasing = skepticalMus.Count(p => p < -cutoff);

            // the hypotheses are equally likely beofre any data
            // if something goes wrong and none of the following "if" statements are triggered, then PEP will evaluate to 1.0
            double nullHypothesisCount = 1.0;
            double alternativeHypothesisCount = 1.0;

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
            else if (numIncreasing == numDecreasing) // this technically doesn't need to be here, but it's here for logical completeness
            {
                nullHypothesisCount = 1.0;
                alternativeHypothesisCount = 1.0;
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
