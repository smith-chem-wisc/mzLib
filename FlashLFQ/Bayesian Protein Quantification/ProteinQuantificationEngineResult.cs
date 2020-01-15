using BayesianEstimation;
using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace FlashLFQ
{
    public abstract class ProteinQuantificationEngineResult
    {
        // protein and peptides
        public readonly ProteinGroup Protein;
        public readonly List<Peptide> Peptides;

        // conditions
        public readonly string ControlCondition;
        public readonly string TreatmentCondition;

        // mu, sigma, nu estimates
        public double FoldChangePointEstimate { get; protected set; }
        public double StandardDeviationPointEstimate { get; protected set; }
        public double NuPointEstimate { get; protected set; }

        // intensities
        public double ControlConditionIntensity { get; protected set; }
        public double TreatmentConditionIntensity { get; protected set; }

        // additional statistics (PEP, FDR, null interval width)
        public double PosteriorErrorProbability { get; private set; }
        public double FalseDiscoveryRate { get; set; }
        public double? NullHypothesisInterval { get; set; }
        public int NMeasurements { get; protected set; }

        protected ProteinQuantificationEngineResult(ProteinGroup protein, List<Peptide> peptides, string controlCondition, string treatmentCondition)
        {
            this.Protein = protein;
            this.Peptides = peptides;
            this.ControlCondition = controlCondition;
            this.TreatmentCondition = treatmentCondition;
        }

        /// <summary>
        /// The posterior error probability (PEP) is the probability that the null hypothesis is true.
        /// 
        /// This is estimated by dividing the number of Markov Chain Monte Carlo iterations that show a fold-change
        /// greater than the cutoff by the number of total MCMC iterations.
        /// 
        /// This means that if the null and alternative hypotheses are equally likely, the PEP equals 0.5.
        /// If the null hypothesis is more likely than the alternative, the PEP is greater than 0.5.
        /// If the alternative hypothesis is more likely than the null hypothesis, then the PEP is less than 0.5.
        /// As the alternative hypothesis becomes increasingly likely, the PEP decreases towards zero.
        /// 
        /// If there are no iterations of the MCMC algorithm that show a protein fold-change less than the cutoff,
        /// then the PEP equals zero. 
        /// 
        /// If there are no iterations of the MCMC algorithm that show a protein fold-change greater than the absolute
        /// value of the cutoff, then the PEP equals 1 and the null hypothesis can be accepted.
        /// 
        /// Additional MCMC iterations can be performed to determine a more precise estimate of the PEP.
        /// </summary>
        protected void CalculatePosteriorErrorProbability(double[] musWithSkepticalPrior)
        {
            if (NullHypothesisInterval == null || musWithSkepticalPrior == null)
            {
                PosteriorErrorProbability = double.NaN;
                return;
            }

            int numIncreasing = musWithSkepticalPrior.Count(p => p > NullHypothesisInterval);
            int numDecreasing = musWithSkepticalPrior.Count(p => p < -NullHypothesisInterval);

            // if something goes wrong and none of the following "if" statements are triggered, then PEP will evaluate to 1.0
            double nullHypothesisCount = musWithSkepticalPrior.Length;
            double alternativeHypothesisCount = 0;

            if (numIncreasing >= numDecreasing && numIncreasing > 0)
            {
                nullHypothesisCount = musWithSkepticalPrior.Count(p => p < NullHypothesisInterval);
                alternativeHypothesisCount = numIncreasing;
            }
            else if (numIncreasing < numDecreasing)
            {
                nullHypothesisCount = musWithSkepticalPrior.Count(p => p > -NullHypothesisInterval);
                alternativeHypothesisCount = numDecreasing;
            }
            // this doesn't need to be here, but it's here for logical completeness
            else if (musWithSkepticalPrior.All(p => Math.Abs(p) <= NullHypothesisInterval))
            {
                nullHypothesisCount = musWithSkepticalPrior.Length;
                alternativeHypothesisCount = 0;
            }

            PosteriorErrorProbability = nullHypothesisCount / musWithSkepticalPrior.Length;
        }
        
        protected static string TabSeparatedHeader()
        {
            return "";
        }
    }
}
