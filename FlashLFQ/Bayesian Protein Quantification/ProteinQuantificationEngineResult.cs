using BayesianEstimation;
using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace FlashLFQ
{
    public class ProteinQuantificationEngineResult
    {
        public readonly ProteinGroup protein;
        public readonly string BaseCondition;
        public readonly string TreatmentCondition;
        public readonly double FoldChangePointEstimate;
        public readonly double ConditionIntensityPointEstimate;
        public readonly double StandardDeviationPointEstimate;
        public readonly double NuPointEstimate;

        // the HDI_95 is the "95% highest density interval". this is the interval where 95% of the 
        // probability density is contained (in this case, for the estimate of the mean). 
        // it can be thought of as analogous to a 95% confidence interval.
        public readonly (double, double) MeanHDI_95;

        public readonly List<(Peptide peptide, List<double> foldChanges)> PeptideFoldChangeMeasurements;
        public double PosteriorErrorProbability { get; private set; }
        public double FalseDiscoveryRate { get; set; }
        public double? NullHypothesisCutoff { get; set; }

        public ProteinQuantificationEngineResult(ProteinGroup protein, string condition1, string condition2, double[] mus, double[] sds, double[] nus,
            List<(Peptide, List<double>)> fcs, double referenceIntensity)
        {
            this.protein = protein;
            this.BaseCondition = condition1;
            this.TreatmentCondition = condition2;
            this.PeptideFoldChangeMeasurements = fcs;
            this.FoldChangePointEstimate = mus.Median();
            this.StandardDeviationPointEstimate = sds.Median();
            this.NuPointEstimate = nus.Median();
            this.MeanHDI_95 = Util.GetHighestDensityInterval(mus, 0.95);
            this.ConditionIntensityPointEstimate = Math.Pow(2, FoldChangePointEstimate) * referenceIntensity;
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
        public void CalculatePosteriorErrorProbability(double[] musWithSkepticalPrior)
        {
            if (NullHypothesisCutoff == null || musWithSkepticalPrior == null)
            {
                PosteriorErrorProbability = double.NaN;
                return;
            }

            int numIncreasing = musWithSkepticalPrior.Count(p => p > NullHypothesisCutoff);
            int numDecreasing = musWithSkepticalPrior.Count(p => p < -NullHypothesisCutoff);

            // if something goes wrong and none of the following "if" statements are triggered, then PEP will evaluate to 1.0
            double nullHypothesisCount = musWithSkepticalPrior.Length;
            double alternativeHypothesisCount = 0;

            if (numIncreasing >= numDecreasing && numIncreasing > 0)
            {
                nullHypothesisCount = musWithSkepticalPrior.Count(p => p < NullHypothesisCutoff);
                alternativeHypothesisCount = numIncreasing;
            }
            else if (numIncreasing < numDecreasing)
            {
                nullHypothesisCount = musWithSkepticalPrior.Count(p => p > -NullHypothesisCutoff);
                alternativeHypothesisCount = numDecreasing;
            }
            // this doesn't need to be here, but it's here for logical completeness
            else if (musWithSkepticalPrior.All(p => Math.Abs(p) <= NullHypothesisCutoff))
            {
                nullHypothesisCount = musWithSkepticalPrior.Length;
                alternativeHypothesisCount = 0;
            }

            PosteriorErrorProbability = nullHypothesisCount / musWithSkepticalPrior.Length;
        }

        public override string ToString()
        {
            int nPeptides = PeptideFoldChangeMeasurements.Count;
            int nMeasurements = PeptideFoldChangeMeasurements.SelectMany(p => p.foldChanges).Count();
            var measurementsString = PeptideFoldChangeMeasurements.Select(p => p.peptide.Sequence + ":" + string.Join(";", p.foldChanges.Select(v => v.ToString("F4"))));

            return
                protein.ProteinGroupName + "\t" +
                protein.GeneName + "\t" +
                protein.Organism + "\t" +
                BaseCondition + "\t" +
                TreatmentCondition + "\t" +
                NullHypothesisCutoff.Value + "\t" +
                FoldChangePointEstimate + "\t" +
                ConditionIntensityPointEstimate + "\t" +
                nPeptides + "\t" +
                nMeasurements + "\t" +
                string.Join(",", measurementsString) + "\t" +
                PosteriorErrorProbability + "\t" +
                FalseDiscoveryRate + "\t\t";
        }

        public static string TabSeparatedHeader()
        {
            return
                "Protein Group" + "\t" +
                "Gene" + "\t" +
                "Organism" + "\t" +
                "Base Condition" + "\t" +
                "Treatment Condition" + "\t" +
                "Log2 Fold-Change Cutoff" + "\t" +
                "Protein Log2 Fold-Change" + "\t" +
                "Protein Intensity for Treatment Condition" + "\t" +
                "Number of Peptides" + "\t" +
                "Number of Fold-Change Measurements" + "\t" +
                "List of Fold-Change Measurements Grouped by Peptide" + "\t" +
                "Posterior Error Probability" + "\t" +
                "False Discovery Rate" + "\t\t";
        }
    }
}
