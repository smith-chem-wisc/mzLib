using BayesianEstimation;
using MathNet.Numerics.Distributions;
using System;
using System.Collections.Generic;
using MzLibUtil;

namespace FlashLFQ
{
    public class ProteinFoldChangeEstimationModel : Model
    {
        public ProteinFoldChangeEstimationModel(double priorMuMean, double priorMuSd, double muInitialGuess,
            IContinuousDistribution sdPrior, double sdInitialGuess, IContinuousDistribution nuPrior, double nuInitialGuess,
            double minimumNu = 2, double minimumSd = 0) : base()
        {
            modelParameters = new Parameter[3];

            if (double.IsNaN(minimumSd) || minimumSd <= 0)
            {
                minimumSd = 0.001;
            }
            if (double.IsNaN(minimumNu) || minimumNu < 1)
            {
                minimumNu = 1;
            }

            // mu (mean)
            // t-distributed prior
            // the protein's fold change is likely around some value (the arithmetic mean, the median, zero, etc. depending on assumptions),
            // but with a possibility that it is a fairly different value, especially given enough data. the Student's t-distribution
            // is a better choice than the normal distribution in this case because it has long tails.
            modelParameters[0] = new Parameter(
                new StudentT(priorMuMean, priorMuSd, 1.0),
                new List<(double, double)> { (double.NegativeInfinity, double.PositiveInfinity) },
                muInitialGuess);

            // sigma (standard deviation)
            if (sdPrior != null)
            {
                modelParameters[1] = new Parameter(
                    sdPrior,
                    new List<(double, double)> { (minimumSd, double.PositiveInfinity) },
                    sdInitialGuess);
            }
            else
            {
                // continuous uniform prior (see Kruschke 2013 J Exper Psych)
                modelParameters[1] = new Parameter(
                    new ContinuousUniform(sdInitialGuess / 1000, sdInitialGuess * 1000),
                    new List<(double, double)> { (0, double.PositiveInfinity) },
                    sdInitialGuess);
            }

            // nu (sometimes referred to as "degrees of freedom")
            if (nuPrior != null)
            {
                modelParameters[2] = new Parameter(
                nuPrior,
                new List<(double, double)> { (minimumNu, double.PositiveInfinity) },
                Math.Max(minimumNu + 1, nuInitialGuess));
            }
            else
            {
                var dist = new Normal(30, 15);
                // truncated normal distribution
                // assumes most data are normally distributed, with relatively conservative outlier rejection
                modelParameters[2] = new Parameter(
                dist,
                new List<(double, double)> { (minimumNu, double.PositiveInfinity) },
                dist.Mode);
            }
        }

        protected override double ProbabilityOfModelGivenADatapoint(double[] paramProposals, Datum datapoint)
        {
            return Math.Pow(StudentT.PDF(paramProposals[0], paramProposals[1], paramProposals[2], datapoint.X), datapoint.Weight);
        }
    }
}
