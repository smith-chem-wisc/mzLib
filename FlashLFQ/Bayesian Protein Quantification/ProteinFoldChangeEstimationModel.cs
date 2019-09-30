using BayesianEstimation;
using MathNet.Numerics.Distributions;
using System.Collections.Generic;

namespace FlashLFQ
{
    public class ProteinFoldChangeEstimationModel : Model
    {
        public ProteinFoldChangeEstimationModel(double priorMuMean, double priorMuSd, double muInitialGuess,
            double priorSdStart, double priorSdEnd, double sdInitialGuess, double priorNuExponent, double nuInitialGuess) : base()
        {
            modelParameters = new Parameter[3];

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
            // continuous uniform prior (see Kruschke 2013 J Exper Psych)
            modelParameters[1] = new Parameter(
                new ContinuousUniform(priorSdStart, priorSdEnd),
                new List<(double, double)> { (0, double.PositiveInfinity) },
                sdInitialGuess);

            // nu (sometimes referred to as "degrees of freedom")
            // exponential prior (see Kruschke 2013 J Exper Psych)
            modelParameters[2] = new Parameter(
                new Exponential(priorNuExponent),
                new List<(double, double)> { (0, double.PositiveInfinity) },
                nuInitialGuess);
        }

        protected override double ProbabilityOfModelGivenADatapoint(double[] paramProposals, double datapoint)
        {
            return StudentT.PDF(paramProposals[0], paramProposals[1], paramProposals[2], datapoint);
        }
    }
}
