using MathNet.Numerics.Distributions;
using System.Collections.Generic;

namespace BayesianEstimation
{
    public class StudentTDistributionModel : Model
    {
        public StudentTDistributionModel(double priorMuMean, double priorMuSd, double muInitialGuess,
            double priorSdStart, double priorSdEnd, double sdInitialGuess, double priorNuExponent, double nuInitialGuess) : base()
        {
            modelParameters = new Parameter[3];

            // mu (mean)
            modelParameters[0] = new Parameter(
                new Normal(priorMuMean, priorMuSd),
                new List<(double, double)> { (double.NegativeInfinity, double.PositiveInfinity) },
                muInitialGuess);

            // sigma (standard deviation)
            modelParameters[1] = new Parameter(
                new ContinuousUniform(priorSdStart, priorSdEnd),
                new List<(double, double)> { (0, double.PositiveInfinity) },
                sdInitialGuess);

            // nu (sometimes referred to as "degrees of freedom")
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
