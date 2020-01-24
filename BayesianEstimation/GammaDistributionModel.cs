using MathNet.Numerics.Distributions;
using System;
using System.Collections.Generic;
using System.Text;

namespace BayesianEstimation
{
    public class GammaDistributionModel : Model
    {
        public GammaDistributionModel(double priorShapeStart, double priorShapeEnd, double shapeInitialGuess,
            double priorRateStart, double priorRateEnd, double rateInitialGuess) : base()
        {
            modelParameters = new Parameter[2];

            // shape
            modelParameters[0] = new Parameter(
                new ContinuousUniform(priorShapeStart, priorShapeEnd),
                new List<(double, double)> { (0, double.PositiveInfinity) },
                shapeInitialGuess);

            // rate
            modelParameters[1] = new Parameter(
                new ContinuousUniform(priorRateStart, priorRateEnd),
                new List<(double, double)> { (0, double.PositiveInfinity) },
                rateInitialGuess);
        }

        protected override double ProbabilityOfModelGivenADatapoint(double[] paramProposals, double datapoint)
        {
            return Gamma.PDF(paramProposals[0], paramProposals[1], datapoint);
        }
    }
}
