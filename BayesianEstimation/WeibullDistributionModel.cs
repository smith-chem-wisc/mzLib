using MathNet.Numerics.Distributions;
using System;
using System.Collections.Generic;
using System.Text;

namespace BayesianEstimation
{
    public class WeibullDistributionModel : Model
    {
        public WeibullDistributionModel(double priorShapeStart, double priorShapeEnd, double shapeInitialGuess,
            double priorScaleStart, double priorScaleEnd, double scaleInitialGuess) : base()
        {
            modelParameters = new Parameter[2];

            // shape
            modelParameters[0] = new Parameter(
                new ContinuousUniform(priorShapeStart, priorShapeEnd),
                new List<(double, double)> { (0, double.PositiveInfinity) },
                shapeInitialGuess);

            // scale
            modelParameters[1] = new Parameter(
                new ContinuousUniform(priorScaleStart, priorScaleEnd),
                new List<(double, double)> { (0, double.PositiveInfinity) },
                scaleInitialGuess);
        }

        protected override double ProbabilityOfModelGivenADatapoint(double[] paramProposals, double datapoint)
        {
            return Weibull.PDF(paramProposals[0], paramProposals[1], datapoint);
        }
    }
}
