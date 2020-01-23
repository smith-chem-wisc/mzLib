using MathNet.Numerics.Distributions;
using System.Collections.Generic;

namespace BayesianEstimation
{
    public class Parameter
    {
        public readonly IContinuousDistribution priorProbabilityDistribution;
        public readonly List<(double, double)> acceptableParameterRanges;
        public readonly double initialGuess;

        public Parameter(IContinuousDistribution priorProbabilityDistribution, List<(double, double)> acceptableParameterRanges, double initialGuess)
        {
            this.priorProbabilityDistribution = priorProbabilityDistribution;
            this.acceptableParameterRanges = acceptableParameterRanges;
            this.initialGuess = initialGuess;
        }
    }
}
