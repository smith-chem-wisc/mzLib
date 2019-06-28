using System;

namespace BayesianEstimation
{
    public abstract class Model
    {
        public Parameter[] modelParameters { get; protected set; }

        public Model()
        {

        }

        /// <summary>
        /// Calculates the log probability of the prior plus the log probability of the model, given the data.
        /// </summary>
        public double LogPosteriorProbability(double[] paramProposals, double[] data)
        {
            return LogPriorProbability(paramProposals) + LogProbabilityOfModelGivenTheData(paramProposals, data);
        }

        /// <summary>
        /// Calculates the probability of the model, given a single data point.
        /// 
        /// This method needs to be overridden by an inheriting class.
        /// </summary>
        protected abstract double ProbabilityOfModelGivenADatapoint(double[] paramProposals, double datapoint);

        /// <summary>
        /// Calculates the log probability of the model with the proposed parameters, given the data.
        /// </summary>
        private double LogProbabilityOfModelGivenTheData(double[] paramProposals, double[] data)
        {
            double log_p = 0;

            for (var i = 0; i < data.Length; i++)
            {
                double datapoint = data[i];
                double probabilityOfModelGivenTheDatapoint = ProbabilityOfModelGivenADatapoint(paramProposals, datapoint);
                log_p += Math.Log(probabilityOfModelGivenTheDatapoint);
            }

            return log_p;
        }

        /// <summary>
        /// Calculates the log prior probability of a set of parameters.
        /// </summary>
        private double LogPriorProbability(double[] paramProposals)
        {
            double log_p = 0;

            for (int p = 0; p < paramProposals.Length; p++)
            {
                double priorProbabilityForParameterProposal = modelParameters[p].priorProbabilityDistribution.Density(paramProposals[p]);
                log_p += Math.Log(priorProbabilityForParameterProposal);
            }

            return log_p;
        }
    }
}
