using MzLibUtil;
using System;

namespace BayesianEstimation
{
    /// <summary>
    /// This class is a generic framework for a model that describes a set of data points. It must be
    /// inherited by a more specific model with defined parameters and a probability function to be used.
    /// 
    /// For example, the Student's t distribution can be used to model a set of data points. See 
    /// StudentTDistributionModel.cs for details.
    /// </summary>
    public abstract class Model
    {
        /// <summary>
        /// These model parameters are null and must be defined by an inheriting class.
        /// </summary>
        public Parameter[] modelParameters { get; protected set; }
        
        protected Model()
        {

        }
        
        /// <summary>
        /// Calculates the probability of the model, given a single data point.
        /// 
        /// This method needs to be overridden by an inheriting class.
        /// </summary>
        protected abstract double ProbabilityOfModelGivenADatapoint(double[] paramProposals, Datum datapoint);

        /// <summary>
        /// Calculates the log probability of the prior plus the log probability of the model, given the data.
        /// </summary>
        public double LogPosteriorProbability(double[] paramProposals, Datum[] data)
        {
            return LogPriorProbability(paramProposals) + LogProbabilityOfModelGivenTheData(paramProposals, data);
        }

        /// <summary>
        /// Calculates the log probability of the model with the proposed parameters, given the data.
        /// </summary>
        private double LogProbabilityOfModelGivenTheData(double[] paramProposals, Datum[] data)
        {
            double log_p = 0;

            for (var i = 0; i < data.Length; i++)
            {
                Datum datapoint = data[i];
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
