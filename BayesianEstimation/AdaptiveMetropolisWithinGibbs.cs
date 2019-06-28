using MathNet.Numerics.Distributions;
using MathNet.Numerics.Random;
using System;
using System.Collections.Generic;
using System.Linq;

namespace BayesianEstimation
{
    /// <summary>
    /// This class uses an adaptive Metropolis within Gibbs sampler to estimate the mean (mu), standard deviation (sigma), and degree of 
    /// normality (nu) of a series of data points. The data points are assumed to come from a t-distribution (similar to a 
    /// normal distribution, but with fatter tails; this method is more tolerant of outliers than fitting a normal distribution). 
    /// The mu, sigma, and nu parameters of the t-distribution are estimated by sampling from normal, uniform, and exponential distributions, 
    /// respectively, using a Markov Chain Monte Carlo (MCMC) method and calculating their probability of describing the data (Bayesian statistics). 
    /// Eventually, after enough iterations, the algorithm hopefully converges on an accurate answer, with some amount of measureable 
    /// uncertainty for each parameter.
    /// 
    /// You can fit one-sample or two-sample data by making the "data2" parameter in the constructor null or not null, respectively.
    /// 
    /// If your data do not conform to a t-distribution, you can inherit this class and then define a new PosteriorProbability
    /// function. The rest of the algorithm remains the same (the Metropolis within Gibbs sampler samples from the distributions 
    /// you define in order to fit the distribution you want to your data). See FlashLfqBayesianEstimator.cs for an example of defining
    /// a different posterior probability function.
    /// 
    /// See:
    /// Bayesian Estimation Supersedes the t Test. Kruschke, J. K., Journal of Experimental Psychology, 2013.
    /// 
    /// This code is adapted from Javascript code taken from http://sumsar.net/best_online/ . This webpage is a useful
    /// interactive example of J. Kruschke's Bayesian equivalent of a t-test.
    /// </summary>
    public class AdaptiveMetropolisWithinGibbs
    {
        public List<double[]> MarkovChain { get; private set; }
        private double[] currentState;
        private double[] proposedState;
        private MersenneTwister random;
        private Model model;
        private double[] data;
        private readonly double[] logSd;
        private readonly double[] acceptanceCount;
        private readonly int batchSize;
        private int batchCount;

        /// <summary>
        /// Construct the adaptive Metropolis within Gibbs sampler. Leave data2 null for one-sample data.
        /// </summary>
        public AdaptiveMetropolisWithinGibbs(double[] data, Model model, int batch_size = 3, int? seed = null)
        {
            if (seed != null)
            {
                random = new MersenneTwister(seed.Value, true);
            }
            else
            {
                random = new MersenneTwister(true);
            }

            this.model = model;
            this.batchSize = batch_size;
            this.data = data;
            logSd = new double[model.modelParameters.Length];
            acceptanceCount = new double[model.modelParameters.Length];
            currentState = new double[model.modelParameters.Length];
            proposedState = new double[model.modelParameters.Length];

            for(int i = 0; i < currentState.Length; i++)
            {
                currentState[i] = model.modelParameters[i].initialGuess;
            }
        }

        /// <summary>
        /// Burn in and then sample the MCMC chain.
        /// </summary>
        public void Run(int burnin = 20000, int n = 20000)
        {
            nSamples(burnin);
            nSamples(n);
        }

        /// <summary>
        /// Samples the parameter distributions via MCMC.
        /// </summary>
        private double[] next_sample()
        {
            MarkovChain.Add(currentState.ToArray());

            for (var p = 0; p < model.modelParameters.Length; p++)
            {
                var parameterProposal = Normal.Sample(random, currentState[p], Math.Exp(logSd[p]));
                
                proposedState = currentState.ToArray();
                proposedState[p] = parameterProposal;

                double acceptProbability = 0;
                var acceptableParamRanges = model.modelParameters[p].acceptableParameterRanges;

                foreach (var acceptableParamRange in acceptableParamRanges)
                {
                    if (parameterProposal < acceptableParamRange.Item1 || parameterProposal > acceptableParamRange.Item2)
                    {
                        acceptProbability = 0;
                    }
                    else
                    {
                        double currentParamPosteriorProbability = model.LogPosteriorProbability(currentState, data);
                        double proposedParamPosteriorProbability = model.LogPosteriorProbability(proposedState, data);

                        if (double.IsNegativeInfinity(currentParamPosteriorProbability) || double.IsNaN(currentParamPosteriorProbability))
                        {
                            // if currentParamPosteriorProbability is as bad as, say, negative infinity or NaN we should always jump
                            acceptProbability = 1;
                            break;
                        }

                        acceptProbability = Math.Exp(proposedParamPosteriorProbability - currentParamPosteriorProbability);
                        break;
                    }
                }

                if (acceptProbability > random.NextDouble())
                {
                    acceptanceCount[p]++;
                    currentState = proposedState;
                } // else do nothing
            }

            if (MarkovChain.Count % batchSize == 0)
            {
                batchCount++;

                for (var param_i = 0; param_i < model.modelParameters.Length; param_i++)
                {
                    if (acceptanceCount[param_i] / batchSize > 0.44)
                    {
                        logSd[param_i] += Math.Min(0.01, 1.0 / Math.Sqrt(batchCount));
                    }
                    else if (acceptanceCount[param_i] / batchSize < 0.44)
                    {
                        logSd[param_i] -= Math.Min(0.01, 1.0 / Math.Sqrt(batchCount));
                    }
                    acceptanceCount[param_i] = 0;
                }
            }

            return currentState;
        }

        /// <summary>
        /// Resets the MCMC chain and samples the parameter distributions n times via MCMC, storing it in a new chain.
        /// </summary>
        private double[] nSamples(int n)
        {
            MarkovChain = new List<double[]>();
            for (int i = 0; i < n - 1; i++)
            {
                next_sample();
            }

            return next_sample();
        }
    }
}
