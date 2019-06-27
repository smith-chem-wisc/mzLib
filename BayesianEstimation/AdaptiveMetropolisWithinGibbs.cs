using MathNet.Numerics.Distributions;
using MathNet.Numerics.Random;
using MathNet.Numerics.Statistics;
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
        public List<double[]> markovChain { get; private set; }
        private double[] logSd;
        private double[] currentState;
        private double[] acceptanceCount;
        private int numParams;
        private int batchCount;
        private int batchSize;
        private MersenneTwister random;
        private double meanMu;
        private double sdMu;
        private double sigmaLow;
        private double sigmaHigh;
        private List<double>[] Data;
        private bool isTwoSample;
        private List<Tuple<double, double>>[] acceptableParameterRanges;

        /// <summary>
        /// Construct the adaptive Metropolis within Gibbs sampler. Leave data2 null for one-sample data.
        /// </summary>
        public AdaptiveMetropolisWithinGibbs(List<double> data1, List<double> data2, int batch_size = 3, int? seed = null)
        {
            if (seed != null)
            {
                random = new MersenneTwister(seed.Value, true);
            }
            else
            {
                random = new MersenneTwister(true);
            }

            batchCount = 0;
            this.batchSize = batch_size;
            isTwoSample = data2 != null;

            if (isTwoSample)
            {
                Data = new List<double>[] { data1, data2 };
            }
            else
            {
                Data = new List<double>[] { data1 };
            }
            
            SetUpAcceptableParameterRanges();
            make_BEST_posterior_func(isTwoSample);
        }

        /// <summary>
        /// Burn in and then sample the MCMC chain
        /// </summary>
        public void Run(int burnin = 20000, int n = 20000)
        {
            n_samples(burnin);
            n_samples(n);
        }

        /// <summary>
        /// Sets up the acceptable range of values for each parameter. For example, standard deviation must be positive.
        /// </summary>
        protected void SetUpAcceptableParameterRanges()
        {
            if (acceptableParameterRanges == null)
            {
                if (isTwoSample)
                {
                    acceptableParameterRanges = new List<Tuple<double, double>>[5];
                    acceptableParameterRanges[0] = new List<Tuple<double, double>> { new Tuple<double, double>(double.NegativeInfinity, double.PositiveInfinity) }; // mu1 (can take any value)
                    acceptableParameterRanges[1] = new List<Tuple<double, double>> { new Tuple<double, double>(double.NegativeInfinity, double.PositiveInfinity) }; // mu2 (can take any value)
                    acceptableParameterRanges[2] = new List<Tuple<double, double>> { new Tuple<double, double>(0, double.PositiveInfinity) }; // sd1 (must be positive)
                    acceptableParameterRanges[3] = new List<Tuple<double, double>> { new Tuple<double, double>(0, double.PositiveInfinity) }; // sd2 (must be positive)
                    acceptableParameterRanges[4] = new List<Tuple<double, double>> { new Tuple<double, double>(0, double.PositiveInfinity) }; // nu (must be positive)
                }
                else
                {
                    acceptableParameterRanges = new List<Tuple<double, double>>[3];
                    acceptableParameterRanges[0] = new List<Tuple<double, double>> { new Tuple<double, double>(double.NegativeInfinity, double.PositiveInfinity) }; // mu (can take any value)
                    acceptableParameterRanges[1] = new List<Tuple<double, double>> { new Tuple<double, double>(0, double.PositiveInfinity) }; // sd (must be positive)
                    acceptableParameterRanges[2] = new List<Tuple<double, double>> { new Tuple<double, double>(0, double.PositiveInfinity) }; // nu (must be positive)
                }
            }
        }
        
        /// <summary>
        /// Sets up the posterior function.
        /// </summary>
        protected void make_BEST_posterior_func(bool isTwoSample)
        {
            // initialize parameters
            var pooled = Data.SelectMany(p => p.Select(v => v)).ToList();

            meanMu = pooled.Mean();
            sdMu = pooled.StandardDeviation();

            sigmaLow = Math.Min(0.1, pooled.StandardDeviation() / 10.0);
            sigmaHigh = Math.Max(10, pooled.StandardDeviation() * 10.0);

            // the initial parameters are set to:
            // mu: the pooled data's mean
            // sd: the pooled data's sd
            // nu: 5
            if (isTwoSample)
            {
                numParams = 5;
                currentState = new double[] { pooled.Average(), pooled.Average(), pooled.StandardDeviation(), pooled.StandardDeviation(), 5 };
            }
            else
            {
                numParams = 3;
                currentState = new double[] { pooled.Average(), pooled.StandardDeviation(), 5 };
            }

            logSd = new double[numParams];
            acceptanceCount = new double[numParams];
        }

        /// <summary>
        /// Calculates the probability that a set of parameters explains the data.
        /// </summary>
        protected double PosteriorProbability(double[] parameters)
        {
            double[] mu = null;
            double[] sigma = null;
            double nu = 0;
            if (isTwoSample)
            {
                mu = new double[] { parameters[0], parameters[1] };
                sigma = new double[] { parameters[2], parameters[3] };
                nu = parameters[4];
            }
            else
            {
                mu = new double[] { parameters[0] };
                sigma = new double[] { parameters[1] };
                nu = parameters[2];
            }

            double log_p = 0;
            log_p += Math.Log(Exponential.PDF(1.0 / 29.0, nu - 1.0)); // estimate nu

            for (var sample = 0; sample < mu.Length; sample++) // estimate mu and sd for each sample
            {
                log_p += Math.Log(ContinuousUniform.PDF(sigmaLow, sigmaHigh, sigma[sample]));

                double normalPriorProb = Normal.PDF(meanMu, sdMu, mu[sample]);
                log_p += Math.Log(normalPriorProb);

                for (var subj_i = 0; subj_i < Data[sample].Count; subj_i++)
                {
                    double mmt = Data[sample][subj_i];
                    double studentTPdf = StudentT.PDF(mu[sample], sigma[sample], nu, mmt);
                    log_p += Math.Log(studentTPdf);
                }
            }

            return log_p;
        }

        /// <summary>
        /// Samples the parameter distributions via MCMC.
        /// </summary>
        private double[] next_sample()
        {
            markovChain.Add(currentState.ToArray());

            for (var param_i = 0; param_i < numParams; param_i++)
            {
                var parameterProposal = Normal.Sample(random, currentState[param_i], Math.Exp(logSd[param_i]));
                var proposed = currentState.ToArray();
                proposed[param_i] = parameterProposal;
                double accept_prob = 0;
                var acceptableParamRanges = acceptableParameterRanges[param_i];

                foreach (var paramRange in acceptableParamRanges)
                {
                    if (parameterProposal < paramRange.Item1 || parameterProposal > paramRange.Item2)
                    {
                        accept_prob = 0;
                    }
                    else
                    {
                        var curr_post_dens = PosteriorProbability(currentState);
                        var prop_post_dens = PosteriorProbability(proposed);
                        if (double.IsNegativeInfinity(curr_post_dens) || double.IsNaN(curr_post_dens))
                        {
                            // if curr_post_dens is as bad as, say, negative infinity or NaN we should always jump
                            accept_prob = 1;
                            break;
                        }

                        accept_prob = Math.Exp(prop_post_dens - curr_post_dens);
                        break;
                    }
                }

                if (accept_prob > random.NextDouble())
                {
                    acceptanceCount[param_i]++;
                    currentState = proposed;
                } // else do nothing
            }

            if (markovChain.Count % batchSize == 0)
            {
                batchCount++;

                for (var param_i = 0; param_i < numParams; param_i++)
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
        private double[] n_samples(int n)
        {
            markovChain = new List<double[]>();
            for (int i = 0; i < n - 1; i++)
            {
                next_sample();
            }

            return next_sample();
        }
    }
}
