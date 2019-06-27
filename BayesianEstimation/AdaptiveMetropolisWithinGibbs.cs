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
    /// respectively, and calculating their probability of describing the data (Bayesian statistics). Eventually, after enough iterations, the 
    /// algorithm hopefully converges on an accurate answer, with some amount of measureable uncertainty for each parameter.
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
        public List<double[]> chain;
        double[] log_sd;
        double[] curr_state;
        double[] acceptance_count;
        int n_params;
        int batch_count;
        int batch_size;
        MersenneTwister random;
        double mean_mu;
        double sd_mu;
        double sigma_low;
        double sigma_high;
        List<double>[] data;
        bool isTwoSample;
        List<Tuple<double, double>>[] acceptableParameterRanges;

        public AdaptiveMetropolisWithinGibbs(List<double> data1, List<double> data2, int batch_size, List<Tuple<double, double>>[] acceptableParameterRanges, int? seed = null)
        {
            if (seed != null)
            {
                random = new MersenneTwister(seed.Value, true);
            }
            else
            {
                random = new MersenneTwister(true);
            }

            batch_count = 0;
            this.batch_size = batch_size;
            isTwoSample = data2 != null;

            if (isTwoSample)
            {
                data = new List<double>[2] { data1, data2 };
            }
            else
            {
                data = new List<double>[1] { data1 };
            }

            make_BEST_posterior_func(isTwoSample);
            this.acceptableParameterRanges = acceptableParameterRanges;
        }

        public double[] next_sample()
        {
            chain.Add(curr_state.ToArray());

            for (var param_i = 0; param_i < n_params; param_i++)
            {
                var parameterProposal = Normal.Sample(random, curr_state[param_i], Math.Exp(log_sd[param_i]));
                var proposed = curr_state.ToArray();
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
                        var curr_post_dens = PosteriorProbability(curr_state);
                        var prop_post_dens = PosteriorProbability(proposed);
                        if (double.IsNegativeInfinity(curr_post_dens) || double.IsNaN(curr_post_dens))
                        {
                            // if curr_post_dens is as bad as, say, negative infinity or NaN we should always jump
                            accept_prob = 1;
                        }

                        accept_prob = Math.Exp(prop_post_dens - curr_post_dens);
                        break;
                    }
                }

                if (accept_prob > random.NextDouble())
                {
                    acceptance_count[param_i]++;
                    curr_state = proposed;
                } // else do nothing
            }

            if (chain.Count % batch_size == 0)
            {
                batch_count++;

                for (var param_i = 0; param_i < n_params; param_i++)
                {
                    if (acceptance_count[param_i] / batch_size > 0.44)
                    {
                        log_sd[param_i] += Math.Min(0.01, 1.0 / Math.Sqrt(batch_count));
                    }
                    else if (acceptance_count[param_i] / batch_size < 0.44)
                    {
                        log_sd[param_i] -= Math.Min(0.01, 1.0 / Math.Sqrt(batch_count));
                    }
                    acceptance_count[param_i] = 0;
                }
            }

            return curr_state;
        }

        public double[] n_samples(int n)
        {
            chain = new List<double[]>();
            for (int i = 0; i < n - 1; i++)
            {
                next_sample();
            }

            return next_sample();
        }

        public void make_BEST_posterior_func(bool isTwoSample)
        {
            // initialize parameters
            var pooled = data.SelectMany(p => p.Select(v => v)).ToList();

            mean_mu = 0;
            sd_mu = 0.5;

            sigma_low = Math.Min(0.1, pooled.StandardDeviation() / 10.0);
            sigma_high = Math.Max(10, pooled.StandardDeviation() * 10.0);

            if (isTwoSample)
            {
                n_params = 5;
                curr_state = new double[] { mean_mu, pooled.Average(), pooled.StandardDeviation(), pooled.StandardDeviation(), 5 };
            }
            else
            {
                n_params = 3;
                curr_state = new double[] { mean_mu, pooled.StandardDeviation(), 5 };
            }

            log_sd = new double[n_params];
            acceptance_count = new double[n_params];
        }

        public double PosteriorProbability(double[] parameters)
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
                log_p += Math.Log(ContinuousUniform.PDF(sigma_low, sigma_high, sigma[sample]));

                double normalPriorProb = Normal.PDF(mean_mu, sd_mu, mu[sample]);
                log_p += Math.Log(normalPriorProb);

                for (var subj_i = 0; subj_i < data[sample].Count; subj_i++)
                {
                    double mmt = data[sample][subj_i];
                    double studentTPdf = StudentT.PDF(mu[sample], sigma[sample], nu, mmt);
                    log_p += Math.Log(studentTPdf);
                }
            }

            return log_p;
        }
    }
}
