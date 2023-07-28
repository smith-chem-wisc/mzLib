using BayesianEstimation;
using MathNet.Numerics.Statistics;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class TestBayesianEstimation
    {
        private static Stopwatch Stopwatch { get; set; }

        [SetUp]
        public static void Setup()
        {
            Stopwatch = new Stopwatch();
            Stopwatch.Start();
        }

        [TearDown]
        public static void TearDown()
        {
            Console.WriteLine($"Analysis time: {Stopwatch.Elapsed.Hours}h {Stopwatch.Elapsed.Minutes}m {Stopwatch.Elapsed.Seconds}s");
        }

        [Test]
        /// <summary>
        /// This test demonstrates the utility of estimating the mean and standard deviation of a set of datapoints using a Bayesian 
        /// method (which uses Monte-Carlo sampling) rather than using the standard equations for the mean and standard deviations. The example data
        /// contains an outlier, which skews the estimate of mean and standard deviation if you compute these using mean/std dev equations.
        /// The Bayesian method fits a series of t-distributions to the data, which is more tolerant of the outlier, providing a better
        /// estimate of the mean and standard deviation.
        /// </summary>
        public static void TestOneSampleBayesianEstimation()
        {
            List<double> data = new List<double> { 1.0, 1.0, 1.1, 0.9, 0.8, 0.9, 1.0, 100.0 };

            // the std dev and mean of the data are very high because the data contain an outlier (100)
            double stdev = data.StandardDeviation(); // sd of the data including outlier is 35.0
            double mean = data.Mean(); // mean of the data including outlier is 13.3

            // let's see what the mean and std dev are without the outlier
            List<double> dataWithoutTheOutlier = new List<double> { 1.0, 1.0, 1.1, 0.9, 0.8, 0.9, 1.0 };
            double stdevWithoutTheOutlier = dataWithoutTheOutlier.StandardDeviation(); // sd of the data excluding the outlier is 0.098
            double meanWithoutTheOutlier = dataWithoutTheOutlier.Mean(); // mean of the data excluding the outlier is 0.96

            // let's do a Bayesian estimate of the mean and standard deviation of the data set that includes the outlier...
            AdaptiveMetropolisWithinGibbs sampler = new AdaptiveMetropolisWithinGibbs(
                data.ToArray(), 
                new StudentTDistributionModel(
                    priorMuMean: data.Average(), priorMuSd: data.StandardDeviation() * 1000000, muInitialGuess: data.Average(), 
                    priorSdStart: data.StandardDeviation() / 1000, priorSdEnd: data.StandardDeviation() * 1000, sdInitialGuess: data.StandardDeviation(), 
                    priorNuExponent: 1.0 / 29.0, nuInitialGuess: 5), 
                seed: 0);

            // burn in and then sample the MCMC chain 
            sampler.Run(1000, 1000);

            var chain = sampler.MarkovChain;

            // each iteration of the MCMC chain produces one t-distribution fit result
            // in this example, 1000 different t-distributions are fit to the data
            // each of these distributions has a mean, sd, and degree of normality
            var mus = chain.Select(p => p[0]).ToList();
            var sigmas = chain.Select(p => p[1]).ToList();
            var nus = chain.Select(p => p[2]).ToList();

            Assert.That(mus.Count == 1000);

            // here are the mean/sd/normality results of the Bayesian estimation:
            // (the point estimate is the average of each parameter over the 1000 iterations)
            double muPointEstimate = mus.Average(); // point estimate of mu (mean)
            double sigmaPointEstimate = sigmas.Average(); // point estimate of sigma (std dev)
            double nuPointEstimate = nus.Average(); // point estimate of nu (degree of normality)

            // the mean is estimated at 0.96
            Assert.That(Math.Round(muPointEstimate, 3) == 0.962);

            // std dev is estimated at 0.10
            Assert.That(Math.Round(sigmaPointEstimate, 3) == 0.104);

            // nu (degree of normality) is estimated at 0.65
            Assert.That(Math.Round(nuPointEstimate, 3) == 0.654);

            // instead of only a point estimate of the mean, the Bayesian method also gives a range of credible values.
            // sort of like a 95% confidence interval, we can construct a 95% highest density interval where 95% of the 
            // probability density for a parameter is contained
            var highestDensityInterval = BayesianEstimation.Util.GetHighestDensityInterval(mus.ToArray());
            Assert.That(Math.Round(highestDensityInterval.hdi_start, 3) == 0.862);
            Assert.That(Math.Round(highestDensityInterval.hdi_end, 3) == 1.062);
        }


        [Test]
        /// <summary>
        /// Bayesian estimation of the difference in means between two samples.
        /// </summary>
        public static void TestTwoSampleBayesianEstimation()
        {
            List<double> data1 = new List<double> { 1.0, 0.9, 1.1, 1.0, 0.9, 1.1, 1.0, 0.9, 1.1 };
            List<double> data2 = new List<double> { 1.0, 0.9, 1.1, 1.0, 0.9, 1.1, 1.0, 0.9, 1.1 };

            double pooledMean = data1.Concat(data2).Mean();
            double pooledSd = data1.Concat(data2).StandardDeviation();

            // construct a t-distribution model fitter for each dataset
            var sampler1 = new AdaptiveMetropolisWithinGibbs(
                data1.ToArray(),
                new StudentTDistributionModel(
                    priorMuMean: pooledMean, priorMuSd: pooledSd * 1000000, muInitialGuess: pooledMean,
                    priorSdStart: pooledSd / 1000, priorSdEnd: pooledSd * 1000, sdInitialGuess: pooledSd,
                    priorNuExponent: 1.0 / 29.0, nuInitialGuess: 5),
                seed: 0);

            var sampler2 = new AdaptiveMetropolisWithinGibbs(
                data2.ToArray(),
                new StudentTDistributionModel(
                    priorMuMean: pooledMean, priorMuSd: pooledSd * 1000000, muInitialGuess: pooledMean,
                    priorSdStart: pooledSd / 1000, priorSdEnd: pooledSd * 1000, sdInitialGuess: pooledSd,
                    priorNuExponent: 1.0 / 29.0, nuInitialGuess: 5),
                seed: 1);

            // burn in and then sample the MCMC chains
            sampler1.Run(1000, 1000);
            sampler2.Run(1000, 1000);

            var chain1 = sampler1.MarkovChain;
            var chain2 = sampler2.MarkovChain;

            // calculate difference in means
            List<double> meanDiffs = new List<double>();
            for(int i = 0; i < chain1.Count; i++)
            {
                meanDiffs.Add(chain1[i][0] - chain2[i][0]);
            }
            
            double avgMeanDiff = meanDiffs.Average(); // point estimate of mean diff
            Assert.That(Math.Round(avgMeanDiff, 3) == -0.011);

            var highestDensityInterval = BayesianEstimation.Util.GetHighestDensityInterval(meanDiffs.ToArray());
            Assert.That(Math.Round(highestDensityInterval.hdi_start, 3) == -0.118);
            Assert.That(Math.Round(highestDensityInterval.hdi_end, 3) == 0.095);
        }
    }
}
