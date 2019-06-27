using BayesianEstimation;
using MathNet.Numerics.Statistics;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace Test
{
    [TestFixture]
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
            AdaptiveMetropolisWithinGibbs sampler = new AdaptiveMetropolisWithinGibbs(data, null, seed: 0);

            // burn in and then sample the MCMC chain
            sampler.Run(1000, 1000);

            var chain = sampler.markovChain;

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
            Assert.That(Math.Round(muPointEstimate, 3) == 0.957);

            // std dev is estimated at 0.17
            Assert.That(Math.Round(sigmaPointEstimate, 3) == 0.167);

            // nu (degree of normality) is estimated at 1.2
            Assert.That(Math.Round(nuPointEstimate, 3) == 1.217);

            // instead of only a point estimate of the mean, the Bayesian method also gives a range of credible values.
            // sort of like a 95% confidence interval, we can construct a 95% highest density interval where 95% of the 
            // probability density for a parameter is contained
            var highestDensityInterval = BayesianEstimation.Util.GetHighestDensityInterval(mus.ToArray());
            Assert.That(Math.Round(highestDensityInterval.hdi_start, 3) == 0.840);
            Assert.That(Math.Round(highestDensityInterval.hdi_end, 3) == 1.086);
        }

        [Test]
        /// <summary>
        /// Bayesian estimation of the difference in means between two samples.
        /// </summary>
        public static void TestTwoSampleBayesianEstimation()
        {
            List<double> data1 = new List<double> { 1.0, 0.9, 1.1 };
            List<double> data2 = new List<double> { 1.0, 0.9, 1.1 };
            
            AdaptiveMetropolisWithinGibbs sampler = new AdaptiveMetropolisWithinGibbs(data1, data2, seed: 0);

            // burn in and then sample the MCMC chain
            sampler.Run(1000, 1000);
            
            var chain = sampler.markovChain;

            var muDiffs = chain.Select(p => p[1] - p[0]).ToList(); // difference in means
            var sigma1s = chain.Select(p => p[2]).ToList(); // std dev estimate for sample 1
            var sigma2s = chain.Select(p => p[3]).ToList(); // std dev estimate for sample 2
            var nus = chain.Select(p => p[4]).ToList(); // nu estimate

            double avgMeanDiff = muDiffs.Average(); // point estimate of mean diff
            Assert.That(Math.Round(avgMeanDiff, 3) == -0.003);

            var highestDensityInterval = BayesianEstimation.Util.GetHighestDensityInterval(muDiffs.ToArray());
            Assert.That(Math.Round(highestDensityInterval.hdi_start, 3) == -0.176);
            Assert.That(Math.Round(highestDensityInterval.hdi_end, 3) == 0.172);
        }
    }
}
