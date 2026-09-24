using MathNet.Numerics.Distributions;
using NUnit.Framework;
using Statistics;
using System;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Test.Statistics;

/// <summary>
/// Closed-form values, textbook formulas, simulations that recover known parameters, and thread-count
/// determinism for logistic regression, the random-intercept mixed model, Spearman correlation and
/// p-value combination. Comparisons with R's glm, nlme and cor.test are in ReferenceComparisonTests.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class GlmMixedCorrelationTests
{
    // ---------------------------------------------------------------- p-value combination

    [Test]
    public void FisherOfOnePValueIsThatPValue()
    {
        var r = PValueCombination.Fisher(new[] { 0.037 });
        Assert.That(r.PValue, Is.EqualTo(0.037).Within(1e-14));
        Assert.That(r.Studies, Is.EqualTo(1));
        Assert.That(r.Method, Is.EqualTo(PValueCombinationMethod.Fisher));
    }

    [Test]
    public void FisherOfTwoPValuesMatchesClosedForm()
    {
        // χ²(4) upper tail at x = −2 ln(p1 p2) is p1 p2 (1 − ln(p1 p2)).
        double p1 = 0.04, p2 = 0.2, prod = p1 * p2;
        var r = PValueCombination.Fisher(new[] { p1, p2 });
        Assert.That(r.Statistic, Is.EqualTo(-2 * Math.Log(prod)).Within(1e-12));
        Assert.That(r.PValue, Is.EqualTo(prod * (1 - Math.Log(prod))).Within(1e-14));
    }

    [Test]
    public void StoufferMatchesClosedForm()
    {
        double z = Normal.InvCDF(0, 1, 0.95);
        var r = PValueCombination.Stouffer(new[] { 0.05, 0.05 });
        Assert.That(r.Statistic, Is.EqualTo(2 * z / Math.Sqrt(2)).Within(1e-12));
        Assert.That(r.PValue, Is.EqualTo(Normal.CDF(0, 1, -Math.Sqrt(2) * z)).Within(1e-14));
        Assert.That(PValueCombination.Stouffer(new[] { 0.3 }).PValue, Is.EqualTo(0.3).Within(1e-14));
    }

    [Test]
    public void StoufferWeightsCountLargerStudiesMore()
    {
        var p = new[] { 0.01, 0.6 };
        var equal = PValueCombination.Stouffer(p, new[] { 2.0, 2.0 });
        Assert.That(equal.PValue, Is.EqualTo(PValueCombination.Stouffer(p).PValue).Within(1e-14), "equal weights are unweighted");
        Assert.That(PValueCombination.Stouffer(p, new[] { 10.0, 1.0 }).PValue, Is.LessThan(equal.PValue));
    }

    [Test]
    public void CombinationOmitsMissingAndReportsNothingAsNaN()
    {
        var r = PValueCombination.Fisher(new[] { 0.1, double.NaN, 0.2 });
        Assert.That(r.Studies, Is.EqualTo(2));
        Assert.That(r.PValue, Is.EqualTo(PValueCombination.Fisher(new[] { 0.1, 0.2 }).PValue));
        Assert.That(PValueCombination.Stouffer(new[] { double.NaN }).PValue, Is.NaN);
        Assert.Throws<ArgumentOutOfRangeException>(() => PValueCombination.Fisher(new[] { 1.5 }));
    }

    [Test]
    public void OneSidedFollowsTheSignOfTheEffect()
    {
        Assert.That(PValueCombination.OneSided(0.1, 2.0), Is.EqualTo(0.05));
        Assert.That(PValueCombination.OneSided(0.1, -2.0), Is.EqualTo(0.95));
        Assert.That(PValueCombination.OneSided(0.1, 0), Is.NaN);
        // Opposite effects cancel under Stouffer instead of reinforcing each other.
        var cancel = PValueCombination.Stouffer(new[] { PValueCombination.OneSided(0.01, 1), PValueCombination.OneSided(0.01, -1) });
        Assert.That(cancel.PValue, Is.EqualTo(0.5).Within(1e-12));
    }

    // ---------------------------------------------------------------- Spearman

    [Test]
    public void SpearmanPerfectOrderIsExactAndOneOverFactorial()
    {
        var r = SpearmanCorrelation.Correlate(new double[] { 1, 2, 3, 4, 5 }, new double[] { 10, 20, 30, 40, 50 });
        Assert.That(r.Rho, Is.EqualTo(1).Within(1e-15));
        Assert.That(r.Method, Is.EqualTo(SpearmanPValueMethod.Exact));
        Assert.That(r.PValue, Is.EqualTo(2.0 / 120).Within(1e-15), "one ordering in 5! is this extreme, both tails");
        var neg = SpearmanCorrelation.Correlate(new double[] { 1, 2, 3, 4, 5 }, new double[] { 5, 4, 3, 2, 1 });
        Assert.That(neg.Rho, Is.EqualTo(-1).Within(1e-15));
        Assert.That(neg.PValue, Is.EqualTo(2.0 / 120).Within(1e-15));
    }

    [Test]
    public void SpearmanIsPearsonOfAverageRanks()
    {
        var x = new double[] { 3.1, 1.0, 2.2, 2.2, 5.0, 4.4, 0.3, 6.1 };
        var y = new double[] { 1.0, 0.5, 2.0, 1.5, 3.0, 3.0, 0.1, 2.5 };
        var rx = SpearmanCorrelation.Ranks(x, out bool tx);
        Assert.That(tx, Is.True);
        Assert.That(rx, Is.EqualTo(new[] { 5, 2, 3.5, 3.5, 7, 6, 1, 8 }));
        var ry = SpearmanCorrelation.Ranks(y, out _);
        double pearson = MathNet.Numerics.Statistics.Correlation.Pearson(rx, ry);
        var r = SpearmanCorrelation.Correlate(x, y);
        Assert.That(r.Rho, Is.EqualTo(pearson).Within(1e-14));
        Assert.That(r.HasTies, Is.True);
        Assert.That(r.Method, Is.EqualTo(SpearmanPValueMethod.Asymptotic), "ties force the t approximation");
        double t = pearson * Math.Sqrt(6 / (1 - pearson * pearson));
        Assert.That(r.PValue, Is.EqualTo(2 * StudentT.CDF(0, 1, 6, -Math.Abs(t))).Within(1e-14));
    }

    [Test]
    public void SpearmanOmitsIncompletePairsAndRefusesTooFew()
    {
        var r = SpearmanCorrelation.Correlate(new[] { 1, 2, double.NaN, 4, 5 }, new[] { 2, 1, 3, double.PositiveInfinity, 5 });
        Assert.That(r.N, Is.EqualTo(3));
        var few = SpearmanCorrelation.Correlate(new[] { 1.0, 2 }, new[] { 1.0, 2 });
        Assert.That(few.Method, Is.EqualTo(SpearmanPValueMethod.NotEstimable));
        Assert.That(few.Rho, Is.NaN);
        var constant = SpearmanCorrelation.Correlate(new[] { 1.0, 1, 1, 1 }, new[] { 1.0, 2, 3, 4 });
        Assert.That(constant.Rho, Is.NaN);
        Assert.That(constant.PValue, Is.NaN);
    }

    [Test]
    public void SpearmanExactPValueIsAProperDistribution()
    {
        // A middling ordering at n = 7 takes the exact path and a valid p-value; n = 10 switches to the approximation.
        var r = SpearmanCorrelation.Correlate(new double[] { 1, 2, 3, 4, 5, 6, 7 }, new double[] { 4, 7, 1, 5, 2, 6, 3 });
        Assert.That(r.Method, Is.EqualTo(SpearmanPValueMethod.Exact));
        Assert.That(r.PValue, Is.InRange(0.0, 1.0));
        var big = SpearmanCorrelation.Correlate(Enumerable.Range(0, 10).Select(i => (double)i).ToArray(),
                                                Enumerable.Range(0, 10).Select(i => (double)(i * i % 7)).ToArray());
        Assert.That(big.Method, Is.EqualTo(SpearmanPValueMethod.Asymptotic), "n above ExactMaxN");
    }

    // ---------------------------------------------------------------- logistic regression

    [Test]
    public void LogisticWithOneBinaryCovariateIsTheTwoByTwoTable()
    {
        // Group 0: 3 events of 10; group 1: 7 of 10. MLE slope = log odds ratio, SE = √(1/a + 1/b + 1/c + 1/d).
        int n = 20;
        var design = new double[n, 2];
        var y = new double[1, n];
        for (int s = 0; s < n; s++)
        {
            design[s, 0] = 1;
            design[s, 1] = s < 10 ? 0 : 1;
            y[0, s] = s < 10 ? (s < 3 ? 1 : 0) : (s < 17 ? 1 : 0);
        }
        var fit = LogisticRegression.Fit(y, design, new[] { "Intercept", "group" });
        Assert.That(fit.Status[0], Is.EqualTo(FeatureFitStatus.Fitted));
        Assert.That(fit.Coefficient(0, 0), Is.EqualTo(Math.Log(3.0 / 7)).Within(1e-10));
        Assert.That(fit.Coefficient(0, 1), Is.EqualTo(Math.Log(7.0 / 3 / (3.0 / 7))).Within(1e-10));
        Assert.That(fit.StandardError(0, 1), Is.EqualTo(Math.Sqrt(1 / 3.0 + 1 / 7.0 + 1 / 7.0 + 1 / 3.0)).Within(1e-9));
        Assert.That(fit.StandardError(0, 0), Is.EqualTo(Math.Sqrt(1 / 3.0 + 1 / 7.0)).Within(1e-9));
        double z = fit.Coefficient(0, 1) / fit.StandardError(0, 1);
        Assert.That(fit.PValue(0, 1), Is.EqualTo(2 * Normal.CDF(0, 1, -Math.Abs(z))).Within(1e-14));
        Assert.That(fit.Events[0], Is.EqualTo(10));
        double dev = -2 * (3 * Math.Log(0.3) + 7 * Math.Log(0.7) + 7 * Math.Log(0.7) + 3 * Math.Log(0.3));
        Assert.That(fit.Deviance[0], Is.EqualTo(dev).Within(1e-9));
    }

    [Test]
    public void LogisticReportsSeparationInsteadOfHugeCoefficients()
    {
        var design = new double[8, 2];
        var y = new double[3, 8];
        for (int s = 0; s < 8; s++)
        {
            design[s, 0] = 1; design[s, 1] = s;
            y[0, s] = s < 4 ? 0 : 1;          // complete separation on x
            y[1, s] = 1;                      // no non-events
            y[2, s] = s < 6 ? double.NaN : s % 2; // two observations for two coefficients
        }
        var fit = LogisticRegression.Fit(y, design);
        Assert.That(fit.Status[0], Is.EqualTo(FeatureFitStatus.Separated));
        Assert.That(fit.Coefficient(0, 1), Is.NaN);
        Assert.That(fit.Status[1], Is.EqualTo(FeatureFitStatus.Separated));
        Assert.That(fit.Status[2], Is.EqualTo(FeatureFitStatus.TooFewObservations));
    }

    [Test]
    public void LogisticRejectsNonBinaryResponses()
    {
        var design = new double[4, 1] { { 1 }, { 1 }, { 1 }, { 1 } };
        Assert.Throws<ArgumentException>(() => LogisticRegression.Fit(new double[1, 4] { { 0, 1, 0.5, 1 } }, design));
    }

    [Test]
    public void LogisticRecoversSimulatedCoefficients()
    {
        var rng = new Random(20260924);
        int n = 4000;
        double b0 = -0.5, b1 = 1.2;
        var design = new double[n, 2];
        var y = new double[1, n];
        for (int s = 0; s < n; s++)
        {
            double x = Normal.Sample(rng, 0, 1);
            design[s, 0] = 1; design[s, 1] = x;
            y[0, s] = rng.NextDouble() < 1 / (1 + Math.Exp(-(b0 + b1 * x))) ? 1 : 0;
        }
        var fit = LogisticRegression.Fit(y, design);
        Assert.That(fit.Coefficient(0, 0), Is.EqualTo(b0).Within(4 * fit.StandardError(0, 0)));
        Assert.That(fit.Coefficient(0, 1), Is.EqualTo(b1).Within(4 * fit.StandardError(0, 1)));
    }

    // ---------------------------------------------------------------- mixed model

    private static (double[,] y, double[,] design, string[] groups) OneWay(double[][] groupValues)
    {
        int n = groupValues.Sum(g => g.Length), s = 0;
        var y = new double[1, n];
        var design = new double[n, 1];
        var groups = new string[n];
        for (int g = 0; g < groupValues.Length; g++)
            foreach (var v in groupValues[g]) { y[0, s] = v; design[s, 0] = 1; groups[s] = $"g{g}"; s++; }
        return (y, design, groups);
    }

    [Test]
    public void BalancedOneWayRemlIsTheAnovaEstimator()
    {
        // Balanced one-way random effects: REML gives σ² = MSW and τ² = (MSB − MSW)/n when that is positive,
        // the grand mean, and SE² = MSB / (G n).
        var data = new[] { new[] { 5.0, 6, 7 }, new[] { 9.0, 10, 12 }, new[] { 2.0, 3, 3 }, new[] { 7.0, 8, 6 } };
        var (y, design, groups) = OneWay(data);
        int G = 4, k = 3;
        double grand = data.SelectMany(v => v).Average();
        double msb = k * data.Sum(g => Math.Pow(g.Average() - grand, 2)) / (G - 1);
        double msw = data.Sum(g => g.Sum(v => Math.Pow(v - g.Average(), 2))) / (G * (k - 1));
        var fit = MixedModel.Fit(y, design, groups);
        Assert.That(fit.Status[0], Is.EqualTo(FeatureFitStatus.Fitted));
        Assert.That(fit.Estimator, Is.EqualTo(VarianceEstimator.Reml));
        Assert.That(fit.ResidualVariance[0], Is.EqualTo(msw).Within(1e-8 * msw));
        Assert.That(fit.GroupVariance[0], Is.EqualTo((msb - msw) / k).Within(1e-8 * msb));
        Assert.That(fit.Coefficient(0, 0), Is.EqualTo(grand).Within(1e-10));
        Assert.That(fit.StandardError(0, 0), Is.EqualTo(Math.Sqrt(msb / (G * k))).Within(1e-8));
        Assert.That(fit.DegreesOfFreedom(0, 0), Is.EqualTo(G * k - G), "the intercept takes the within-group df");
    }

    [Test]
    public void NegligibleGroupVarianceGivesTheBoundaryAndOrdinaryLeastSquares()
    {
        // Group means closer together than chance allows: MSB < MSW, so τ̂² = 0 and the fit is OLS.
        var data = new[] { new[] { 1.0, 9, 5 }, new[] { 2.0, 8, 5.2 }, new[] { 0.5, 9.5, 4.9 } };
        var (y, design, groups) = OneWay(data);
        var fit = MixedModel.Fit(y, design, groups);
        var ols = LinearModel.Fit(y, design);
        Assert.That(fit.GroupVariance[0], Is.EqualTo(0));
        Assert.That(fit.Coefficient(0, 0), Is.EqualTo(ols.Coefficient(0, 0)).Within(1e-12));
        Assert.That(fit.ResidualVariance[0], Is.EqualTo(ols.Sigma[0] * ols.Sigma[0]).Within(1e-10));
    }

    [Test]
    public void ContainmentDegreesOfFreedomSeparateWithinAndBetweenCovariates()
    {
        // 6 groups of 4. 'age' varies within groups, 'site' is constant within each group.
        int G = 6, k = 4, n = G * k;
        var rng = new Random(7);
        var design = new double[n, 3];
        var y = new double[1, n];
        var groups = new string[n];
        for (int s = 0; s < n; s++)
        {
            int g = s / k;
            design[s, 0] = 1; design[s, 1] = s % k; design[s, 2] = g % 2;
            groups[s] = $"d{g}";
            y[0, s] = 1 + 0.5 * design[s, 1] + 2 * design[s, 2] + g * 0.3 + Normal.Sample(rng, 0, 0.5);
        }
        var fit = MixedModel.Fit(y, design, groups, new[] { "Intercept", "age", "site" });
        Assert.That(fit.DegreesOfFreedom(0, 1), Is.EqualTo(n - G - 1), "within-group covariate");
        Assert.That(fit.DegreesOfFreedom(0, 2), Is.EqualTo(G - 1 - 1), "group-level covariate");
        Assert.That(fit.DegreesOfFreedom(0, 0), Is.EqualTo(n - G - 1));
        Assert.That(fit.Groups[0], Is.EqualTo(G));
    }

    [Test]
    public void MixedModelRecoversSimulatedParameters()
    {
        var rng = new Random(20260924);
        int G = 60, k = 20, n = G * k;
        double beta = 0.8, tau = 1.5, sigma = 0.7;
        var design = new double[n, 2];
        var y = new double[1, n];
        var groups = new string[n];
        for (int g = 0; g < G; g++)
        {
            double u = Normal.Sample(rng, 0, tau);
            for (int i = 0; i < k; i++)
            {
                int s = g * k + i;
                double x = Normal.Sample(rng, 0, 1);
                design[s, 0] = 1; design[s, 1] = x; groups[s] = $"g{g}";
                y[0, s] = 2 + beta * x + u + Normal.Sample(rng, 0, sigma);
            }
        }
        var fit = MixedModel.Fit(y, design, groups);
        Assert.That(fit.Coefficient(0, 1), Is.EqualTo(beta).Within(4 * fit.StandardError(0, 1)));
        Assert.That(Math.Sqrt(fit.GroupVariance[0]), Is.EqualTo(tau).Within(0.35 * tau));
        Assert.That(Math.Sqrt(fit.ResidualVariance[0]), Is.EqualTo(sigma).Within(0.05 * sigma));
        var ml = MixedModel.Fit(y, design, groups, estimator: VarianceEstimator.MaximumLikelihood);
        Assert.That(ml.Estimator, Is.EqualTo(VarianceEstimator.MaximumLikelihood));
        Assert.That(ml.GroupVariance[0], Is.LessThan(fit.GroupVariance[0]), "ML shrinks variances relative to REML");
    }

    [Test]
    public void MixedModelRefusesWhatItCannotEstimate()
    {
        var design = new double[4, 1] { { 1 }, { 1 }, { 1 }, { 1 } };
        var y = new double[2, 4] { { 1, 2, 3, 4 }, { 1, 2, 3, 4 } };
        var oneGroup = MixedModel.Fit(y, design, new[] { "a", "a", "a", "a" });
        Assert.That(oneGroup.Status[0], Is.EqualTo(FeatureFitStatus.TooFewGroups));
        var singletons = MixedModel.Fit(y, design, new[] { "a", "b", "c", "d" });
        Assert.That(singletons.Status[0], Is.EqualTo(FeatureFitStatus.TooFewGroups));
        Assert.That(singletons.Coefficient(0, 0), Is.NaN);
        Assert.Throws<ArgumentException>(() => MixedModel.Fit(y, design, new[] { "a", "b", "", "d" }));
    }

    // ---------------------------------------------------------------- determinism

    [Test]
    public void OutputIsIdenticalForEveryThreadCount()
    {
        var rng = new Random(11);
        int features = 50, n = 24;
        var design = new double[n, 2];
        var cont = new double[features, n];
        var bin = new double[features, n];
        var groups = new string[n];
        for (int s = 0; s < n; s++) { design[s, 0] = 1; design[s, 1] = s / 3.0; groups[s] = $"g{s % 4}"; }
        for (int f = 0; f < features; f++)
            for (int s = 0; s < n; s++)
            {
                cont[f, s] = rng.NextDouble() < 0.1 ? double.NaN : Normal.Sample(rng, s * 0.01 * f, 1);
                bin[f, s] = rng.NextDouble() < 0.1 ? double.NaN : (rng.NextDouble() < 0.5 ? 1 : 0);
            }
        var m1 = MixedModel.Fit(cont, design, groups, maxThreads: 1);
        var m8 = MixedModel.Fit(cont, design, groups, maxThreads: 8);
        var l1 = LogisticRegression.Fit(bin, design, maxThreads: 1);
        var l8 = LogisticRegression.Fit(bin, design, maxThreads: 8);
        for (int f = 0; f < features; f++)
            for (int j = 0; j < 2; j++)
            {
                Assert.That(BitConverter.DoubleToInt64Bits(m8.Coefficient(f, j)), Is.EqualTo(BitConverter.DoubleToInt64Bits(m1.Coefficient(f, j))));
                Assert.That(BitConverter.DoubleToInt64Bits(m8.StandardError(f, j)), Is.EqualTo(BitConverter.DoubleToInt64Bits(m1.StandardError(f, j))));
                Assert.That(BitConverter.DoubleToInt64Bits(l8.Coefficient(f, j)), Is.EqualTo(BitConverter.DoubleToInt64Bits(l1.Coefficient(f, j))));
                Assert.That(BitConverter.DoubleToInt64Bits(l8.StandardError(f, j)), Is.EqualTo(BitConverter.DoubleToInt64Bits(l1.StandardError(f, j))));
            }
        Assert.That(m8.GroupVariance, Is.EqualTo(m1.GroupVariance));
        Assert.That(l8.Status, Is.EqualTo(l1.Status));
    }
}
