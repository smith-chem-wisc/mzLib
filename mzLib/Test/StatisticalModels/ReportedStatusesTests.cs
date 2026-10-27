using MathNet.Numerics.Distributions;
using NUnit.Framework;
using StatisticalModels;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Test.StatisticalModels;

/// <summary>
/// Each model's <c>ReportedStatuses</c>: what it lists, that it cannot be changed, and that every status the model
/// reports on inputs built to reach each of its outcomes is in its list, and that the list names nothing those inputs
/// fail to reach except <see cref="FeatureFitStatus.NotConverged"/>, which no small input reaches reliably.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class ReportedStatusesTests
{
    private static readonly FeatureFitStatus[] Base =
        { FeatureFitStatus.Fitted, FeatureFitStatus.TooFewObservations, FeatureFitStatus.RankDeficient };

    private static IEnumerable<TestCaseData> Lists()
    {
        yield return new TestCaseData(LinearModel.ReportedStatuses, Base).SetName("LinearModel");
        yield return new TestCaseData(LogisticRegression.ReportedStatuses,
            Base.Append(FeatureFitStatus.Separated).Append(FeatureFitStatus.NotConverged).ToArray()).SetName("LogisticRegression");
        yield return new TestCaseData(MixedModel.ReportedStatuses,
            Base.Append(FeatureFitStatus.NotConverged).Append(FeatureFitStatus.TooFewGroups)
                .Append(FeatureFitStatus.VarianceRatioAtLimit).ToArray()).SetName("MixedModel");
        yield return new TestCaseData(LogitNormalMixedModel.ReportedStatuses,
            Base.Append(FeatureFitStatus.Separated).Append(FeatureFitStatus.NotConverged)
                .Append(FeatureFitStatus.TooFewGroups).ToArray()).SetName("LogitNormalMixedModel");
        yield return new TestCaseData(NestedMixedModel.ReportedStatuses,
            Base.Append(FeatureFitStatus.NotConverged).Append(FeatureFitStatus.TooFewGroups).ToArray()).SetName("NestedMixedModel");
    }

    [TestCaseSource(nameof(Lists))]
    public void ListsItsStatusesInEnumOrderAndCannotBeChanged(IReadOnlyList<FeatureFitStatus> reported, FeatureFitStatus[] expected)
    {
        Assert.That(reported, Is.EqualTo(expected));
        Assert.That(reported, Is.Ordered.And.Unique);
        Assert.That(reported.All(s => Enum.IsDefined(s)));
        Assert.Throws<NotSupportedException>(() => ((IList<FeatureFitStatus>)reported).Add(FeatureFitStatus.Fitted));
        Assert.Throws<NotSupportedException>(() => ((IList<FeatureFitStatus>)reported)[0] = FeatureFitStatus.NotConverged);
    }

    [Test]
    public void LinearModelReportsOnlyAndAllOfItsList()
    {
        // Columns: intercept, age, male. Samples 0-3 are female, 4-7 male.
        var design = new double[8, 3];
        for (int s = 0; s < 8; s++) { design[s, 0] = 1; design[s, 1] = 30 + 5 * s; design[s, 2] = s >= 4 ? 1 : 0; }
        var y = new double[3, 8];
        for (int s = 0; s < 8; s++)
        {
            y[0, s] = s < 3 ? s : double.NaN;                  // 3 observations, 3 coefficients
            y[1, s] = s >= 4 ? 0.5 * s + s % 2 * 0.1 : double.NaN; // all male: sex not estimable
            y[2, s] = 0.1 * s + (s % 3) * 0.2;                 // fitted
        }
        var fit = LinearModel.Fit(y, design);
        AssertReportsWithin(fit.Status, LinearModel.ReportedStatuses, unreached: Array.Empty<FeatureFitStatus>());
    }

    [Test]
    public void LogisticRegressionReportsOnlyAndAllOfItsList()
    {
        var design = new double[8, 2];
        var y = new double[5, 8];
        for (int s = 0; s < 8; s++)
        {
            design[s, 0] = 1; design[s, 1] = s < 4 ? 0 : 1;
            y[0, s] = s % 3 == 0 ? 1 : 0;                       // fitted (not converged in one iteration)
            y[1, s] = s < 4 ? 0 : 1;                            // complete separation
            y[2, s] = s < 6 ? double.NaN : s % 2;               // two observations for two coefficients
            y[3, s] = s < 4 ? s % 2 : double.NaN;               // observed only where group = 0
            y[4, s] = 1;                                        // no non-events
        }
        var statuses = LogisticRegression.Fit(y, design).Status
            .Concat(LogisticRegression.Fit(y, design, maxIterations: 1).Status).ToList();
        AssertReportsWithin(statuses, LogisticRegression.ReportedStatuses, unreached: Array.Empty<FeatureFitStatus>());
    }

    [Test]
    public void MixedModelReportsOnlyItsList()
    {
        var statuses = new List<FeatureFitStatus>();
        foreach (var (y, design, groups) in MixedInputs())
            statuses.AddRange(MixedModel.Fit(y, design, groups).Status);
        statuses.AddRange(MixedModel.Fit(AtLimit(out var design1, out var groups1), design1, groups1).Status);
        AssertReportsWithin(statuses, MixedModel.ReportedStatuses, unreached: new[] { FeatureFitStatus.NotConverged });
    }

    [Test]
    public void LogitNormalMixedModelReportsOnlyItsList()
    {
        var statuses = new List<FeatureFitStatus>();
        foreach (var (y, design, groups) in MixedInputs())
            statuses.AddRange(LogitNormalMixedModel.Fit(y, new BinomialCount?[y.GetLength(0), y.GetLength(1)], design, groups,
                DetectionLink.Identity).Status);

        // Animals with trait 1 show every trial succeed, animals with trait 0 none: the trait runs to infinity.
        int g = 6, k = 3, n = g * k;
        var design2 = new double[n, 2];
        var groups2 = new string[n];
        var y2 = new double[1, n];
        var counts = new BinomialCount?[1, n];
        for (int s = 0; s < n; s++)
        {
            design2[s, 0] = 1; design2[s, 1] = s / k % 2; groups2[s] = $"a{s / k}";
            y2[0, s] = double.NaN;
            counts[0, s] = new BinomialCount(s / k % 2 == 1 ? 4 : 0, 4);
        }
        statuses.AddRange(LogitNormalMixedModel.Fit(y2, counts, design2, groups2, DetectionLink.Identity,
            residualVariance: new[] { 0.25 }).Status);
        AssertReportsWithin(statuses, LogitNormalMixedModel.ReportedStatuses, unreached: new[] { FeatureFitStatus.NotConverged });
    }

    [Test]
    public void NestedMixedModelReportsOnlyItsList()
    {
        var problems = new[]
        {
            Nested(peptides: 2),                                         // fitted
            Nested(peptides: 1),                                         // one value per sample: too few groups
            Nested(peptides: 2, observed: i => i % 8 == 0),              // 2 observations, 3 coefficients
            Nested(peptides: 2, observed: i => i % 2 == 0),              // x is 0 wherever observed
        };
        var fit = NestedMixedModel.Fit(problems, new[] { "intercept", "x" });
        AssertReportsWithin(fit.Status, NestedMixedModel.ReportedStatuses, unreached: new[] { FeatureFitStatus.NotConverged });
    }

    /// <summary>
    /// Every reported status is in the list, and every listed status but <paramref name="unreached"/> was reported, so
    /// the list neither misses nor invents an outcome these inputs reach.
    /// </summary>
    private static void AssertReportsWithin(IEnumerable<FeatureFitStatus> statuses, IReadOnlyList<FeatureFitStatus> reported,
        IReadOnlyList<FeatureFitStatus> unreached)
    {
        var seen = statuses.Distinct().OrderBy(s => s).ToArray();
        Assert.That(reported, Is.SupersetOf(seen), "a status the model reported is missing from its list");
        Assert.That(seen, Is.EquivalentTo(reported.Except(unreached)), "the inputs should reach every other listed status");
    }

    /// <summary>One group; singleton groups; one observation for one coefficient; a covariate constant where observed; fitted.</summary>
    private static IEnumerable<(double[,] Y, double[,] Design, string[] Groups)> MixedInputs()
    {
        var design = new double[4, 1] { { 1 }, { 1 }, { 1 }, { 1 } };
        var y = new double[1, 4] { { 1, 2, 3, 4 } };
        yield return (y, design, new[] { "a", "a", "a", "a" });
        yield return (y, design, new[] { "a", "b", "c", "d" });

        var design2 = new double[6, 2] { { 1, 0 }, { 1, 0 }, { 1, 0 }, { 1, 0 }, { 1, 1 }, { 1, 1 } };
        var y2 = new double[2, 6]
        {
            { 1, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN },
            { 1, 2, 3, 4, double.NaN, double.NaN },
        };
        yield return (y2, design2, new[] { "a", "a", "b", "b", "c", "c" });

        var rng = new Random(20261027);
        int groups = 6, k = 4, n = groups * k;
        var design3 = new double[n, 2];
        var y3 = new double[1, n];
        var groups3 = new string[n];
        for (int s = 0; s < n; s++)
        {
            int g = s / k;
            design3[s, 0] = 1; design3[s, 1] = s % 2; groups3[s] = $"g{g}";
            y3[0, s] = 2 + 0.5 * (s % 2) + 0.8 * g % 3 + Normal.Sample(rng, 0, 0.5);
        }
        yield return (y3, design3, groups3);
    }

    /// <summary>Two near-identical replicates per group: τ²/σ² ≈ 1e10, far beyond the search's e¹².</summary>
    private static double[,] AtLimit(out double[,] design, out string[] groups)
    {
        int g = 8, k = 2, n = g * k;
        var rng = new Random(3);
        design = new double[n, 1];
        groups = new string[n];
        var y = new double[1, n];
        for (int s = 0; s < n; s++)
        {
            design[s, 0] = 1; groups[s] = $"g{s / k}";
            y[0, s] = 10 * (s / k) + Normal.Sample(rng, 0, 1e-4);
        }
        return y;
    }

    /// <summary>Peptides in 8 samples, x = 0/1 by sample, deterministic values; NaN where <paramref name="observed"/> is false.</summary>
    private static NestedMixedModelProblem Nested(int peptides, Func<int, bool>? observed = null)
    {
        int samples = 8, n = samples * peptides;
        var y = new double[n];
        var common = new double[n, 2];
        var nuisance = peptides > 1 ? new double[n, peptides - 1] : null;
        var sample = new string[n];
        double[] offset = { 0.3, -0.2, 0.1, -0.35, 0.25, -0.05, 0.15, -0.2 };
        for (int p = 0; p < peptides; p++)
            for (int j = 0; j < samples; j++)
            {
                int i = p * samples + j;
                common[i, 0] = 1;
                common[i, 1] = j % 2;
                if (p > 0) nuisance![i, p - 1] = 1;
                sample[i] = $"s{j}";
                y[i] = observed is null || observed(i) ? 20 + 0.5 * (j % 2) + 0.4 * p + offset[j] + 0.01 * ((i * 7) % 5) : double.NaN;
            }
        return new NestedMixedModelProblem(y, common, nuisance, sample, null);
    }
}
