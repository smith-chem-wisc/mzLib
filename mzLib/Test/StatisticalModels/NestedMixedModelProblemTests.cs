using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using StatisticalModels;

namespace Test.StatisticalModels;

/// <summary>
/// <see cref="NestedMixedModelProblem.FromNuisanceLabels"/>: peptide labels coded as the lme4 reference fits code them,
/// so a caller (a language binding, a pipeline) passes labels instead of building treatment-coded columns by hand.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class NestedMixedModelProblemTests
{
    private static readonly string[] Names = { "intercept", "x" };

    [TestCase("nested")]
    [TestCase("sample_only")]
    public void LabelsGiveTheSameFitAsTheReferenceCoding(string scenario)
    {
        // The reference coding (LmmData) drops missing values, then codes the observed peptides against the first in
        // ordinal order, as R's lme4 fixtures did. Built from the raw rows, missing values included, the labels must
        // give the same problems, so the same fit to the last bit.
        var reference = LmmData.Read(scenario);
        var samples = LmmData.Rows($"lmm_{scenario}_samples.tsv").ToDictionary(s => s["sample"]);
        var problems = LmmData.Rows($"lmm_{scenario}_values.tsv").GroupBy(v => v["protein"]).Select(protein =>
        {
            var rows = protein.ToList();
            var common = new double[rows.Count, reference.CommonNames.Count];
            for (int i = 0; i < rows.Count; i++)
                for (int j = 0; j < reference.CommonNames.Count; j++)
                    common[i, j] = LmmData.Parse(samples[rows[i]["sample"]][reference.CommonNames[j]]);
            return NestedMixedModelProblem.FromNuisanceLabels(
                rows.Select(r => LmmData.Parse(r["y"])).ToArray(), common, rows.Select(r => r["peptide"]).ToArray(),
                rows.Select(r => r["sample"]).ToArray(),
                rows.Select(r => samples[r["sample"]]["individual"] is "NA" or "" ? null : samples[r["sample"]]["individual"]).ToArray());
        }).ToList();

        var expected = NestedMixedModel.Fit(reference.Problems, reference.CommonNames);
        var actual = NestedMixedModel.Fit(problems, reference.CommonNames);

        Assert.That(actual.Status, Is.EqualTo(expected.Status));
        Assert.That(actual.Structure, Is.EqualTo(expected.Structure));
        Assert.That(actual.Observed, Is.EqualTo(expected.Observed));
        Assert.That(actual.RemlCriterion, Is.EqualTo(expected.RemlCriterion));
        for (int f = 0; f < actual.FeatureCount; f++)
            for (int c = 0; c < reference.CommonNames.Count; c++)
                Assert.That(actual.Coefficient(f, c), Is.EqualTo(expected.Coefficient(f, c)), $"protein {f}, coefficient {c}");
    }

    [Test]
    public void APeptideSeenOnlyInMissingValuesGetsNoColumn()
    {
        // pep2 and pep3 are observed in 8 samples; pep1 only as missing values. Coding pep1 too would add an all-zero
        // column and make the protein rank deficient; the labels drop it, and pep2 becomes the reference.
        var (y, common, labels, sample) = TwoPeptidesAndAMissingOne();

        var fromLabels = NestedMixedModelProblem.FromNuisanceLabels(y, common, labels, sample, null);
        Assert.That(fromLabels.NuisanceLevels, Is.EqualTo(new[] { "pep2", "pep3" }));
        Assert.That(NestedMixedModel.Fit(new[] { fromLabels }, Names).Status[0], Is.EqualTo(FeatureFitStatus.Fitted));

        // Coding every label instead makes pep1 the reference, so pep2 and pep3 each get a column. Over the observed
        // rows those two columns sum to the intercept, so the design is collinear.
        var observed = Enumerable.Range(0, y.Length).Where(i => double.IsFinite(y[i])).ToArray();
        var everyLabel = new double[observed.Length, 2];
        var x = new double[observed.Length, 2];
        for (int r = 0; r < observed.Length; r++)
        {
            int i = observed[r];
            x[r, 0] = common[i, 0];
            x[r, 1] = common[i, 1];
            everyLabel[r, labels[i] == "pep2" ? 0 : 1] = 1;
        }
        var naive = new NestedMixedModelProblem(observed.Select(i => y[i]).ToArray(), x, everyLabel,
            observed.Select(i => sample[i]).ToArray(), null);
        Assert.That(NestedMixedModel.Fit(new[] { naive }, Names).Status[0], Is.EqualTo(FeatureFitStatus.RankDeficient));
    }

    [Test]
    public void OneLabelGivesNoNuisanceColumnAndNoIndividualMeansNone()
    {
        // One peptide measured twice in each of 8 samples (e.g. two charge states).
        double[] offset = { 0.3, -0.2, 0.1, -0.35, 0.25, -0.05, 0.15, -0.2 };
        double[] noise = { 0.02, -0.05, 0.04, 0.01, -0.03, 0.05, -0.01, 0.03, -0.04, 0.02, 0.0, -0.02, 0.03, -0.01, 0.04, -0.03 };
        var y = Enumerable.Range(0, 16).Select(i => 20 + 0.5 * (i % 8 % 2) + offset[i % 8] + noise[i]).ToArray();
        var common = new double[16, 2];
        for (int i = 0; i < 16; i++) { common[i, 0] = 1; common[i, 1] = i % 8 % 2; }
        var samples = Enumerable.Range(0, 16).Select(i => $"s{i % 8}").ToArray();

        var problem = NestedMixedModelProblem.FromNuisanceLabels(y, common, Enumerable.Repeat("pep", 16).ToArray(), samples,
            new string?[16]);

        Assert.That(problem.NuisanceLevels, Is.EqualTo(new[] { "pep" }));
        var fit = NestedMixedModel.Fit(new[] { problem }, Names);
        Assert.That(fit.Status[0], Is.EqualTo(FeatureFitStatus.Fitted));
        Assert.That(fit.Structure[0], Is.EqualTo(RandomStructure.Sample), "no individual given on any kept row");
    }

    [Test]
    public void LevelsAreOrdinalNotFirstSeen()
    {
        var (y, common, labels, sample) = TwoPeptidesAndAMissingOne();
        var reversed = labels.Select(l => l == "pep2" ? "b" : l == "pep3" ? "B" : l).ToArray();

        var problem = NestedMixedModelProblem.FromNuisanceLabels(y, common, reversed, sample, null);

        Assert.That(problem.NuisanceLevels, Is.EqualTo(new[] { "B", "b" }), "ordinal: upper case sorts first");
    }

    [Test]
    public void ProblemsBuiltDirectlyHaveNoLevels()
    {
        var problem = new NestedMixedModelProblem(new[] { 1.0 }, new double[1, 1] { { 1 } }, null, new[] { "s" }, null);
        Assert.That(problem.NuisanceLevels, Is.Null);
    }

    [Test]
    public void RefusesMismatchedOrEmptyLabels()
    {
        var common = new double[2, 1] { { 1 }, { 1 } };
        Assert.Throws<ArgumentException>(() =>
            NestedMixedModelProblem.FromNuisanceLabels(new[] { 1.0, 2.0 }, common, new[] { "a" }, new[] { "s1", "s2" }, null));
        Assert.Throws<ArgumentException>(() =>
            NestedMixedModelProblem.FromNuisanceLabels(new[] { 1.0, 2.0 }, common, new[] { "a", "" }, new[] { "s1", "s2" }, null));
        Assert.Throws<ArgumentException>(() =>
            NestedMixedModelProblem.FromNuisanceLabels(new[] { 1.0, 2.0 }, common, new[] { "a", "b" }, new[] { "s1" }, null));
        Assert.Throws<ArgumentException>(() =>
            NestedMixedModelProblem.FromNuisanceLabels(new[] { 1.0, 2.0 }, new double[1, 1], new[] { "a", "b" }, new[] { "s1", "s2" }, null));
        Assert.Throws<ArgumentException>(() =>
            NestedMixedModelProblem.FromNuisanceLabels(new[] { 1.0, 2.0 }, common, new[] { "a", "b" }, new[] { "s1", "s2" }, new string?[] { "i" }));
        Assert.Throws<ArgumentNullException>(() =>
            NestedMixedModelProblem.FromNuisanceLabels(new[] { 1.0 }, new double[1, 1], null!, new[] { "s1" }, null));
    }

    /// <summary>Two observed peptides in 8 samples, x = 0/1 by sample, and a third peptide that is missing everywhere.</summary>
    private static (double[] y, double[,] common, string[] labels, string[] sample) TwoPeptidesAndAMissingOne()
    {
        const int s = 8;
        double[] offset = { 0.3, -0.2, 0.1, -0.35, 0.25, -0.05, 0.15, -0.2 };
        double[] noise = { 0.02, -0.05, 0.04, 0.01, -0.03, 0.05, -0.01, 0.03, -0.04, 0.02, 0.0, -0.02, 0.03, -0.01, 0.04, -0.03 };
        string[] peptides = { "pep1", "pep2", "pep3" };
        int n = s * peptides.Length;
        var y = new double[n];
        var common = new double[n, 2];
        var labels = new string[n];
        var sample = new string[n];
        for (int k = 0; k < peptides.Length; k++)
            for (int j = 0; j < s; j++)
            {
                int i = k * s + j;
                common[i, 0] = 1;
                common[i, 1] = j % 2;
                labels[i] = peptides[k];
                sample[i] = $"s{j}";
                y[i] = k == 0 ? double.NaN : 20 + 0.5 * (j % 2) + 0.4 * k + offset[j] + noise[i % noise.Length];
            }
        return (y, common, labels, sample);
    }
}
