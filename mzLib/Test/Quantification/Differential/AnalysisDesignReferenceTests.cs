using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Quantification.Differential;

namespace Test.Quantification.Differential;

/// <summary>
/// <see cref="AnalysisDesign"/>'s design matrix equals R's <c>stats::model.matrix</c> for the same design, value for value,
/// on fourteen designs (balanced, unbalanced, a lost sample, paired, within and between individuals, batch, an age slope
/// with and without a stated centre, two factors, an interaction, a non-alphabetical reference, a confounded design, no
/// reference, and a factor with batch and a covariate). The fixtures were generated once by
/// <c>ReferenceData/make_design_fixtures.R</c> (base R only); R is never run by these tests.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class AnalysisDesignReferenceTests
{
    private const double Tolerance = 1e-12;

    private static readonly CovariateSpec AgePerDecadeAt50 = new("age", "years", Scale: 10, ScaledUnit: "decade", Centre: 50);

    private static readonly Dictionary<string, DesignSpec> Specs = new()
    {
        ["balanced"] = Factors(("condition", "A")),
        ["unbalanced"] = Factors(("condition", "A")),
        ["lost_sample"] = Factors(("condition", "ctrl")),
        ["paired"] = Factors(("time", "pre")),
        ["within_between"] = Factors(("genotype", "wt"), ("time", "d0")),
        ["batch"] = Factors(("condition", "A")),
        ["age_slope"] = new DesignSpec { Covariates = new[] { AgePerDecadeAt50 } },
        ["age_default_centre"] = Factors(("sex", "female")) with
        {
            Covariates = new[] { new CovariateSpec("age", "years", Scale: 10, ScaledUnit: "decade") }
        },
        ["two_factors"] = Factors(("age_group", "young"), ("treatment", "vehicle")),
        ["interaction"] = Factors(("age_group", "young"), ("treatment", "vehicle")) with
        {
            Interactions = new[] { ("age_group", "treatment") }
        },
        ["nonalpha_reference"] = Factors(("age_group", "young")),
        ["confounded"] = Factors(("condition", "A")),
        ["no_reference"] = Factors(("condition", null)),
        ["factor_batch_covariate"] = Factors(("condition", "ctrl")) with { Covariates = new[] { AgePerDecadeAt50 } },
    };

    private static IEnumerable<string> DesignNames => Specs.Keys;

    [TestCaseSource(nameof(DesignNames))]
    public void DesignMatrixEqualsModelMatrix(string name)
    {
        var (samples, _) = ReadSamples(name);
        var (rColumns, expected) = ReadMatrix(name);
        var design = AnalysisDesign.Create(samples, Specs[name]);

        Assert.That(design.Samples.Select(s => s.SampleId), Is.EqualTo(samples.Select(s => s.SampleId)),
            "rows are in ordinal sample order, as the fixture is");
        var actual = design.Matrix;
        Assert.That(actual.GetLength(0), Is.EqualTo(expected.Length), "rows");
        Assert.That(actual.GetLength(1), Is.EqualTo(rColumns.Length),
            $"columns: R has {string.Join(", ", rColumns)}; ours {string.Join(", ", design.ColumnNames)}");
        for (int i = 0; i < expected.Length; i++)
            for (int j = 0; j < rColumns.Length; j++)
                Assert.That(actual[i, j], Is.EqualTo(expected[i][j]).Within(Tolerance),
                    $"{name}: sample {samples[i].SampleId}, column {design.ColumnNames[j]} (R: {rColumns[j]})");
    }

    [TestCaseSource(nameof(DesignNames))]
    public void RankEqualsR(string name)
    {
        var (samples, _) = ReadSamples(name);
        var design = AnalysisDesign.Create(samples, Specs[name]);
        int expectedRank = name == "confounded" ? 2 : design.ColumnNames.Count;
        Assert.That(design.Rank, Is.EqualTo(expectedRank), "R's qr()$rank, printed when the fixtures were made");
    }

    private static DesignSpec Factors(params (string Name, string? Reference)[] factors) =>
        new() { Factors = factors.Select(f => new FactorSpec(f.Name, f.Reference)).ToArray() };

    internal static (List<AnalysisSample> Samples, string[] Header) ReadSamples(string name)
    {
        var lines = ReadTsv($"design_{name}_samples.tsv");
        var header = lines[0];
        var samples = new List<AnalysisSample>();
        foreach (var row in lines.Skip(1))
        {
            var factors = new Dictionary<string, string>();
            var covariates = new Dictionary<string, double>();
            string? individual = null, batch = null;
            for (int c = 1; c < header.Length; c++)
            {
                switch (header[c])
                {
                    case "individual": individual = row[c]; break;
                    case "batch": batch = row[c]; break;
                    case "age": covariates["age"] = double.Parse(row[c], CultureInfo.InvariantCulture); break;
                    default: factors[header[c]] = row[c]; break;
                }
            }
            samples.Add(new AnalysisSample(row[0])
            {
                Factors = factors, Covariates = covariates, Individual = individual, Batch = batch
            });
        }
        return (samples, header);
    }

    private static (string[] Columns, double[][] Rows) ReadMatrix(string name)
    {
        var lines = ReadTsv($"design_{name}_matrix.tsv");
        return (lines[0], lines.Skip(1)
            .Select(r => r.Select(v => double.Parse(v, CultureInfo.InvariantCulture)).ToArray()).ToArray());
    }

    private static List<string[]> ReadTsv(string file) =>
        File.ReadAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, "Quantification", "Differential",
                "ReferenceData", file))
            .Where(l => l.Length > 0).Select(l => l.Split('\t')).ToList();
}
