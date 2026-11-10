using System;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using Quantification.Differential;
using StatisticalModels;

namespace Test.Quantification.Differential;

/// <summary>
/// A feature's own design (STAT1 M5b ruling 1): reduced to what its observations estimate, every contrast judged by what
/// is estimable (with the nuisance columns), never by which column survived.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class FeatureDesignTests
{
    private static readonly string[] GroupTerms = { "intercept", "group", "group" };

    /// <summary>Stratum rows for groups A (reference), B, C: [intercept, group=B, group=C].</summary>
    private static double[] Row(char group) => new[] { 1.0, group == 'B' ? 1 : 0, group == 'C' ? 1 : 0 };

    private static double[,] Rows(params char[] groups)
    {
        var m = new double[groups.Length, 3];
        for (int i = 0; i < groups.Length; i++)
            for (int j = 0; j < 3; j++) m[i, j] = Row(groups[i])[j];
        return m;
    }

    [Test]
    public void AFeatureSeenInEveryGroupKeepsEveryColumn()
    {
        var design = FeatureDesign.Build(Rows('A', 'A', 'B', 'B', 'C', 'C'), null, GroupTerms);

        Assert.That(design.IsFull, Is.True);
        Assert.That(design.CommonKey, Is.EqualTo("0,1,2"));
        Assert.That(design.Estimability(new[] { 0.0, 0, 1 }), Is.EqualTo(AnalysisDesign.Estimable));
        Assert.That(design.VarianceDfRows(new[] { 1, 2 }), Is.EqualTo(new[] { new[] { 0.0, 1, 0 }, new[] { 0.0, 0, 1 } }));
    }

    /// <summary>
    /// The absent-level probe's trap (results/stat1/m5b): with the REFERENCE group A absent, no column is zero; B + C equals
    /// the intercept, C is left out, and the surviving B column means B against C. C against B must come out as C minus B
    /// (here exactly 1), and neither "B against A" nor "C against A" may be estimable.
    /// </summary>
    [Test]
    public void WithTheReferenceGroupAbsentOnlyDifferencesOfThePresentGroupsAreEstimable()
    {
        var design = FeatureDesign.Build(Rows('B', 'B', 'C', 'C'), null, GroupTerms);
        double[] cVsA = { 0, 0, 1 }, bVsA = { 0, 1, 0 }, cVsB = { 0, -1, 1 };

        Assert.That(design.KeptCommon, Is.EqualTo(new[] { 0, 1 }));
        Assert.That(design.Estimability(cVsA), Is.Not.EqualTo(AnalysisDesign.Estimable));
        Assert.That(design.Estimability(bVsA), Is.Not.EqualTo(AnalysisDesign.Estimable));
        Assert.That(design.Estimability(cVsB), Is.EqualTo(AnalysisDesign.Estimable));

        // B samples at 21, C samples at 22: C minus B is 1.
        var fit = LinearModel.Fit(new double[,] { { 21, 21, 22, 22 } }, new double[,] { { 1, 1 }, { 1, 1 }, { 1, 0 }, { 1, 0 } },
            new[] { "intercept", "group=B" });
        var weights = design.KeptWeights(cVsB);
        Assert.That(fit.ContrastEstimate(0, weights), Is.EqualTo(1).Within(1e-12));

        var rows = design.VarianceDfRows(new[] { 1, 2 });
        Assert.That(rows, Has.Count.EqualTo(1), "one estimable direction: C against B");
        Assert.That(rows[0][0], Is.EqualTo(0).Within(1e-12));
        Assert.That(Math.Abs(rows[0][1]), Is.GreaterThan(0.1));
    }

    [Test]
    public void APeptideSeenInOneGroupOnlyConfoundsThatGroupWithThePeptide()
    {
        // pep1 only in A, pep2 only in B: the pep2 column equals the B column.
        var common = Rows('A', 'A', 'B', 'B');
        var design = FeatureDesign.Build(new[,] { { 1.0, 0 }, { 1, 0 }, { 1, 1 }, { 1, 1 } }, new double[,] { { 0 }, { 0 }, { 1 }, { 1 } },
            new[] { "intercept", "group" });

        Assert.That(design.KeptNuisance, Is.Empty, "the peptide column is the later of the two aliased columns");
        Assert.That(design.Estimability(new[] { 0.0, 1 }), Is.EqualTo(AnalysisDesign.NotEstimableConfounded));

        // pep2 also in A: separable.
        var shared = FeatureDesign.Build(new[,] { { 1.0, 0 }, { 1, 0 }, { 1, 1 }, { 1, 1 } }, new double[,] { { 0 }, { 1 }, { 1 }, { 1 } },
            new[] { "intercept", "group" });
        Assert.That(shared.IsFull, Is.True);
        Assert.That(shared.Estimability(new[] { 0.0, 1 }), Is.EqualTo(AnalysisDesign.Estimable));
        Assert.That(common.GetLength(0), Is.EqualTo(4));
    }

    [Test]
    public void ABatchTheFeatureIsMissingFromIsDroppedAndTheGroupsStillCompare()
    {
        // [intercept, group=B, batch=b2]; the feature is seen in batch b1 only.
        var design = FeatureDesign.Build(new[,] { { 1.0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 1, 1, 0 } }, null,
            new[] { "intercept", "group", "batch" });

        Assert.That(design.KeptCommon, Is.EqualTo(new[] { 0, 1 }));
        Assert.That(design.Estimability(new[] { 0.0, 1, 0 }), Is.EqualTo(AnalysisDesign.Estimable));
        Assert.That(design.KeptWeights(new[] { 0.0, 1, 0 }), Is.EqualTo(new[] { 0.0, 1 }));
        Assert.That(design.VarianceDfRows(new[] { 1 }), Is.EqualTo(new[] { new[] { 0.0, 1 } }));
    }

    /// <summary>
    /// A model mean at the stratum's proportions (ruling 3) needs every column it puts weight on: with batch b2 never seen,
    /// the group contrast is estimable but a mean that averages over both batches is not; over batch b1 alone it is.
    /// </summary>
    [Test]
    public void AMeanOverABatchTheFeatureNeverSawCannotBeEstimated()
    {
        var design = FeatureDesign.Build(new[,] { { 1.0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 1, 1, 0 } }, null,
            new[] { "intercept", "group", "batch" });

        Assert.That(design.IsEstimable(new[] { 0.0, 1, 0 }, Array.Empty<double>()), Is.True);
        Assert.That(design.IsEstimable(new[] { 1.0, 1, 0.5 }, Array.Empty<double>()), Is.False);
        Assert.That(design.IsEstimable(new[] { 1.0, 1, 0 }, Array.Empty<double>()), Is.True);
    }

    /// <summary>
    /// A mean averages the peptides' offsets with equal weight. pep2 seen only in group B is aliased with B: B's mean over
    /// both peptides cannot be estimated, nor B against A; with pep2 also seen in A, both can.
    /// </summary>
    [Test]
    public void AMeanOverPeptidesNeedsEachPeptideSeparableFromTheGroups()
    {
        var common = new[,] { { 1.0, 0 }, { 1, 0 }, { 1, 1 }, { 1, 1 } };
        var confounded = FeatureDesign.Build(common, new double[,] { { 0 }, { 0 }, { 1 }, { 1 } }, new[] { "intercept", "group" });
        var shared = FeatureDesign.Build(common, new double[,] { { 0 }, { 1 }, { 1 }, { 1 } }, new[] { "intercept", "group" });

        Assert.That(confounded.IsEstimable(new[] { 1.0, 1 }, new[] { 0.5 }), Is.False);
        // What pep2 measures in B (its whole weight on pep2) is B's own observations, so it is estimable; the same
        // without the peptide's weight is not.
        Assert.That(confounded.IsEstimable(new[] { 1.0, 1 }, new[] { 1.0 }), Is.True);
        Assert.That(confounded.IsEstimable(new[] { 1.0, 1 }, new[] { 0.0 }), Is.False);
        Assert.That(shared.IsEstimable(new[] { 1.0, 1 }, new[] { 0.5 }), Is.True);
        Assert.That(shared.IsEstimable(new[] { 1.0, 0 }, new[] { 0.5 }), Is.True);
        Assert.Throws<ArgumentException>(() => shared.IsEstimable(new[] { 1.0, 0 }, Array.Empty<double>()));
    }

    [Test]
    public void ATwoLevelFactorWithItsReferenceAbsentHasNoEstimableDirection()
    {
        var design = FeatureDesign.Build(new[,] { { 1.0, 1 }, { 1, 1 }, { 1, 1 } }, null, new[] { "intercept", "group" });

        Assert.That(design.Estimability(new[] { 0.0, 1 }), Is.Not.EqualTo(AnalysisDesign.Estimable));
        Assert.That(design.VarianceDfRows(new[] { 1 }), Is.Empty);
    }

    [Test]
    public void AContrastOfTheWrongLengthIsRefused()
    {
        var design = FeatureDesign.Build(Rows('A', 'B', 'C', 'A'), null, GroupTerms);
        Assert.Throws<ArgumentException>(() => design.Estimability(new[] { 0.0, 1 }));
        Assert.Throws<ArgumentException>(() => FeatureDesign.Build(Rows('A', 'B'), null, new[] { "intercept" }));
    }
}
