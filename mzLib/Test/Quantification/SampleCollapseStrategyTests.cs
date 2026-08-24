using MassSpectrometry;
using NUnit.Framework;
using Quantification;
using Quantification.Strategies;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Test.Quantification;

/// <summary>
/// Covers <see cref="SampleCollapseStrategy"/>: the parameterized collapse that reduces one named
/// sample dimension at a time, in contrast to <see cref="SumCollapse"/> and <see cref="MeanCollapse"/>
/// which group by Condition + BiologicalReplicate and collapse everything else at once.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class SampleCollapseStrategyTests
{
    #region Fixture

    /// <summary>A row key that is just a name, so the tests are about columns only.</summary>
    private sealed record Row(string Name) : IEquatable<Row>;

    private static SpectraFileInfo Lfq(string condition, int biorep, int techrep, int fraction) =>
        new($"{condition}_b{biorep}_t{techrep}_f{fraction}.raw", condition, biorep, techrep, fraction);

    private static QuantMatrix<Row> Matrix(IList<ISampleInfo> columns, params double[] values)
    {
        var row = new Row("r1");
        var matrix = new QuantMatrix<Row>(new List<Row> { row }, columns, null);
        matrix.SetRow(row, values);
        return matrix;
    }

    private static double[] SingleRow<T>(QuantMatrix<T> matrix) where T : IEquatable<T> =>
        matrix.GetRow(0);

    #endregion

    /// <summary>
    /// Two fractions of one sample collapse to one column; the technical replicate is preserved.
    /// </summary>
    [Test]
    public void CollapseFraction_SumsAcrossFractionsOnly()
    {
        var columns = new List<ISampleInfo>
        {
            Lfq("Control", 1, 1, 1),
            Lfq("Control", 1, 1, 2),
            Lfq("Control", 1, 2, 1),
            Lfq("Control", 1, 2, 2)
        };
        var matrix = Matrix(columns, 10, 20, 30, 40);

        var collapsed = new SampleCollapseStrategy(CollapseDimension.Fraction, AggregationType.Sum)
            .CollapseSamples(matrix);

        Assert.That(collapsed.ColumnCount, Is.EqualTo(2), "two technical replicates survive");
        Assert.That(SingleRow(collapsed), Is.EqualTo(new[] { 30.0, 70.0 }));

        var survivors = collapsed.ColumnKeys.ToList();
        Assert.That(survivors.Select(s => s.TechnicalReplicate), Is.EqualTo(new[] { 1, 2 }));
        Assert.That(survivors.Select(s => s.Fraction), Is.EqualTo(new[] { 0, 0 }), "fraction is zeroed");
    }

    /// <summary>
    /// The complementary case: collapsing technical replicates must leave the fractions alone.
    /// This is the distinction master's SumCollapse/MeanCollapse cannot express.
    /// </summary>
    [Test]
    public void CollapseTechnicalReplicate_LeavesFractionsIntact()
    {
        var columns = new List<ISampleInfo>
        {
            Lfq("Control", 1, 1, 1),
            Lfq("Control", 1, 2, 1),
            Lfq("Control", 1, 1, 2),
            Lfq("Control", 1, 2, 2)
        };
        var matrix = Matrix(columns, 10, 30, 20, 40);

        var collapsed = new SampleCollapseStrategy(CollapseDimension.TechnicalReplicate, AggregationType.Median)
            .CollapseSamples(matrix);

        Assert.That(collapsed.ColumnCount, Is.EqualTo(2), "two fractions survive");
        Assert.That(SingleRow(collapsed), Is.EqualTo(new[] { 20.0, 30.0 }));
        Assert.That(collapsed.ColumnKeys.Select(s => s.Fraction), Is.EqualTo(new[] { 1, 2 }));
        Assert.That(collapsed.ColumnKeys.Select(s => s.TechnicalReplicate), Is.EqualTo(new[] { 0, 0 }));
    }

    [Test]
    public void CollapseBiologicalReplicate_GroupsByConditionTechAndFraction()
    {
        var columns = new List<ISampleInfo>
        {
            Lfq("Control", 1, 1, 1),
            Lfq("Control", 2, 1, 1),
            Lfq("Treated", 1, 1, 1),
            Lfq("Treated", 2, 1, 1)
        };
        var matrix = Matrix(columns, 10, 30, 100, 200);

        var collapsed = new SampleCollapseStrategy(CollapseDimension.BiologicalReplicate, AggregationType.Average)
            .CollapseSamples(matrix);

        Assert.That(collapsed.ColumnCount, Is.EqualTo(2));
        Assert.That(SingleRow(collapsed), Is.EqualTo(new[] { 20.0, 150.0 }));
        Assert.That(collapsed.ColumnKeys.Select(s => s.Condition), Is.EqualTo(new[] { "Control", "Treated" }));
        Assert.That(collapsed.ColumnKeys.Select(s => s.BiologicalReplicate), Is.EqualTo(new[] { 0, 0 }));
    }

    [TestCase(AggregationType.Sum, 60.0)]
    [TestCase(AggregationType.Average, 20.0)]
    [TestCase(AggregationType.Median, 20.0)]
    [TestCase(AggregationType.Max, 30.0)]
    public void Aggregation_OddCount(AggregationType aggregation, double expected)
    {
        var columns = new List<ISampleInfo>
        {
            Lfq("C", 1, 1, 1), Lfq("C", 1, 1, 2), Lfq("C", 1, 1, 3)
        };
        var matrix = Matrix(columns, 10, 20, 30);

        var collapsed = new SampleCollapseStrategy(CollapseDimension.Fraction, aggregation).CollapseSamples(matrix);

        Assert.That(collapsed.ColumnCount, Is.EqualTo(1));
        Assert.That(SingleRow(collapsed)[0], Is.EqualTo(expected).Within(1e-9));
    }

    /// <summary>An even-sized group must average the two middle values, not pick one.</summary>
    [Test]
    public void Median_EvenCount_AveragesTheTwoMiddleValues()
    {
        var columns = new List<ISampleInfo>
        {
            Lfq("C", 1, 1, 1), Lfq("C", 1, 1, 2), Lfq("C", 1, 1, 3), Lfq("C", 1, 1, 4)
        };
        var matrix = Matrix(columns, 40, 10, 30, 20);

        var collapsed = new SampleCollapseStrategy(CollapseDimension.Fraction, AggregationType.Median)
            .CollapseSamples(matrix);

        Assert.That(SingleRow(collapsed)[0], Is.EqualTo(25.0).Within(1e-9));
    }

    /// <summary>
    /// Collapsing an isobaric column must yield an isobaric column. Losing the channel identity would
    /// silently disable ReferenceChannelNormalization for every step downstream.
    /// </summary>
    [Test]
    public void CollapsingIsobaricSamples_PreservesChannelIdentity()
    {
        var columns = new List<ISampleInfo>
        {
            new IsobaricQuantSampleInfo("f1.raw", "Reference", 1, 1, 1, 3, "126", 126.127, true),
            new IsobaricQuantSampleInfo("f2.raw", "Reference", 1, 2, 1, 3, "126", 126.127, true)
        };
        var matrix = Matrix(columns, 100, 300);

        var collapsed = new SampleCollapseStrategy(CollapseDimension.TechnicalReplicate, AggregationType.Sum)
            .CollapseSamples(matrix);

        Assert.That(collapsed.ColumnCount, Is.EqualTo(1));
        var survivor = collapsed.ColumnKeys[0];

        Assert.That(survivor, Is.InstanceOf<IsobaricQuantSampleInfo>());
        var isobaric = (IsobaricQuantSampleInfo)survivor;
        Assert.Multiple(() =>
        {
            Assert.That(isobaric.ChannelLabel, Is.EqualTo("126"));
            Assert.That(isobaric.ReporterIonMz, Is.EqualTo(126.127).Within(1e-9));
            Assert.That(isobaric.PlexId, Is.EqualTo(3));
            Assert.That(isobaric.IsReferenceChannel, Is.True);
            Assert.That(isobaric.TechnicalReplicate, Is.EqualTo(0), "the collapsed dimension is zeroed");
            Assert.That(isobaric.Fraction, Is.EqualTo(1), "other dimensions survive");
        });
        Assert.That(SingleRow(collapsed)[0], Is.EqualTo(400.0).Within(1e-9));
    }

    /// <summary>A label-free column stays label-free.</summary>
    [Test]
    public void CollapsingLabelFreeSamples_YieldsSpectraFileInfo()
    {
        var columns = new List<ISampleInfo> { Lfq("C", 1, 1, 1), Lfq("C", 1, 1, 2) };
        var collapsed = new SampleCollapseStrategy(CollapseDimension.Fraction, AggregationType.Sum)
            .CollapseSamples(Matrix(columns, 1, 2));

        Assert.That(collapsed.ColumnKeys[0], Is.InstanceOf<SpectraFileInfo>());
    }

    /// <summary>
    /// Applying the strategy twice is the standard bottom-up recipe: sum fractions, then take the
    /// median across technical replicates.
    /// </summary>
    [Test]
    public void Composed_SumFractionsThenMedianTechnicalReplicates()
    {
        var columns = new List<ISampleInfo>
        {
            Lfq("C", 1, 1, 1), Lfq("C", 1, 1, 2),   // techrep 1 -> 10 + 20 = 30
            Lfq("C", 1, 2, 1), Lfq("C", 1, 2, 2),   // techrep 2 -> 40 + 60 = 100
            Lfq("C", 1, 3, 1), Lfq("C", 1, 3, 2)    // techrep 3 -> 25 + 25 = 50
        };
        var matrix = Matrix(columns, 10, 20, 40, 60, 25, 25);

        var summed = new SampleCollapseStrategy(CollapseDimension.Fraction, AggregationType.Sum)
            .CollapseSamples(matrix);
        Assert.That(SingleRow(summed), Is.EqualTo(new[] { 30.0, 100.0, 50.0 }));

        var final = new SampleCollapseStrategy(CollapseDimension.TechnicalReplicate, AggregationType.Median)
            .CollapseSamples(summed);

        Assert.That(final.ColumnCount, Is.EqualTo(1));
        Assert.That(SingleRow(final)[0], Is.EqualTo(50.0).Within(1e-9), "median of 30, 100, 50");
    }

    /// <summary>Column order in equals column order out — output must be deterministic.</summary>
    [Test]
    public void CollapsedColumns_FollowFirstAppearanceOrder()
    {
        var columns = new List<ISampleInfo>
        {
            Lfq("Zebra", 1, 1, 1),
            Lfq("Apple", 1, 1, 1),
            Lfq("Zebra", 1, 1, 2)
        };
        var collapsed = new SampleCollapseStrategy(CollapseDimension.Fraction, AggregationType.Sum)
            .CollapseSamples(Matrix(columns, 1, 2, 3));

        Assert.That(collapsed.ColumnKeys.Select(s => s.Condition), Is.EqualTo(new[] { "Zebra", "Apple" }));
        Assert.That(SingleRow(collapsed), Is.EqualTo(new[] { 4.0, 2.0 }));
    }

    /// <summary>Nothing to collapse is not an error; the matrix passes through unchanged in value.</summary>
    [Test]
    public void NoGroupHasMoreThanOneMember_ValuesAreUnchanged()
    {
        var columns = new List<ISampleInfo> { Lfq("A", 1, 1, 1), Lfq("B", 1, 1, 1) };
        var collapsed = new SampleCollapseStrategy(CollapseDimension.Fraction, AggregationType.Sum)
            .CollapseSamples(Matrix(columns, 7, 9));

        Assert.That(collapsed.ColumnCount, Is.EqualTo(2));
        Assert.That(SingleRow(collapsed), Is.EqualTo(new[] { 7.0, 9.0 }));
    }

    [Test]
    public void Name_DescribesDimensionAndAggregation()
    {
        Assert.That(new SampleCollapseStrategy(CollapseDimension.Fraction, AggregationType.Median).Name,
            Is.EqualTo("Collapse_Fraction_Median"));
    }
}
