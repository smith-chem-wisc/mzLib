using MassSpectrometry;
using NUnit.Framework;
using Quantification;
using Quantification.Interfaces;
using Quantification.Strategies;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Test.Quantification;

/// <summary>
/// Covers the <see cref="ICollapseStrategy"/> implementations, which decide *which* samples are
/// combined. How the grouped values are reduced to one is <see cref="IAggregationStrategy"/>'s job
/// and is covered by <see cref="AggregationStrategyTests"/>.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class CollapseStrategyTests
{
    #region Fixture

    /// <summary>A row key that is just a name, so these tests are about columns only.</summary>
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

    private static double[] SingleRow<T>(QuantMatrix<T> matrix) where T : IEquatable<T> => matrix.GetRow(0);

    #endregion

    /// <summary>
    /// Collapsing fractions must leave technical replicates standing. This is the distinction the
    /// former SumCollapse and MeanCollapse could not express — they grouped by condition and
    /// biological replicate only, so they collapsed everything below that in one step.
    /// </summary>
    [Test]
    public void CollapseFractions_KeepsTechnicalReplicatesApart()
    {
        var columns = new List<ISampleInfo>
        {
            Lfq("Control", 1, 1, 1),
            Lfq("Control", 1, 1, 2),
            Lfq("Control", 1, 2, 1),
            Lfq("Control", 1, 2, 2)
        };

        var collapsed = new CollapseFractions().CollapseSamples(Matrix(columns, 10, 20, 30, 40), new SumAggregation());

        Assert.That(collapsed.ColumnCount, Is.EqualTo(2), "two technical replicates survive");
        Assert.That(SingleRow(collapsed), Is.EqualTo(new[] { 30.0, 70.0 }));
        Assert.That(collapsed.ColumnKeys.Select(s => s.TechnicalReplicate), Is.EqualTo(new[] { 1, 2 }));
        Assert.That(collapsed.ColumnKeys.Select(s => s.Fraction), Is.EqualTo(new[] { 0, 0 }), "fraction is zeroed");
    }

    /// <summary>The complementary case: collapsing technical replicates leaves the fractions alone.</summary>
    [Test]
    public void CollapseTechnicalReplicates_KeepsFractionsApart()
    {
        var columns = new List<ISampleInfo>
        {
            Lfq("Control", 1, 1, 1),
            Lfq("Control", 1, 2, 1),
            Lfq("Control", 1, 1, 2),
            Lfq("Control", 1, 2, 2)
        };

        var collapsed = new CollapseTechnicalReplicates()
            .CollapseSamples(Matrix(columns, 10, 30, 20, 40), new MedianAggregation());

        Assert.That(collapsed.ColumnCount, Is.EqualTo(2), "two fractions survive");
        Assert.That(SingleRow(collapsed), Is.EqualTo(new[] { 20.0, 30.0 }));
        Assert.That(collapsed.ColumnKeys.Select(s => s.Fraction), Is.EqualTo(new[] { 1, 2 }));
        Assert.That(collapsed.ColumnKeys.Select(s => s.TechnicalReplicate), Is.EqualTo(new[] { 0, 0 }));
    }

    [Test]
    public void CollapseBiologicalReplicates_GroupsByConditionFractionAndTechnicalReplicate()
    {
        var columns = new List<ISampleInfo>
        {
            Lfq("Control", 1, 1, 1),
            Lfq("Control", 2, 1, 1),
            Lfq("Treated", 1, 1, 1),
            Lfq("Treated", 2, 1, 1)
        };

        var collapsed = new CollapseBiologicalReplicates()
            .CollapseSamples(Matrix(columns, 10, 30, 100, 200), new MeanAggregation());

        Assert.That(collapsed.ColumnCount, Is.EqualTo(2));
        Assert.That(SingleRow(collapsed), Is.EqualTo(new[] { 20.0, 150.0 }));
        Assert.That(collapsed.ColumnKeys.Select(s => s.Condition), Is.EqualTo(new[] { "Control", "Treated" }));
        Assert.That(collapsed.ColumnKeys.Select(s => s.BiologicalReplicate), Is.EqualTo(new[] { 0, 0 }));
    }

    /// <summary>
    /// The grouping the former SumCollapse and MeanCollapse both performed — one column per condition
    /// and biological replicate. Paired with a Sum or Mean aggregation it reproduces them exactly.
    /// </summary>
    [Test]
    public void CollapseFractionsAndTechnicalReplicates_LeavesOneColumnPerConditionAndBiorep()
    {
        var columns = new List<ISampleInfo>
        {
            Lfq("Control", 1, 1, 1),
            Lfq("Control", 1, 2, 1),
            Lfq("Control", 1, 1, 2),
            Lfq("Control", 2, 1, 1)
        };
        var strategy = new CollapseFractionsAndTechnicalReplicates();

        var summed = strategy.CollapseSamples(Matrix(columns, 10, 20, 30, 500), new SumAggregation());
        Assert.That(summed.ColumnCount, Is.EqualTo(2));
        Assert.That(SingleRow(summed), Is.EqualTo(new[] { 60.0, 500.0 }));

        var meaned = strategy.CollapseSamples(Matrix(columns, 10, 20, 30, 500), new MeanAggregation());
        Assert.That(SingleRow(meaned), Is.EqualTo(new[] { 20.0, 500.0 }));

        Assert.That(summed.ColumnKeys.Select(s => s.BiologicalReplicate), Is.EqualTo(new[] { 1, 2 }),
            "biological replicate is kept");
        Assert.That(summed.ColumnKeys.Select(s => s.Fraction), Is.EqualTo(new[] { 0, 0 }));
        Assert.That(summed.ColumnKeys.Select(s => s.TechnicalReplicate), Is.EqualTo(new[] { 0, 0 }));
    }

    /// <summary>
    /// The same grouping under four aggregations — the point of splitting the two concerns is that
    /// this table is expressible at all.
    /// </summary>
    [Test]
    public void OneGrouping_ComposesWithEveryAggregation()
    {
        var columns = new List<ISampleInfo> { Lfq("C", 1, 1, 1), Lfq("C", 1, 1, 2), Lfq("C", 1, 1, 3) };
        var strategy = new CollapseFractions();

        Assert.Multiple(() =>
        {
            Assert.That(SingleRow(strategy.CollapseSamples(Matrix(columns, 10, 20, 30), new SumAggregation()))[0],
                Is.EqualTo(60.0).Within(1e-9));
            Assert.That(SingleRow(strategy.CollapseSamples(Matrix(columns, 10, 20, 30), new MeanAggregation()))[0],
                Is.EqualTo(20.0).Within(1e-9));
            Assert.That(SingleRow(strategy.CollapseSamples(Matrix(columns, 10, 20, 30), new MedianAggregation()))[0],
                Is.EqualTo(20.0).Within(1e-9));
            Assert.That(SingleRow(strategy.CollapseSamples(Matrix(columns, 10, 20, 30), new MaxAggregation()))[0],
                Is.EqualTo(30.0).Within(1e-9));
        });
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

        var collapsed = new CollapseTechnicalReplicates()
            .CollapseSamples(Matrix(columns, 100, 300), new SumAggregation());

        Assert.That(collapsed.ColumnCount, Is.EqualTo(1));
        Assert.That(collapsed.ColumnKeys[0], Is.InstanceOf<IsobaricQuantSampleInfo>());

        var isobaric = (IsobaricQuantSampleInfo)collapsed.ColumnKeys[0];
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

    [Test]
    public void CollapsingLabelFreeSamples_YieldsSpectraFileInfo()
    {
        var columns = new List<ISampleInfo> { Lfq("C", 1, 1, 1), Lfq("C", 1, 1, 2) };
        var collapsed = new CollapseFractions().CollapseSamples(Matrix(columns, 1, 2), new SumAggregation());

        Assert.That(collapsed.ColumnKeys[0], Is.InstanceOf<SpectraFileInfo>());
    }

    /// <summary>
    /// Strategies compose in sequence — the standard bottom-up recipe is sum the fractions, then take
    /// the median across technical replicates.
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

        var summed = new CollapseFractions()
            .CollapseSamples(Matrix(columns, 10, 20, 40, 60, 25, 25), new SumAggregation());
        Assert.That(SingleRow(summed), Is.EqualTo(new[] { 30.0, 100.0, 50.0 }));

        var final = new CollapseTechnicalReplicates().CollapseSamples(summed, new MedianAggregation());

        Assert.That(final.ColumnCount, Is.EqualTo(1));
        Assert.That(SingleRow(final)[0], Is.EqualTo(50.0).Within(1e-9), "median of 30, 100, 50");
    }

    /// <summary>Column order in equals column order out, so output is deterministic.</summary>
    [Test]
    public void CollapsedColumns_FollowFirstAppearanceOrder()
    {
        var columns = new List<ISampleInfo>
        {
            Lfq("Zebra", 1, 1, 1),
            Lfq("Apple", 1, 1, 1),
            Lfq("Zebra", 1, 1, 2)
        };

        var collapsed = new CollapseFractions().CollapseSamples(Matrix(columns, 1, 2, 3), new SumAggregation());

        Assert.That(collapsed.ColumnKeys.Select(s => s.Condition), Is.EqualTo(new[] { "Zebra", "Apple" }));
        Assert.That(SingleRow(collapsed), Is.EqualTo(new[] { 4.0, 2.0 }));
    }

    [Test]
    public void NoGroupHasMoreThanOneMember_ValuesAreUnchanged()
    {
        var columns = new List<ISampleInfo> { Lfq("A", 1, 1, 1), Lfq("B", 1, 1, 1) };
        var collapsed = new CollapseFractions().CollapseSamples(Matrix(columns, 7, 9), new SumAggregation());

        Assert.That(collapsed.ColumnCount, Is.EqualTo(2));
        Assert.That(SingleRow(collapsed), Is.EqualTo(new[] { 7.0, 9.0 }));
    }

    /// <summary>NoCollapse passes the matrix straight through and never consults the aggregation.</summary>
    [Test]
    public void NoCollapse_ReturnsTheSameMatrix()
    {
        var columns = new List<ISampleInfo> { Lfq("C", 1, 1, 1), Lfq("C", 1, 1, 2) };
        var matrix = Matrix(columns, 1, 2);

        Assert.That(new NoCollapse().CollapseSamples(matrix, null), Is.SameAs(matrix));
    }

    /// <summary>
    /// A grouping strategy cannot do anything without an aggregation, so it must say so rather than
    /// failing with a null reference somewhere inside the row loop.
    /// </summary>
    [Test]
    public void CollapseSamples_WithoutAnAggregation_Throws()
    {
        var columns = new List<ISampleInfo> { Lfq("C", 1, 1, 1) };

        Assert.That(() => new CollapseFractions().CollapseSamples(Matrix(columns, 1), null),
            Throws.ArgumentNullException);
    }

    [Test]
    public void Names_DescribeTheGroupingAndAreDistinct()
    {
        var strategies = new ICollapseStrategy[]
        {
            new NoCollapse(),
            new CollapseFractions(),
            new CollapseTechnicalReplicates(),
            new CollapseBiologicalReplicates(),
            new CollapseFractionsAndTechnicalReplicates()
        };

        Assert.That(strategies.Select(s => s.Name), Is.Unique);
        Assert.That(new CollapseFractions().Name, Is.EqualTo("Collapse Fractions"));
    }
}
