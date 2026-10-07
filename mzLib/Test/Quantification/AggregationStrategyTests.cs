using NUnit.Framework;
using Quantification.Interfaces;
using Quantification.Strategies;
using System;
using System.Diagnostics.CodeAnalysis;

namespace Test.Quantification;

/// <summary>
/// Covers the <see cref="IAggregationStrategy"/> implementations, which reduce several intensity
/// values to one. Separated from <see cref="ICollapseStrategy"/>, which decides which values are
/// grouped in the first place.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class AggregationStrategyTests
{
    private static readonly double[] Unsorted = { 40, 10, 30, 20 };
    private static readonly double[] OddCount = { 30, 10, 20 };

    [Test]
    public void Sum_AddsEveryValue() =>
        Assert.That(new SumAggregation().Aggregate(Unsorted), Is.EqualTo(100.0).Within(1e-9));

    [Test]
    public void Mean_DividesByCount() =>
        Assert.That(new MeanAggregation().Aggregate(Unsorted), Is.EqualTo(25.0).Within(1e-9));

    [Test]
    public void Max_TakesTheLargest() =>
        Assert.That(new MaxAggregation().Aggregate(Unsorted), Is.EqualTo(40.0).Within(1e-9));

    /// <summary>An even count averages the two middle values rather than picking one.</summary>
    [Test]
    public void Median_EvenCount_AveragesTheTwoMiddleValues() =>
        Assert.That(new MedianAggregation().Aggregate(Unsorted), Is.EqualTo(25.0).Within(1e-9));

    [Test]
    public void Median_OddCount_TakesTheMiddleValue() =>
        Assert.That(new MedianAggregation().Aggregate(OddCount), Is.EqualTo(20.0).Within(1e-9));

    /// <summary>
    /// The median copies to the stack below a threshold and rents an array above it. Both paths must
    /// agree, so this crosses the boundary.
    /// </summary>
    [TestCase(2)]
    [TestCase(63)]
    [TestCase(64)]
    [TestCase(65)]
    [TestCase(500)]
    public void Median_AgreesAcrossTheStackAllocBoundary(int count)
    {
        // Values 1..count in reverse, so the median is the mean of the two central integers
        // for an even count and the central integer for an odd one.
        var values = new double[count];
        for (int i = 0; i < count; i++) values[i] = count - i;

        double expected = count % 2 == 0 ? (count / 2) + 0.5 : (count + 1) / 2.0;

        Assert.That(new MedianAggregation().Aggregate(values), Is.EqualTo(expected).Within(1e-9));
    }

    /// <summary>Aggregating an empty group yields 0 rather than throwing or producing NaN.</summary>
    [Test]
    public void EveryStrategy_OnAnEmptySpan_ReturnsZero()
    {
        foreach (IAggregationStrategy strategy in AllStrategies())
            Assert.That(strategy.Aggregate(ReadOnlySpan<double>.Empty), Is.EqualTo(0.0),
                strategy.Name + " should aggregate an empty group to zero");
    }

    /// <summary>A one-value group aggregates to that value, whatever the strategy.</summary>
    [Test]
    public void EveryStrategy_OnASingleValue_ReturnsThatValue()
    {
        foreach (IAggregationStrategy strategy in AllStrategies())
            Assert.That(strategy.Aggregate(new double[] { 7.5 }), Is.EqualTo(7.5).Within(1e-9),
                strategy.Name + " should pass a single value through");
    }

    /// <summary>Callers reuse the buffer, so no strategy may reorder or overwrite it.</summary>
    [Test]
    public void EveryStrategy_DoesNotMutateTheInput()
    {
        foreach (IAggregationStrategy strategy in AllStrategies())
        {
            var buffer = new double[] { 40, 10, 30, 20 };
            strategy.Aggregate(buffer);
            Assert.That(buffer, Is.EqualTo(new double[] { 40, 10, 30, 20 }),
                strategy.Name + " must not mutate the values it is given");
        }
    }

    [Test]
    public void Names_AreDistinctAndDescriptive()
    {
        var names = Array.ConvertAll(AllStrategies(), s => s.Name);
        Assert.That(names, Is.Unique);
        Assert.That(names, Is.EquivalentTo(new[] { "Sum", "Mean", "Median", "Max" }));
    }

    private static IAggregationStrategy[] AllStrategies() =>
        new IAggregationStrategy[]
        {
            new SumAggregation(), new MeanAggregation(), new MedianAggregation(), new MaxAggregation()
        };
}
