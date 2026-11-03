using MathNet.Numerics.Statistics;
using NUnit.Framework;
using StatisticalModels;
using System.Diagnostics.CodeAnalysis;

namespace Test.StatisticalModels;

/// <summary>
/// Guards the project's name. mzLib ships MathNet.Numerics, and downstream code calls its static class as
/// <c>Statistics.Median(...)</c> after <c>using MathNet.Numerics.Statistics;</c>. When this project's namespace was
/// named <c>Statistics</c>, that call bound to the namespace and MetaMorpheus stopped compiling (CS0234 in
/// BinTreeStructure.cs). This file imports both, exactly as a consumer would, so a colliding name fails the build here.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class NamespaceCollisionTests
{
    [Test]
    public void MathNetStatisticsStillResolvesBesideThisProject()
    {
        Assert.That(Statistics.Median(new[] { 3.0, 1.0, 2.0 }), Is.EqualTo(2.0));
        Assert.That(MultipleTesting.BenjaminiHochberg(new[] { 0.01 })[0], Is.EqualTo(0.01));
    }
}
