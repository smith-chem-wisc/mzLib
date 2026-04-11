using NUnit.Framework;

namespace Test.TopDownEngine.Shared;

internal static class TopDownAssertions
{
    internal static void AssertWithinTolerance(double expected, double actual, double tolerance, string message)
    {
        Assert.That(actual, Is.EqualTo(expected).Within(tolerance), message);
    }
}
