using NUnit.Framework;
using TopDownSimulator.Extraction;

namespace Test.TopDownSimulator;

[TestFixture]
public class CoeluterFinderTests
{
    [Test]
    public void FindsNearbyRecordsInSameFile()
    {
        var anchor = new MmResultRecord("sample", 100, 101, 10, 10000.0, 20.0, 50.0, "A", "P1", "anchor");
        var results = new[]
        {
            anchor,
            new MmResultRecord("sample", 105, 106, 11, 10007.0, 20.2, 40.0, "B", "P2", "nearby"),
            new MmResultRecord("sample", 110, 111, 12, 10030.0, 20.1, 40.0, "C", "P3", "far-mass"),
            new MmResultRecord("other", 115, 116, 12, 10005.0, 20.1, 40.0, "D", "P4", "other-file"),
        };

        var finder = new CoeluterFinder();
        var coeluters = finder.FindCoeluters(results, anchor, rtHalfWidth: 0.5, massHalfWidthDa: 10.0);

        Assert.That(coeluters, Has.Count.EqualTo(1));
        Assert.That(coeluters[0].Identifier, Is.EqualTo("nearby"));
    }
}
