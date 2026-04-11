using System.Diagnostics.CodeAnalysis;
using NUnit.Framework;
using Test.TopDownEngine.Shared;
using System.IO;
using Readers;

namespace Test.TopDownEngine.Harness;

[TestFixture]
[ExcludeFromCodeCoverage]
public class TopDownRawDataSmokeTests
{
    [Test]
    public void RawDataFolder_IsDiscoverable()
    {
        var raw1 = TopDownTestPaths.FindExistingFile("02-18-20_jurkat_td_rep1_fract7.raw");
        var raw2 = TopDownTestPaths.FindExistingFile("02-18-20_jurkat_td_rep2_fract7.raw");
        if (raw1 is null || raw2 is null)
        {
            Assert.Ignore("Top-down RAW fixtures are not present in local test data.");
        }

        Assert.That(raw1, Is.Not.Null);
        Assert.That(raw2, Is.Not.Null);
    }

    [Test]
    public void MetaMorpheusResultFolder_IsDiscoverable()
    {
        var allPsms = TopDownTestPaths.FindExistingFile("MM_1p1p4_GPTMD_Search", "AllPSMs.psmtsv");
        if (allPsms is null)
        {
            Assert.Ignore("MetaMorpheus top-down fixtures are not present in local test data.");
        }

        Assert.That(allPsms, Is.Not.Null);
    }

    [Test]
    public void RawFile_IsReadableByPath()
    {
        var rawPath = TopDownTestPaths.FindExistingFile("02-18-20_jurkat_td_rep1_fract7.raw");
        if (rawPath is null)
        {
            Assert.Ignore("Top-down RAW fixture is not present in local test data.");
        }

        Assert.That(rawPath, Is.Not.Null);
        Assert.That(File.Exists(rawPath), Is.True);
    }
}
