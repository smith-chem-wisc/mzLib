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
        Assert.That(TopDownTestPaths.FindExistingFile("02-18-20_jurkat_td_rep1_fract7.raw"), Is.Not.Null);
        Assert.That(TopDownTestPaths.FindExistingFile("02-18-20_jurkat_td_rep2_fract7.raw"), Is.Not.Null);
    }

    [Test]
    public void MetaMorpheusResultFolder_IsDiscoverable()
    {
        Assert.That(TopDownTestPaths.FindExistingFile("MM_1p1p4_GPTMD_Search", "AllPSMs.psmtsv"), Is.Not.Null);
    }

    [Test]
    public void RawFile_IsReadableByPath()
    {
        var rawPath = TopDownTestPaths.FindExistingFile("02-18-20_jurkat_td_rep1_fract7.raw");
        Assert.That(rawPath, Is.Not.Null);
        Assert.That(File.Exists(rawPath), Is.True);
    }
}
