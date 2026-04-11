using System.Diagnostics.CodeAnalysis;
using NUnit.Framework;
using Readers;
using Test.TopDownEngine.Shared;

namespace Test.TopDownEngine.Parsers.MetaMorpheus;

[TestFixture]
[ExcludeFromCodeCoverage]
public class MetaMorpheusSpectrumMatchTests
{
    [Test]
    public void AllPsms_File_Parses()
    {
        var path = TopDownTestPaths.FindExistingFile("MM_1p1p4_GPTMD_Search", "AllPSMs.psmtsv");
        Assert.That(path, Is.Not.Null);

        var results = SpectrumMatchTsvReader.ReadPsmTsv(path!, out var warnings);

        Assert.That(results, Is.Not.Empty);
        Assert.That(warnings, Is.Not.Null);
        Assert.That(results[0].FileName, Is.Not.Empty);
        Assert.That(results[0].PrecursorCharge, Is.GreaterThan(0));
        Assert.That(results[0].RetentionTime, Is.GreaterThanOrEqualTo(-1));
        Assert.That(results[0].MonoisotopicMass, Is.GreaterThan(0));
    }

    [Test]
    public void AllProteoforms_File_Parses()
    {
        var path = TopDownTestPaths.FindExistingFile("MM_1p1p4_GPTMD_Search", "AllProteoforms.psmtsv");
        Assert.That(path, Is.Not.Null);

        var results = SpectrumMatchTsvReader.ReadTsv(path!, out var warnings);

        Assert.That(results, Is.Not.Empty);
        Assert.That(warnings, Is.Not.Null);
        Assert.That(results[0].FileName, Is.Not.Empty);
        Assert.That(results[0].PrecursorCharge, Is.GreaterThan(0));
    }
}
