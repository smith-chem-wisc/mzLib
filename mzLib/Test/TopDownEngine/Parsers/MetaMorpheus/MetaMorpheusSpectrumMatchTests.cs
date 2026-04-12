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
        var path = TopDownTestPaths.FindExistingFile(
            "Frac7_GPTMD_Search/Task2-TopDownSearch/AllPSMs.psmtsv",
            "Frac7_GPTMD_Search\\Task2-TopDownSearch\\AllPSMs.psmtsv",
            "MM_1p1p4_GPTMD_Search/AllPSMs.psmtsv",
            "MM_1p1p4_GPTMD_Search\\AllPSMs.psmtsv",
            "AllPSMs.psmtsv");
        if (path is null)
        {
            Assert.Ignore("MetaMorpheus parser fixtures are not present in local test data.");
        }

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
        var path = TopDownTestPaths.FindExistingFile(
            "Frac7_GPTMD_Search/Task2-TopDownSearch/AllProteoforms.psmtsv",
            "Frac7_GPTMD_Search\\Task2-TopDownSearch\\AllProteoforms.psmtsv",
            "MM_1p1p4_GPTMD_Search/AllProteoforms.psmtsv",
            "MM_1p1p4_GPTMD_Search\\AllProteoforms.psmtsv",
            "AllProteoforms.psmtsv");
        if (path is null)
        {
            Assert.Ignore("MetaMorpheus parser fixtures are not present in local test data.");
        }

        Assert.That(path, Is.Not.Null);

        var results = SpectrumMatchTsvReader.ReadTsv(path!, out var warnings);

        Assert.That(results, Is.Not.Empty);
        Assert.That(warnings, Is.Not.Null);
        Assert.That(results[0].FileName, Is.Not.Empty);
        Assert.That(results[0].PrecursorCharge, Is.GreaterThan(0));
    }
}
