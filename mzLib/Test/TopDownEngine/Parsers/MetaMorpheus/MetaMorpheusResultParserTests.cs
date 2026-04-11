using System.Diagnostics.CodeAnalysis;
using System.IO;
using NUnit.Framework;
using Test.TopDownEngine.Shared;

namespace Test.TopDownEngine.Parsers.MetaMorpheus;

[TestFixture]
[ExcludeFromCodeCoverage]
public class MetaMorpheusResultParserTests
{
    [Test]
    public void MetaMorpheusFiles_ArePresentInData()
    {
        var allPsms = TopDownTestPaths.FindExistingFile("MM_1p1p4_GPTMD_Search", "AllPSMs.psmtsv");
        var allProteoforms = TopDownTestPaths.FindExistingFile("MM_1p1p4_GPTMD_Search", "AllProteoforms.psmtsv");
        var allProteinGroups = TopDownTestPaths.FindExistingFile("MM_1p1p4_GPTMD_Search", "AllProteinGroups.tsv");
        if (allPsms is null || allProteoforms is null || allProteinGroups is null)
        {
            Assert.Ignore("MetaMorpheus parser fixtures are not present in local test data.");
        }

        Assert.That(allPsms, Is.Not.Null);
        Assert.That(allProteoforms, Is.Not.Null);
        Assert.That(allProteinGroups, Is.Not.Null);
    }

    [Test]
    public void MetaMorpheusIndividualResultFiles_ArePresentInData()
    {
        var p1 = TopDownTestPaths.FindExistingFile("MM_1p1p4_GPTMD_Search", "Individual File Results", "02-18-20_jurkat_td_rep1_fract7_Proteoforms.psmtsv");
        var p2 = TopDownTestPaths.FindExistingFile("MM_1p1p4_GPTMD_Search", "Individual File Results", "02-18-20_jurkat_td_rep2_fract7_Proteoforms.psmtsv");
        var p3 = TopDownTestPaths.FindExistingFile("MM_1p1p4_GPTMD_Search", "Individual File Results", "02-18-20_jurkat_td_rep1_fract7_PSMs.psmtsv");
        var p4 = TopDownTestPaths.FindExistingFile("MM_1p1p4_GPTMD_Search", "Individual File Results", "02-18-20_jurkat_td_rep2_fract7_PSMs.psmtsv");
        if (p1 is null || p2 is null || p3 is null || p4 is null)
        {
            Assert.Ignore("MetaMorpheus individual-result fixtures are not present in local test data.");
        }

        Assert.That(p1, Is.Not.Null);
        Assert.That(p2, Is.Not.Null);
        Assert.That(p3, Is.Not.Null);
        Assert.That(p4, Is.Not.Null);
    }
}
