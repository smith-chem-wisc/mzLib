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
        Assert.That(TopDownTestPaths.FindExistingFile("MM_1p1p4_GPTMD_Search", "AllPSMs.psmtsv"), Is.Not.Null);
        Assert.That(TopDownTestPaths.FindExistingFile("MM_1p1p4_GPTMD_Search", "AllProteoforms.psmtsv"), Is.Not.Null);
        Assert.That(TopDownTestPaths.FindExistingFile("MM_1p1p4_GPTMD_Search", "AllProteinGroups.tsv"), Is.Not.Null);
    }

    [Test]
    public void MetaMorpheusIndividualResultFiles_ArePresentInData()
    {
        Assert.That(TopDownTestPaths.FindExistingFile("MM_1p1p4_GPTMD_Search", "Individual File Results", "02-18-20_jurkat_td_rep1_fract7_Proteoforms.psmtsv"), Is.Not.Null);
        Assert.That(TopDownTestPaths.FindExistingFile("MM_1p1p4_GPTMD_Search", "Individual File Results", "02-18-20_jurkat_td_rep2_fract7_Proteoforms.psmtsv"), Is.Not.Null);
        Assert.That(TopDownTestPaths.FindExistingFile("MM_1p1p4_GPTMD_Search", "Individual File Results", "02-18-20_jurkat_td_rep1_fract7_PSMs.psmtsv"), Is.Not.Null);
        Assert.That(TopDownTestPaths.FindExistingFile("MM_1p1p4_GPTMD_Search", "Individual File Results", "02-18-20_jurkat_td_rep2_fract7_PSMs.psmtsv"), Is.Not.Null);
    }
}
