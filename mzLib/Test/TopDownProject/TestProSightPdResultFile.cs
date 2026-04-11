using System.Diagnostics.CodeAnalysis;
using NUnit.Framework;
using Readers;

namespace Test.TopDownProject;

[TestFixture]
[ExcludeFromCodeCoverage]
internal class TestProSightPdResultFile
{
    [Test]
    public void ReadsTargetPsmsFromSqliteDatabase()
    {
        var path = @"DataFiles/02-18-20_jurkat_td_rep2_fract7.pdResult";
        var file = new ProSightPdResultFile(path);

        Assert.That(file.Results, Is.Not.Empty);
        Assert.That(file.Results[0].Sequence, Is.Not.Null.Or.Empty);
        Assert.That(file.Results[0].Qvalue, Is.GreaterThanOrEqualTo(0));

        var proteins = file.LoadProteins();
        Assert.That(proteins, Is.Not.Empty);
        Assert.That(proteins[0].Accession, Is.Not.Null.Or.Empty);
        Assert.That(proteins[0].Qvalue, Is.GreaterThanOrEqualTo(0));
    }
}
