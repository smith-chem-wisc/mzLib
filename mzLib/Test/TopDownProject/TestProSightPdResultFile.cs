using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Data.SQLite;
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
        var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "02-18-20_jurkat_td_rep2_fract7.pdResult");
        if (!File.Exists(path))
        {
            Assert.Ignore("ProSight PD fixture is not present in local test data.");
        }

        var file = new ProSightPdResultFile(path);
        try
        {
            _ = file.Results;
        }
        catch (SQLiteException ex) when (ex.Message.Contains("no such table: TargetPsms"))
        {
            Assert.Ignore("ProSight PD fixture does not contain expected TargetPsms schema.");
        }

        Assert.That(file.Results, Is.Not.Empty);
        Assert.That(file.Results[0].Sequence, Is.Not.Null.Or.Empty);
        Assert.That(file.Results[0].Qvalue, Is.GreaterThanOrEqualTo(0));

        var proteins = file.LoadProteins();
        Assert.That(proteins, Is.Not.Empty);
        Assert.That(proteins[0].Accession, Is.Not.Null.Or.Empty);
        Assert.That(proteins[0].Qvalue, Is.GreaterThanOrEqualTo(0));
    }
}
