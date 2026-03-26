using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests;

[TestFixture]
public class TestFdrBenchResultFiles
{
    private static string GetTestFile(string relativePath) =>
        Path.Combine(TestContext.CurrentContext.TestDirectory, relativePath);

    [Test]
    public void CanReadFdrBenchPeptideFile()
    {
        var filePath = GetTestFile("FileReadingTests/ExternalFileTypes/peptides_fdrbench_peptide.tsv");
        var peptideFile = new FdrBenchPeptideFile(filePath);
        var results = peptideFile.Results;

        Assert.That(results, Has.Count.EqualTo(9));

        var first = results.First();
        Assert.Multiple(() =>
        {
            Assert.That(first.Run, Is.EqualTo("20230406_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_I_1"));
            Assert.That(first.Peptide, Is.EqualTo("AAAPAPEEEMDECEQALAAEPK"));
            Assert.That(first.ModifiedPeptide, Is.EqualTo("AAAPAPEEEMDEC(UniMod:4)EQALAAEPK"));
            Assert.That(first.Charge, Is.EqualTo(3));
            Assert.That(first.QValue, Is.EqualTo(1e-9));
            Assert.That(first.Pep, Is.EqualTo(0));
            Assert.That(first.Protein, Is.EqualTo("AAAPAPEEEMDECEQALAAEPK_target"));
            Assert.That(first.Score, Is.EqualTo(1));
        });
    }

    [Test]
    public void CanReadFdrBenchProteinFile()
    {
        var filePath = GetTestFile("FileReadingTests/ExternalFileTypes/peptides_fdrbench_protein.tsv");
        var proteinFile = new FdrBenchProteinFile(filePath);
        var results = proteinFile.Results;

        Assert.That(results, Has.Count.EqualTo(9));

        var first = results.First();
        Assert.Multiple(() =>
        {
            Assert.That(first.Protein, Is.EqualTo("P62857"));
            Assert.That(first.QValue, Is.EqualTo(1.64447e-4).Within(1e-10));
            Assert.That(first.Score, Is.EqualTo(1));
        });
    }

    [Test]
    public void WritingPeptideRecordsProducesEquivalentFile()
    {
        var sourcePath = GetTestFile("FileReadingTests/ExternalFileTypes/peptides_fdrbench_peptide.tsv");
        var sourceRecords = new FdrBenchPeptideFile(sourcePath).Results
            .Select(ClonePeptide)
            .ToList();

        var tempPath = Path.Combine(Path.GetTempPath(), Guid.NewGuid().ToString() + "_fdrbench_peptide.tsv");
        try
        {
            var writer = new FdrBenchPeptideFile()
            {
                Results = sourceRecords
            };
            writer.WriteResults(tempPath);

            var roundTripRecords = new FdrBenchPeptideFile(tempPath).Results;
            AssertPeptideRecordsEqual(sourceRecords, roundTripRecords);
        }
        finally
        {
            if (File.Exists(tempPath))
            {
                File.Delete(tempPath);
            }
        }
    }

    [Test]
    public void WritingProteinRecordsProducesEquivalentFile()
    {
        var sourcePath = GetTestFile("FileReadingTests/ExternalFileTypes/peptides_fdrbench_protein.tsv");
        var sourceRecords = new FdrBenchProteinFile(sourcePath).Results
            .Select(CloneProtein)
            .ToList();

        var tempPath = Path.Combine(Path.GetTempPath(), Guid.NewGuid().ToString() + "_fdrbench_protein.tsv");
        try
        {
            var writer = new FdrBenchProteinFile()
            {
                Results = sourceRecords
            };
            writer.WriteResults(tempPath);

            var roundTripRecords = new FdrBenchProteinFile(tempPath).Results;
            AssertProteinRecordsEqual(sourceRecords, roundTripRecords);
        }
        finally
        {
            if (File.Exists(tempPath))
            {
                File.Delete(tempPath);
            }
        }
    }

    private static void AssertPeptideRecordsEqual(IReadOnlyList<FdrBenchPeptide> expected, IReadOnlyList<FdrBenchPeptide> actual)
    {
        Assert.That(actual, Has.Count.EqualTo(expected.Count));
        for (int i = 0; i < expected.Count; i++)
        {
            Assert.Multiple(() =>
            {
                Assert.That(actual[i].Run, Is.EqualTo(expected[i].Run));
                Assert.That(actual[i].Peptide, Is.EqualTo(expected[i].Peptide));
                Assert.That(actual[i].ModifiedPeptide, Is.EqualTo(expected[i].ModifiedPeptide));
                Assert.That(actual[i].Charge, Is.EqualTo(expected[i].Charge));
                Assert.That(actual[i].QValue, Is.EqualTo(expected[i].QValue));
                Assert.That(actual[i].Pep, Is.EqualTo(expected[i].Pep));
                Assert.That(actual[i].Protein, Is.EqualTo(expected[i].Protein));
                Assert.That(actual[i].Score, Is.EqualTo(expected[i].Score));
            });
        }
    }

    private static void AssertProteinRecordsEqual(IReadOnlyList<FdrBenchProtein> expected, IReadOnlyList<FdrBenchProtein> actual)
    {
        Assert.That(actual, Has.Count.EqualTo(expected.Count));
        for (int i = 0; i < expected.Count; i++)
        {
            Assert.Multiple(() =>
            {
                Assert.That(actual[i].Protein, Is.EqualTo(expected[i].Protein));
                Assert.That(actual[i].QValue, Is.EqualTo(expected[i].QValue));
                Assert.That(actual[i].Score, Is.EqualTo(expected[i].Score));
            });
        }
    }

    private static FdrBenchPeptide ClonePeptide(FdrBenchPeptide peptide) => new()
    {
        Run = peptide.Run,
        Peptide = peptide.Peptide,
        ModifiedPeptide = peptide.ModifiedPeptide,
        Charge = peptide.Charge,
        QValue = peptide.QValue,
        Pep = peptide.Pep,
        Protein = peptide.Protein,
        Score = peptide.Score
    };

    private static FdrBenchProtein CloneProtein(FdrBenchProtein protein) => new()
    {
        Protein = protein.Protein,
        QValue = protein.QValue,
        Score = protein.Score
    };
}
