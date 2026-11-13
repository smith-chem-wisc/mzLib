using NUnit.Framework;
using Readers;
using Readers.DeNovoMapping;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;

namespace Test.FileReadingTests.DeNovoMapping
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestDeNovoMappingResultFile
    {
        private static string FixturePath => Path.Combine(TestContext.CurrentContext.TestDirectory,
            @"FileReadingTests\DeNovoMapping\TestData\DeNovo_H1N1_C1_3_MappingSummary.tsv");

        private string _directory = null!;

        [SetUp]
        public void SetUp()
        {
            _directory = Path.Combine(TestContext.CurrentContext.TestDirectory, "DeNovoMappingResultFileTests", Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(_directory);
        }

        [TearDown]
        public void TearDown()
        {
            if (Directory.Exists(_directory))
                Directory.Delete(_directory, true);
        }

        [Test]
        public void LoadsFixtureWrittenFromRealData()
        {
            // H1N1_Virus_Spike_in_C1_3 Casanovo predictions mapped onto three viral proteomes.
            var file = FileReader.ReadFile<DeNovoMappingResultFile>(FixturePath);

            Assert.That(file.Software, Is.EqualTo(Software.Mzlib));
            Assert.That(file.FileType, Is.EqualTo(SupportedFileType.DenovoMappingResult));
            Assert.That(file.Select(r => r.DatabaseIdentifier), Is.EqualTo(new[] { "UP000009255", "UP000240399", "UP000464024" }));

            var h1n1 = file.First();
            Assert.Multiple(() =>
            {
                Assert.That(h1n1.DatabaseProteinCount, Is.EqualTo(13));
                Assert.That(h1n1.TargetPredictions, Is.EqualTo(10));
                Assert.That(h1n1.DecoyPredictions, Is.EqualTo(0));
                Assert.That(h1n1.UniqueProteinsMapped, Is.EqualTo(2));
                Assert.That(h1n1.TargetPredictionScores, Has.Count.EqualTo(10));
                Assert.That(h1n1.RetentionTimeErrors, Has.Count.EqualTo(10));
                // Calibrated minutes; ten errors, so the median is the mean of the middle two (0.3476 and 0.7126).
                Assert.That(h1n1.MedianRtError, Is.EqualTo(0.5301).Within(1e-3));
                Assert.That(h1n1.PredictionCountsByFile["H1N1_Virus_Spike_in_C1_3_DefaultParams"], Is.EqualTo(10));
            });

            var empty = file.Last();
            Assert.That(empty.TotalPredictions, Is.EqualTo(0));
            Assert.That(empty.MeanRtError, Is.NaN);
            Assert.That(empty.TargetPredictionScores, Is.Empty);
        }

        [Test]
        public void GenericReaderReturnsMappingSummary()
        {
            IResultFile file = FileReader.ReadResultFile(FixturePath);

            Assert.That(file, Is.InstanceOf<DeNovoMappingResultFile>());
            Assert.That(file.FileType, Is.EqualTo(SupportedFileType.DenovoMappingResult));
            Assert.That(file.Software, Is.EqualTo(Software.Mzlib));
        }

        [Test]
        public void RenamedSummaryIsRecognizedByHeader()
        {
            string renamed = Path.Combine(_directory, "viral_denovo_results.tsv");
            File.Copy(FixturePath, renamed);

            Assert.That(renamed.ParseFileType(), Is.EqualTo(SupportedFileType.DenovoMappingResult));
            Assert.That(FileReader.ReadResultFile(renamed), Is.InstanceOf<DeNovoMappingResultFile>());
        }

        [Test]
        public void WriteAppendsTheExtensionWhenMissing()
        {
            var file = new DeNovoMappingResultFile { Results = [new DeNovoMappingResult { DatabaseIdentifier = "Db" }] };
            string requested = Path.Combine(_directory, "summary");

            file.WriteResults(requested);

            Assert.That(file.FilePath, Is.EqualTo(requested + "MappingSummary.tsv"));
            Assert.That(File.Exists(file.FilePath));
            Assert.That(file.FilePath.ParseFileType(), Is.EqualTo(SupportedFileType.DenovoMappingResult));
        }

        [Test]
        public void PerFileCountsRoundTrip()
        {
            string path = Path.Combine(_directory, "CountsMappingSummary.tsv");
            var original = new DeNovoMappingResultFile(path)
            {
                Results =
                [
                    new DeNovoMappingResult
                    {
                        DatabaseIdentifier = "DbA",
                        TargetPredictions = 3,
                        TargetPredictionScores = [0.9, 0.8, 0.7],
                        PredictionCountsByFile = new Dictionary<string, int> { ["RunOne"] = 2, ["RunTwo"] = 1 },
                    },
                    new DeNovoMappingResult
                    {
                        DatabaseIdentifier = "DbB",
                        TargetPredictions = 1,
                        TargetPredictionScores = [0.5],
                        PredictionCountsByFile = new Dictionary<string, int> { ["RunThree"] = 1 },
                    },
                ]
            };
            original.WriteResults(path);

            var loaded = FileReader.ReadFile<DeNovoMappingResultFile>(path);

            var dbA = loaded.Single(r => r.DatabaseIdentifier == "DbA");
            var dbB = loaded.Single(r => r.DatabaseIdentifier == "DbB");
            Assert.That(dbA.PredictionCountsByFile, Is.EqualTo(new Dictionary<string, int> { ["RunOne"] = 2, ["RunThree"] = 0, ["RunTwo"] = 1 }));
            Assert.That(dbB.PredictionCountsByFile["RunThree"], Is.EqualTo(1));
            Assert.That(dbA.TotalPredictions, Is.EqualTo(3));
            Assert.That(dbA.MeanPredictionScore, Is.EqualTo(0.8).Within(1e-9));
        }

        [Test]
        public void LegacyEntrapmentColumnsAreNotTreatedAsFiles()
        {
            // MetaCurator's summaries carried entrapment columns that mzLib's mapper does not produce.
            string path = Path.Combine(_directory, "LegacyMappingSummary.tsv");
            File.WriteAllLines(path,
            [
                "DatabaseIdentifier\tDatabaseProteinCount\tDatabasePeptideCount\tTotalPredictions\tTargetPredictions\tEntrapmentTargetPredictions\tDecoyPredictions\tEntrapmentDecoyPredictions\tUniquePeptidesMapped\tUniqueProteinsMapped\tMeanRtError\tMedianRtError\tStdDevRtError\tMeanPredictionScore\tMedianPredictionScore\tStdDevPredictionScore\tRetentionTimeErrors\tPredictionScores\tTargetPredictionScores\tDecoyPredictionScores\tRunOne",
                "DB001\t100\t500\t50\t45\t0\t5\t0\t40\t30\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\t\t\t\t\t7",
            ]);

            var result = FileReader.ReadFile<DeNovoMappingResultFile>(path).Single();

            Assert.That(result.TargetPredictions, Is.EqualTo(45));
            Assert.That(result.DecoyPredictions, Is.EqualTo(5));
            Assert.That(result.PredictionCountsByFile, Is.EqualTo(new Dictionary<string, int> { ["RunOne"] = 7 }));
        }

        [Test]
        public void EmptyFileLoadsNoResults()
        {
            string path = Path.Combine(_directory, "BlankMappingSummary.tsv");
            File.WriteAllText(path, string.Empty);

            Assert.That(FileReader.ReadFile<DeNovoMappingResultFile>(path).Results, Is.Empty);
        }
    }
}
