using NUnit.Framework;
using Omics;
using Proteomics;
using Readers.DeNovoMapping;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using UsefulProteomicsDatabases;

namespace Test.FileReadingTests.DeNovoMapping
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestDeNovoMapper
    {
        // Tryptic peptides (7+ residues): VLSAADK, VGAHAGEYGAEALER, MFLSFPTTK
        private const string ProteinA = "MKVLSAADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTK";
        // Tryptic peptides: AGLIDEFLSK, TPEPTIDER
        private const string ProteinB = "AGLIDEFLSKTPEPTIDER";
        private const string UnrelatedProtein = "GGGGWWWWGGGGKHHHHHHHHHHR";

        private string _directory = null!;

        [SetUp]
        public void SetUp()
        {
            _directory = Path.Combine(TestContext.CurrentContext.TestDirectory, "DeNovoMappingTests", Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(_directory);
        }

        [TearDown]
        public void TearDown()
        {
            if (Directory.Exists(_directory))
                Directory.Delete(_directory, true);
        }

        [Test]
        public void MapsTargetsDecoysAndUnmappedPredictions()
        {
            string fasta = WriteFasta("TargetDb.fasta", ("PA", ProteinA), ("PB", ProteinB));
            string decoyPeptide = FirstDecoyOnlyPeptide(fasta);

            string mzTab = WriteMzTab("RunA.mztab", "file:///Z:/missing/RunA.mgf",
                ("VLSAADK", 0.9, 10.0),            // target, protein A
                ("AGLLDEFLSK", -0.25, 20.0),        // target, protein B, I written as L; score shifted to 0.75
                ("TPEPTLDER", 0.5, null),           // target, protein B; no retention time
                (decoyPeptide, 0.1, 30.0),          // decoy
                ("QQQQQQQK", 0.95, 40.0));          // not in the database

            var mapper = new DeNovoMapper([mzTab], new DeNovoMappingParameters { MaxThreads = 2 });
            var result = mapper.MapDatabases([fasta]).Single();

            Assert.Multiple(() =>
            {
                Assert.That(result.DatabaseIdentifier, Is.EqualTo("TargetDb"));
                Assert.That(result.DatabaseProteinCount, Is.EqualTo(2));
                Assert.That(result.TargetPredictions, Is.EqualTo(3));
                Assert.That(result.DecoyPredictions, Is.EqualTo(1));
                Assert.That(result.TotalPredictions, Is.EqualTo(4));
                Assert.That(result.UniquePeptidesMapped, Is.EqualTo(3));
                Assert.That(result.UniqueProteinsMapped, Is.EqualTo(2));
                Assert.That(result.TargetPredictionScores.OrderBy(x => x), Is.EqualTo(new[] { 0.5, 0.75, 0.9 }).Within(1e-9));
                Assert.That(result.DecoyPredictionScores, Is.EqualTo(new[] { 0.1 }).Within(1e-9));
                Assert.That(result.PredictionScores, Has.Count.EqualTo(4));
                Assert.That(result.PredictionCountsByFile["RunA"], Is.EqualTo(4));

                // Two target hits have an observed time; the third had none and its spectra file is missing.
                Assert.That(result.RetentionTimeErrors, Has.Count.EqualTo(2));
                Assert.That(result.MeanRtError, Is.Not.NaN);
                Assert.That(mapper.Warnings, Has.Count.EqualTo(1));
            });
        }

        [Test]
        public void IsoleucineIsOnlyEquivalentToLeucineWhenRequested()
        {
            string fasta = WriteFasta("IlDb.fasta", ("PB", ProteinB));
            string mzTab = WriteMzTab("RunIl.mztab", "file:///Z:/RunIl.mgf", ("AGLLDEFLSK", 0.9, 1.0));

            var normalized = new DeNovoMapper([mzTab], new DeNovoMappingParameters { MapRetentionTimes = false });
            var literal = new DeNovoMapper([mzTab], new DeNovoMappingParameters { MapRetentionTimes = false, ReplaceIWithL = false });

            Assert.That(normalized.MapDatabases([fasta]).Single().TargetPredictions, Is.EqualTo(1));
            Assert.That(literal.MapDatabases([fasta]).Single().TargetPredictions, Is.EqualTo(0));
        }

        [Test]
        public void SharedPeptideCountsOncePerPrediction()
        {
            // Both proteins contain VLSAADK.
            string fasta = WriteFasta("SharedDb.fasta", ("PA", ProteinA), ("PA2", "VLSAADKGGGGGGGR"));
            string mzTab = WriteMzTab("RunShared.mztab", "file:///Z:/RunShared.mgf", ("VLSAADK", 0.9, 1.0));

            var mapper = new DeNovoMapper([mzTab], new DeNovoMappingParameters { MapRetentionTimes = false });
            var result = mapper.MapDatabases([fasta]).Single();

            Assert.That(result.TargetPredictions, Is.EqualTo(1));
            Assert.That(result.UniqueProteinsMapped, Is.EqualTo(2));
            Assert.That(result.UniquePeptidesMapped, Is.EqualTo(1));
        }

        [Test]
        public void EachDatabaseGetsItsOwnRow()
        {
            string targetDb = WriteFasta("DbOne.fasta", ("PA", ProteinA));
            string otherDb = WriteFasta("DbTwo.fasta", ("PX", UnrelatedProtein));
            string mzTab = WriteMzTab("RunDbs.mztab", "file:///Z:/RunDbs.mgf", ("VLSAADK", 0.9, 1.0));

            var mapper = new DeNovoMapper([mzTab], new DeNovoMappingParameters { MapRetentionTimes = false });
            var results = mapper.MapDatabases([targetDb, otherDb]);

            Assert.That(results.Select(r => r.DatabaseIdentifier), Is.EqualTo(new[] { "DbOne", "DbTwo" }));
            Assert.That(results[0].TargetPredictions, Is.EqualTo(1));
            Assert.That(results[1].TotalPredictions, Is.EqualTo(0));
            Assert.That(results[1].MeanPredictionScore, Is.NaN);
        }

        [Test]
        public void MapsProteinsAlreadyInMemory()
        {
            string mzTab = WriteMzTab("RunMemory.mztab", "file:///Z:/RunMemory.mgf", ("VLSAADK", 0.9, 1.0));
            var proteins = new List<IBioPolymer>
            {
                new Protein(ProteinA, "PA"),
                new Protein("GGGGGGGRVLSAADK", "DECOY_PA", isDecoy: true),
            };

            var mapper = new DeNovoMapper([mzTab], new DeNovoMappingParameters { MapRetentionTimes = false });
            var result = mapper.MapDatabase("InMemory", proteins);

            Assert.That(result.TargetPredictions, Is.EqualTo(1));
            Assert.That(result.DecoyPredictions, Is.EqualTo(1));
            Assert.That(result.DatabaseProteinCount, Is.EqualTo(1));
        }

        [Test]
        public void ReadsRetentionTimesFromSpectraWhenMzTabHasNone()
        {
            // SCANS numbers deliberately differ from the positions Casanovo's index= refers to.
            string mgf = WriteMgf("RunSpectra.mgf", (101, 60), (105, 120), (110, 180));
            string mzTab = WriteMzTab("RunSpectra.mztab", "file:///Z:/elsewhere/RunSpectra.mgf",
                ("VLSAADK", 0.9, null), ("AGLLDEFLSK", 0.8, null), ("TPEPTLDER", 0.7, null));

            var mapper = new DeNovoMapper([mzTab], new DeNovoMappingParameters { SpectraFilePaths = [mgf] });

            Assert.That(mapper.Warnings, Is.Empty);
            Assert.That(mapper.Predictions.Select(p => p.ObservedRetentionTime),
                Is.EqualTo(new double?[] { 1.0, 2.0, 3.0 }).Within(1e-9));
            Assert.That(mapper.Predictions.All(p => p.PredictedRetentionTime.HasValue));
        }

        [Test]
        public void TrustsRetentionTimesAlreadyInMzTab()
        {
            // The spectra file would give different times; it must not be read.
            string mgf = WriteMgf("RunTrusted.mgf", (1, 600));
            string mzTab = WriteMzTab("RunTrusted.mztab", "file:///Z:/RunTrusted.mgf", ("VLSAADK", 0.9, 12.5));

            var mapper = new DeNovoMapper([mzTab], new DeNovoMappingParameters { SpectraFilePaths = [mgf] });

            Assert.That(mapper.Predictions.Single().ObservedRetentionTime, Is.EqualTo(12.5));
        }

        [Test]
        public void SummaryFileRoundTrips()
        {
            string fasta = WriteFasta("RoundTripDb.fasta", ("PA", ProteinA), ("PB", ProteinB));
            string first = WriteMzTab("RunOne.mztab", "file:///Z:/RunOne.mgf", ("VLSAADK", 0.9, 10.0), ("TPEPTLDER", 0.4, 11.0));
            string second = WriteMzTab("RunTwo.mztab", "file:///Z:/RunTwo.mgf", ("AGLLDEFLSK", 0.6, 12.0));

            var mapper = new DeNovoMapper([first, second], new DeNovoMappingParameters { MaxThreads = 1 });
            var results = mapper.MapDatabases([fasta]);
            string summaryPath = Path.Combine(_directory, "MappingSummary.tsv");
            DeNovoMapper.WriteSummary(results, summaryPath);

            var header = File.ReadLines(summaryPath).First().Split('\t');
            Assert.That(header.Take(18), Is.EqualTo(new[]
            {
                "DatabaseIdentifier", "DatabaseProteinCount", "DatabasePeptideCount", "TotalPredictions",
                "TargetPredictions", "DecoyPredictions", "UniquePeptidesMapped", "UniqueProteinsMapped",
                "MeanRtError", "MedianRtError", "StdDevRtError", "MeanPredictionScore", "MedianPredictionScore",
                "StdDevPredictionScore", "RetentionTimeErrors", "PredictionScores", "TargetPredictionScores",
                "DecoyPredictionScores",
            }));
            Assert.That(header.Skip(18), Is.EqualTo(new[] { "RunOne", "RunTwo" }));

            var file = new DeNovoMappingResultFile(summaryPath);
            file.LoadResults();
            var loaded = file.Results.Single();
            var original = results.Single();

            Assert.Multiple(() =>
            {
                Assert.That(loaded.DatabaseIdentifier, Is.EqualTo("RoundTripDb"));
                Assert.That(loaded.TargetPredictions, Is.EqualTo(original.TargetPredictions));
                Assert.That(loaded.DatabasePeptideCount, Is.EqualTo(original.DatabasePeptideCount));
                Assert.That(loaded.TargetPredictionScores.OrderBy(x => x), Is.EqualTo(original.TargetPredictionScores.OrderBy(x => x)).Within(1e-4));
                Assert.That(loaded.RetentionTimeErrors.OrderBy(x => x), Is.EqualTo(original.RetentionTimeErrors.OrderBy(x => x)).Within(1e-4));
                Assert.That(loaded.MeanRtError, Is.EqualTo(original.MeanRtError).Within(1e-9));
            });
        }

        [Test]
        public void SummaryFileReadsEmptyLists()
        {
            string path = Path.Combine(_directory, "EmptyMappingSummary.tsv");
            DeNovoMapper.WriteSummary([new DeNovoMappingResult { DatabaseIdentifier = "Empty" }], path);

            var file = new DeNovoMappingResultFile(path);
            file.LoadResults();

            Assert.That(file.Results.Single().RetentionTimeErrors, Is.Empty);
            Assert.That(file.Results.Single().MeanPredictionScore, Is.NaN);
            Assert.That(file.Results.Single().DatabaseProteinCount, Is.EqualTo(-1));
        }

        #region Helpers

        private string WriteFasta(string name, params (string Accession, string Sequence)[] proteins)
        {
            string path = Path.Combine(_directory, name);
            var builder = new StringBuilder();
            foreach (var (accession, sequence) in proteins)
                builder.AppendLine($">sp|{accession}|{accession}_TEST Test protein OS=Test OX=1 GN={accession} PE=1 SV=1").AppendLine(sequence);
            File.WriteAllText(path, builder.ToString());
            return path;
        }

        /// <summary>
        /// A peptide that digesting the reversed decoys produces but the targets do not.
        /// </summary>
        private static string FirstDecoyOnlyPeptide(string fasta)
        {
            var parameters = new DeNovoMappingParameters();
            var proteins = ProteinDbLoader.LoadProteinFasta(fasta, true, DecoyType.Reverse, false, out _);
            List<string> Peptides(bool decoy) => proteins.Where(p => p.IsDecoy == decoy)
                .SelectMany(p => p.Digest(parameters.DigestionParams, [], []))
                .Select(p => p.BaseSequence.Replace('I', 'L'))
                .ToList();

            var targets = Peptides(false).ToHashSet();
            return Peptides(true).First(p => !targets.Contains(p));
        }

        private string WriteMzTab(string name, string msRunLocation, params (string Sequence, double Score, double? RetentionTime)[] psms)
        {
            string path = Path.Combine(_directory, name);
            var builder = new StringBuilder();
            builder.AppendLine("MTD\tmzTab-version\t1.0.0");
            builder.AppendLine("MTD\tsoftware[1]\t[MS, MS:1003281, Casanovo, 5.0.0]");
            builder.AppendLine($"MTD\tms_run[1]-location\t{msRunLocation}");
            builder.AppendLine("PSH\tsequence\tPSM_ID\taccession\tunique\tdatabase\tdatabase_version\tsearch_engine\tsearch_engine_score[1]\tmodifications\tretention_time\tcharge\texp_mass_to_charge\tcalc_mass_to_charge\tspectra_ref\tpre\tpost\tstart\tend\topt_ms_run[1]_aa_scores");

            for (int i = 0; i < psms.Length; i++)
            {
                var (sequence, score, rt) = psms[i];
                string rtText = rt?.ToString(CultureInfo.InvariantCulture) ?? "null";
                string aaScores = string.Join(',', Enumerable.Repeat("0.90000", sequence.Length));
                builder.AppendLine($"PSM\t{sequence}\t{i + 1}\tnull\tnull\tnull\tnull\t[MS, MS:1003281, Casanovo, 5.0.0]\t" +
                                   $"{score.ToString(CultureInfo.InvariantCulture)}\tnull\t{rtText}\t2\t500.0\t500.0\tms_run[1]:index={i}\tnull\tnull\tnull\tnull\t{aaScores}");
            }

            File.WriteAllText(path, builder.ToString());
            return path;
        }

        private string WriteMgf(string name, params (int Scan, double RtSeconds)[] spectra)
        {
            string path = Path.Combine(_directory, name);
            var builder = new StringBuilder();
            foreach (var (scan, rtSeconds) in spectra)
            {
                builder.AppendLine("BEGIN IONS");
                builder.AppendLine($"TITLE={Path.GetFileNameWithoutExtension(name)}.{scan}.{scan}.2");
                builder.AppendLine("PEPMASS=500.0");
                builder.AppendLine("CHARGE=2+");
                builder.AppendLine($"SCANS={scan}");
                builder.AppendLine($"RTINSECONDS={rtSeconds.ToString(CultureInfo.InvariantCulture)}");
                builder.AppendLine("200.0 100.0");
                builder.AppendLine("300.0 200.0");
                builder.AppendLine("END IONS");
            }
            File.WriteAllText(path, builder.ToString());
            return path;
        }

        #endregion
    }
}
