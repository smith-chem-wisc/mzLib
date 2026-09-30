using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using FlashLFQ;
using MzLibUtil;
using NUnit.Framework;
using Parquet;
using Parquet.Schema;
using Readers;

namespace Test.FileReadingTests.ExternalFileReading;

/// <summary>
/// The test file is a real DIA-NN 2.3.2 report.parquet. It holds 559 rows (532 targets and 27 decoys)
/// from one run, and it was made in two steps:
/// <list type="number">
/// <item>An in-silico library was predicted from a 202-protein human FASTA.</item>
/// <item>The PXD005573 HeLa run Fig2HeLa-0-5h_MHRM_R01_T0 was searched against that library with
/// <c>--qvalue 0.05 --report-decoys</c>.</item>
/// </list>
/// DIA-NN 2.x writes its main report only as parquet, and several of the 1.8/1.9 columns
/// <see cref="DiaNnPrecursor"/> maps are absent from it.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
internal class TestDiaNnParquetReportFile
{
    private const string ReportPath = @"FileReadingTests\ExternalFileTypes\DiaNn2_report.parquet";
    private const int RowsInFile = 559;
    private const int TargetRows = 532;

    private static string TestFilePath => Path.Combine(TestContext.CurrentContext.TestDirectory, ReportPath);

    private string _outputDirectory = "";

    [OneTimeSetUp]
    public void SetUp()
    {
        _outputDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\DiaNnParquetTests");
        Directory.CreateDirectory(_outputDirectory);
    }

    [OneTimeTearDown]
    public void TearDown()
    {
        Directory.Delete(_outputDirectory, true);
    }

    [Test]
    public void LoadsEveryTargetRow()
    {
        var file = new DiaNnParquetReportFile(TestFilePath);

        Assert.That(file.Count(), Is.EqualTo(TargetRows));
        Assert.That(file.CanRead(TestFilePath));
        Assert.That(file.FileType, Is.EqualTo(SupportedFileType.DiaNnReportParquet));
        Assert.That(file.Software, Is.EqualTo(Software.DiaNn));
    }

    /// <summary>
    /// <see cref="DiaNnPrecursor.IsDecoy"/> is always false, so a decoy row that was loaded would be
    /// quantified as a target by every consumer. The fixture really does contain decoys; the reader
    /// must drop them.
    /// </summary>
    [Test]
    public void DecoyRowsAreNotLoaded()
    {
        long rowsInFile;
        using (var reader = ParquetReader.CreateAsync(TestFilePath).GetAwaiter().GetResult())
            rowsInFile = Enumerable.Range(0, reader.RowGroupCount).Sum(i => reader.OpenRowGroupReader(i).RowCount);
        Assert.That(rowsInFile, Is.EqualTo(RowsInFile), "the fixture must still hold its decoys for this test to mean anything");

        var file = new DiaNnParquetReportFile(TestFilePath);

        Assert.That(file.Count(), Is.EqualTo(TargetRows));
        Assert.That(file.Any(precursor => precursor.PrecursorId == "AGEFLQLHNLR2"), Is.False, "AGEFLQLHNLR2 is a decoy row in the fixture");
        Assert.That(file.All(precursor => !precursor.IsDecoy));
    }

    [Test]
    public void FirstRecordMatchesTheFile()
    {
        DiaNnPrecursor first = new DiaNnParquetReportFile(TestFilePath).First();

        Assert.That(first.PrecursorId, Is.EqualTo("AAAEVAGQFVIK2"));
        Assert.That(first.Run, Is.EqualTo("Fig2HeLa-0-5h_MHRM_R01_T0"));
        Assert.That(first.FileName, Is.EqualTo("Fig2HeLa-0-5h_MHRM_R01_T0"));
        Assert.That(first.ProteinGroup, Is.EqualTo("P02786"));
        Assert.That(first.ProteinIds, Is.EqualTo("P02786"));
        Assert.That(first.ProteinNames, Is.EqualTo("TFR1_HUMAN"));
        Assert.That(first.Genes, Is.EqualTo("TFRC"));
        Assert.That(first.ModifiedSequence, Is.EqualTo("AAAEVAGQFVIK"));
        Assert.That(first.BaseSequence, Is.EqualTo("AAAEVAGQFVIK"));
        Assert.That(first.PrecursorCharge, Is.EqualTo(2));
        Assert.That(first.PrecursorMz, Is.EqualTo(602.3402).Within(1e-4));
        Assert.That(first.Proteotypic, Is.True);

        // DIA-NN stores these as 32-bit floats, so compare at float precision
        Assert.That(first.QValue, Is.EqualTo(0.00013613857).Within(1e-9));
        Assert.That(first.PosteriorErrorProbability, Is.EqualTo(0.001142032).Within(1e-8));
        Assert.That(first.GlobalQValue, Is.EqualTo(0.00013613857).Within(1e-9));
        Assert.That(first.ProteinGroupQValue, Is.EqualTo(0.012974026).Within(1e-8));
        Assert.That(first.PrecursorQuantity, Is.EqualTo(1.3388264e8).Within(10));
        Assert.That(first.ProteinGroupMaxLfq, Is.EqualTo(4.613557e8).Within(100));
        Assert.That(first.GenesMaxLfq, Is.EqualTo(4.613557e8).Within(100));

        Assert.That(first.RetentionTime, Is.EqualTo(28.714846).Within(1e-5));
        Assert.That(first.RetentionTimeStart, Is.EqualTo(28.629341).Within(1e-5));
        Assert.That(first.RetentionTimeStop, Is.EqualTo(28.879116).Within(1e-5));
        Assert.That(first.PredictedRetentionTime, Is.EqualTo(28.801697).Within(1e-5));
        Assert.That(first.IndexedRetentionTime, Is.EqualTo(33.976067).Within(1e-5));
        Assert.That(first.PredictedIndexedRetentionTime, Is.EqualTo(36.654892).Within(1e-5));
        Assert.That(first.IndexedIonMobility, Is.EqualTo(0.9424178).Within(1e-6));
    }

    /// <summary>
    /// DIA-NN 2.x renamed Lib.Index to Precursor.Lib.Index.
    /// </summary>
    [Test]
    public void LibraryIndexIsReadFromItsRenamedColumn()
    {
        DiaNnPrecursor first = new DiaNnParquetReportFile(TestFilePath).First();

        Assert.That(first.LibraryIndex, Is.EqualTo(2));
    }

    /// <summary>
    /// A column DIA-NN 2.x no longer writes is "not reported", not zero. Numbers become NaN or null,
    /// never 0.
    /// </summary>
    [Test]
    public void ColumnsAbsentFromDiaNn2AreNotReportedRatherThanZero()
    {
        DiaNnPrecursor first = new DiaNnParquetReportFile(TestFilePath).First();

        Assert.That(first.SpectraFilePath, Is.Null, "File.Name");
        Assert.That(first.ProteinGroupQuantity, Is.NaN, "PG.Quantity");
        Assert.That(first.ProteinGroupNormalized, Is.NaN, "PG.Normalised");
        Assert.That(first.GenesQuantity, Is.Null, "Genes.Quantity");
        Assert.That(first.GenesNormalized, Is.Null, "Genes.Normalised");
    }

    [Test]
    public void LastRecordMatchesTheFile()
    {
        Assert.That(new DiaNnParquetReportFile(TestFilePath).Last().PrecursorId, Is.EqualTo("YVLPNFEVK2"));
    }

    /// <summary>
    /// DIA-NN lets the user name the report anything, so a parquet is recognized by its columns,
    /// just as the TSV report is recognized by its header.
    /// </summary>
    [Test]
    public void IsDetectedByItsColumnsRegardlessOfFileName()
    {
        string renamedPath = Path.Combine(_outputDirectory, "some_unrelated_name.parquet");
        File.Copy(TestFilePath, renamedPath, true);

        Assert.That(renamedPath.ParseFileType(), Is.EqualTo(SupportedFileType.DiaNnReportParquet));

        IQuantifiableResultFile quantifiable = FileReader.ReadQuantifiableResultFile(renamedPath);
        Assert.That(quantifiable, Is.TypeOf<DiaNnParquetReportFile>());
        Assert.That(quantifiable.GetQuantifiableResults().Count(), Is.EqualTo(TargetRows));
    }

    [Test]
    public void AParquetFromAnotherToolIsRefused()
    {
        string otherPath = Path.Combine(_outputDirectory, "not_diann.parquet");
        var schema = new ParquetSchema(new DataField<int>("id"), new DataField<string>("name"));
        using (var stream = File.Create(otherPath))
        using (var writer = ParquetWriter.CreateAsync(schema, stream).GetAwaiter().GetResult())
        using (var rowGroup = writer.CreateRowGroup())
        {
            rowGroup.WriteAsync(schema.DataFields[0], new[] { 1 }).GetAwaiter().GetResult();
            rowGroup.WriteAsync(schema.DataFields[1], new[] { "a" }).GetAwaiter().GetResult();
        }

        var e = Assert.Throws<MzLibException>(() => otherPath.ParseFileType());
        Assert.That(e!.Message, Is.EqualTo("Parquet file type not supported"));
    }

    [Test]
    public void AFileThatIsNotParquetIsRefusedByName()
    {
        string fakePath = Path.Combine(_outputDirectory, "fake_report.parquet");
        File.WriteAllText(fakePath, "Precursor.Id\tStripped.Sequence\tRun\n");

        var e = Assert.Throws<MzLibException>(() => new DiaNnParquetReportFile(fakePath).LoadResults());
        Assert.That(e!.Message, Does.Contain(fakePath));
    }

    [Test]
    public void LoadingAMissingFileThrowsFileNotFound()
    {
        string missing = Path.Combine(_outputDirectory, "missing_report.parquet");

        Assert.Throws<FileNotFoundException>(() => new DiaNnParquetReportFile(missing).LoadResults());
    }

    [Test]
    public void WritingIsNotSupported()
    {
        var file = new DiaNnParquetReportFile(TestFilePath);

        Assert.That(() => file.WriteResults(Path.Combine(_outputDirectory, "out.parquet")),
            Throws.TypeOf<NotSupportedException>());
    }

    /// <summary>
    /// The factory builds result files with the parameterless constructor.
    /// </summary>
    [Test]
    public void TheFactoryConstructorStillReportsDiaNn()
    {
        var file = new DiaNnParquetReportFile { FilePath = TestFilePath };

        Assert.That(file.Software, Is.EqualTo(Software.DiaNn));
        Assert.That(file.Count(), Is.EqualTo(TargetRows));
    }

    /// <summary>
    /// The full path through FlashLFQ's adapter. Each target row becomes one Identification on the
    /// spectra file its Run names.
    /// </summary>
    [Test]
    public void MakeIdentificationsProducesOnePerTargetRow()
    {
        var file = new DiaNnParquetReportFile(TestFilePath);
        List<SpectraFileInfo> spectraFiles = new()
        {
            new SpectraFileInfo(@"D:\Data\Fig2HeLa-0-5h_MHRM_R01_T0.raw", "HeLa", 0, 0, 0),
        };

        List<Identification> identifications = MzLibExtensions.MakeIdentifications(file, spectraFiles);

        Assert.That(identifications.Count, Is.EqualTo(TargetRows));
        Assert.That(identifications[0].BaseSequence, Is.EqualTo("AAAEVAGQFVIK"));
        Assert.That(identifications[0].Ms2RetentionTimeInMinutes, Is.EqualTo(28.714846).Within(1e-5));
        Assert.That(identifications[0].FileInfo.FullFilePathWithExtension, Is.EqualTo(@"D:\Data\Fig2HeLa-0-5h_MHRM_R01_T0.raw"));
    }
}
