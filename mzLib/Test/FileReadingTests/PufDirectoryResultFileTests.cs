using System;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Readers;
using Readers.Puf;

namespace Test.FileReadingTests;

[TestFixture]
public class PufDirectoryResultFileTests
{
    private static readonly string TestDir = Path.Combine(
        TestContext.CurrentContext.TestDirectory,
        @"FileReadingTests\ExternalFileTypes\Puf"
    );

    [Test]
    public void CanReadAllPufFilesInDirectory_IResultFile()
    {
        Assert.That(Directory.Exists(TestDir), $"Test directory not found: {TestDir}");

        var dirFile = new PufDirectoryResultFile(TestDir);
        dirFile.LoadResults();

        // Should find all three files
        Assert.That(dirFile.PufFiles.Count, Is.EqualTo(3));

        // Should aggregate all experiments
        var allExperimentIds = dirFile.Results.Select(e => e.Id).ToList();
        Assert.That(allExperimentIds, Does.Contain("953"));
        Assert.That(allExperimentIds, Does.Contain("955"));
        Assert.That(allExperimentIds, Does.Contain("956"));
    }

    [Test]
    public void CanReadAllPufFilesInDirectory_MsDataFile()
    {
        Assert.That(Directory.Exists(TestDir), $"Test directory not found: {TestDir}");

        var dirFile = MsDataFileReader.GetDataFile(TestDir) as PufMsDataFile;
        dirFile!.LoadAllStaticData();

        // Should find all three files
        Assert.That(dirFile.PufFiles.Count, Is.EqualTo(3));

        // Should aggregate all experiments
        var allExperimentIds = dirFile.Scans.Select(e => e.OneBasedScanNumber).ToList();
        Assert.That(allExperimentIds, Does.Contain(953));
        Assert.That(allExperimentIds, Does.Contain(955));
        Assert.That(allExperimentIds, Does.Contain(956));
    }

    [Test]
    public void CanWriteAndReadBackPufDirectory()
    {
        Assert.That(Directory.Exists(TestDir), $"Test directory not found: {TestDir}");

        var dirFile = new PufDirectoryResultFile(TestDir);
        dirFile.LoadResults();

        var tempDir = Path.Combine(Path.GetTempPath(), Path.GetRandomFileName()).Split('.')[0];
        Directory.CreateDirectory(tempDir);

        try
        {
            // Write all PUF files to a new directory
            dirFile.WriteResults(tempDir);

            // Check that all files were written
            var writtenFiles = Directory.GetFiles(tempDir, "*.puf");
            Assert.That(writtenFiles.Length, Is.EqualTo(3));

            // Read back using a new PufDirectoryResultFile
            var roundTrip = new PufDirectoryResultFile(tempDir);
            roundTrip.LoadResults();

            // Should aggregate all experiments again
            var allExperimentIds = roundTrip.Results.Select(e => e.Id).ToList();
            Assert.That(allExperimentIds, Does.Contain("953"));
            Assert.That(allExperimentIds, Does.Contain("955"));
            Assert.That(allExperimentIds, Does.Contain("956"));
        }
        finally
        {
            if (Directory.Exists(tempDir))
                Directory.Delete(tempDir, true);
        }
    }

    [Test]
    public void CanWritePufDirectoryToMzmlAndReadBack()
    {
        Assert.That(Directory.Exists(TestDir), $"Test directory not found: {TestDir}");

        // Load PUF directory as MsDataFile
        var pufMsDataFile = new PufMsDataFile(TestDir);
        pufMsDataFile.LoadAllStaticData();

        // Write to mzML
        var tempMzml = Path.Combine(Path.GetTempPath(), Path.GetRandomFileName() + ".mzML");
        try
        {
            // Use your mzML writer (replace with your actual method if different)
            pufMsDataFile.ExportAsMzML(tempMzml, true);

            // Read back mzML
            var mzml = MsDataFileReader.GetDataFile(tempMzml);
            mzml.LoadAllStaticData();

            // Assert number of scans matches
            Assert.That(mzml.Scans.Length, Is.EqualTo(pufMsDataFile.Scans.Length));
            Assert.That(mzml.NumSpectra, Is.EqualTo(pufMsDataFile.Scans.Length));

            // Optionally, check a few scan properties
            for (int i = 0; i < Math.Min(3, mzml.Scans.Length); i++)
            {
                var orig = pufMsDataFile.Scans[i];
                var round = mzml.Scans[i];
                Assert.That(round.MassSpectrum.Size, Is.EqualTo(orig.MassSpectrum.Size));
                Assert.That(round.MsnOrder, Is.EqualTo(orig.MsnOrder));
            }
        }
        finally
        {
            if (File.Exists(tempMzml))
                File.Delete(tempMzml);
        }
    }
}
