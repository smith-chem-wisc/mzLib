using System.IO;
using System.Linq;
using NUnit.Framework;
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
    public void CanReadAllPufFilesInDirectory()
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
}
