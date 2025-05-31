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
    public void PufDirectoryResultFile_Scans_AggregateAndOrderCorrectly()
    {
        // Arrange
        var testDir = Path.Combine(TestContext.CurrentContext.TestDirectory, "FileReadingTests", "ExternalFileTypes", "Puf");
        var dirFile = new PufDirectoryResultFile(testDir);
        dirFile.LoadResults();

        // Act
        var scans = dirFile.Scans;

        // Assert
        Assert.That(scans, Is.Not.Null.And.Count.EqualTo(3));
        var scanNumbers = scans.Select(s => s.OneBasedScanNumber).ToList();
        Assert.That(scanNumbers, Is.EquivalentTo(new[] { 953, 955, 956 }));

        // Check precursor m/z and intensity for one scan (target953.puf)
        var scan953 = scans.First(s => s.OneBasedScanNumber == 953);
        Assert.That(scan953.SelectedIonMZ, Is.EqualTo(0).Within(1e-6));
        Assert.That(scan953.SelectedIonIntensity, Is.EqualTo(102053.91).Within(1e-2));
        Assert.That(scan953.ScanDescription, Does.Contain("Characterization score: 365.02"));
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

    [Test]
    public void PufDirectoryResultFile_AggregatesIdentificationsAcrossMultipleFiles()
    {
        // Arrange: directory containing the three PUF files
        var testDir = Path.Combine(TestContext.CurrentContext.TestDirectory, "FileReadingTests", "ExternalFileTypes", "Puf");
        Assert.That(Directory.Exists(testDir), $"Test directory not found: {testDir}");

        // Act: load the directory result file
        var dirResultFile = new Readers.Puf.PufDirectoryResultFile(testDir);
        dirResultFile.LoadResults();

        // Assert: should aggregate all identifications from all files
        var allIds = dirResultFile.Identifications;
        Assert.That(allIds, Is.Not.Null.And.Count.EqualTo(3), "Should have one identification per PUF file");

        // Check that each peptide sequence is as expected (all three files have the same sequence in this test set)
        foreach (var (peptide, matchedIons) in allIds)
        {
            Assert.That(peptide.FullSequence, Does.StartWith("MRHYEIVFLVHPDQSEQV"));
            Assert.That(peptide.Length, Is.EqualTo(139));
            Assert.That(matchedIons, Is.Not.Null.And.Not.Empty);
        }

        // Optionally, check that the union of all matched ion names covers the expected set
        var allIonNames = allIds.SelectMany(t => t.Item2)
            .Select(m => $"{m.NeutralTheoreticalProduct.ProductType}{m.NeutralTheoreticalProduct.FragmentNumber}")
            .ToList();

        // Example: check that B11 and Y138 are present in at least one file
        Assert.That(allIonNames, Does.Contain("b11"));
        Assert.That(allIonNames, Does.Contain("y138"));
    }
}
