using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Security.Cryptography;
using System.Threading.Tasks;
using MassSpectrometry;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests.SpectraFileReading
{
    /// <summary>
    /// The Thermo reader hashes each raw file once per process and serves later loads of the same, unchanged
    /// file from that hash. These tests hold the checksum to the one the reader computed before, byte for
    /// byte, and show that a changed file, a failed hash and concurrent loads are all handled.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public sealed class TestRawFileChecksumCache
    {
        private string _scratchDirectory;

        [SetUp]
        public void SetUp()
        {
            _scratchDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "RawChecksumCache_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(_scratchDirectory);
        }

        [TearDown]
        public void TearDown()
        {
            Directory.Delete(_scratchDirectory, recursive: true);
        }

        /// <summary>
        /// The checksum as ThermoRawFileReader.GetSourceFile computed it before the cache, verbatim.
        /// </summary>
        private static string ChecksumAsBefore(string path)
        {
            string sendCheckSum;
            using (FileStream stream = File.OpenRead(path))
            {
                SHA1 sha = SHA1.Create();
                byte[] checksum = sha.ComputeHash(stream);
                sendCheckSum = BitConverter.ToString(checksum)
                    .Replace("-", string.Empty);
            }
            return sendCheckSum;
        }

        private string PrivateCopy(string fileName)
        {
            string copy = Path.Combine(_scratchDirectory, fileName);
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", fileName), copy);
            return copy;
        }

        [TestCase("small.RAW")]
        [TestCase("testFileWMS2.raw")]
        [TestCase("05-13-16_cali_MS_60K-res_MS.raw")]
        [TestCase("ScanDescriptionTestData.raw")]
        [TestCase("sliced_ethcd.raw")]
        public void ChecksumIsTheSameStringAsBefore(string fileName)
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", fileName);
            string expected = ChecksumAsBefore(path);

            Assert.That(ThermoRawFileReader.ComputeSha1Checksum(path), Is.EqualTo(expected));
            Assert.That(MsDataFileReader.GetDataFile(path).GetSourceFile().CheckSum, Is.EqualTo(expected));
            Assert.That(MsDataFileReader.GetDataFile(path).LoadAllStaticData().SourceFile.CheckSum, Is.EqualTo(expected),
                "the load, which hashes on its background task, reports the same checksum");
        }

        /// <summary>
        /// The second load of an unchanged file must not read it again. Holding the file open with no sharing
        /// makes any read throw, so a checksum coming back proves it came from the cache.
        /// </summary>
        [Test]
        public void SecondChecksumOfAnUnchangedFileDoesNotReadIt()
        {
            string path = PrivateCopy("small.RAW");
            string expected = ChecksumAsBefore(path);
            Assert.That(MsDataFileReader.GetDataFile(path).GetSourceFile().CheckSum, Is.EqualTo(expected));

            using (new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.None))
            {
                Assert.That(() => ThermoRawFileReader.ComputeSha1Checksum(path), Throws.InstanceOf<IOException>(),
                    "the lock does not stop a read, so this test proves nothing");
                Assert.That(MsDataFileReader.GetDataFile(path).GetSourceFile().CheckSum, Is.EqualTo(expected));
            }
        }

        [Test]
        public void FileWhoseLengthChangedIsHashedAgain()
        {
            string path = PrivateCopy("small.RAW");
            string before = ThermoRawFileReader.GetSha1Checksum(path);

            using (var stream = new FileStream(path, FileMode.Append))
            {
                stream.WriteByte(0x42);
            }

            string after = ThermoRawFileReader.GetSha1Checksum(path);
            Assert.That(after, Is.Not.EqualTo(before));
            Assert.That(after, Is.EqualTo(ChecksumAsBefore(path)));
        }

        [Test]
        public void FileRewrittenAtTheSameLengthIsHashedAgainWhenItsWriteTimeChanged()
        {
            string path = PrivateCopy("small.RAW");
            DateTime written = File.GetLastWriteTimeUtc(path);
            string before = ThermoRawFileReader.GetSha1Checksum(path);

            byte[] bytes = File.ReadAllBytes(path);
            bytes[bytes.Length / 2] ^= 0xFF;
            File.WriteAllBytes(path, bytes);
            File.SetLastWriteTimeUtc(path, written.AddMinutes(1));

            string after = ThermoRawFileReader.GetSha1Checksum(path);
            Assert.That(after, Is.Not.EqualTo(before));
            Assert.That(after, Is.EqualTo(ChecksumAsBefore(path)));
        }

        [Test]
        public void FailedHashIsNotKept()
        {
            string path = PrivateCopy("small.RAW");

            using (new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.None))
            {
                Assert.That(() => ThermoRawFileReader.GetSha1Checksum(path), Throws.InstanceOf<IOException>());
            }

            Assert.That(ThermoRawFileReader.GetSha1Checksum(path), Is.EqualTo(ChecksumAsBefore(path)),
                "the failure from the locked read was cached");
        }

        [Test]
        public void ConcurrentChecksumsOfOneFileAgree()
        {
            string path = PrivateCopy("testFileWMS2.raw");
            string expected = ChecksumAsBefore(path);

            string[] checksums = new string[16];
            Parallel.For(0, checksums.Length, new ParallelOptions { MaxDegreeOfParallelism = 8 },
                i => checksums[i] = ThermoRawFileReader.GetSha1Checksum(path));

            Assert.That(checksums.Distinct(), Is.EquivalentTo(new[] { expected }));
        }
    }
}
