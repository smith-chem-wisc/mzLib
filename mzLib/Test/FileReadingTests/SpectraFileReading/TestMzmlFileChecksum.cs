using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Security.Cryptography;
using System.Text;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests.SpectraFileReading
{
    /// <summary>
    /// The indexedmzML schema defines fileChecksum as "SHA-1 checksum from beginning of file to end
    /// of 'fileChecksum' open tag" -- from byte 0, so the UTF-8 byte-order mark XmlWriter emits counts.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public sealed class TestMzmlFileChecksum
    {
        private const string OpenTag = "<fileChecksum>";
        private const string CloseTag = "</fileChecksum>";

        [Test]
        public static void WrittenFileChecksumCoversTheWholeFileIncludingTheByteOrderMark()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "indexed_checksum.mzML");
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(TwoScanFile(), path, writeIndexed: true);

            AssertWrittenChecksumIsCorrect(path);
        }

        /// <summary>
        /// Same assertion over a file bigger than the hash's read buffer, so the &lt;fileChecksum&gt;
        /// tag is found in a later block than the one the file starts in.
        /// </summary>
        [Test]
        public static void WrittenFileChecksumIsCorrectForAFileLargerThanTheReadBuffer()
        {
            MsDataFile source = MsDataFileReader.GetDataFile(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "SmallCalibratibleYeast.mzml"));
            source.LoadAllStaticData();

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "indexed_checksum_large.mzML");
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(source, path, writeIndexed: true);

            Assert.That(new FileInfo(path).Length, Is.GreaterThan(65536), "file is too small to span read buffers");
            AssertWrittenChecksumIsCorrect(path);
        }

        private static void AssertWrittenChecksumIsCorrect(string path)
        {
            byte[] bytes = File.ReadAllBytes(path);

            // The BOM is what the bug turned on, so fail loudly if a future writer stops emitting one
            // and this test silently stops covering the case.
            Assert.That(bytes.Length, Is.GreaterThan(3));
            Assert.That(new[] { bytes[0], bytes[1], bytes[2] }, Is.EqualTo(Encoding.UTF8.GetPreamble()),
                "the written mzML has no UTF-8 byte-order mark, so this test no longer covers the BOM");

            int openTagStart = IndexOfAscii(bytes, OpenTag);
            Assert.That(openTagStart, Is.GreaterThanOrEqualTo(0), $"no {OpenTag} in the written file");
            int bytesCovered = openTagStart + OpenTag.Length;

            int closeTagStart = IndexOfAscii(bytes, CloseTag);
            Assert.That(closeTagStart, Is.GreaterThan(bytesCovered));
            string written = Encoding.UTF8.GetString(bytes, bytesCovered, closeTagStart - bytesCovered).Trim();

            string expected = Hex(SHA1.HashData(bytes.AsSpan(0, bytesCovered)));
            string withoutBom = Hex(SHA1.HashData(bytes.AsSpan(3, bytesCovered - 3)));

            Assert.That(withoutBom, Is.Not.EqualTo(expected), "SHA-1 collision, or the BOM was not actually skipped");
            Assert.That(written, Is.EqualTo(expected),
                written == withoutBom
                    ? "the written checksum skips the byte-order mark"
                    : "the written checksum covers the wrong byte range");
        }

        /// <summary>
        /// The index offsets are absolute from byte 0 and must stay that way -- they are the reason the
        /// checksum was the only thing wrong about the file.
        /// </summary>
        [Test]
        public static void WrittenIndexOffsetsAreAbsoluteFromByteZero()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "indexed_offsets.mzML");
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(TwoScanFile(), path, writeIndexed: true);

            byte[] bytes = File.ReadAllBytes(path);
            string text = Encoding.UTF8.GetString(bytes);

            long indexListOffset = long.Parse(Between(text, "<indexListOffset>", "</indexListOffset>"));
            Assert.That(AsciiAt(bytes, indexListOffset), Does.StartWith("<indexList "));

            int firstOffsetElement = text.IndexOf("<offset idRef=", StringComparison.Ordinal);
            Assert.That(firstOffsetElement, Is.GreaterThanOrEqualTo(0));
            long firstSpectrumOffset = long.Parse(Between(text.Substring(firstOffsetElement), "\">", "</offset>"));
            Assert.That(AsciiAt(bytes, firstSpectrumOffset), Does.StartWith("<spectrum "));
        }

        private static MsDataFile TwoScanFile()
        {
            MzSpectrum ms1 = new MzSpectrum(new double[] { 400.1, 500.2 }, new double[] { 10, 20 }, false);
            MzSpectrum ms2 = new MzSpectrum(new double[] { 100.3, 200.4 }, new double[] { 30, 40 }, false);

            MsDataScan[] scans =
            {
                new MsDataScan(ms1, 1, 1, true, Polarity.Positive, 1.0, new MzRange(100, 1000), "f",
                    MZAnalyzerType.Orbitrap, ms1.SumOfAllY, null, null, "scan=1"),
                new MsDataScan(ms2, 2, 2, true, Polarity.Positive, 2.0, new MzRange(100, 1000), "f",
                    MZAnalyzerType.Orbitrap, ms2.SumOfAllY, null, null, "scan=2", 500.2, 1, null, 500.2,
                    2.0, DissociationType.HCD, 1, null)
            };

            return new FakeMsDataFile(scans);
        }

        private static string Hex(byte[] hash) => Convert.ToHexString(hash).ToLowerInvariant();

        private static int IndexOfAscii(byte[] haystack, string needle)
        {
            byte[] pattern = Encoding.ASCII.GetBytes(needle);

            for (int i = 0; i <= haystack.Length - pattern.Length; i++)
            {
                int j = 0;
                while (j < pattern.Length && haystack[i + j] == pattern[j])
                {
                    j++;
                }

                if (j == pattern.Length)
                {
                    return i;
                }
            }

            return -1;
        }

        private static string AsciiAt(byte[] bytes, long offset) =>
            Encoding.UTF8.GetString(bytes, (int)offset, Math.Min(40, bytes.Length - (int)offset));

        private static string Between(string text, string start, string end)
        {
            int from = text.IndexOf(start, StringComparison.Ordinal) + start.Length;
            return text.Substring(from, text.IndexOf(end, from, StringComparison.Ordinal) - from);
        }
    }
}
