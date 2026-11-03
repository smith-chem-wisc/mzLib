using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Security.Cryptography;
using System.Text;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests.SpectraFileReading
{
    /// <summary>
    /// An indexed mzML used to be completed by reading the written file back line by line, serializing the
    /// whole document twice more and hashing it. It is now completed in one read with only the tail rewritten.
    /// These tests hold the new files to the old ones byte for byte, and the one-read scan to the old
    /// line-by-line read, including with read boundaries falling inside lines and line endings.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public sealed class TestMzmlIndexCompletion
    {
        private string _scratchDirectory;

        [SetUp]
        public void SetUp()
        {
            _scratchDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "MzmlIndexCompletion_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(Path.Combine(_scratchDirectory, "onePass"));
            Directory.CreateDirectory(Path.Combine(_scratchDirectory, "threePass"));
        }

        [TearDown]
        public void TearDown()
        {
            Directory.Delete(_scratchDirectory, recursive: true);
        }

        /// <summary>
        /// Writes <paramref name="file"/> both ways under the same file name (the name is written into the
        /// file) and returns the two files' bytes and whether the new way took its one-pass route.
        /// </summary>
        private (byte[] OnePass, byte[] ThreePass, bool TookOnePass) WriteBothWays(MsDataFile file, string fileName)
        {
            string onePass = Path.Combine(_scratchDirectory, "onePass", fileName);
            string threePass = Path.Combine(_scratchDirectory, "threePass", fileName);

            bool tookOnePass = MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(file, onePass, writeIndexed: true,
                completeIndexInOnePass: true);
            bool threePassTookOnePass = MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(file, threePass, writeIndexed: true,
                completeIndexInOnePass: false);
            Assert.That(threePassTookOnePass, Is.False);

            return (File.ReadAllBytes(onePass), File.ReadAllBytes(threePass), tookOnePass);
        }

        // input files with and without chromatograms; the writer always writes a TIC chromatogram
        [TestCase("SmallCalibratibleYeast.mzml")]
        [TestCase("SmallCalibratibleYeast_mixed_bits.mzml")]
        [TestCase("sliced_ethcd.mzML")]
        [TestCase("MS3_TMT11_Mouse_snip.mzML")]
        [TestCase("12-18-17_frac7_calib_ms1_663_665.mzML")]
        [TestCase("FXN11_tr1_032017-calib_ms1_scans716_718.mzML")]
        [TestCase("GUACUG_NegativeMode_Sliced.mzML")]
        [TestCase("badScan7192.mzML")]
        [TestCase("tester.mzML")]
        [TestCase("small.RAW")]
        [TestCase("sliced_ethcd.raw")]
        public void OnePassFileIsByteIdenticalToThreePassFile(string inputFile)
        {
            MsDataFile file = MsDataFileReader.GetDataFile(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", inputFile));
            file.LoadAllStaticData();

            var (onePass, threePass, tookOnePass) = WriteBothWays(file, Path.GetFileNameWithoutExtension(inputFile) + "-calib.mzML");

            Assert.That(tookOnePass, Is.True, "the one-pass route declined an ordinary file");
            Assert.That(onePass, Is.EqualTo(threePass));
        }

        /// <summary>
        /// A file whose scans repeat a native id cannot be indexed: the line-by-line read throws on the repeated
        /// id. The one-pass route must leave it to throw exactly that, and leave the same file behind.
        /// </summary>
        [Test]
        public void RepeatedSpectrumIdFailsTheSameWayBothRoutes()
        {
            MsDataFile file = MsDataFileReader.GetDataFile(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "BinGenerationTest.mzML"));
            file.LoadAllStaticData();
            string onePass = Path.Combine(_scratchDirectory, "onePass", "repeated.mzML");
            string threePass = Path.Combine(_scratchDirectory, "threePass", "repeated.mzML");

            var onePassError = Assert.Throws<ArgumentException>(() => MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(
                file, onePass, writeIndexed: true, completeIndexInOnePass: true));
            var threePassError = Assert.Throws<ArgumentException>(() => MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(
                file, threePass, writeIndexed: true, completeIndexInOnePass: false));

            Assert.That(onePassError.Message, Is.EqualTo(threePassError.Message));
            Assert.That(File.ReadAllBytes(onePass), Is.EqualTo(File.ReadAllBytes(threePass)));
        }

        /// <summary>
        /// A non-ASCII character in the file - here in the source file's name - is where the line reader's
        /// offsets could differ from byte positions, so the one-pass route must decline it and the file must
        /// come out as the three-pass code writes it.
        /// </summary>
        [Test]
        public void FileWithANonAsciiCharacterIsCompletedTheOldWay()
        {
            MsDataFile source = MsDataFileReader.GetDataFile(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "sliced_ethcd.mzML"));
            source.LoadAllStaticData();
            SourceFile original = source.SourceFile;
            var renamed = new SourceFile(original.NativeIdFormat, original.MassSpectrometerFileFormat, original.CheckSum,
                original.FileChecksumType, original.Uri, original.Id, "José_µLC.raw");
            MsDataFile file = new GenericMsDataFile(source.GetAllScansList().ToArray(), renamed);

            var (onePass, threePass, tookOnePass) = WriteBothWays(file, "nonAscii.mzML");

            Assert.That(Encoding.UTF8.GetString(threePass), Does.Contain("José"), "the name did not reach the file");
            Assert.That(tookOnePass, Is.False);
            Assert.That(onePass, Is.EqualTo(threePass));
        }

        /// <summary>
        /// The one-read scan against a verbatim copy of the old line-by-line read, at read sizes that put
        /// boundaries everywhere in the file, on a written mzML and on a synthetic file that mixes "\n", "\r\n"
        /// and "\r" line endings, deep indentation, a line with only whitespace and upper-case tags.
        /// </summary>
        [Test]
        public void ScanFindsTheSameOffsetsAsTheLineReaderAtEveryReadSize()
        {
            MsDataFile source = MsDataFileReader.GetDataFile(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "sliced_ethcd.mzML"));
            source.LoadAllStaticData();
            string written = Path.Combine(_scratchDirectory, "written.mzML");
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(source, written, writeIndexed: true);

            string synthetic = Path.Combine(_scratchDirectory, "synthetic.txt");
            var text = new StringBuilder("﻿<?xml version=\"1.0\"?>\n<indexedmzML>\r\n");
            for (int i = 0; i < 60; i++)
            {
                string indent = new string(' ', (i * 7) % 41);
                string ending = (i % 3) switch { 0 => "\n", 1 => "\r\n", _ => "\r" };
                text.Append(indent).Append(i % 5 == 0 ? "<SPECTRUM " : "<spectrum ")
                    .Append("index=\"").Append(i).Append("\" id=\"scan=").Append(i).Append("\" defaultArrayLength=\"3\">").Append(ending);
                text.Append(indent).Append("  <binary>").Append(new string('A', 50 + 37 * i)).Append("</binary>").Append(ending);
                text.Append(indent).Append(" \t ").Append(ending);
                text.Append(indent).Append("</spectrum>").Append(ending);
            }
            text.Append("  <chromatogram index=\"0\" id=\"TIC\" defaultArrayLength=\"2\">\r\n");
            text.Append("  <spectrumList count=\"60\">\n  <chromatogramList count=\"1\">\n");
            text.Append("  <indexList count=\"2\">\n    <index name=\"spectrum\" />\n  </indexList>\n");
            text.Append("  <indexListOffset xsi:nil=\"true\" />\n  <fileChecksum>0</fileChecksum>\n</indexedmzML>");
            File.WriteAllText(synthetic, text.ToString(), new UTF8Encoding(false));

            foreach (string path in new[] { written, synthetic })
            {
                var (spectra, chromatograms, indexListOffset) = OffsetsAsTheLineReaderTookThem(path);
                Assert.That(indexListOffset, Is.Not.Null);
                byte[] bytes = File.ReadAllBytes(path);
                byte[] prefixHash = SHA1.HashData(bytes.AsSpan(0, (int)indexListOffset.Value));

                IEnumerable<int> readSizes = Enumerable.Range(1, 120).Select(n => 300 + 7 * n)
                    .Concat(new[] { 4096, 4097, 65536, MzmlMethods.ScanChunkBytes });
                foreach (int readSize in readSizes)
                {
                    Assert.That(MzmlMethods.TryScanIndexedFileBody(path, out var scannedSpectra, out var scannedChromatograms,
                        out long scannedIndexListOffset, out IncrementalHash hash, readSize), Is.True, $"{path} at {readSize}");
                    using (hash)
                    {
                        Assert.That(scannedSpectra.ToList(), Is.EqualTo(spectra.ToList()), $"{path} at {readSize}");
                        Assert.That(scannedChromatograms.ToList(), Is.EqualTo(chromatograms.ToList()), $"{path} at {readSize}");
                        Assert.That(scannedIndexListOffset, Is.EqualTo(indexListOffset.Value), $"{path} at {readSize}");
                        Assert.That(hash.GetHashAndReset(), Is.EqualTo(prefixHash), $"{path} at {readSize}");
                    }
                }
            }
        }

        [Test]
        public void ScanDeclinesAFileWithTwoIndexLists()
        {
            string path = Path.Combine(_scratchDirectory, "twoIndexLists.txt");
            File.WriteAllText(path, "<a>\n  <indexList count=\"2\">\n  <indexList count=\"2\">\n</a>");

            Assert.That(MzmlMethods.TryScanIndexedFileBody(path, out _, out _, out _, out IncrementalHash hash,
                MzmlMethods.ScanChunkBytes), Is.False);
            hash.Dispose();
        }

        /// <summary>
        /// The old writer's offset read, verbatim, from "open the initial file and compute byte offset".
        /// </summary>
        private static (Dictionary<string, long> Spectra, Dictionary<string, long> Chromatograms, long? IndexListOffset)
            OffsetsAsTheLineReaderTookThem(string outputFile)
        {
            Dictionary<string, long> spectrumNativeIdToByteOffset = new Dictionary<string, long>();
            Dictionary<string, long> chromatogramIdToByteOffset = new Dictionary<string, long>();
            long? indexListOffset = null;

            using (var stream = new FileStream(outputFile, FileMode.Open, FileAccess.Read, FileShare.Read, 4096, FileOptions.SequentialScan))
            {
                using (StreamReader reader = new StreamReader(stream, bufferSize: 4096))
                {
                    reader.BaseStream.Position = 0;
                    reader.DiscardBufferedData();

                    while (reader.Peek() > 0)
                    {
                        long byteOffset = TextFileReading.GetByteOffsetAtCurrentPosition(reader);
                        var line = reader.ReadLine();

                        for (int i = 0; i < line.Length; i++)
                        {
                            char c = line[i];

                            if (c == ' ')
                            {
                                // spaces are one byte in utf8 encoding
                                byteOffset += 1;
                            }
                            else
                            {
                                break;
                            }
                        }

                        line = line.Trim();

                        if (line.StartsWith("<spectrum ", StringComparison.InvariantCultureIgnoreCase))
                        {
                            string scanId = GetIdFromLine(line);
                            spectrumNativeIdToByteOffset.Add(scanId, byteOffset);
                        }
                        else if (line.StartsWith("<chromatogram ", StringComparison.InvariantCultureIgnoreCase))
                        {
                            string chromatogramId = GetIdFromLine(line);
                            chromatogramIdToByteOffset.Add(chromatogramId, byteOffset);
                        }
                        else if (line.StartsWith("<indexList ", StringComparison.InvariantCultureIgnoreCase))
                        {
                            indexListOffset = byteOffset;
                        }
                    }
                }
            }

            return (spectrumNativeIdToByteOffset, chromatogramIdToByteOffset, indexListOffset);
        }

        /// <summary> MzmlMethods.GetIdFromLine, verbatim. </summary>
        private static string GetIdFromLine(string line)
        {
            int ind = line.IndexOf("id=\"") + 4;
            int len = 0;

            for (int i = ind; i < line.Length; i++)
            {
                if (line[i] == '\"')
                {
                    break;
                }

                len++;
            }

            return line.Substring(ind, len);
        }
    }
}
