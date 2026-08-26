using System;
using System.Collections.Generic;
using System.IO;
using System.Xml;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests.SpectraFileReading
{
    /// <summary>
    /// run/@id and sourceFile/@id are xs:ID in the mzML 1.1 schema, so their values must be
    /// NCNames. Instrument file names routinely are not: dated names lead with a digit, and
    /// spaces and parentheses are common. mzLib issue 259 is a calibrated mzML rejected by
    /// ProteomeXchange for exactly this.
    ///
    /// VerifyNCName is the same production rule xs:ID validation applies, so asserting on it
    /// checks the reported failure without vendoring the mzML XSD into the test data.
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestMzmlXmlIdValidity
    {
        [Test]
        // The name from issue 259. The hyphens are legal mid-name; the leading digit is what breaks it.
        [TestCase("12-10-16_A17A_yeast_BU_fract5_rep1", "12-10-16_A17A_yeast_BU_fract5_rep1.raw")]
        [TestCase("yeast rep1 (2)-Calibrated", "yeast rep1 (2).raw")]
        // outputBaseName becomes a real path, so it stays legal on Windows; sourceFileName is only
        // metadata and can carry characters (':') no file system would accept.
        [TestCase("-leading-dash", "colon:in:name.raw")]
        public static void WrittenIdsAreValidNCNames(string outputBaseName, string sourceFileName)
        {
            string outPath = Path.Combine(TestContext.CurrentContext.TestDirectory, outputBaseName + ".mzML");
            WriteOneScanMzml(outPath, sourceFileName);

            try
            {
                var ids = ReadIdAttributes(outPath);
                Assert.That(ids.ContainsKey("run"), Is.True, "no run element was written");
                Assert.That(ids.ContainsKey("sourceFile"), Is.True, "no sourceFile element was written");

                foreach (var (element, id) in ids)
                {
                    Assert.That(() => XmlConvert.VerifyNCName(id), Throws.Nothing,
                        $"{element}/@id = '{id}' is not a valid NCName, so the mzML fails xs:ID validation");
                }

                // The sanitised id must not be the only record of the file name.
                Assert.That(ReadSourceFileName(outPath), Is.EqualTo(sourceFileName));

                var reread = MsDataFileReader.GetDataFile(outPath);
                reread.LoadAllStaticData();
                Assert.That(reread.NumSpectra, Is.EqualTo(1));
            }
            finally
            {
                if (File.Exists(outPath)) File.Delete(outPath);
            }
        }

        [Test]
        public static void IdsAreLeftAloneWhenTheFileNameIsAlreadyAnNCName()
        {
            const string outputBaseName = "A17A_yeast_BU_fract5_rep1-Calibrated";
            const string sourceFileName = "A17A_yeast_BU_fract5_rep1.raw";

            string outPath = Path.Combine(TestContext.CurrentContext.TestDirectory, outputBaseName + ".mzML");
            WriteOneScanMzml(outPath, sourceFileName);

            try
            {
                var ids = ReadIdAttributes(outPath);
                Assert.That(ids["run"], Is.EqualTo(outputBaseName));
                Assert.That(ids["sourceFile"], Is.EqualTo(sourceFileName));
            }
            finally
            {
                if (File.Exists(outPath)) File.Delete(outPath);
            }
        }

        [Test]
        public static void NoSourceFileListIsWrittenWhenItCannotBePopulated()
        {
            // fileDescription/sourceFileList is optional, but once written the schema requires a
            // count attribute and at least one sourceFile child. The writer only has the values for
            // those when all three of NativeIdFormat, MassSpectrometerFileFormat and
            // FileChecksumType are known, so with any of them missing the element must be omitted
            // rather than emitted empty.
            var sourceFile = new SourceFile(
                nativeIdFormat: "scan number only nativeID format",
                massSpectrometerFileFormat: null,
                checkSum: null,
                fileChecksumType: null,
                uri: new Uri("file:///test"),
                id: "test",
                fileName: "test.mzML");

            string outPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "no_source_file_list.mzML");
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(
                new GenericMsDataFile(new[] { OneScan() }, sourceFile), outPath, false);

            try
            {
                using var reader = XmlReader.Create(outPath);
                while (reader.Read())
                {
                    if (reader.NodeType == XmlNodeType.Element && reader.LocalName == "sourceFileList")
                    {
                        Assert.Fail("an empty sourceFileList was written, which is schema-invalid");
                    }
                }
            }
            finally
            {
                if (File.Exists(outPath)) File.Delete(outPath);
            }
        }

        [Test]
        public static void ASourceFileWithNoUriStillWrites()
        {
            // The SourceFile constructor that takes neither a path nor a Uri leaves Uri null, and
            // the writer has to put something in the required location attribute.
            var sourceFile = new SourceFile(
                nativeIdFormat: "scan number only nativeID format",
                massSpectrometerFileFormat: "Thermo RAW format",
                checkSum: null,
                fileChecksumType: "SHA-1",
                id: "test");
            Assert.That(sourceFile.Uri, Is.Null, "this test is pointless if the constructor sets a Uri");

            string outPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "no_source_uri.mzML");

            try
            {
                Assert.That(() => MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(
                    new GenericMsDataFile(new[] { OneScan() }, sourceFile), outPath, false), Throws.Nothing);

                // Read the attributes back. Asserting only that a sourceFile id EXISTS let this test
                // pass against a bare Uri?.ToString(): a null string attribute is omitted rather than
                // written, and Mzml.GetSourceFile treats a missing location as absent and substitutes
                // the opened path -- so NumSpectra was still 1 and the half-fix shipped green.
                var sourceFileElement = ReadSourceFileAttributes(outPath);

                Assert.Multiple(() =>
                {
                    // Both are use="required" in the mzML schema, so neither may be omitted.
                    Assert.That(sourceFileElement.TryGetValue("location", out string location), Is.True,
                        "sourceFile/@location was omitted, which is schema-invalid");
                    Assert.That(location, Is.EqualTo(SourceFile.UnknownLocation));

                    Assert.That(sourceFileElement.TryGetValue("name", out string name), Is.True,
                        "sourceFile/@name was omitted, which is schema-invalid");
                    Assert.That(name, Is.EqualTo(SourceFile.UnknownName));

                    // The only test that exercises the null-FileName branch of ToValidXmlId.
                    Assert.That(sourceFileElement.TryGetValue("id", out string id), Is.True);
                    Assert.That(() => XmlConvert.VerifyNCName(id), Throws.Nothing,
                        $"sourceFile/@id '{id}' is not a valid NCName");
                });

                var reread = MsDataFileReader.GetDataFile(outPath);
                reread.LoadAllStaticData();
                Assert.That(reread.NumSpectra, Is.EqualTo(1));
            }
            finally
            {
                if (File.Exists(outPath)) File.Delete(outPath);
            }
        }

        /// <summary>
        /// Every attribute of the written <c>sourceFile</c> element. Distinct from
        /// <c>ReadIdAttributes</c>, which only collects <c>@id</c> and so cannot see an omitted
        /// required attribute.
        /// </summary>
        private static Dictionary<string, string> ReadSourceFileAttributes(string mzmlPath)
        {
            using var reader = XmlReader.Create(mzmlPath);
            while (reader.Read())
            {
                if (reader.NodeType != XmlNodeType.Element || reader.LocalName != "sourceFile")
                    continue;

                var attributes = new Dictionary<string, string>();
                if (reader.MoveToFirstAttribute())
                {
                    do { attributes[reader.LocalName] = reader.Value; } while (reader.MoveToNextAttribute());
                }
                return attributes;
            }

            return new Dictionary<string, string>();
        }

        private static void WriteOneScanMzml(string outPath, string sourceFileName)
        {
            // All three of these must be non-null for the writer to emit a sourceFileList at all.
            var sourceFile = new SourceFile(
                nativeIdFormat: "scan number only nativeID format",
                massSpectrometerFileFormat: "Thermo RAW format",
                checkSum: null,
                fileChecksumType: "SHA-1",
                uri: new Uri("file:///test"),
                id: "test",
                fileName: sourceFileName);

            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(
                new GenericMsDataFile(new[] { OneScan() }, sourceFile), outPath, false);
        }

        private static MsDataScan OneScan()
        {
            return new MsDataScan(
                new MzSpectrum(new[] { 100.0, 200.0 }, new[] { 1000.0, 2000.0 }, shouldCopy: false),
                oneBasedScanNumber: 1,
                msnOrder: 1,
                isCentroid: true,
                polarity: Polarity.Positive,
                retentionTime: 1.0,
                scanWindowRange: new MzRange(50, 1000),
                scanFilter: "FTMS + p NSI Full ms",
                mzAnalyzer: MZAnalyzerType.Orbitrap,
                totalIonCurrent: 3000,
                injectionTime: null,
                noiseData: null,
                nativeId: "scan=1");
        }

        private static Dictionary<string, string> ReadIdAttributes(string path)
        {
            var ids = new Dictionary<string, string>();
            using var reader = XmlReader.Create(path);
            while (reader.Read())
            {
                if (reader.NodeType != XmlNodeType.Element) continue;
                if (reader.LocalName is not ("run" or "sourceFile")) continue;

                string id = reader.GetAttribute("id");
                if (id != null) ids[reader.LocalName] = id;
            }
            return ids;
        }

        private static string ReadSourceFileName(string path)
        {
            using var reader = XmlReader.Create(path);
            while (reader.Read())
            {
                if (reader.NodeType == XmlNodeType.Element && reader.LocalName == "sourceFile")
                {
                    return reader.GetAttribute("name");
                }
            }
            return null;
        }
    }
}
