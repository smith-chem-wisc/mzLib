using System;
using System.IO;
using MassSpectrometry;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests.SpectraFileReading
{
    /// <summary>
    /// sourceFile/@location is xsd:anyURI, but converters put arbitrary strings there, and
    /// new Uri(string) rejects them with two different messages:
    ///
    ///   "the hostname could not be parsed"          -- issue 300, location="file://."
    ///   "the format of the URI could not be determined" -- issue 526, any relative string
    ///
    /// Both are the same unguarded new Uri(simpler.location) call. Mzml.GetSourceFile has parsed
    /// defensively since #1060; nothing pinned that, which is what this fixture is for.
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestMzmlSourceFileLocation
    {
        // The one attribute value in the fixture that GetSourceFile actually reads.
        private const string FixtureLocation =
            "location=\"file:///C:/Users/o310362/Desktop/ok/small_calibratible_yeast\"";

        [Test]
        [TestCase("file://.", TestName = "issue 300: leading dot authority")]
        [TestCase("Company", TestName = "issue 526: bare relative string")]
        [TestCase("../relative/dir", TestName = "issue 526: relative path")]
        [TestCase("file://")]
        [TestCase(@"Z:\some\windows\path")]
        [TestCase("")]
        public static void AMalformedSourceFileLocationDoesNotStopTheFileLoading(string location)
        {
            string template = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "DataFiles", "SmallCalibratibleYeast.mzml");
            string contents = File.ReadAllText(template);
            Assert.That(contents, Does.Contain(FixtureLocation),
                "fixture no longer carries the location attribute this test rewrites");

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "malformed_location_" + Math.Abs(location.GetHashCode()) + ".mzML");
            File.WriteAllText(path, contents.Replace(FixtureLocation, $"location=\"{location}\""));

            try
            {
                MsDataFile data = MsDataFileReader.GetDataFile(path);
                Assert.That(() => data.LoadAllStaticData(), Throws.Nothing);
                Assert.That(data.NumSpectra, Is.GreaterThan(0));

                // Metadata only, but it must still be a usable absolute Uri rather than null.
                Assert.That(data.SourceFile?.Uri, Is.Not.Null);
                Assert.That(data.SourceFile.Uri.IsAbsoluteUri, Is.True);
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [Test]
        [TestCase("file://.")]
        [TestCase("Company")]
        [TestCase("../relative/dir")]
        public static void TheReportedValuesAreStillUnparseableSoTheGuardIsStillNeeded(string location)
        {
            // If a future runtime starts accepting these, the guard becomes belt-and-braces rather
            // than load-bearing, and this test says so out loud instead of quietly passing.
            Assert.That(() => new Uri(location), Throws.TypeOf<UriFormatException>());
        }
    }
}
