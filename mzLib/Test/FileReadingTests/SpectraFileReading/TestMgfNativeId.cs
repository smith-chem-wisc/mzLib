using System.Collections.Generic;
using System.IO;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests.SpectraFileReading
{
    /// <summary>
    /// An MGF spectrum's native ID is index=N, N being its block's zero-based position in the file (PSI-MS
    /// MS:1000774, "multiple peak list nativeID format"). That is how mzIdentML and PRIDE name a spectrum in an MGF, and
    /// the scan number cannot stand in for it: a block's scan number is its SCANS value when it has one.
    /// </summary>
    [TestFixture]
    public static class TestMgfNativeId
    {
        // tester.mgf's blocks carry SCANS=14, SCANS=19, none, none and SCANS=2, so the reader numbers them 14, 19, 20, 21
        // and 2. Their TITLEs, written by the tool that made the file, name them NativeID:"index=0" to "index=4".
        private static readonly Dictionary<int, string> TesterNativeIdByScanNumber = new()
        {
            { 14, "index=0" },
            { 19, "index=1" },
            { 20, "index=2" },
            { 21, "index=3" },
            { 2, "index=4" },
        };

        private static string TesterMgf => Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "tester.mgf");

        [Test]
        public static void StaticReadNamesEachSpectrumByItsBlockPosition()
        {
            var file = MsDataFileReader.GetDataFile(TesterMgf);
            file.LoadAllStaticData();

            var nativeIds = file.GetAllScansList().ToDictionary(s => s.OneBasedScanNumber, s => s.NativeId);
            Assert.That(nativeIds, Is.EquivalentTo(TesterNativeIdByScanNumber));

            // the file's own TITLEs name the same spectra the same way, in file order
            var titled = File.ReadLines(TesterMgf).Where(l => l.StartsWith("TITLE="))
                .Select(l => l.Split("NativeID:\"")[1].TrimEnd('"')).ToList();
            Assert.That(titled, Is.EqualTo(new[] { "index=0", "index=1", "index=2", "index=3", "index=4" }));
        }

        [Test]
        public static void DynamicReadNamesEachSpectrumAsTheStaticReadDoes()
        {
            var file = MsDataFileReader.GetDataFile(TesterMgf);
            file.InitiateDynamicConnection();
            try
            {
                foreach (var (scanNumber, nativeId) in TesterNativeIdByScanNumber)
                {
                    Assert.That(file.GetOneBasedScanFromDynamicConnection(scanNumber).NativeId, Is.EqualTo(nativeId));
                }
            }
            finally
            {
                file.CloseDynamicConnection();
            }
        }

        [Test]
        public static void BlockWithNoPeaksStillTakesItsPosition()
        {
            // The reader drops a block with no peaks, but the block is still in the file, so the next one is index=2.
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestMgfNativeId_emptyBlock.mgf");
            File.WriteAllLines(path, new[]
            {
                "BEGIN IONS", "PEPMASS=500.0", "CHARGE=2+", "100.0 10.0", "200.0 20.0", "END IONS",
                "BEGIN IONS", "PEPMASS=600.0", "CHARGE=2+", "END IONS",
                "BEGIN IONS", "PEPMASS=700.0", "CHARGE=2+", "300.0 30.0", "400.0 40.0", "END IONS",
            });
            try
            {
                var staticFile = MsDataFileReader.GetDataFile(path);
                staticFile.LoadAllStaticData();
                Assert.That(staticFile.GetAllScansList().Select(s => s.NativeId), Is.EqualTo(new[] { "index=0", "index=2" }));

                // the dynamic path numbers the scans its own way; their native IDs are the same
                var dynamicFile = MsDataFileReader.GetDataFile(path);
                dynamicFile.InitiateDynamicConnection();
                try
                {
                    var dynamicIds = Enumerable.Range(1, 3)
                        .Select(n => dynamicFile.GetOneBasedScanFromDynamicConnection(n))
                        .Where(s => s != null)
                        .Select(s => s.NativeId);
                    Assert.That(dynamicIds, Is.EqualTo(new[] { "index=0", "index=2" }));
                }
                finally
                {
                    dynamicFile.CloseDynamicConnection();
                }
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public static void SourceFileDescribesTheNativeIds()
        {
            var sourceFile = MsDataFileReader.GetDataFile(TesterMgf).GetSourceFile();
            Assert.That(sourceFile.NativeIdFormat, Is.EqualTo("multiple peak list nativeID format"));
            Assert.That(sourceFile.MassSpectrometerFileFormat, Is.EqualTo("Mascot MGF format"));
        }

        [Test]
        public static void MzmlWrittenFromAnMgfGivesEverySpectrumAnId()
        {
            // mzML's spectrum id is required, and the writer takes it from the native ID, which an MGF spectrum lacked.
            // The writer needs scans numbered without gaps (tester.mgf's are not) and an isolation width on each MS2
            // (mzLib's own MGF writer gives one), so this file has no SCANS lines and does have ISOLATIONWIDTH.
            string mgfPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestMgfNativeId_toMzml.mgf");
            string mzmlPath = Path.ChangeExtension(mgfPath, ".mzML");
            File.WriteAllLines(mgfPath, new[]
            {
                "BEGIN IONS", "PEPMASS=500.0", "CHARGE=2+", "RTINSECONDS=60", "ISOLATIONWIDTH=2.0", "100.0 10.0", "200.0 20.0", "END IONS",
                "BEGIN IONS", "PEPMASS=600.0", "CHARGE=2+", "RTINSECONDS=120", "ISOLATIONWIDTH=2.0", "300.0 30.0", "400.0 40.0", "END IONS",
            });
            try
            {
                var mgf = MsDataFileReader.GetDataFile(mgfPath);
                mgf.LoadAllStaticData();
                mgf.ExportAsMzML(mzmlPath, writeIndexed: false);

                var mzml = MsDataFileReader.GetDataFile(mzmlPath);
                mzml.LoadAllStaticData();
                Assert.That(mzml.GetAllScansList().Select(s => s.NativeId), Is.EqualTo(new[] { "index=0", "index=1" }));
            }
            finally
            {
                File.Delete(mgfPath);
                if (File.Exists(mzmlPath)) File.Delete(mzmlPath);
            }
        }
    }
}
