using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Readers;
using Assert = NUnit.Framework.Legacy.ClassicAssert;

namespace Test.FileReadingTests.SpectraFileReading
{
    /// <summary>
    /// Guards the Thermo CommonCore RawFileReader dependency upgrade (5.0.0.7 -> 8.0.37, net8).
    ///
    /// These tests deliberately exercise every ThermoFisher.CommonCore API surface that mzLib's
    /// <see cref="Readers.ThermoRawFileReader"/> relies on, and assert version-independent invariants,
    /// so that a behavioral regression introduced by this (or any future) CommonCore bump fails loudly
    /// here instead of silently corrupting parsed spectra. API surfaces covered:
    ///   RawFileReaderAdapter.FileFactory, RawFileReaderFactory.CreateThreadManager/CreateThreadAccessor,
    ///   RunHeaderEx, GetFilterForScanNumber (+ MSOrder/Polarity/MassAnalyzer enums),
    ///   GetScanStatsForScanNumber, GetCentroidStream, Scan.FromFile (PreferredMasses fallback),
    ///   GetTrailerExtraInformation, GetScanEventForScanNumber + GetReaction, RetentionTimeFromScanNumber.
    /// </summary>
    [TestFixture]
    public sealed class TestRawFileReaderVersionCompat
    {
        private static string DataPath(string fileName) =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", fileName);

        // Representative Thermo .raw files already committed to the Test project, spanning
        // MS1-only, MS1/MS2, FAIMS (compensation voltage), EThcD, and scan-description data.
        private static readonly string[] RawFiles =
        {
            "small.raw",
            "testFileWMS2.raw",
            "05-13-16_cali_MS_60K-res_MS.raw",
            "sliced_ethcd.raw",
            "ScanDescriptionTestData.raw",
            "TestCompensationVoltageReading.raw",
        };

        /// <summary>
        /// Canary: the CommonCore assemblies actually deployed to the test output must be the
        /// upgraded 8.x build. If the DLL swap did not propagate (stale bin, bad HintPath/copy),
        /// this fails immediately rather than letting the rest of the suite validate the old reader.
        /// Read via the deployed file version so the Test project needs no direct Thermo reference.
        /// </summary>
        [Test]
        [TestCase("ThermoFisher.CommonCore.RawFileReader.dll")]
        [TestCase("ThermoFisher.CommonCore.Data.dll")]
        public void ThermoCommonCoreAssembliesAreUpgradedToMajor8(string dllName)
        {
            string dllPath = Directory
                .GetFiles(TestContext.CurrentContext.TestDirectory, dllName, SearchOption.AllDirectories)
                .FirstOrDefault();

            Assert.IsNotNull(dllPath, $"{dllName} was not deployed to the test output directory.");

            var fileVersion = System.Diagnostics.FileVersionInfo.GetVersionInfo(dllPath);
            Assert.GreaterOrEqual(fileVersion.FileMajorPart, 8,
                $"{dllName} deployed version is {fileVersion.FileVersion}, expected >= 8.0. " +
                "The upgraded Thermo DLL did not reach the output directory.");
        }

        /// <summary>
        /// Static load path: every scan must parse into a coherent MsDataScan. Exercises the
        /// thread-manager accessor, filter-enum mappings (MS order / polarity / mass analyzer),
        /// centroid/PreferredMasses spectrum extraction, and MS2 precursor resolution.
        /// </summary>
        [Test]
        [TestCaseSource(nameof(RawFiles))]
        public void StaticLoad_ProducesCoherentScans(string fileName)
        {
            var reader = MsDataFileReader.GetDataFile(DataPath(fileName));
            reader.LoadAllStaticData(null, maxThreads: 1);
            List<MsDataScan> scans = reader.GetAllScansList();

            Assert.IsTrue(scans.Count > 0, $"{fileName}: no scans loaded");

            foreach (MsDataScan scan in scans)
            {
                // MS order comes straight from the scan filter enum; out-of-range means a broken mapping.
                Assert.IsTrue(scan.MsnOrder >= 1 && scan.MsnOrder <= 10,
                    $"{fileName} scan {scan.OneBasedScanNumber}: implausible MS order {scan.MsnOrder}");

                // Filter string is the backbone of nearly all downstream parsing.
                Assert.IsFalse(string.IsNullOrWhiteSpace(scan.ScanFilter),
                    $"{fileName} scan {scan.OneBasedScanNumber}: empty scan filter");

                // MassAnalyzerType -> MZAnalyzerType: Unknown means the FilterEnums mapping regressed.
                Assert.AreNotEqual(MZAnalyzerType.Unknown, scan.MzAnalyzer,
                    $"{fileName} scan {scan.OneBasedScanNumber}: mass analyzer resolved to Unknown");

                // GetPolarity throws on anything but +/-, so reaching here already proves resolution;
                // assert explicitly to document the invariant.
                Assert.IsTrue(scan.Polarity == Polarity.Positive || scan.Polarity == Polarity.Negative,
                    $"{fileName} scan {scan.OneBasedScanNumber}: unresolved polarity");

                Assert.IsTrue(scan.RetentionTime >= 0,
                    $"{fileName} scan {scan.OneBasedScanNumber}: negative retention time");

                // Centroided reads must never carry zero-intensity peaks (GetSpectrum strips them).
                Assert.IsFalse(scan.MassSpectrum.YArray.Any(y => y == 0),
                    $"{fileName} scan {scan.OneBasedScanNumber}: zero-intensity peak survived centroiding");

                // TIC is the sum of intensities; a populated spectrum must have positive TIC.
                if (scan.MassSpectrum.Size > 0)
                    Assert.Greater(scan.TotalIonCurrent, 0,
                        $"{fileName} scan {scan.OneBasedScanNumber}: non-empty spectrum with non-positive TIC");

                if (scan.MsnOrder > 1)
                {
                    // These come from GetScanEventForScanNumber().GetReaction(0) + the trailer, and from
                    // the precursor-scan back-walk; a null here means MS2 metadata resolution broke.
                    Assert.IsTrue(scan.IsolationMz.HasValue,
                        $"{fileName} scan {scan.OneBasedScanNumber}: MS2 without isolation m/z");
                    Assert.IsTrue(scan.OneBasedPrecursorScanNumber.HasValue,
                        $"{fileName} scan {scan.OneBasedScanNumber}: MS2 without precursor scan linkage");
                    Assert.IsTrue(scan.OneBasedPrecursorScanNumber.Value < scan.OneBasedScanNumber,
                        $"{fileName} scan {scan.OneBasedScanNumber}: precursor scan is not earlier in the run");
                }
            }
        }

        /// <summary>
        /// Trailer-extra parsing (GetTrailerExtraInformation Labels/Values) is the most format-sensitive
        /// surface across CommonCore versions. Assert that it is actually being read: across the corpus,
        /// at least one scan must expose an ion injection time, and at least one MS2 scan must expose a
        /// charge-state or monoisotopic-guess pulled from the trailer.
        /// </summary>
        [Test]
        public void TrailerExtraInformation_IsParsedAcrossCorpus()
        {
            int scansWithInjectionTime = 0;
            int ms2WithTrailerDerivedPrecursorInfo = 0;

            foreach (string fileName in RawFiles)
            {
                var reader = MsDataFileReader.GetDataFile(DataPath(fileName));
                reader.LoadAllStaticData(null, maxThreads: 1);

                foreach (MsDataScan scan in reader.GetAllScansList())
                {
                    if (scan.InjectionTime.HasValue)
                        scansWithInjectionTime++;

                    if (scan.MsnOrder > 1 &&
                        (scan.SelectedIonChargeStateGuess.HasValue || scan.SelectedIonMonoisotopicGuessMz.HasValue))
                        ms2WithTrailerDerivedPrecursorInfo++;
                }
            }

            Assert.Greater(scansWithInjectionTime, 0,
                "No scan across the corpus reported an ion injection time - trailer parsing likely broke.");
            Assert.Greater(ms2WithTrailerDerivedPrecursorInfo, 0,
                "No MS2 scan across the corpus reported a trailer-derived charge state or monoisotopic m/z.");
        }

        /// <summary>
        /// Static vs. dynamic (RawFileReaderAdapter.FileFactory) reads of the same file must agree
        /// scan-for-scan. This is a strong version-independent oracle: the two independent CommonCore
        /// code paths can only match if both parse the file identically.
        /// </summary>
        [Test]
        [TestCase("testFileWMS2.raw")]
        [TestCase("sliced_ethcd.raw")]
        public void StaticAndDynamicReadsAgree(string fileName)
        {
            var staticReader = MsDataFileReader.GetDataFile(DataPath(fileName));
            staticReader.LoadAllStaticData(null, maxThreads: 1);
            staticReader.InitiateDynamicConnection();

            try
            {
                foreach (MsDataScan staticScan in staticReader.GetAllScansList())
                {
                    MsDataScan dynamicScan =
                        staticReader.GetOneBasedScanFromDynamicConnection(staticScan.OneBasedScanNumber);

                    Assert.IsNotNull(dynamicScan,
                        $"{fileName}: dynamic read returned null for scan {staticScan.OneBasedScanNumber}");
                    Assert.AreEqual(staticScan.MsnOrder, dynamicScan.MsnOrder);
                    Assert.AreEqual(staticScan.ScanFilter, dynamicScan.ScanFilter);
                    Assert.AreEqual(staticScan.DissociationType, dynamicScan.DissociationType);
                    Assert.AreEqual(staticScan.MassSpectrum.XArray.Length, dynamicScan.MassSpectrum.XArray.Length,
                        $"{fileName} scan {staticScan.OneBasedScanNumber}: peak count differs static vs dynamic");
                }
            }
            finally
            {
                staticReader.CloseDynamicConnection();
            }
        }
    }
}
