using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    /// <summary>
    /// The instrument model must survive being written back out.
    ///
    /// Before MzmlMethods carried it through, every mzML mzLib wrote declared only the bare parent
    /// term MS:1000031 "instrument model" and nothing else. Calibration and spectral averaging both
    /// rewrite the mzML, so by the time a downstream step looked at the file it no longer said which
    /// mass spectrometer produced the data -- and the loss was invisible, because MS:1000031 looks
    /// like a real declaration rather than an absent one.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestInstrumentModelRoundTrip
    {
        private static string Data(string relativePath) =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, relativePath);

        /// <summary>
        /// A Thermo RAW records its model as free text, so the writer receives it named but not accessioned
        /// (ThermoRawFileReader.BuildInstrumentModel). Written as the bare parent term, every calibrated or averaged
        /// copy of a RAW file lost its instrument: a MetaMorpheus Calibrate -> Search of PXD034564 (TMT10, 32 files)
        /// read its 16 calibrated files back with no instrument, and its SDRF said "not available" for them. The
        /// writer now resolves the name to the PSI-MS instrument model on the way out.
        /// </summary>
        [Test]
        [TestCase(@"DataFiles\testFileWMS2.raw", "MS:1001742", "LTQ Orbitrap Velos")]
        [TestCase(@"DataFiles\sliced_ethcd.raw", "MS:1002416", "Orbitrap Fusion")]
        public void ThermoRawInstrumentModel_SurvivesAWriteAsAnAccessionedTerm(string relativePath,
            string expectedAccession, string expectedName)
        {
            var original = MsDataFileReader.GetDataFile(Data(relativePath));
            original.LoadAllStaticData();
            var source = original.GetSourceFile().InstrumentModel;

            Assert.That(source?.Name, Is.EqualTo(expectedName), "precondition: the RAW file names its instrument");
            Assert.That(source.Accession, Is.Empty, "precondition: a RAW file carries no CV term");

            string outputPath = Path.Combine(
                TestContext.CurrentContext.TestDirectory, $"instrument_raw_{Guid.NewGuid():N}.mzML");
            try
            {
                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(original, outputPath, false);

                var rewritten = MsDataFileReader.GetDataFile(outputPath);
                rewritten.LoadAllStaticData();
                var model = rewritten.GetSourceFile().InstrumentModel;

                Assert.That(model, Is.Not.Null, "the rewritten file lost its instrument model");
                Assert.That(model.Accession, Is.EqualTo(expectedAccession));
                Assert.That(model.Name, Is.EqualTo(expectedName));
            }
            finally
            {
                if (File.Exists(outputPath)) File.Delete(outputPath);
            }
        }

        [Test]
        [TestCase(@"AveragingTests\TestData\TDYeastFractionMS1.mzML", "MS:1002732", "Orbitrap Fusion Lumos")]
        [TestCase(@"DatabaseTests\sliced_b6.mzML", "MS:1002416", "Orbitrap Fusion")]
        public void InstrumentModel_SurvivesAWrite(string relativePath, string expectedAccession, string expectedName)
        {
            var original = MsDataFileReader.GetDataFile(Data(relativePath));
            original.LoadAllStaticData();

            Assert.That(original.GetSourceFile().InstrumentModel?.Accession, Is.EqualTo(expectedAccession),
                "precondition: the source file declares an instrument");

            string outputPath = Path.Combine(
                TestContext.CurrentContext.TestDirectory, $"instrument_rt_{Guid.NewGuid():N}.mzML");
            try
            {
                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(original, outputPath, false);

                var rewritten = MsDataFileReader.GetDataFile(outputPath);
                rewritten.LoadAllStaticData();
                var model = rewritten.GetSourceFile().InstrumentModel;

                Assert.That(model, Is.Not.Null, "the rewritten file lost its instrument model");
                Assert.That(model.Accession, Is.EqualTo(expectedAccession));
                Assert.That(model.Name, Is.EqualTo(expectedName));
            }
            finally
            {
                if (File.Exists(outputPath)) File.Delete(outputPath);
            }
        }

        [Test]
        public void UnknownInstrument_StillWritesTheParentTerm()
        {
            // A file whose instrument we never knew must keep writing MS:1000031, both because that
            // is what mzML expects and because the accession attribute is mandatory -- a name-only
            // cvParam would not validate.
            var original = MsDataFileReader.GetDataFile(Data(@"DataFiles\SmallCalibratibleYeast.mzml"));
            original.LoadAllStaticData();

            Assert.That(original.GetSourceFile().InstrumentModel, Is.Null,
                "precondition: this file declares only the bare parent term");

            string outputPath = Path.Combine(
                TestContext.CurrentContext.TestDirectory, $"instrument_none_{Guid.NewGuid():N}.mzML");
            try
            {
                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(original, outputPath, false);

                string text = File.ReadAllText(outputPath);
                Assert.That(text, Does.Contain("MS:1000031"));

                var rewritten = MsDataFileReader.GetDataFile(outputPath);
                rewritten.LoadAllStaticData();
                Assert.That(rewritten.GetSourceFile().InstrumentModel, Is.Null,
                    "the bare parent term must still not be reported as a real model");
            }
            finally
            {
                if (File.Exists(outputPath)) File.Delete(outputPath);
            }
        }
    }
}
