using System.Diagnostics.CodeAnalysis;
using System.IO;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Tests for SourceFile.InstrumentModel -- the instrument that produced a data file.
    ///
    /// Nothing in mzLib exposed this before, even though mzML has always carried it and
    /// Mzml.GetMsDataOneBasedScanFromConnection was already walking the very element it lives on
    /// (it kept componentList.analyzer[0] for the mass analyzer and discarded the rest).
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestInstrumentModel
    {
        private static string Data(params string[] parts) =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, Path.Combine(parts));

        [Test]
        [TestCase(@"AveragingTests\TestData\TDYeastFractionMS1.mzML", "MS:1002732", "Orbitrap Fusion Lumos")]
        [TestCase(@"DatabaseTests\sliced_b6.mzML", "MS:1002416", "Orbitrap Fusion")]
        [TestCase(@"FlashLFQ\TestData\20100614_Velos1_TaGe_SA_K562_3.mzML", "MS:1001742", "LTQ Orbitrap Velos")]
        public void MzML_ReadsInstrumentModelFromReferenceableParamGroup(
            string relativePath, string expectedAccession, string expectedName)
        {
            // The ProteoWizard layout: the model is hoisted into a referenceableParamGroup and the
            // instrumentConfiguration only points at it. Reading the element's own cvParams alone
            // finds nothing on any real converted file.
            var file = MsDataFileReader.GetDataFile(Data(relativePath));
            var sourceFile = file.GetSourceFile();

            Assert.That(sourceFile.InstrumentModel, Is.Not.Null);
            Assert.That(sourceFile.InstrumentModel.Accession, Is.EqualTo(expectedAccession));
            Assert.That(sourceFile.InstrumentModel.Name, Is.EqualTo(expectedName));
            Assert.That(sourceFile.InstrumentModel.CvLabel, Is.EqualTo("MS"));
        }

        [Test]
        public void MzML_SerialNumberIsNotMistakenForTheModel()
        {
            // MS:1000529 (instrument serial number) sits immediately beside the model in the same
            // group. It is excluded both by name and by the fact that it carries a value.
            var file = MsDataFileReader.GetDataFile(Data(@"AveragingTests\TestData\TDYeastFractionMS1.mzML"));
            var model = file.GetSourceFile().InstrumentModel;

            Assert.That(model.Accession, Is.Not.EqualTo("MS:1000529"));
            Assert.That(model.Value, Is.Empty, "an instrument model is a bare presence flag");
        }

        [Test]
        public void MzML_BareParentTermIsNotReportedAsAModel()
        {
            // mzLib's own writer emits only MS:1000031 ("instrument model"), the PARENT term, with
            // no specific instrument. Reporting that as the model would be worse than reporting
            // nothing: downstream it would look like a real, joinable instrument identity when it
            // says nothing at all.
            var file = MsDataFileReader.GetDataFile(Data(@"DataFiles\SmallCalibratibleYeast.mzml"));
            var sourceFile = file.GetSourceFile();

            Assert.That(sourceFile.InstrumentModel, Is.Null);
        }

        [Test]
        public void MzML_VendorBranchNodeIsNotReportedAsAModel()
        {
            // badScan7192.mzML declares MS:1000483 "Thermo Fisher Scientific instrument model" --
            // the abstract vendor node, valueless, sitting exactly where a real model would sit.
            // The value test alone lets it through, and since the writer change it would then be
            // written into every mzML as though it identified an instrument.
            var file = MsDataFileReader.GetDataFile(Data(@"DataFiles\badScan7192.mzML"));
            var model = file.GetSourceFile().InstrumentModel;

            Assert.That(model, Is.Null,
                "a vendor branch node says which company made the instrument, not which instrument");
        }

        [Test]
        public void InstrumentModel_DefaultsToNullWhenUnset()
        {
            // Init-only and unset: every pre-existing SourceFile construction site is unaffected.
            var sourceFile = new MassSpectrometry.SourceFile(
                "no nativeID format", "mzML format", "abc", "SHA-1", "id");

            Assert.That(sourceFile.InstrumentModel, Is.Null);
        }
    }
}
