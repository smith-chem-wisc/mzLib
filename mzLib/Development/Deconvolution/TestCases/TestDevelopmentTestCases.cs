using MzLibUtil;
using NUnit.Framework;
using MassSpectrometry;
using System.Diagnostics.CodeAnalysis;
using Readers;

namespace Development.Deconvolution

{
    [TestFixture]
    [Ignore("Only needed when developing deconvolution methods")]
    [ExcludeFromCodeCoverage]
    public static class TestDevelopmentTestCases
    {
        [Test]
        public static void TestSinglePeakDeconvolutionTestCase()
        {
            Deconvoluter classicTopDownDeconvoluter = new Deconvoluter(DeconvolutionType.ClassicDeconvolution,
                new ClassicDeconvolutionParameters(1, 60, 4, 3));
            const string sampleInformation = "Direct Injection Cytochrome C, Averaged";
            var pathToDataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, "Deconvolution", "TestData",
                "Averaged_221110_CytoOnly.mzML");
            const int scanNumber = 1;
            const double expectedMostAbundantObservedIsotopicMass = 12367.44;
            const int expectedIonChargeState = 9;
            const double selectedIonMz = 1374.16;
            const int precursorPpmTolerance = 20;
            var range = new MzRange(selectedIonMz - 8.5, selectedIonMz + 8.5);
            var spectrum = MsDataFileReader.GetDataFile(pathToDataFile)
                .LoadAllStaticData()
                .GetAllScansList()
                .First(p => p.OneBasedScanNumber == scanNumber).MassSpectrum;

            var testCase = new SinglePeakDeconvolutionTestCase(classicTopDownDeconvoluter, sampleInformation, pathToDataFile,
                scanNumber, expectedMostAbundantObservedIsotopicMass, expectedIonChargeState, selectedIonMz,
                precursorPpmTolerance);

            Assert.IsNotNull(testCase);
            Assert.That(testCase.SampleInformation, Is.EqualTo(sampleInformation));
            Assert.That(testCase.ExpectedMostAbundantObservedIsotopicMass, Is.EqualTo(expectedMostAbundantObservedIsotopicMass));
            Assert.That(testCase.ExpectedIonChargeState, Is.EqualTo(expectedIonChargeState));
            Assert.That(testCase.SelectedIonMz, Is.EqualTo(selectedIonMz));
            Assert.That(testCase.RangeToDeconvolute.Width, Is.EqualTo(range.Width));
            Assert.That(testCase.RangeToDeconvolute.Maximum, Is.EqualTo(range.Maximum));
            Assert.That(testCase.RangeToDeconvolute.Minimum, Is.EqualTo(range.Minimum));
            Assert.That(testCase.SpectrumToDeconvolute, Is.EqualTo(spectrum));
            Assert.That(testCase.DeconvolutionPPmTolerance.Value, Is.EqualTo(precursorPpmTolerance));
            Assert.That(testCase.ToString(), Is.EqualTo($"{DeconvolutionType.ClassicDeconvolution}: {sampleInformation} Charge: {expectedIonChargeState}"));
        }

        [Test]
        public static void TestWholeSpectrumDeconvolutionTestCase()
        {
            Deconvoluter classicTopDownDeconvoluter = new Deconvoluter(DeconvolutionType.ClassicDeconvolution,
                new ClassicDeconvolutionParameters(1, 60, 4, 3));
            const string sampleInformation = "Direct Injection Cytochrome C, Averaged";
            var pathToDataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, "Deconvolution", "TestData",
                "Averaged_221110_CytoOnly.mzML");
            const int scanNumber = 1;
            double[] expectedMostAbundantObservedIsotopicMass = new[] { 12367.44, };
            int[] expectedIonChargeState = new[] { 9 };
            double[] selectedIonMz = new[] { 1374.16 };
            const int precursorPpmTolerance = 20;
            var spectrum = MsDataFileReader.GetDataFile(pathToDataFile)
                .LoadAllStaticData()
                .GetAllScansList()
                .First(p => p.OneBasedScanNumber == scanNumber).MassSpectrum;

            var testCase = new WholeSpectrumDeconvolutionTestCase(classicTopDownDeconvoluter, sampleInformation, pathToDataFile, scanNumber,
                precursorPpmTolerance, expectedMostAbundantObservedIsotopicMass, expectedIonChargeState, selectedIonMz);

            Assert.IsNotNull(testCase);
            Assert.That(testCase.SampleInformation, Is.EqualTo(sampleInformation));
            Assert.That(testCase.ExpectedMostAbundantObservedIsotopicMasses, Is.EqualTo(expectedMostAbundantObservedIsotopicMass));
            Assert.That(testCase.ExpectedIonChargeStates, Is.EqualTo(expectedIonChargeState));
            Assert.That(testCase.SelectedIonMzs, Is.EqualTo(selectedIonMz));
            Assert.That(testCase.SpectrumToDeconvolute, Is.EqualTo(spectrum));
            Assert.That(testCase.DeconvolutionPPmTolerance.Value, Is.EqualTo(precursorPpmTolerance));
            Assert.That(testCase.Count, Is.EqualTo(1));
            Assert.That(testCase.ToString(), Is.EqualTo($"{DeconvolutionType.ClassicDeconvolution}: {sampleInformation}"));

            selectedIonMz = new[] { 2.0, 15.3 };
            MzLibException exception = Assert.Throws<MzLibException>(() => new WholeSpectrumDeconvolutionTestCase(classicTopDownDeconvoluter, sampleInformation,
                pathToDataFile, scanNumber,
                precursorPpmTolerance, expectedMostAbundantObservedIsotopicMass, expectedIonChargeState,
                selectedIonMz));
            Assert.That(exception!.Message, Is.EqualTo("Must have same number of masses, charges, and mzs"));
        }
    }
}
