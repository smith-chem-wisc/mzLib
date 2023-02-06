using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;

namespace Development.Deconvolution
{
    public enum SampleType
    {
        TopDown,
        BottomUp,
    }

    [ExcludeFromCodeCoverage]
    public class DeconvolutionTestCase
    {
        public SampleType SampleType { get; set; }
        public string SampleInformation { get; init; }
        public string TestName { get; init; }
        public double MostAbundantMass { get; init; }
        public double SelectedIonChargeState { get; init; }
        public double SelectedIonMz { get; init; }
        public MzRange RangeToDeconvolute { get; init; }
        public MzSpectrum SpectrumToDeconvolute { get; init; }
        public DeconvolutionTestCase(SampleType sampleType, string sampleInformation, string testName, double mostAbundantMass, double selectedIonChargeState,
            double selectedIonMz, MzSpectrum spectrumToDeconvolute)
        {
            SampleType = sampleType;
            SampleInformation = sampleInformation;
            TestName = testName;
            MostAbundantMass = mostAbundantMass;
            SelectedIonChargeState = selectedIonChargeState;
            SelectedIonMz = selectedIonMz;
            SpectrumToDeconvolute = spectrumToDeconvolute;

            // 8.5 was selected as this is the magic number found in Classic Deconvolution
            RangeToDeconvolute = new MzRange(selectedIonMz - 8.5, selectedIonMz + 8.5);
        }

        public override string ToString()
        {
            return SampleInformation + " " + TestName + $"Charge: {SelectedIonChargeState}";
        }

        /// <summary>
        /// Generates a test case for each charge state and mz inputted
        /// </summary>
        /// <param name="sampleType">top down or bottom up?</param>
        /// <param name="sampleInformation">How was this data acquired?</param>
        /// <param name="testName">Sample specific information</param>
        /// <param name="monoIsotopicMass">mono mass of deconvoluted protein/peptide</param>
        /// <param name="selectedIonChargeStates">charge state of selected ion</param>
        /// <param name="selectedIonMzs">mz of selected on</param>
        /// <param name="fileToTest">mass spec file to test</param>
        /// <returns></returns>
        public static IEnumerable<DeconvolutionTestCase> GenerateTestCases(SampleType sampleType, string sampleInformation, string testName, double mostAbundantMass, double[] selectedIonChargeStates,
            double[] selectedIonMzs, string fileToTest)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Deconvolution", "TestData", fileToTest);
            MzSpectrum spectrum = Mzml.LoadAllStaticData(filePath).GetAllScansList().First().MassSpectrum;
            for (int i = 0; i < selectedIonMzs.Length; i++)
            {
                yield return new DeconvolutionTestCase(sampleType, sampleInformation, testName, mostAbundantMass, selectedIonChargeStates[i], selectedIonMzs[i], spectrum);
            }
        }
    }
}
