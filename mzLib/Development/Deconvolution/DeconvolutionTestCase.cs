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
        public int SelectedIonChargeState { get; init; }
        public double SelectedIonMz { get; init; }
        public MzRange RangeToDeconvolute { get; init; }
        public MzSpectrum SpectrumToDeconvolute { get; init; }
        public PpmTolerance DeconvolutionPPmTolerance { get; init; }
        public DeconvolutionTestCase(SampleType sampleType, string sampleInformation, string testName, string spectrumPath, 
            double mostAbundantMass, int selectedIonChargeState, double selectedIonMz, double precursorPpmMassTolerance)
        {
            SampleType = sampleType;
            SampleInformation = sampleInformation;
            TestName = testName;
            MostAbundantMass = mostAbundantMass;
            SelectedIonChargeState = selectedIonChargeState;
            SelectedIonMz = selectedIonMz;
            SpectrumToDeconvolute = Mzml.LoadAllStaticData(spectrumPath).GetAllScansList().First().MassSpectrum;
            DeconvolutionPPmTolerance = new PpmTolerance(precursorPpmMassTolerance);

            // 8.5 was selected as this is the magic number found in Classic Deconvolution
            RangeToDeconvolute = new MzRange(selectedIonMz - 8.5, selectedIonMz + 8.5);
        }

        public override string ToString()
        {
            return SampleInformation + " " + TestName + $" Charge: {SelectedIonChargeState}";
        }
    }
}
