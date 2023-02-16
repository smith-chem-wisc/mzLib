using System.Diagnostics.CodeAnalysis;
using MassSpectrometry;
using MzLibUtil;
using SpectralAveraging;

namespace Development.Deconvolution
{
    public enum SampleType
    {
        TopDown,
        BottomUp
    }

    /// <summary>
    /// Contains the expected results from deconvolution a single isotopic peak as well as information on the sample
    /// </summary>
    public class SinglePeakDeconvolutionTestCase
    {
        /// <summary>
        /// Instantiate a SinglePeakDeconvolutionTestCase
        /// </summary>
        /// <param name="sampleType">Type of sample being deconvoluted, can be used to set deconvolution parameters</param>
        /// <param name="sampleInformation">Quick information relevant to the sample, will be visible on test failing. Use this field to allow other to quickly identify which tests are failing</param>
        /// <param name="spectrumPath">path to the spectrum of interest, current implementation expects the file to </param>
        /// <param name="scanNumber">One based scan number of spectrum to deconvolute</param>
        /// <param name="expectedMostAbundantObservedIsotopicMass">Expected mass from deconvolution result</param>
        /// <param name="expectedIonChargeState">Expected charge state from deconvolution result</param>
        /// <param name="selectedIonMz">M/z of peak to deconvolute from spectrum</param>
        /// <param name="precursorPpmMassTolerance">Tolerance which deconvolution results must match expected value</param>
        public SinglePeakDeconvolutionTestCase(SampleType sampleType, string sampleInformation, string spectrumPath, int scanNumber,
            double expectedMostAbundantObservedIsotopicMass, int expectedIonChargeState, double selectedIonMz, double precursorPpmMassTolerance)
        {
            SampleType = sampleType;
            SampleInformation = sampleInformation;
            ExpectedMostAbundantObservedIsotopicMass = expectedMostAbundantObservedIsotopicMass;
            ExpectedIonChargeState = expectedIonChargeState;
            SelectedIonMz = selectedIonMz;
            DeconvolutionPPmTolerance = new PpmTolerance(precursorPpmMassTolerance);
            SpectrumToDeconvolute = SpectraFileHandler.LoadAllScansFromFile(spectrumPath)
                .First(p => p.OneBasedScanNumber == scanNumber).MassSpectrum;

            // 8.5 was selected as this is the magic number found in Classic Deconvolution
            RangeToDeconvolute = new MzRange(selectedIonMz - 8.5, selectedIonMz + 8.5);
        }

        /// <summary>
        /// Type of sample being deconvoluted, can be used to set deconvolution parameters
        /// </summary>
        public SampleType SampleType { get; set; }

        /// <summary>
        /// Quick information relevant to the sample, will be visible on test failing
        /// Use this field to allow other to quickly identify which tests are failing
        /// </summary>
        public string SampleInformation { get; init; }

        /// <summary>
        /// Expected mass from deconvolution result
        /// </summary>
        public double ExpectedMostAbundantObservedIsotopicMass { get; init; }

        /// <summary>
        /// Expected charge state from deconvolution result
        /// </summary>
        public int ExpectedIonChargeState { get; init; }

        /// <summary>
        /// M/z of peak to deconvolute from spectrum
        /// </summary>
        public double SelectedIonMz { get; init; }

        /// <summary>
        /// Range within the spectrum to deconvolute
        /// Currently set to +- 8.5 mz from SelectedIonMz per MetaMorpheus default heuristic
        /// </summary>
        public MzRange RangeToDeconvolute { get; init; }

        /// <summary>
        /// Spectrum to Deconvolute
        /// </summary>
        public MzSpectrum SpectrumToDeconvolute { get; init; }

        /// <summary>
        /// Tolerance which deconvolution results must match expected value
        /// </summary>
        public PpmTolerance DeconvolutionPPmTolerance { get; init; }

        public override string ToString()
        {
            return SampleInformation + $" Charge: {ExpectedIonChargeState}";
        }
    }
}

