using MassSpectrometry;
using MzLibUtil;
using Readers;

namespace Development.Deconvolution
{
    /// <summary>
    /// Contains expected results from deconvolution of a whole spectrum as well as information on the sample
    /// </summary>
    public class WholeSpectrumDeconvolutionTestCase
    {
        /// <summary>
        /// Instantiate a WholeSpectrumDeconvolutionTestCase
        /// </summary>
        /// <param name="deconvoluter">The object which will be performing the deconvolution when tested</param>
        /// <param name="sampleInformation">Quick information relevant to the sample, will be visible on test failing. Use this field to allow other to quickly identify which tests are failing</param>
        /// <param name="spectrumPath">path to the spectrum of interest, current implementation expects the file to </param>
        /// <param name="scanNumber">One based scan number of spectrum to deconvolute</param>
        /// <param name="precursorPpmMassTolerance">Tolerance which deconvolution results must match expected value</param>
        /// <param name="expectedMostAbundantObservedIsotopicMasses">Expected masses from deconvolution result</param>
        /// <param name="expectedIonChargeStates">Expected charge states from deconvolution result</param>
        /// <param name="selectedIonMzs">M/z of peaks to deconvolute from spectrum</param>
        public WholeSpectrumDeconvolutionTestCase(Deconvoluter deconvoluter, string sampleInformation, string spectrumPath, int scanNumber,
            double precursorPpmMassTolerance, double[] expectedMostAbundantObservedIsotopicMasses, int[] expectedIonChargeStates, double[] selectedIonMzs)
        {
            if (!new[] { expectedMostAbundantObservedIsotopicMasses.Length, expectedIonChargeStates.Length, selectedIonMzs.Length }.AllSame())
                throw new MzLibException("Must have same number of masses, charges, and mzs");

            Deconvoluter = deconvoluter;
            SampleInformation = sampleInformation;
            ExpectedMostAbundantObservedIsotopicMasses = expectedMostAbundantObservedIsotopicMasses;
            ExpectedIonChargeStates = expectedIonChargeStates;
            SelectedIonMzs = selectedIonMzs;
            DeconvolutionPPmTolerance = new PpmTolerance(precursorPpmMassTolerance);
            SpectrumToDeconvolute = MsDataFileReader.GetDataFile(spectrumPath)
                .LoadAllStaticData()
                .GetAllScansList()
                .First(p => p.OneBasedScanNumber == scanNumber).MassSpectrum;
            Count = selectedIonMzs.Length;
        }

        /// <summary>
        /// The object which will be performing the deconvolution when tested
        /// </summary>
        public Deconvoluter Deconvoluter { get; set; }

        /// <summary>
        /// Quick information relevant to the sample, will be visible on test failing
        /// Use this field to allow other to quickly identify which tests are failing
        /// </summary>
        public string SampleInformation { get; init; }

        /// <summary>
        /// Expected masses from deconvolution result
        /// </summary>
        public double[] ExpectedMostAbundantObservedIsotopicMasses { get; init; }

        /// <summary>
        /// Expected charge states from deconvolution result
        /// </summary>
        public int[] ExpectedIonChargeStates { get; init; }

        /// <summary>
        /// M/z of peaks to deconvolute from spectrum
        /// </summary>
        public double[] SelectedIonMzs { get; init; }

        /// <summary>
        /// Spectrum to Deconvolute
        /// </summary>
        public MzSpectrum SpectrumToDeconvolute { get; init; }

        /// <summary>
        /// Tolerance which deconvolution results must match expected value
        /// </summary>
        public PpmTolerance DeconvolutionPPmTolerance { get; init; }

        /// <summary>
        /// Number of deconvoluted masses that are to be assessed
        /// </summary>
        public double Count { get; init; }

        public override string ToString()
        {
            return $"{Deconvoluter.DeconvolutionType}: {SampleInformation}";
        }
    }
}
