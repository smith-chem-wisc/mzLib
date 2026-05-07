using MassSpectrometry;
using MassSpectrometry.MzSpectra;
using Omics.SpectrumMatch;

namespace Omics.SpectralLibrary
{
    public static class LibraryExtensions
    {
        /// <summary>
        /// Computes spectral similarity between two library spectra.
        /// </summary>
        public static SpectralSimilarity ComputeSpectralSimilarity(
            this LibrarySpectrum experimental,
            LibrarySpectrum theoretical,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return new SpectralSimilarity(experimental, theoretical, scheme, toleranceInPpm, allPeaks);
        }

        /// <summary>
        /// Computes spectral similarity between a library spectrum (experimental) and an MzSpectrum (theoretical).
        /// </summary>
        public static SpectralSimilarity ComputeSpectralSimilarity(
            this LibrarySpectrum experimental,
            MzSpectrum theoretical,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return new SpectralSimilarity(experimental, theoretical, scheme, toleranceInPpm, allPeaks);
        }

        /// <summary>
        /// Computes spectral similarity between an MzSpectrum (experimental) and a library spectrum (theoretical).
        /// </summary>
        public static SpectralSimilarity ComputeSpectralSimilarity(
            this MzSpectrum experimental,
            LibrarySpectrum theoretical,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return new SpectralSimilarity(experimental, theoretical, scheme, toleranceInPpm, allPeaks);
        }

        /// <summary>
        /// Computes spectral similarity between a library spectrum (experimental) and explicit m/z and intensity arrays (theoretical).
        /// </summary>
        public static SpectralSimilarity ComputeSpectralSimilarity(
            this LibrarySpectrum experimental,
            double[] theoreticalX,
            double[] theoreticalY,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return new SpectralSimilarity(
                experimental.XArray, experimental.YArray,
                theoreticalX, theoreticalY,
                scheme, toleranceInPpm, allPeaks);
        }

        /// <summary>
        /// Computes spectral similarity between explicit m/z and intensity arrays (experimental) and a library spectrum (theoretical).
        /// </summary>
        public static SpectralSimilarity ComputeSpectralSimilarity(
            this double[] experimentalX,
            double[] experimentalY,
            LibrarySpectrum theoretical,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return new SpectralSimilarity(
                experimentalX, experimentalY,
                theoretical.XArray, theoretical.YArray,
                scheme, toleranceInPpm, allPeaks);
        }

        #region MsDataScan overloads

        /// <summary>
        /// Computes spectral similarity between an MsDataScan (experimental) and a library spectrum (theoretical).
        /// </summary>
        public static SpectralSimilarity ComputeSpectralSimilarity(
            this MsDataScan experimental,
            LibrarySpectrum theoretical,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return new SpectralSimilarity(experimental.MassSpectrum, theoretical, scheme, toleranceInPpm, allPeaks);
        }

        /// <summary>
        /// Computes spectral similarity between an MsDataScan (experimental) and an MzSpectrum (theoretical).
        /// </summary>
        public static SpectralSimilarity ComputeSpectralSimilarity(
            this MsDataScan experimental,
            MzSpectrum theoretical,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return new SpectralSimilarity(experimental.MassSpectrum, theoretical, scheme, toleranceInPpm, allPeaks);
        }

        /// <summary>
        /// Computes spectral similarity between a library spectrum (experimental) and an MsDataScan (theoretical).
        /// </summary>
        public static SpectralSimilarity ComputeSpectralSimilarity(
            this LibrarySpectrum experimental,
            MsDataScan theoretical,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return new SpectralSimilarity(experimental, theoretical.MassSpectrum, scheme, toleranceInPpm, allPeaks);
        }

        /// <summary>
        /// Computes spectral similarity between an MzSpectrum (experimental) and an MsDataScan (theoretical).
        /// </summary>
        public static SpectralSimilarity ComputeSpectralSimilarity(
            this MzSpectrum experimental,
            MsDataScan theoretical,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return new SpectralSimilarity(experimental, theoretical.MassSpectrum, scheme, toleranceInPpm, allPeaks);
        }

        /// <summary>
        /// Computes spectral similarity between an MsDataScan (experimental) and explicit m/z and intensity arrays (theoretical).
        /// </summary>
        public static SpectralSimilarity ComputeSpectralSimilarity(
            this MsDataScan experimental,
            double[] theoreticalX,
            double[] theoreticalY,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return new SpectralSimilarity(
                experimental.MassSpectrum.XArray, experimental.MassSpectrum.YArray,
                theoreticalX, theoreticalY,
                scheme, toleranceInPpm, allPeaks);
        }

        /// <summary>
        /// Computes spectral similarity between explicit m/z and intensity arrays (experimental) and an MsDataScan (theoretical).
        /// </summary>
        public static SpectralSimilarity ComputeSpectralSimilarity(
            this double[] experimentalX,
            double[] experimentalY,
            MsDataScan theoretical,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return new SpectralSimilarity(
                experimentalX, experimentalY,
                theoretical.MassSpectrum.XArray, theoretical.MassSpectrum.YArray,
                scheme, toleranceInPpm, allPeaks);
        }

        #endregion
    }
}
