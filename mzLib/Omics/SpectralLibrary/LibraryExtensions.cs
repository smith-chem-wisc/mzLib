using MassSpectrometry;
using MassSpectrometry.MzSpectra;
using Omics.SpectrumMatch;
using System.Collections.Generic;

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

        #region GetSimilarityMeasure overloads

        /// <summary>
        /// Computes spectral similarity between two library spectra and returns the requested measure.
        /// </summary>
        public static double? GetSimilarityMeasure(
            this LibrarySpectrum experimental,
            LibrarySpectrum theoretical,
            SpectralSimilarity.SimilarityMeasures measure,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoretical, scheme, toleranceInPpm, allPeaks)
                .GetSimilarityMeasure(measure);
        }

        /// <summary>
        /// Computes spectral similarity between a library spectrum and an MzSpectrum and returns the requested measure.
        /// </summary>
        public static double? GetSimilarityMeasure(
            this LibrarySpectrum experimental,
            MzSpectrum theoretical,
            SpectralSimilarity.SimilarityMeasures measure,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoretical, scheme, toleranceInPpm, allPeaks)
                .GetSimilarityMeasure(measure);
        }

        /// <summary>
        /// Computes spectral similarity between an MzSpectrum and a library spectrum and returns the requested measure.
        /// </summary>
        public static double? GetSimilarityMeasure(
            this MzSpectrum experimental,
            LibrarySpectrum theoretical,
            SpectralSimilarity.SimilarityMeasures measure,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoretical, scheme, toleranceInPpm, allPeaks)
                .GetSimilarityMeasure(measure);
        }

        /// <summary>
        /// Computes spectral similarity between a library spectrum and explicit m/z/intensity arrays and returns the requested measure.
        /// </summary>
        public static double? GetSimilarityMeasure(
            this LibrarySpectrum experimental,
            double[] theoreticalX,
            double[] theoreticalY,
            SpectralSimilarity.SimilarityMeasures measure,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoreticalX, theoreticalY, scheme, toleranceInPpm, allPeaks)
                .GetSimilarityMeasure(measure);
        }

        /// <summary>
        /// Computes spectral similarity between explicit m/z/intensity arrays and a library spectrum and returns the requested measure.
        /// </summary>
        public static double? GetSimilarityMeasure(
            this double[] experimentalX,
            double[] experimentalY,
            LibrarySpectrum theoretical,
            SpectralSimilarity.SimilarityMeasures measure,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimentalX.ComputeSpectralSimilarity(experimentalY, theoretical, scheme, toleranceInPpm, allPeaks)
                .GetSimilarityMeasure(measure);
        }

        /// <summary>
        /// Computes spectral similarity between an MsDataScan and a library spectrum and returns the requested measure.
        /// </summary>
        public static double? GetSimilarityMeasure(
            this MsDataScan experimental,
            LibrarySpectrum theoretical,
            SpectralSimilarity.SimilarityMeasures measure,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoretical, scheme, toleranceInPpm, allPeaks)
                .GetSimilarityMeasure(measure);
        }

        /// <summary>
        /// Computes spectral similarity between an MsDataScan and an MzSpectrum and returns the requested measure.
        /// </summary>
        public static double? GetSimilarityMeasure(
            this MsDataScan experimental,
            MzSpectrum theoretical,
            SpectralSimilarity.SimilarityMeasures measure,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoretical, scheme, toleranceInPpm, allPeaks)
                .GetSimilarityMeasure(measure);
        }

        /// <summary>
        /// Computes spectral similarity between a library spectrum and an MsDataScan and returns the requested measure.
        /// </summary>
        public static double? GetSimilarityMeasure(
            this LibrarySpectrum experimental,
            MsDataScan theoretical,
            SpectralSimilarity.SimilarityMeasures measure,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoretical, scheme, toleranceInPpm, allPeaks)
                .GetSimilarityMeasure(measure);
        }

        /// <summary>
        /// Computes spectral similarity between an MzSpectrum and an MsDataScan and returns the requested measure.
        /// </summary>
        public static double? GetSimilarityMeasure(
            this MzSpectrum experimental,
            MsDataScan theoretical,
            SpectralSimilarity.SimilarityMeasures measure,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoretical, scheme, toleranceInPpm, allPeaks)
                .GetSimilarityMeasure(measure);
        }

        /// <summary>
        /// Computes spectral similarity between an MsDataScan and explicit m/z/intensity arrays and returns the requested measure.
        /// </summary>
        public static double? GetSimilarityMeasure(
            this MsDataScan experimental,
            double[] theoreticalX,
            double[] theoreticalY,
            SpectralSimilarity.SimilarityMeasures measure,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoreticalX, theoreticalY, scheme, toleranceInPpm, allPeaks)
                .GetSimilarityMeasure(measure);
        }

        /// <summary>
        /// Computes spectral similarity between explicit m/z/intensity arrays and an MsDataScan and returns the requested measure.
        /// </summary>
        public static double? GetSimilarityMeasure(
            this double[] experimentalX,
            double[] experimentalY,
            MsDataScan theoretical,
            SpectralSimilarity.SimilarityMeasures measure,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimentalX.ComputeSpectralSimilarity(experimentalY, theoretical, scheme, toleranceInPpm, allPeaks)
                .GetSimilarityMeasure(measure);
        }

        #endregion

        #region GetAllSimilarityMeasures overloads

        /// <summary>
        /// Computes spectral similarity between two library spectra and returns all similarity measures.
        /// </summary>
        public static IEnumerable<(SpectralSimilarity.SimilarityMeasures, double?)> GetAllSimilarityMeasures(
            this LibrarySpectrum experimental,
            LibrarySpectrum theoretical,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoretical, scheme, toleranceInPpm, allPeaks)
                .GetAllSimilarityMeasures();
        }

        /// <summary>
        /// Computes spectral similarity between a library spectrum and an MzSpectrum and returns all similarity measures.
        /// </summary>
        public static IEnumerable<(SpectralSimilarity.SimilarityMeasures, double?)> GetAllSimilarityMeasures(
            this LibrarySpectrum experimental,
            MzSpectrum theoretical,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoretical, scheme, toleranceInPpm, allPeaks)
                .GetAllSimilarityMeasures();
        }

        /// <summary>
        /// Computes spectral similarity between an MzSpectrum and a library spectrum and returns all similarity measures.
        /// </summary>
        public static IEnumerable<(SpectralSimilarity.SimilarityMeasures, double?)> GetAllSimilarityMeasures(
            this MzSpectrum experimental,
            LibrarySpectrum theoretical,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoretical, scheme, toleranceInPpm, allPeaks)
                .GetAllSimilarityMeasures();
        }

        /// <summary>
        /// Computes spectral similarity between a library spectrum and explicit m/z/intensity arrays and returns all similarity measures.
        /// </summary>
        public static IEnumerable<(SpectralSimilarity.SimilarityMeasures, double?)> GetAllSimilarityMeasures(
            this LibrarySpectrum experimental,
            double[] theoreticalX,
            double[] theoreticalY,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoreticalX, theoreticalY, scheme, toleranceInPpm, allPeaks)
                .GetAllSimilarityMeasures();
        }

        /// <summary>
        /// Computes spectral similarity between explicit m/z/intensity arrays and a library spectrum and returns all similarity measures.
        /// </summary>
        public static IEnumerable<(SpectralSimilarity.SimilarityMeasures, double?)> GetAllSimilarityMeasures(
            this double[] experimentalX,
            double[] experimentalY,
            LibrarySpectrum theoretical,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimentalX.ComputeSpectralSimilarity(experimentalY, theoretical, scheme, toleranceInPpm, allPeaks)
                .GetAllSimilarityMeasures();
        }

        /// <summary>
        /// Computes spectral similarity between an MsDataScan and a library spectrum and returns all similarity measures.
        /// </summary>
        public static IEnumerable<(SpectralSimilarity.SimilarityMeasures, double?)> GetAllSimilarityMeasures(
            this MsDataScan experimental,
            LibrarySpectrum theoretical,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoretical, scheme, toleranceInPpm, allPeaks)
                .GetAllSimilarityMeasures();
        }

        /// <summary>
        /// Computes spectral similarity between an MsDataScan and an MzSpectrum and returns all similarity measures.
        /// </summary>
        public static IEnumerable<(SpectralSimilarity.SimilarityMeasures, double?)> GetAllSimilarityMeasures(
            this MsDataScan experimental,
            MzSpectrum theoretical,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoretical, scheme, toleranceInPpm, allPeaks)
                .GetAllSimilarityMeasures();
        }

        /// <summary>
        /// Computes spectral similarity between a library spectrum and an MsDataScan and returns all similarity measures.
        /// </summary>
        public static IEnumerable<(SpectralSimilarity.SimilarityMeasures, double?)> GetAllSimilarityMeasures(
            this LibrarySpectrum experimental,
            MsDataScan theoretical,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoretical, scheme, toleranceInPpm, allPeaks)
                .GetAllSimilarityMeasures();
        }

        /// <summary>
        /// Computes spectral similarity between an MzSpectrum and an MsDataScan and returns all similarity measures.
        /// </summary>
        public static IEnumerable<(SpectralSimilarity.SimilarityMeasures, double?)> GetAllSimilarityMeasures(
            this MzSpectrum experimental,
            MsDataScan theoretical,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoretical, scheme, toleranceInPpm, allPeaks)
                .GetAllSimilarityMeasures();
        }

        /// <summary>
        /// Computes spectral similarity between an MsDataScan and explicit m/z/intensity arrays and returns all similarity measures.
        /// </summary>
        public static IEnumerable<(SpectralSimilarity.SimilarityMeasures, double?)> GetAllSimilarityMeasures(
            this MsDataScan experimental,
            double[] theoreticalX,
            double[] theoreticalY,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoreticalX, theoreticalY, scheme, toleranceInPpm, allPeaks)
                .GetAllSimilarityMeasures();
        }

        /// <summary>
        /// Computes spectral similarity between explicit m/z/intensity arrays and an MsDataScan and returns all similarity measures.
        /// </summary>
        public static IEnumerable<(SpectralSimilarity.SimilarityMeasures, double?)> GetAllSimilarityMeasures(
            this double[] experimentalX,
            double[] experimentalY,
            MsDataScan theoretical,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimentalX.ComputeSpectralSimilarity(experimentalY, theoretical, scheme, toleranceInPpm, allPeaks)
                .GetAllSimilarityMeasures();
        }

        #endregion

        #region GetSelectedSimilarityMeasures overloads

        /// <summary>
        /// Computes spectral similarity between two library spectra and returns selected similarity measures.
        /// </summary>
        public static IEnumerable<(SpectralSimilarity.SimilarityMeasures, double?)> GetSelectedSimilarityMeasures(
            this LibrarySpectrum experimental,
            LibrarySpectrum theoretical,
            SpectralSimilarity.SimilarityMeasures[] measures,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoretical, scheme, toleranceInPpm, allPeaks)
                .GetSelectedSimilarityMeasures(measures);
        }

        /// <summary>
        /// Computes spectral similarity between a library spectrum and an MzSpectrum and returns selected similarity measures.
        /// </summary>
        public static IEnumerable<(SpectralSimilarity.SimilarityMeasures, double?)> GetSelectedSimilarityMeasures(
            this LibrarySpectrum experimental,
            MzSpectrum theoretical,
            SpectralSimilarity.SimilarityMeasures[] measures,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoretical, scheme, toleranceInPpm, allPeaks)
                .GetSelectedSimilarityMeasures(measures);
        }

        /// <summary>
        /// Computes spectral similarity between an MzSpectrum and a library spectrum and returns selected similarity measures.
        /// </summary>
        public static IEnumerable<(SpectralSimilarity.SimilarityMeasures, double?)> GetSelectedSimilarityMeasures(
            this MzSpectrum experimental,
            LibrarySpectrum theoretical,
            SpectralSimilarity.SimilarityMeasures[] measures,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoretical, scheme, toleranceInPpm, allPeaks)
                .GetSelectedSimilarityMeasures(measures);
        }

        /// <summary>
        /// Computes spectral similarity between a library spectrum and explicit m/z/intensity arrays and returns selected similarity measures.
        /// </summary>
        public static IEnumerable<(SpectralSimilarity.SimilarityMeasures, double?)> GetSelectedSimilarityMeasures(
            this LibrarySpectrum experimental,
            double[] theoreticalX,
            double[] theoreticalY,
            SpectralSimilarity.SimilarityMeasures[] measures,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoreticalX, theoreticalY, scheme, toleranceInPpm, allPeaks)
                .GetSelectedSimilarityMeasures(measures);
        }

        /// <summary>
        /// Computes spectral similarity between explicit m/z/intensity arrays and a library spectrum and returns selected similarity measures.
        /// </summary>
        public static IEnumerable<(SpectralSimilarity.SimilarityMeasures, double?)> GetSelectedSimilarityMeasures(
            this double[] experimentalX,
            double[] experimentalY,
            LibrarySpectrum theoretical,
            SpectralSimilarity.SimilarityMeasures[] measures,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimentalX.ComputeSpectralSimilarity(experimentalY, theoretical, scheme, toleranceInPpm, allPeaks)
                .GetSelectedSimilarityMeasures(measures);
        }

        /// <summary>
        /// Computes spectral similarity between an MsDataScan and a library spectrum and returns selected similarity measures.
        /// </summary>
        public static IEnumerable<(SpectralSimilarity.SimilarityMeasures, double?)> GetSelectedSimilarityMeasures(
            this MsDataScan experimental,
            LibrarySpectrum theoretical,
            SpectralSimilarity.SimilarityMeasures[] measures,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoretical, scheme, toleranceInPpm, allPeaks)
                .GetSelectedSimilarityMeasures(measures);
        }

        /// <summary>
        /// Computes spectral similarity between an MsDataScan and an MzSpectrum and returns selected similarity measures.
        /// </summary>
        public static IEnumerable<(SpectralSimilarity.SimilarityMeasures, double?)> GetSelectedSimilarityMeasures(
            this MsDataScan experimental,
            MzSpectrum theoretical,
            SpectralSimilarity.SimilarityMeasures[] measures,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoretical, scheme, toleranceInPpm, allPeaks)
                .GetSelectedSimilarityMeasures(measures);
        }

        /// <summary>
        /// Computes spectral similarity between a library spectrum and an MsDataScan and returns selected similarity measures.
        /// </summary>
        public static IEnumerable<(SpectralSimilarity.SimilarityMeasures, double?)> GetSelectedSimilarityMeasures(
            this LibrarySpectrum experimental,
            MsDataScan theoretical,
            SpectralSimilarity.SimilarityMeasures[] measures,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoretical, scheme, toleranceInPpm, allPeaks)
                .GetSelectedSimilarityMeasures(measures);
        }

        /// <summary>
        /// Computes spectral similarity between an MzSpectrum and an MsDataScan and returns selected similarity measures.
        /// </summary>
        public static IEnumerable<(SpectralSimilarity.SimilarityMeasures, double?)> GetSelectedSimilarityMeasures(
            this MzSpectrum experimental,
            MsDataScan theoretical,
            SpectralSimilarity.SimilarityMeasures[] measures,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoretical, scheme, toleranceInPpm, allPeaks)
                .GetSelectedSimilarityMeasures(measures);
        }

        /// <summary>
        /// Computes spectral similarity between an MsDataScan and explicit m/z/intensity arrays and returns selected similarity measures.
        /// </summary>
        public static IEnumerable<(SpectralSimilarity.SimilarityMeasures, double?)> GetSelectedSimilarityMeasures(
            this MsDataScan experimental,
            double[] theoreticalX,
            double[] theoreticalY,
            SpectralSimilarity.SimilarityMeasures[] measures,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimental.ComputeSpectralSimilarity(theoreticalX, theoreticalY, scheme, toleranceInPpm, allPeaks)
                .GetSelectedSimilarityMeasures(measures);
        }

        /// <summary>
        /// Computes spectral similarity between explicit m/z/intensity arrays and an MsDataScan and returns selected similarity measures.
        /// </summary>
        public static IEnumerable<(SpectralSimilarity.SimilarityMeasures, double?)> GetSelectedSimilarityMeasures(
            this double[] experimentalX,
            double[] experimentalY,
            MsDataScan theoretical,
            SpectralSimilarity.SimilarityMeasures[] measures,
            SpectralSimilarity.SpectrumNormalizationScheme scheme = SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
            double toleranceInPpm = 20,
            bool allPeaks = true)
        {
            return experimentalX.ComputeSpectralSimilarity(experimentalY, theoretical, scheme, toleranceInPpm, allPeaks)
                .GetSelectedSimilarityMeasures(measures);
        }

        #endregion
    }
}
