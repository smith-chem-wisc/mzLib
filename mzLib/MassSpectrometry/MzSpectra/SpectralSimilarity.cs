using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace MassSpectrometry.MzSpectra
{
    public class SpectralSimilarity
    {
        public SpectralSimilarity(MzSpectrum experimentalSpectrum, MzSpectrum theoreticalSpectrum, SpectrumNormalizationScheme scheme, double toleranceInPpm, bool allPeaks, double filterOutBelowThisMz = 300)
        {
            ExperimentalYArray = Normalize(FilterOutIonsBelowThisMz(experimentalSpectrum.XArray,experimentalSpectrum.YArray, filterOutBelowThisMz).Select(p=>p.Item2).ToArray(),scheme);
            ExperimentalXArray = FilterOutIonsBelowThisMz(experimentalSpectrum.XArray, experimentalSpectrum.YArray, filterOutBelowThisMz).Select(p => p.Item1).ToArray();
            TheoreticalYArray = Normalize(FilterOutIonsBelowThisMz(theoreticalSpectrum.XArray, theoreticalSpectrum.YArray, filterOutBelowThisMz).Select(p => p.Item2).ToArray(), scheme);
            TheoreticalXArray = FilterOutIonsBelowThisMz(theoreticalSpectrum.XArray, theoreticalSpectrum.YArray, filterOutBelowThisMz).Select(p => p.Item1).ToArray();
            _localPpmTolerance = toleranceInPpm;
            IntensityPairs = GetIntensityPairs(allPeaks);
        }

        public SpectralSimilarity(MzSpectrum experimentalSpectrum, double[] theoreticalX, double[] theoreticalY, SpectrumNormalizationScheme scheme, double toleranceInPpm, bool allPeaks, double filterOutBelowThisMz = 300)
        {
            ExperimentalYArray = Normalize(FilterOutIonsBelowThisMz(experimentalSpectrum.XArray, experimentalSpectrum.YArray, filterOutBelowThisMz).Select(p => p.Item2).ToArray(), scheme);
            ExperimentalXArray = FilterOutIonsBelowThisMz(experimentalSpectrum.XArray, experimentalSpectrum.YArray, filterOutBelowThisMz).Select(p => p.Item1).ToArray();
            TheoreticalYArray = Normalize(FilterOutIonsBelowThisMz(theoreticalX, theoreticalY, filterOutBelowThisMz).Select(p => p.Item2).ToArray(), scheme);
            TheoreticalXArray = FilterOutIonsBelowThisMz(theoreticalX, theoreticalY, filterOutBelowThisMz).Select(p => p.Item1).ToArray();
            _localPpmTolerance = toleranceInPpm;
            IntensityPairs = GetIntensityPairs(allPeaks);
        }

        public SpectralSimilarity(double[] P_XArray, double[] P_YArray, double[] Q_XArray, double[] Q_YArray, SpectrumNormalizationScheme scheme, double toleranceInPpm, bool allPeaks, double filterOutBelowThisMz = 300)
        {
            ExperimentalYArray = Normalize(FilterOutIonsBelowThisMz(P_XArray, P_YArray, filterOutBelowThisMz).Select(p => p.Item2).ToArray(), scheme);
            ExperimentalXArray = FilterOutIonsBelowThisMz(P_XArray, P_YArray, filterOutBelowThisMz).Select(p => p.Item1).ToArray();
            TheoreticalYArray = Normalize(FilterOutIonsBelowThisMz(Q_XArray, Q_YArray, filterOutBelowThisMz).Select(p => p.Item2).ToArray(), scheme);
            TheoreticalXArray = FilterOutIonsBelowThisMz(Q_XArray, Q_YArray, filterOutBelowThisMz).Select(p => p.Item1).ToArray();
            _localPpmTolerance = toleranceInPpm;
            IntensityPairs = GetIntensityPairs(allPeaks);
        }
        public double[] ExperimentalYArray { get; private set; }
        public double[] ExperimentalXArray { get; private set; }
        public double[] TheoreticalYArray { get; private set; }
        public double[] TheoreticalXArray { get; private set; }

        private readonly double _localPpmTolerance;

        public List<(double, double)> IntensityPairs { get; } = new();

        /// All peaks with mz less than the cutOff will be filtered out. default to zero to remove an mz values that are accidentally negative. this is an unexpected error. 
        private static IEnumerable<(double, double)> FilterOutIonsBelowThisMz(IReadOnlyList<double> spectrumX, IReadOnlyList<double> spectrumY,double filterOutBelowThisMz = 0)
        {
            if (spectrumY.Count == 0)
            {
                throw new MzLibException(string.Format(CultureInfo.InvariantCulture, "Empty YArray in spectrum."));
            }
            if (spectrumY.Sum() == 0)
            {
                throw new MzLibException(string.Format(CultureInfo.InvariantCulture, "Spectrum has no intensity."));
            }

            List<(double, double)> spectrumWithMzCutoff = new List<(double, double)>();
            for (int i = 0; i < spectrumX.Count; i++)
            {
                //second conditional to avoid getting an accidental negative intensities
                if (spectrumX[i] >= filterOutBelowThisMz && spectrumY[i] >= 0)
                {
                    spectrumWithMzCutoff.Add((spectrumX[i], spectrumY[i]));
                }
            }
            return spectrumWithMzCutoff;
        }

        /// <summary>
        /// Every spectrum gets normalized when the SpectralSimilarity object gets created. This methods sends the spectra to the appropriate normalization.
        /// </summary>
        /// <param name="spectrum"></param>
        /// <param name="scheme"></param>
        /// <returns></returns>
        private static double[] Normalize(double[] spectrum, SpectrumNormalizationScheme scheme)
        {
            if (spectrum.Length == 0)
            {
                return null; 
            }
           
            return scheme switch
            {
                SpectrumNormalizationScheme.mostAbundantPeak => NormalizeMostAbundantPeak(spectrum),
                SpectrumNormalizationScheme.spectrumSum => NormalizeSpectrumSum(spectrum),
                SpectrumNormalizationScheme.squareRootSpectrumSum => NormalizeSquareRootSpectrumSum(spectrum),
                _ => spectrum,
            };
        }

        /// Intensity Pairs a computed immediately upon creation of the SpectralSimilarity object. That way they can be used in all the methods without being recomputed.
        private List<(double,double)> GetIntensityPairs(bool allPeaks, double[] experimentalYArray = null, double[] theoreticalYArray = null)
        {
            if (experimentalYArray == null) experimentalYArray = ExperimentalYArray;
            if (theoreticalYArray == null) theoreticalYArray = TheoreticalYArray;

            if (experimentalYArray == null || theoreticalYArray == null)
            {
                //when all mz of theoretical peaks or experimental peaks are less than mz cut off , it is treated as no corresponding library spectrum is found and later the similarity score will be assigned as null.
                return new List<(double, double)> { (-1, -1) };
            }

            List<(double, double)> intensityPairs = new();
            int expIndex = 0;
            int theoIndex = 0;
            do
            {
                if (Within(ExperimentalXArray[expIndex], TheoreticalXArray[theoIndex]))
                {
                    intensityPairs.Add((ExperimentalYArray[expIndex], TheoreticalYArray[theoIndex]));
                    expIndex++;
                    theoIndex++;
                }
                else if(ExperimentalXArray[expIndex] < TheoreticalXArray[theoIndex])
                {
                    if (allPeaks)
                    {
                        intensityPairs.Add((ExperimentalYArray[expIndex], 0));
                    }
                    expIndex++;
                }
                else
                {
                    if (allPeaks)
                    {
                        intensityPairs.Add((0, TheoreticalYArray[theoIndex]));
                    }
                    theoIndex++;
                }
            }
            while (expIndex < ExperimentalXArray.Length && theoIndex < TheoreticalXArray.Length);

            //if the theoretical peak count is different than the experimental peak count, and the bool allPeaks = TRUE then
            //we need to add zero intensity pairs for each peak that does not have a corresponding experimental peak
            if (allPeaks)
            {
                while (expIndex < ExperimentalXArray.Length)
                {
                    intensityPairs.Add((ExperimentalYArray[expIndex], 0));
                    expIndex++;
                }

                while (theoIndex < TheoreticalXArray.Length)
                {
                    intensityPairs.Add((0, TheoreticalYArray[theoIndex]));
                    theoIndex++;
                }
            }

            //if there are no intensity pairs, then we are required to return a single pair of (-1,-1) to indicate that no peaks were found
            if (intensityPairs.Count == 0)
            {
                intensityPairs.Add((-1, -1));
            }
            return intensityPairs;
        }

        #region normalization

        public static double[] NormalizeSquareRootSpectrumSum(double[] spectrum)
        {
            double sqrtSum = spectrum.Select(y => Math.Sqrt(y)).Sum();
            double[] normalizedSpectrum = new double[spectrum.Length];

            for (int i = 0; i < spectrum.Length; i++)
            {
                normalizedSpectrum[i] = Math.Sqrt(spectrum[i]) / sqrtSum;
            }
            return normalizedSpectrum;
        }

        public static double[] NormalizeMostAbundantPeak(double[] spectrum)
        {
            double max = spectrum.Max();
            double[] normalizedSpectrum = new double[spectrum.Length];

            for (int i = 0; i < spectrum.Length; i++)
            {
                normalizedSpectrum[i] = spectrum[i] / max;
            }
            return normalizedSpectrum;
        }

        public static double[] NormalizeSpectrumSum(double[] spectrum)
        {
            double sum = spectrum.Sum();
            double[] normalizedSpectrum = new double[spectrum.Length];

            for (int i = 0; i < spectrum.Length; i++)
            {
                normalizedSpectrum[i] = spectrum[i] / sum;
            }
            return normalizedSpectrum;
        }

        #endregion normalization

        #region similarityMethods

        //The cosine similarity returns values between 1 and -1 with 1 being closes and -1 being opposite and 0 being orthoganal
        public double? CosineSimilarity()
        {
            if (IntensityPairs.First().Item1 < 0)//if the first pair is (-1,-1) then there are no peaks to compare
            {
                return null;
            }
            double numerator = 0;
            double denominatorValue1 = 0;
            double denominatorValue2 = 0;
            foreach ((double, double) pair in IntensityPairs)
            {
                numerator += pair.Item1 * pair.Item2;
                denominatorValue1 += Math.Pow(pair.Item1, 2);
                denominatorValue2 += Math.Pow(pair.Item2, 2);
            }
            double denominatorProduct = denominatorValue1 * denominatorValue2;
            
            //because we keep all secondary spectrum peaks, denominatorValue1 can equal zero
            if(denominatorProduct == 0)
            {
                return 0;
            }
            return numerator / Math.Sqrt(denominatorProduct);
        }

        //Spectral contrast angle should expect values between 1 and -1;
        public double? SpectralContrastAngle()
        {
            if (IntensityPairs.First().Item1 < 0)//if the first pair is (-1,-1) then there are no peaks to compare
            {
                return null;
            }
            return (1 - 2 * Math.Acos((double)CosineSimilarity()) / Math.PI);

        }

        public double? EuclideanDistance()
        {
            if (IntensityPairs.First().Item1 < 0)
            {
                return null;
            }
            double sum = 0;
            foreach ((double, double) pair in IntensityPairs)
            {
                sum += Math.Pow(pair.Item1 - pair.Item2, 2);
            }
            return 1 - Math.Sqrt(sum);
        }

        public double? BrayCurtis()
        {
            if (IntensityPairs.First().Item1 < 0)
            {
                return null;
            }
            double numerator = 0;
            double denominator = 0;
            foreach ((double, double) pair in IntensityPairs)
            {
                numerator += Math.Abs(pair.Item1 - pair.Item2);
                denominator += pair.Item1 + pair.Item2;
            }
            return (1 - numerator / denominator);
        }

        public double? PearsonsCorrelation()
        {
            if (IntensityPairs.First().Item1 < 0)
            {
                return null;
            }
            double numerator = 0;
            double denominator = 0;
            double denominator1 = 0;
            double denominator2 = 0;
            double averagePrimaryIntensity = IntensityPairs.Select(a => a.Item1).Average();
            double averageSecondaryIntensity = IntensityPairs.Select(a => a.Item2).Average();
            foreach ((double, double) pair in IntensityPairs)
            {
                numerator += (pair.Item1 - averagePrimaryIntensity) * (pair.Item2 - averageSecondaryIntensity);
                denominator1 += Math.Pow((pair.Item1 - averagePrimaryIntensity), 2);
                denominator2 += Math.Pow((pair.Item2 - averageSecondaryIntensity), 2);
            }
            denominator = denominator1 * denominator2;
            if(denominator > 0)
            {
                return numerator / Math.Sqrt(denominator);
            }
            return -1;
        }

        public double? DotProduct()
        {
            if (IntensityPairs.First().Item1 < 0)//if the first pair is (-1,-1) then there are no peaks to compare
            {
                return null;
            }
            double sum = 0;
            foreach ((double, double) pair in IntensityPairs)
            {
                sum += pair.Item1 * pair.Item2;
            }
            return sum;
        }

        // This method should only be used with the SpectrumSum normalization method
        // This method should only be used when allPeaks is set to true
        public double? SpectralEntropy()
        {
            double theoreticalEntropy = 0;
            foreach (double intensity in TheoreticalYArray)
            {
                theoreticalEntropy += -1 * intensity * Math.Log(intensity);
            }
            double experimentalEntropy = 0;
            foreach (double intensity in ExperimentalYArray)
            {
                experimentalEntropy += -1 * intensity * Math.Log(intensity);
            }

            double combinedEntropy = 0;
            foreach ( (double,double) intensityPair in IntensityPairs)
            {
                double combinedIntensity = intensityPair.Item1 / 2 + intensityPair.Item2 / 2;
                combinedEntropy += -1* combinedIntensity * Math.Log(combinedIntensity);
            }

            double similarityScore = 1 - (2 * combinedEntropy - theoreticalEntropy - experimentalEntropy) / Math.Log(4);
            return similarityScore;
        }

        /// <summary>
        /// This is used to calculate the KullbackLeibler divergence between a theoretical and experimental isotopic envelope
        /// When using, the allPeaks argument in the SpectralSimilarity constructor should be set to "true"
        /// </summary>
        /// <param name="correctionConstant"> A constant value added to each intensity pair in order to penalize zero peaks.
        /// Default is 1e-9 (one OoM smaller than the theoretical isotopic intensity minimum)</param>
        /// <returns> A nullable double between 0 and positive infinity. More similar envelopes score lower </returns>
        public double? KullbackLeiblerDivergence_P_Q(double correctionConstant = 1e-9)
        {
            // counts the number of zero values. If zeroCount == intensityPairs.Count, then there are no shared peaks
            int zeroCount = IntensityPairs.Count(p => p.Item1 == 0 | p.Item2 == 0);
            if (zeroCount == IntensityPairs.Count) return null;
            double divergence = 0;

            if (zeroCount == 0) // | correctionConstant == 0)
            {
                foreach (var pair in IntensityPairs)
                {
                    if (pair.Item1 != 0 && pair.Item2 != 0)
                    {
                        divergence += pair.Item1 * Math.Log(pair.Item1 / pair.Item2);
                    }
                }
            }
            else
            {
                // Add correctionConstant and renormalize
                // Need to use temp variables to avoid modifying the Y array fields
                double[] tempExperimentalYArray = Normalize(ExperimentalYArray.Select(i => i + correctionConstant).ToArray(), SpectrumNormalizationScheme.spectrumSum);
                double[] tempTheoreticalYArray = Normalize(TheoreticalYArray.Select(i => i + correctionConstant).ToArray(), SpectrumNormalizationScheme.spectrumSum);
                List<(double, double)> correctedIntensityPairs = GetIntensityPairs(
                    allPeaks: true,
                    experimentalYArray: tempExperimentalYArray,
                    theoreticalYArray: tempTheoreticalYArray);

                foreach (var pair in correctedIntensityPairs)
                {
                    if (pair.Item1 != 0 && pair.Item2 != 0)
                    {
                        divergence += pair.Item1 * Math.Log(pair.Item1 / pair.Item2);
                    }
                }
            }

            return divergence;
        }

        //This formula created by Brian Searle and reported at ASMS 2022 in poster "Scribe: next generation library searching for DDA experiments"
        //Please you the square root normalization with this calculation
        public double? SearleSimilarity()
        {
            double squaredSumDifferences = 0;

            //there must be some legitimate pairs to enter this function so no need to test if pairs exist
            foreach ((double, double) pair in IntensityPairs)
            {
                squaredSumDifferences += Math.Pow((pair.Item1 - pair.Item2),2);
            }
            return squaredSumDifferences > 0 ? Math.Log(Math.Pow(squaredSumDifferences, -1)) : double.MaxValue;
        }

        #endregion similarityMethods

        //use Math.Max() in the denominator for consistency
        private bool Within(double mz1, double mz2)
        {
            return ((Math.Abs(mz1 - mz2) / Math.Max(mz1,mz2) * 1000000.0) < _localPpmTolerance);
        }

        public enum SpectrumNormalizationScheme
        {
            squareRootSpectrumSum,
            spectrumSum,
            mostAbundantPeak,
            unnormalized
        }
    }
}