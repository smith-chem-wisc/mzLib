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
            LocalPpmTolerance = toleranceInPpm;
            normalizationScheme = scheme;
            _intensityPairs = IntensityPairs(allPeaks);
        }

        public SpectralSimilarity(MzSpectrum experimentalSpectrum, double[] theoreticalX, double[] theoreticalY, SpectrumNormalizationScheme scheme, double toleranceInPpm, bool allPeaks, double filterOutBelowThisMz = 300)
        {
            ExperimentalYArray = Normalize(FilterOutIonsBelowThisMz(experimentalSpectrum.XArray, experimentalSpectrum.YArray, filterOutBelowThisMz).Select(p => p.Item2).ToArray(), scheme);
            ExperimentalXArray = FilterOutIonsBelowThisMz(experimentalSpectrum.XArray, experimentalSpectrum.YArray, filterOutBelowThisMz).Select(p => p.Item1).ToArray();
            TheoreticalYArray = Normalize(FilterOutIonsBelowThisMz(theoreticalX, theoreticalY, filterOutBelowThisMz).Select(p => p.Item2).ToArray(), scheme);
            TheoreticalXArray = FilterOutIonsBelowThisMz(theoreticalX, theoreticalY, filterOutBelowThisMz).Select(p => p.Item1).ToArray();
            LocalPpmTolerance = toleranceInPpm;
            normalizationScheme = scheme;
            _intensityPairs = IntensityPairs(allPeaks);
        }

        public SpectralSimilarity(double[] P_XArray, double[] P_YArray, double[] Q_XArray, double[] Q_YArray, SpectrumNormalizationScheme scheme, double toleranceInPpm, bool allPeaks, double filterOutBelowThisMz = 300)
        {
            ExperimentalYArray = Normalize(FilterOutIonsBelowThisMz(P_XArray, P_YArray, filterOutBelowThisMz).Select(p => p.Item2).ToArray(), scheme);
            ExperimentalXArray = FilterOutIonsBelowThisMz(P_XArray, P_YArray, filterOutBelowThisMz).Select(p => p.Item1).ToArray();
            TheoreticalYArray = Normalize(FilterOutIonsBelowThisMz(Q_XArray, Q_YArray, filterOutBelowThisMz).Select(p => p.Item2).ToArray(), scheme);
            TheoreticalXArray = FilterOutIonsBelowThisMz(Q_XArray, Q_YArray, filterOutBelowThisMz).Select(p => p.Item1).ToArray();
            LocalPpmTolerance = toleranceInPpm;
            _intensityPairs = IntensityPairs(allPeaks);
        }
        public double[] ExperimentalYArray { get; private set; }
        public double[] ExperimentalXArray { get; private set; }
        public double[] TheoreticalYArray { get; private set; }
        public double[] TheoreticalXArray { get; private set; }

        private SpectrumNormalizationScheme normalizationScheme;

        private double LocalPpmTolerance;

        private readonly List<(double, double)> _intensityPairs = new();
        
        public List<(double, double)> intensityPairs
        { get { return _intensityPairs; } }


        /// <summary>
        /// All peaks with mz less than the cutOff will be filtered out. default to zero to remove an mz values that are accidentally negative. this is an unexpected error. 
        private static List<(double, double)> FilterOutIonsBelowThisMz(double[] spectrumX, double[] spectrumY,double filterOutBelowThisMz = 0)
        {
            if (spectrumY.Length == 0)
            {
                throw new MzLibException(string.Format(CultureInfo.InvariantCulture, "Empty YArray in spectrum."));
            }
            if (spectrumY.Sum() == 0)
            {
                throw new MzLibException(string.Format(CultureInfo.InvariantCulture, "Spectrum has no intensity."));
            }

            List<(double, double)> spectrumWithMzCutoff = new List<(double, double)>();
            for (int i = 0; i < spectrumX.Length; i++)
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
        private double[] Normalize(double[] spectrum, SpectrumNormalizationScheme scheme)
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

        /// <summary>
        /// Intensity Pairs a computed immediately upon creation of the SpectralSimilarity object. That way they can be used in all the methods without being recomputed.
        /// We loop throught the secondaryXArray under the assumption that it is the shorter of the two arrays (i.e. typically the theoretical spectrum).
        /// Experimental spectrum defaults to 200 peaks and is therefore usually longer.
        /// We sort intensities in descending order so that when we make peak pairs, we're choosing pairs with the highest intensity so long as they are with mz range. 
        /// Sometimes you could have two peaks in mz range and I don't think you want to pair the lesser intensity peak first just because it is closer in mass.
        /// 
        /// Intensity Pair: (Experimental Intensity , Theoretical Intensity)
        /// 
        /// </summary>


        private List<(double, double)> IntensityPairs(bool allPeaks)
        {
            if (ExperimentalYArray==null || TheoreticalYArray == null)
            {
                //when all mz of theoretical peaks or experimental peaks are less than mz cut off , it is treated as no corresponding library spectrum is found and later the similarity score will be assigned as null.
                return new List<(double, double)> { (-1, -1) };
            }

            List<(double, double)> intensityPairs = new();
            List<(double, double)> experimental = new();
            List<(double, double)> theoretical = new();

            for (int i = 0; i < ExperimentalXArray.Length; i++)
            {
                experimental.Add((ExperimentalXArray[i], ExperimentalYArray[i]));
            }
            for (int i = 0; i < TheoreticalXArray.Length; i++)
            {
                theoretical.Add((TheoreticalXArray[i], TheoreticalYArray[i]));
            }
            
            experimental = experimental.OrderByDescending(i => i.Item2).ToList();
            theoretical = theoretical.OrderByDescending(i => i.Item2).ToList();

            foreach ((double,double) xyPair in theoretical)
            {
                int index = 0;
                while(experimental.Count >0 && index < experimental.Count)
                {
                    if (Within(experimental[index].Item1, xyPair.Item1))
                    {
                        intensityPairs.Add((experimental[index].Item2, xyPair.Item2));
                        experimental.RemoveAt(index);
                        index = -1;
                        break;
                    }
                    index++;
                }
                if (experimental.Count == 0)
                {
                    index++;
                }
                if (index > 0)
                {
                    //didn't find a experimental mz in range
                    intensityPairs.Add((0, xyPair.Item2));
                }
            }

            //If we're keeping all experimental and theoretical peaks, then we add intensity pairs for all unpaired experimental peaks here.
            if(experimental.Count > 0 && allPeaks)
            {
                foreach ((double, double) xyPair in experimental)
                {
                    intensityPairs.Add((xyPair.Item2, 0));
                }
            }
            return intensityPairs;
        }

        #region normalization

        private double[] NormalizeSquareRootSpectrumSum(double[] spectrum)
        {
            double sqrtSum = spectrum.Select(y => Math.Sqrt(y)).Sum();

            for (int i = 0; i < spectrum.Length; i++)
            {
                spectrum[i] = Math.Sqrt(spectrum[i]) / sqrtSum;
            }
            return spectrum;
        }

        private double[] NormalizeMostAbundantPeak(double[] spectrum)
        {
            double max = spectrum.Max();

            for (int i = 0; i < spectrum.Length; i++)
            {
                spectrum[i] = spectrum[i] / max;
            }
            return spectrum;
        }

        private double[] NormalizeSpectrumSum(double[] spectrum)
        {
            double sum = spectrum.Sum();

            for (int i = 0; i < spectrum.Length; i++)
            {
                spectrum[i] = spectrum[i] / sum;
            }
            return spectrum;
        }

        #endregion normalization

        #region similarityMethods

        //The cosine similarity returns values between 1 and -1 with 1 being closes and -1 being opposite and 0 being orthoganal
        public double? CosineSimilarity()
        {
            if (_intensityPairs.First().Item1==-1)
            {
                return null;
            }
            double numerator = 0;
            double denominatorValue1 = 0;
            double denominatorValue2 = 0;
            foreach ((double, double) pair in _intensityPairs)
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
            if (_intensityPairs.First().Item1 == -1)
            {
                return null;
            }
            return (1 - 2 * Math.Acos((double)CosineSimilarity()) / Math.PI);

        }

        public double? EuclideanDistance()
        {
            if (_intensityPairs.First().Item1 == -1)
            {
                return null;
            }
            double sum = 0;
            foreach ((double, double) pair in _intensityPairs)
            {
                sum += Math.Pow(pair.Item1 - pair.Item2, 2);
            }
            return 1 - Math.Sqrt(sum);
        }

        public double? BrayCurtis()
        {
            if (_intensityPairs.First().Item1 == -1)
            {
                return null;
            }
            double numerator = 0;
            double denominator = 0;
            foreach ((double, double) pair in _intensityPairs)
            {
                numerator += Math.Abs(pair.Item1 - pair.Item2);
                denominator += pair.Item1 + pair.Item2;
            }
            return (1 - numerator / denominator);
        }

        public double? PearsonsCorrelation()
        {
            if (_intensityPairs.First().Item1 == -1)
            {
                return null;
            }
            double numerator = 0;
            double denominator = 0;
            double denominator1 = 0;
            double denominator2 = 0;
            double averagePrimaryIntensity = intensityPairs.Select(a => a.Item1).Sum() / intensityPairs.Count;
            double averageSecondaryIntensity = intensityPairs.Select(a => a.Item2).Sum() / intensityPairs.Count;
            foreach ((double, double) pair in _intensityPairs)
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
            if (_intensityPairs.First().Item1 == -1)
            {
                return null;
            }
            double sum = 0;
            foreach ((double, double) pair in _intensityPairs)
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
            foreach ( (double,double) intensityPair in _intensityPairs)
            {
                double combinedIntensity = intensityPair.Item1 / 2 + intensityPair.Item2 / 2;
                combinedEntropy += -1* combinedIntensity * Math.Log(combinedIntensity);
            }

            double similarityScore = 1 - (2 * combinedEntropy - theoreticalEntropy - experimentalEntropy) / Math.Log(4);
            return similarityScore;
        }

        // should only be used with SpectrumSum normalization scheme.
        public double? KullbackLeiblerDivergence_P_Q()
        {
            double divergence = 0;
            bool correctionNeeded = false;
            foreach (var pair in intensityPairs)
            {
                if (pair.Item1 != 0 && pair.Item2 != 0)
                {
                    divergence += pair.Item1 * Math.Log(pair.Item1 / pair.Item2);
                }
                else
                {
                    correctionNeeded = true;
                }
            }
            if (correctionNeeded)
            {
                // correct by adding 1e-8 to every intensity value (minimum probability for isotope distribution calculation) and renormalize
                TheoreticalYArray = Normalize(TheoreticalYArray.Select(i => i + 1e-8).ToArray(), SpectrumNormalizationScheme.spectrumSum);
                ExperimentalYArray = Normalize(ExperimentalYArray.Select(i => i + 1e-8).ToArray(), SpectrumNormalizationScheme.spectrumSum);
                List<(double, double)> correctedIntensityPairs = IntensityPairs(allPeaks: true);

                divergence = 0;
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
            foreach ((double, double) pair in _intensityPairs)
            {
                squaredSumDifferences += Math.Pow((pair.Item1 - pair.Item2),2);
            }
            return squaredSumDifferences > 0 ? Math.Log(Math.Pow(squaredSumDifferences, -1)) : double.MaxValue;
        }

        #endregion similarityMethods

        //use Math.Max() in the denominator for consistancy
        private bool Within(double mz1, double mz2)
        {
            return ((Math.Abs(mz1 - mz2) / Math.Max(mz1,mz2) * 1000000.0) < LocalPpmTolerance);
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