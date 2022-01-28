﻿using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace MassSpectrometry.MzSpectra
{
    public class SpectralSimilarity
    {
        public SpectralSimilarity(MzSpectrum experimentalSpectrum, MzSpectrum theoreticalSpectrum, SpectrumNormalizationScheme scheme, double toleranceInPpm, bool allPeaks)
        {
            experimentalYArray = Normalize(experimentalSpectrum.YArray, scheme);
            experimentalXArray = experimentalSpectrum.XArray;
            theoreticalYArray = Normalize(theoreticalSpectrum.YArray, scheme);
            theoreticalXArray = theoreticalSpectrum.XArray;
            localPpmTolerance = toleranceInPpm;
            _intensityPairs = IntensityPairs(allPeaks);
        }

        public SpectralSimilarity(MzSpectrum experimentalSpectrum, double[] theoreticalX, double[] theoreticalY, SpectrumNormalizationScheme scheme, double toleranceInPpm, bool allPeaks)
        {
            experimentalYArray = Normalize(experimentalSpectrum.YArray, scheme);
            experimentalXArray = experimentalSpectrum.XArray;
            theoreticalYArray = Normalize(theoreticalY, scheme);
            theoreticalXArray = theoreticalX;
            localPpmTolerance = toleranceInPpm;
            _intensityPairs = IntensityPairs(allPeaks);
        }

        public double[] experimentalYArray { get; private set; }
        public double[] experimentalXArray { get; private set; }
        public double[] theoreticalYArray { get; private set; }
        public double[] theoreticalXArray { get; private set; }

        private double localPpmTolerance;

        private List<(double, double)> _intensityPairs = new List<(double, double)>();
        
        public List<(double, double)> intensityPairs
        { get { return _intensityPairs; } }

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
                throw new MzLibException(string.Format(CultureInfo.InvariantCulture, "Empty YArray in spectrum."));
            }
            if (spectrum.Sum() == 0)
            {
                throw new MzLibException(string.Format(CultureInfo.InvariantCulture, "Spectrum has no intensity."));
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
        /// </summary>
        /// <returns></returns>

        private List<(double, double)> IntensityPairs(bool allPeaks)
        {
            List<(double, double)> intensityPairs = new List<(double, double)>();
            List<(double, double)> experimental = new List<(double, double)>();
            List<(double, double)> theoretical = new List<(double, double)>();

            for (int i = 0; i < experimentalXArray.Length; i++)
            {
                experimental.Add((experimentalXArray[i], experimentalYArray[i]));
            }
            for (int i = 0; i < theoreticalXArray.Length; i++)
            {
                theoretical.Add((theoreticalXArray[i], theoreticalYArray[i]));
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
                if(index > 0)
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
        public double CosineSimilarity()
        {
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
        public double SpectralContrastAngle()
        {
            return (1 - 2 * Math.Acos(CosineSimilarity()) / Math.PI);

        }

        public double EuclideanDistance()
        {
            double sum = 0;
            foreach ((double, double) pair in _intensityPairs)
            {
                sum += Math.Pow(pair.Item1 - pair.Item2, 2);
            }
            return 1 - Math.Sqrt(sum);
        }

        public double BrayCurtis()
        {
            double numerator = 0;
            double denominator = 0;
            foreach ((double, double) pair in _intensityPairs)
            {
                numerator += Math.Abs(pair.Item1 - pair.Item2);
                denominator += pair.Item1 + pair.Item2;
            }
            return (1 - numerator / denominator);
        }

        public double PearsonsCorrelation()
        {
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

        public double DotProduct()
        {
            double sum = 0;
            foreach ((double, double) pair in _intensityPairs)
            {
                sum += pair.Item1 * pair.Item2;
            }
            return sum;
        }

        #endregion similarityMethods

        //use Math.Max() in the denominator for consistancy
        private bool Within(double mz1, double mz2)
        {
            return ((Math.Abs(mz1 - mz2) / Math.Max(mz1,mz2) * 1000000.0) < localPpmTolerance);
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