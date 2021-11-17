using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace MassSpectrometry.MzSpectra
{
    public class SpectralSimilarity
    {
        public SpectralSimilarity(MzSpectrum primary, MzSpectrum secondary, SpectrumNormalizationScheme scheme, double toleranceInPpm)
        {
            primaryYArray = Normalize(primary.YArray, scheme);
            primaryXArray = primary.XArray;
            secondaryYarray = Normalize(secondary.YArray, scheme);
            secondaryXArray = secondary.XArray;
            localTolerance = toleranceInPpm / 1000000.0;
            _intensityPairs = IntensityPairs();
        }

        public double[] primaryYArray { get; private set; }
        public double[] primaryXArray { get; private set; }
        public double[] secondaryYarray { get; private set; }
        public double[] secondaryXArray { get; private set; }
        private double localTolerance;
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
        /// </summary>
        /// <returns></returns>

        private List<(double, double)> IntensityPairs()
        {
            List<(double, double)> intensityPairs = new List<(double, double)>();
            double[] localPrimaryXarray = primaryXArray.ToArray();
            for (int i = 0; i < secondaryXArray.Length; i++)
            {
                double secondaryMz = secondaryXArray[i];
                double nearestMz = FindClosest(primaryXArray, secondaryMz);
                if (Within(secondaryMz, nearestMz))
                {
                    int theIndex = Array.IndexOf(localPrimaryXarray, nearestMz);
                    intensityPairs.Add((primaryYArray[theIndex], secondaryYarray[i]));
                    localPrimaryXarray[theIndex] = -1;
                }
                else
                {
                    intensityPairs.Add((0, secondaryYarray[i]));
                }
            }
            for (int i = 0; i < localPrimaryXarray.Length; i++)
            {
                if (localPrimaryXarray[i] > 0)
                {
                    intensityPairs.Add((primaryYArray[i], 0));
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
            return numerator / Math.Sqrt(denominatorProduct);
        }

        public double SpectralContrastAngle()
        {
            return (1.0 - ((2 * Math.Acos(CosineSimilarity())) / Math.PI));
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
            denominator += denominator1 * denominator2;
            return numerator / Math.Sqrt(denominator);
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

        #region binarySearch

        private static double FindClosest(double[] arr,
                              double target)
        {
            int n = arr.Length;

            // Corner cases
            if (target <= arr[0])
                return arr[0];
            if (target >= arr[n - 1])
                return arr[n - 1];

            // Doing binary search
            int i = 0, j = n, mid = 0;
            while (i < j)
            {
                mid = (i + j) / 2;

                if (arr[mid] == target)
                    return arr[mid];

                /* If target is less
                than array element,
                then search in left */
                if (target < arr[mid])
                {
                    // If target is greater
                    // than previous to mid,
                    // return closest of two
                    if (mid > 0 && target > arr[mid - 1])
                        return getClosest(arr[mid - 1],
                                     arr[mid], target);

                    /* Repeat for left half */
                    j = mid;
                }

                // If target is
                // greater than mid
                else
                {
                    if (mid < n - 1 && target < arr[mid + 1])
                        return getClosest(arr[mid],
                             arr[mid + 1], target);
                    i = mid + 1; // update i
                }
            }

            // Only single element
            // left after search
            return arr[mid];
        }

        // Method to compare which one
        // is the more close We find the
        // closest by taking the difference
        // between the target and both
        // values. It assumes that val2 is
        // greater than val1 and target
        // lies between these two.
        private static double getClosest(double val1, double val2,
                                     double target)
        {
            if (target - val1 >= val2 - target)
                return val2;
            else
                return val1;
        }

        #endregion binarySearch

        private bool Within(double mz1, double mz2)
        {
            return (Math.Abs(mz1 - mz2) < localTolerance);
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