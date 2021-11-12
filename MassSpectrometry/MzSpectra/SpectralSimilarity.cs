using System;
using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry.MzSpectra
{
    public class SpectralSimilarity
    {
        public SpectralSimilarity(MzSpectrum primary, MzSpectrum secondary, double toleranceInPpm)
        {
            primarySpectrum = primary;
            secondarySpectrum = secondary;
            localTolerance = toleranceInPpm / 1000000.0;
            _intensityPairs = IntensityPairs();
        }

        public MzSpectrum primarySpectrum { get; private set; }
        public MzSpectrum secondarySpectrum { get; private set; }
        private double localTolerance { get; }
        private List<(double, double)> _intensityPairs = new List<(double, double)>();
        public List<(double, double)> intensityPairs
        { get { return _intensityPairs; } }

        public double CosineSimilarity()
        {
            if (intensityPairs.Count > 0)
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
                if (!denominatorProduct.Equals(0))
                {
                    return numerator / Math.Sqrt(denominatorProduct);
                }
                else
                {
                    return 0;
                }
            }
            return 0;
        }

        public double SpectralContrastAngle()
        {
            return (1.0 - ((2 * Math.Acos(CosineSimilarity())) / Math.PI));
        }

        public double EuclideanDistance()
        {
            double sum = 0;
            if (intensityPairs.Count > 0)
            {
                foreach ((double, double) pair in _intensityPairs)
                {
                    sum += Math.Pow(pair.Item1 - pair.Item2, 2);
                }
                return 1 - Math.Sqrt(sum);
            }
            else
            {
                return 1-Math.Sqrt(2);
            }
        }



        public double BrayCurtis()
        {
            if (intensityPairs.Count > 0)
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
            return 0;
        }

        public double PearsonsCorrelation()
        {
            if (intensityPairs.Count > 0)
            {
                double numerator = 0;
                double denominator = 0;
                double averagePrimaryIntensity = intensityPairs.Select(a => a.Item1).Sum() / intensityPairs.Count;
                double averageSecondaryIntensity = intensityPairs.Select(a => a.Item2).Sum() / intensityPairs.Count;
                foreach ((double, double) pair in _intensityPairs)
                {
                    numerator += (pair.Item1 - averagePrimaryIntensity) * (pair.Item2 - averageSecondaryIntensity);
                    double denominator1 = 0;
                    denominator1 = Math.Pow((pair.Item1 - averagePrimaryIntensity), 2);
                    double denominator2 = 0;
                    denominator2 = Math.Pow((pair.Item2 - averageSecondaryIntensity), 2);
                    denominator += denominator1 * denominator2;
                }
                if(!denominator.Equals(0))
                {
                    return numerator / Math.Sqrt(denominator);
                }
                return -1;
            }
            return -1;
        }

        public double DotProduct()
        {
            if(intensityPairs.Count > 0)
            {
                double sum = 0;
                foreach ((double, double) pair in _intensityPairs)
                {
                    sum += pair.Item1 * pair.Item2;
                }
                return sum;
            }
            return 0;
        }
        private List<(double, double)> IntensityPairs()
        {
            List<(double, double)> intensityPairs = new List<(double, double)>();
            double[] normalizedPrimaryIntensities = NormalizeSpectrumByTotalIntensity(primarySpectrum.YArray);
            double[] normalizedSecondaryIntensities = NormalizeSpectrumByTotalIntensity(secondarySpectrum.YArray);
            for (int i = 0; i < primarySpectrum.XArray.Length; i++)
            {
                double primaryMz = primarySpectrum.XArray[i];
                var nearestMz = secondarySpectrum.XArray.OrderBy(x => Math.Abs((long)x - primaryMz)).First();
                if (Within(primaryMz, nearestMz))
                {
                    intensityPairs.Add((normalizedPrimaryIntensities[i], normalizedSecondaryIntensities[Array.IndexOf(normalizedSecondaryIntensities, nearestMz)]));
                }
            }
            return intensityPairs;
        }

        private bool Within(double mz1, double mz2)
        {
            return (Math.Abs(mz1 - mz2) < localTolerance);
        }

        public double[] NormalizeSpectrumByTotalIntensity(double[] spectrum)
        {
            double[] normalizedArray = new double[spectrum.Length];
            double intensitySum = spectrum.Sum();
            for (int i = 0; i < spectrum.Length; i++)
            {
                normalizedArray[i] = spectrum[i] / intensitySum;
            }
            return normalizedArray;
        }
    }
}