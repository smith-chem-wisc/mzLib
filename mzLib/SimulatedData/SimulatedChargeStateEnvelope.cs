using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using Easy.Common.Extensions;
using MathNet.Numerics;
using MathNet.Numerics.Distributions;

namespace SimulatedData
{
    public class SimulatedChargeStateEnvelope : SimulatedData
    {
        private double _mzLow { get; set; }
        private double _mzHigh { get; set; }
        private int _chargeStateLow { get; set; }
        private int _chargeStateHigh { get; set; }
        private IsotopicDistribution _parentDistribution { get; set; }
        private double _adductMass { get; set; }
        private int _minPresentChargeState { get; set; }
        private int _maxPresentChargeState { get; set; }
        public SimulatedChargeStateEnvelope(double mzLow, double mzHigh, double stepSize, int length,
            int chargeStateLow, int chargeStateHigh, IsotopicDistribution parentIsotopicDistribution,
            (double mu, double sigma)? envelopeDistribution = null, double adductMass = 1.007276) 
        : base(length, mzLow, stepSize)
        {
            _mzLow = mzLow; 
            _mzHigh = mzHigh;
            _chargeStateHigh = chargeStateHigh; 
            _chargeStateLow = chargeStateLow;
            _parentDistribution = parentIsotopicDistribution;
            _adductMass = adductMass;
            var chargeStatesList = GenerateChargeStates();
            var binIndices = GetBinIndex(chargeStatesList);

            if (envelopeDistribution.HasValue)
            {
                var intensityMultiples = CreateIntensityMultiplesLogNorm(envelopeDistribution.Value.mu, 
                    envelopeDistribution.Value.sigma);
                CreateChargeStateEnvelope(binIndices, intensityMultiples);
                return; 
            }
            CreateChargeStateEnvelope(binIndices);
            
        }

        private List<double[]> GenerateChargeStates()
        {
            List<double[]> chargeStateTransformedMzs = new();
            
            for (int chargeState = _chargeStateLow; chargeState <= _chargeStateHigh; chargeState++)
            {
                double[] mzArray = new double[_parentDistribution.Masses.Length]; 
                _parentDistribution.Masses.CopyTo(mzArray, 0);
                for (int j = 0; j < mzArray.Length; j++)
                {
                    mzArray[j] = (mzArray[j] + (double)chargeState * _adductMass) / (double)chargeState;
                }
                chargeStateTransformedMzs.Add(mzArray);
            }

            return chargeStateTransformedMzs; 
        }

        //private void FindMinAndMaxPresentChargeState(List<int[]> arrayList)
        //{
        //    bool[] chargeStatePresentBool = new bool[arrayList.Count];
        //    for (int i = 0; i < arrayList.Count; i++)
        //    {
        //         chargeStatePresentBool[i] = arrayList[i].Any(k => k > 0); 
        //    }

        //    int z = 0;
        //    int indexFirstVal = -1;
        //    int indexLastVal = -1; 
        //    while (z < chargeStatePresentBool.Length)
        //    {
        //        if (chargeStatePresentBool[z] && indexFirstVal == -1)
        //        {
        //            indexFirstVal = z; 
        //        }

        //        if (chargeStatePresentBool[^(z+1)] 
        //            && indexLastVal == -1)
        //        {
        //            indexLastVal = chargeStatePresentBool.Length - z; 
        //        }

        //        if (indexFirstVal != -1 && indexLastVal != -1) break; 
        //        z++; 
        //    }
            
        //    _minPresentChargeState = indexFirstVal + _chargeStateLow;
        //    _maxPresentChargeState = indexLastVal + _chargeStateLow - 1;

        //}

        private double[] CreateIntensityMultiplesLogNorm(double mu, double sigma)
        {
            double[] intensityMultiplesArray = new double[_chargeStateHigh - _chargeStateLow];
            // divide charges by max charge
            for (double chargeState = (double)_chargeStateLow; chargeState < (double)_chargeStateHigh; chargeState++)
            {
                double normChargeState = chargeState / _chargeStateHigh;
                intensityMultiplesArray[(int)chargeState - _chargeStateLow] = LogNormal.PDF(mu, sigma, normChargeState);
            }
            return intensityMultiplesArray;
        }

        private List<int[]> GetBinIndex(List<double[]> arrayList)
        {
            List<int[]> indexList = new();
            foreach (double[] array in arrayList)
            {
                int[] binIndex = new int[array.Length];

                for (int i = 0; i < array.Length; i++)
                {
                    if (array[i] < _mzLow || array[i] > _mzHigh)
                    {
                        binIndex[i] = -1;
                        continue;
                    }

                    double approxIndex = (array[i] - _mzLow) / _stepSize;
                    // need to figure out how to round to the correct bin
                    binIndex[i] = (int)approxIndex.Round(0);
                }
                indexList.Add(binIndex);
            }
            return indexList;
        }

        private void CreateChargeStateEnvelope(List<int[]> binIndices, double[] intensityMultiples)
        {
            int intensityMultiplesIndexer = 0;
            int chargeStateIndexer = 0; 
            ApplyElementwise(i => i = 0, Yarray);
            foreach (int[] binArray in binIndices)
            {
                for (int i = 0; i < binArray.Length; i++)
                {
                    if (binArray[i] == -1) continue; 
                    Yarray[binArray[i]] += _parentDistribution.Intensities[i] * intensityMultiples[intensityMultiplesIndexer];
                }

                intensityMultiplesIndexer++;
            }
        }
        private void CreateChargeStateEnvelope(List<int[]> binIndices)
        {
            int intensityMultiplesIndexer = 0;
            int chargeStateIndexer = 0;
            ApplyElementwise(i => i = 0, Yarray);
            foreach (int[] binArray in binIndices)
            {
                for (int i = 0; i < binArray.Length; i++)
                {
                    if (binArray[i] == -1) continue;
                    Yarray[binArray[i]] += _parentDistribution.Intensities[i];
                }

                intensityMultiplesIndexer++;
            }
        }

        public void MakeChargeStateDistributionGaussian(double mean, double stddev)
        {
            for (int i = 0; i < Yarray.Length; i++)
            {
                double gaussianVal = GaussianFunc(Xarray[i], mean, stddev);
                Yarray[i] *= gaussianVal;
            }
        }
        public void MakeChargeStateDistributionWeibull(double shape, double scale)
        {
            for (int i = 0; i < Yarray.Length; i++)
            {
                double val = WeibullFunc(Xarray[i], shape, scale);
                Yarray[i] *= val;
            }
        }

        public void MakeLogNormDistribution(double mu, double sigma)
        {
            for (int i = 0; i < Yarray.Length; i++)
            {
                double val = LogNormFunc(Xarray[i], mu, sigma);
                Yarray[i] *= val;
            }
        }
        protected double GaussianFunc(double d, double mean, double stddev) => Normal.PDF(mean, stddev, d);
        protected double ChiSquareFunc(double d, double x) => ChiSquared.PDF(x, d);
        protected double WeibullFunc(double d, double shape, double scale) => Weibull.PDF(shape, scale, d); 
        protected double LogNormFunc(double d, double mu, double sigma) => LogNormal.PDF(mu, sigma, d);
    }
}
