using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using Easy.Common.Extensions;
using MathNet.Numerics;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.IntegralTransforms;

namespace SimulatedData
{
    public class SimulatedChargeStateEnvelope : SimulatedData
    {
        private double MzLow { get; }
        private double MzHigh { get; }
        private int ChargeStateLow { get; }
        private int ChargeStateHigh { get; }
        private IsotopicDistribution ParentDistribution { get; }
        private double AdductMass { get; }
        public SimulatedChargeStateEnvelope(double mzLow, double mzHigh, double stepSize, int length,
            int chargeStateLow, int chargeStateHigh, IsotopicDistribution parentIsotopicDistribution,
            (double mu, double sigma)? envelopeDistribution = null, double adductMass = 1.007276) 
        : base(length, mzLow, stepSize)
        {
            MzLow = mzLow; 
            MzHigh = mzHigh;
            ChargeStateHigh = chargeStateHigh; 
            ChargeStateLow = chargeStateLow;
            ParentDistribution = parentIsotopicDistribution;
            AdductMass = adductMass;
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
            
            for (int chargeState = ChargeStateLow; chargeState <= ChargeStateHigh; chargeState++)
            {
                double[] mzArray = new double[ParentDistribution.Masses.Length]; 
                ParentDistribution.Masses.CopyTo(mzArray, 0);
                for (int j = 0; j < mzArray.Length; j++)
                {
                    mzArray[j] = (mzArray[j] + (double)chargeState * AdductMass) / (double)chargeState;
                }
                chargeStateTransformedMzs.Add(mzArray);
            }

            return chargeStateTransformedMzs; 
        }
        private double[] CreateIntensityMultiplesLogNorm(double mu, double sigma)
        {
            double[] intensityMultiplesArray = new double[ChargeStateHigh - ChargeStateLow + 1];
            // divide charges by max charge
            for (double chargeState = (double)ChargeStateLow; chargeState <= (double)ChargeStateHigh; chargeState++)
            {
                double normChargeState = chargeState / (double)ChargeStateHigh;
                intensityMultiplesArray[(int)chargeState - ChargeStateLow] = LogNormal.PDF(mu, sigma, normChargeState);
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
                    if (array[i] < MzLow || array[i] > MzHigh)
                    {
                        binIndex[i] = -1;
                        continue;
                    }

                    double approxIndex = (array[i] - MzLow) / _stepSize;
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
                    Yarray[binArray[i]] += ParentDistribution.Intensities[i] * intensityMultiples[intensityMultiplesIndexer];
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
                    if (binArray[i] > Yarray.Length - 1) continue;
                    Yarray[binArray[i]] += ParentDistribution.Intensities[i];
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
		/// <summary>
		/// Creates a lower resolution spectra than the original by convoluting a Gaussian kernel with the Y array. This method implements convolution by first applying a FFT
		/// to both the kernel and the Yarray, performs elementwise multiplication of the kernel and signal, and finally perform the inverse FFT to get the signal back to the
		/// frequency domain. 
		/// </summary>
		/// <param name="width">THIS VALUE MUST BE ODD. The width of the Gaussian kernel. Increasing this value decreases the resolution of the spectra. The width must be odd because the Gaussian
		/// blur method shifts the Yarray spectra to the ri</param>
		/// <param name="sigma">Standard deviation of the Gaussian used in the blur. Increasing value decreases resolution of the blurred spectra.</param>
		/// <param name="iterations">Number of convolutions with the Gaussian kernel to perform. Using more iterations with smaller sigma values and kernel widths increases the
		/// time it takes to run the method, but improves the realism of the final spectra.</param>
		/// <remarks>Note that the blurred spectra is shifted as a result of the convolution, so the end result is shifted back to the correct position by a combination of
		/// zero padding and empirical selection of the final, blurred array indices to correct for the shift.</remarks>
		public void Blur(int width, double sigma, int iterations)
        {
	        if (width % 2 != 1) throw new ArgumentException("Width must be odd."); 
            for (int i = 0; i < iterations; i++)
            {
                var window = Window.Gauss(width, sigma);
                var outputArray = ConvolveWithWindow(Yarray, window);
                // convert back to real by getting the magnitude
                // because the circular convolution causes some shift, I needed to pad the arrays to the right length. 
                // The array selection returns the shift-corrected Y array. I determined it empirically because I was tired. 
                Yarray = outputArray.Select(x => (double)x.Magnitude).ToArray()[(width/2)..(Yarray.Length + (width - 1)/2)];
            }
        }

        public void NormalizeToMaxIntensity()
        {
            double max = Yarray.Max(); 
            for (int i = 0; i < Yarray.Length; i++)
            {
                Yarray[i] /= max;
            }
        }
        // convolution by fft and then elementwise multiplication followed by 
        // fft back to original domain
        private Complex32[] ConvolveWithWindow(double[] signal, double[] window)
        {
            // zero-pad window to length 
            // remember that in C#, arrays are initialized with zeroes, so 
            // all I need to do is copy in the original data
            double[] paddedWindow = new double[signal.Length + window.Length];
            double[] paddedSignal = new double[signal.Length + window.Length];
            Buffer.BlockCopy(window, 0, paddedWindow, 0, sizeof(double) * window.Length);
            Buffer.BlockCopy(signal, 0, paddedSignal, 0, sizeof(double) * signal.Length);
            // convert to complex
            Complex32[] signalComplex = paddedSignal.Select(i => new Complex32((float)i, 0)).ToArray();
            Complex32[] windowComplex = paddedWindow.Select(i => new Complex32((float)i, 0)).ToArray(); 
            // perform fourier transforms of the arrays
            Fourier.Forward(signalComplex);
            Fourier.Forward(windowComplex);
            // perform the elementwise mutliplcation
            for (int i = 0; i < signalComplex.Length; i++)
            {
                signalComplex[i] *= windowComplex[i];
            }

            Fourier.Inverse(signalComplex); 
            return signalComplex; 
        }
        protected double GaussianFunc(double d, double mean, double stddev) => Normal.PDF(mean, stddev, d);
        protected double WeibullFunc(double d, double shape, double scale) => Weibull.PDF(shape, scale, d); 
        protected double LogNormFunc(double d, double mu, double sigma) => LogNormal.PDF(mu, sigma, d);
    }
}
