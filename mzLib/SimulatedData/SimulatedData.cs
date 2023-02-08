using MathNet.Numerics.Distributions;
using System.Collections.Generic;
using Easy.Common.Extensions;
using MassSpectrometry;
using Plotly.NET.CSharp;
using SpectralAveraging;

namespace SimulatedData
{
    /// <summary>
    /// Abstract class from which all simulated data objects can be derived. Provides methods for adding high frequency noise, low frequency noise, and operators
    /// to add, subtract, or multiply classes derived from SimulatedData and other SImulatedData objects or scalar values. 
    /// </summary>
    public abstract class SimulatedData
    {
        /// <summary>
        /// m/z array. 
        /// </summary>
        public double[] Xarray { get; protected set; }
        /// <summary>
        /// Array of intensity values.
        /// </summary>
        public double[] Yarray { get; protected set; }
        protected int Length { get; }
        protected double _stepSize { get; }
        private double _startValue { get; }

        protected SimulatedData(int length, double startValue, double stepSize)
        {
            Length = length;
            _stepSize = stepSize;
            _startValue = startValue;
            Xarray = new double[length];
            Yarray = new double[length];
            FillArray(Xarray);
            // when I had a derived class that included a single Gaussian peak, the way it was implemented 
            // required that the Yarray be filled with the same Xarray values. This isn't required for SimulatedData, and the 
            // implementation of my Gaussian class was pretty quick and dirty. SimulatedChargeStateEnvelope is a much more clean and 
            // well-thought out chunk of code. 
            //FillArray(Yarray);
        }
        /// <summary>
        /// Fills m/z values based on the step size. Currently only does a linear step, however to truly replicate orbitrap m/z,
        /// you need to have m/z values that progress from some minimum value to some maximum values as a function of Sqrt(m/z). 
        /// </summary>
        /// <param name="array"></param>
        protected void FillArray(double[] array)
        {
            //TODO: Implement option for m/z axis that progress as a function of the square root of m/z
            for (int i = 0; i < Length; i++)
            {
                array[i] = _startValue + _stepSize * i;
            }
        }

        /// <summary>
        /// Applies a function to each element in an array. 
        /// </summary>
        /// <param name="function">Function that takes a double as input and returns a double as an output.</param>
        /// <param name ="array">Array of doubles that will have the function in the first parameter applied to each element individually.</param>
        protected void ApplyElementwise(Func<double, double> function, double[] array)
        {
            for (int i = 0; i < Length; i++)
            {
                array[i] = function(array[i]);
            }
        }

		#region Noise Addition 
        /// <summary>
        /// Adds noise to every intensity value based on a normal distribution.
        /// </summary>
        /// <param name="noiseDistribution"></param>
        /// <param name="excludedZones">Optional parameter that will prevent noise values from being added to regions of signal.</param>>
        public void AddHighFrequencyNoise(Normal noiseDistribution, 
            List<(int Min, int Max)>? excludedZones = null)
        {
            if (excludedZones == null)
            { 
                AddHighFrequencyNoise(noiseDistribution);
                return; 
            }

            excludedZones.Sort();
            int excludedZonesIndex = 0; 
            for (int i = 0; i < Yarray.Length; i++)
            {
                if (i >= excludedZones[excludedZonesIndex].Min 
                    && i <= excludedZones[excludedZonesIndex].Max)
                {
                    continue;
                }

                Yarray[i] += noiseDistribution.Sample(); 
            }

        }
        /// <summary>
        /// Adds a random number of noise peaks to random locations in the spectra with random widths and heights based on the parameters set in the LF noise parameters class.
        /// </summary>
        /// <param name="noiseParams"></param>
        public void AddLowFrequencyNoise(LowFrequencyNoiseParameters noiseParams)
        {
            double range = noiseParams.PeakLocationLimitHigh - noiseParams.PeakLocationLimitLow;
            Random random = new Random();

            int numberOfRandomPeaks = random.Next(noiseParams.PeakNumberLimitLow, noiseParams.PeakNumberLimitHigh);
            for (int i = 0; i < numberOfRandomPeaks; i++)
            {
                Random randomInnerLoop = new();
                double peakLocation = randomInnerLoop.NextDouble() * range + noiseParams.PeakLocationLimitLow;
                int peakIndex = Xarray.IndexOf(peakLocation);
                // NOTE: Removed this option as I didn't think it was fair to prevent noise being added to signal. 
                // because we have the possibility of adding many peaks to the actual signal peak, 
                // we may want to exclude noise values being added to the true signal. 
                // the while loop generates new peak location values if the original value falls 
                // within an excluded zone. 
                //if (noiseParams.ExcludedZone.HasValue)
                //{
                //    while (peakLocation > noiseParams.ExcludedZone.Value.Min
                //           && peakLocation < noiseParams.ExcludedZone.Value.Max)
                //    {
                //        peakLocation = randomInnerLoop.Next((int)noiseParams.PeakLocationLimitLow,
                //            (int)noiseParams.PeakLocationLimitHigh);
                //    }
                //}

                double peakWidth =
                    randomInnerLoop.NextDouble(noiseParams.PeakWidthLimitLow, noiseParams.PeakWidthLimitHigh);
                double peakIntensity = randomInnerLoop.NextDouble(noiseParams.PeakIntensityLimitLow,
                    noiseParams.PeakIntensityLimitHigh);

                double[] yArrayToAdd = Enumerable.Range(0, Length)
                    .Select(z => z * _stepSize + _startValue)
                    .ToArray();
                switch (noiseParams.PeakShapeOptions)
                {
                    // Lorentzian peaks are possible, and in the future, we may want to add more peak shape types for added realism.
                    case PeakShapeOptions.Gaussian:
                        {
                            for (int k = 0; k < yArrayToAdd.Length; k++)
                            {
                                yArrayToAdd[k] = GaussianFunc(yArrayToAdd[k], peakLocation, peakWidth, peakIntensity);
                            }
                        }
                        break;
                    case PeakShapeOptions.None:
                        break;
                }

                Yarray = (DoubleArray)Yarray + yArrayToAdd;
            }
        }
        #endregion

        private double GaussianFunc(double d, double mean, double stddev, double intensity) => intensity * Normal.PDF(mean, stddev, d);
        /// <summary>
        /// Calculates the standard deviation of the noise for derived classes of SimulatedData. 
        /// </summary>
        /// <param name="epsilon">The level of convergergence required before standard deviation estimates stops iterating. </param>
        /// <param name="noiseEstimate">Output noise estimate.</param>
        /// <param name="maxIterations">Maximum number of iterations to perform.</param>
        /// <returns></returns>
        protected bool TryMrsNoiseEstimation(double epsilon, out double noiseEstimate,
            int maxIterations = 25)
        {
            bool mrsNoiseEstimationSuccess = MRSNoiseEstimator.MRSNoiseEstimation(Yarray,
                epsilon, out double tempNoiseEstimate, maxIterations);
            if (mrsNoiseEstimationSuccess)
            {
                noiseEstimate = tempNoiseEstimate;
            }
            else
            {
                noiseEstimate = BasicStatistics.CalculateStandardDeviation(Yarray);
            }

            return mrsNoiseEstimationSuccess;
        }
        /// <summary>
        /// Normalizes derived classes' Yarrays to the sum of the entire Yarray. 
        /// </summary>
        public void NormalizeByTic()
        {
            double sumOfYArray = Yarray.Sum();
            ApplyElementwise(i => i / sumOfYArray, Yarray);
        }
        /// <summary>
        /// Converts derived objects to an MzSpectrum object. 
        /// </summary>
        /// <param name="shouldCopy"></param>
        /// <returns></returns>
        public MzSpectrum ToMzSpectrum(bool shouldCopy)
        {
	        return new MzSpectrum(Xarray, Yarray, shouldCopy); 
        }
        /// <summary>
        /// Converts derived objects to an MsDataScan object. 
        /// </summary>
        /// <returns></returns>
        /// <exception cref="NotImplementedException"></exception>
        public MsDataScan ToMsDataScan()
        {
            throw new NotImplementedException();
        }

        // Operators allow you to easily add multiple simulated data objects so you can easily replicate spectra with multiple, co-eluting charge state envelopes.
        // Note that the SimulatedData objects must be of the same length and step size.
		#region Operators
		public static SimulatedData operator +(SimulatedData a, SimulatedData b)
		{
			if (a.Length != b.Length) throw new ArgumentException("Invalid lengths for addition.");
			if (Math.Abs(a._stepSize - b._stepSize) > 0.00001)
				throw new ArgumentException("Incompatible spacing in data");

			for (int i = 0; i < a.Length; i++)
			{
				a.Yarray[i] += b.Yarray[i];
			}

			return a;
		}

		public static SimulatedData operator +(SimulatedData a, double[] b)
		{
			if (a.Length != b.Length) throw new ArgumentException("Invalid lengths for addition.");

			for (int i = 0; i < a.Length; i++)
			{
				a.Yarray[i] += b[i];
			}

			return a;
		}
        /// <summary>
        /// Allows multiplication of two simulated data objects.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException">Thrown in the event that either the lengths of simulated data objects are invalid or
        /// the step size is not the same.</exception>
		public static SimulatedData operator *(SimulatedData a, SimulatedData b)
		{
			if (a.Length != b.Length) throw new ArgumentException("Invalid lengths for addition.");
			if (Math.Abs(a._stepSize - b._stepSize) > 0.001)
				throw new ArgumentException("Incompatible spacing in data");

			for (int i = 0; i < a.Length; i++)
			{
				a.Yarray[i] *= b.Yarray[i];
			}

			return a;
		}

        public static SimulatedData operator -(SimulatedData a, SimulatedData b)
		{
			if (a.Length != b.Length) throw new ArgumentException("Invalid lengths for addition.");
			if (Math.Abs(a._stepSize - b._stepSize) > 0.001)
				throw new ArgumentException("Incompatible spacing in data");

			for (int i = 0; i < a.Length; i++)
			{
				a.Yarray[i] -= b.Yarray[i];
			}

			return a;
		}

		public static SimulatedData operator *(SimulatedData a, double scalar)
		{
			for (int i = 0; i < a.Length; i++)
			{
				a.Yarray[i] *= scalar;
			}

			return a;
		}

		public static SimulatedData operator /(SimulatedData a, double scalar)
		{
			for (int i = 0; i < a.Length; i++)
			{
				a.Yarray[i] /= scalar;
			}

			return a;
		}
		#endregion
	}
}