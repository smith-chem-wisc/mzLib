using MathNet.Numerics.Distributions;
using Plotly.NET;

namespace SimulatedData
{
    public abstract class SimulatedData
    {
        public double[] Xarray { get; protected set; }
        public double[] Yarray { get; protected set; }
        public int Length { get; }
        private double _stepSize { get; }
        private double _startValue { get; }

        protected SimulatedData(int length, double startValue, double stepSize)
        {
            Length = length;
            _stepSize = stepSize;
            _startValue = startValue;
            Xarray = new double[length];
            Yarray = new double[length]; 
            FillArray(Xarray);
            FillArray(Yarray);
        }

        protected SimulatedData()
        {

        }

        private void FillArray(double[] array)
        {
            for (int i = 0; i < Length; i++)
            {
                array[i] = _startValue + _stepSize * i;
            }
        }

        /// <summary>
        /// Applies a function to each element in an array. 
        /// </summary>
        /// <param name="function"></param>
        public void ApplyElementwise(Func<double, double> function, double[] array)
        {
            for (int i = 0; i < Length; i++)
            {
                array[i] = function(array[i]);
            }
        }

        public void ElementwiseArrayAddition(SimulatedData peak)
        {
            for (int i = 0; i < peak.Length; i++)
            {
                Yarray[i] += peak.Yarray[i];
            }
        }

        public void ElementwiseArrayAddition(double[] yArray)
        {
            if (yArray.Length != Length)
                throw new ArgumentException("Incompatible lengths for addition");
            for (int i = 0; i < yArray.Length; i++)
            {
                Yarray[i] += yArray[i];
            }
        }

        public static SimulatedData operator +(SimulatedData a, SimulatedData b)
        {
            if (a.Length != b.Length) throw new ArgumentException("Invalid lengths for addition.");
            if (Math.Abs(a._stepSize - b._stepSize) > 0.001)
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

        public void AddHighFrequencyNoise(Normal noiseDistribution)
        {
            for (int i = 0; i < Yarray.Length; i++)
            {
                Yarray[i] += noiseDistribution.Sample(); 
            }
        }

        public void AddHighFrequencyNoise(Normal noiseDistribution, 
            List<(int Min, int Max)>? excludedZones)
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

        public void AddLowFrequencyNoise(LowFrequencyNoiseParameters noiseParams)
        {
            Random random = new Random();

            int numberOfRandomPeaks = random.Next(noiseParams.PeakNumberLimitLow, noiseParams.PeakNumberLimitHigh);
            for (int i = 0; i < numberOfRandomPeaks; i++)
            {
                Random randomInnerLoop = new();
                int peakLocation = randomInnerLoop.Next((int)noiseParams.PeakLocationLimitLow,
                    (int)noiseParams.PeakLocationLimitHigh);
                // because we have the possibility of adding many peaks to the actual signal peak, 
                // we may want to exclude noise values being added to the true signal. 
                // the while loop generates new peak location values if the original value falls 
                // within an excluded zone. 
                if (noiseParams.ExcludedZone.HasValue)
                {
                    while (peakLocation > noiseParams.ExcludedZone.Value.Min
                           && peakLocation < noiseParams.ExcludedZone.Value.Max)
                    {
                        peakLocation = randomInnerLoop.Next((int)noiseParams.PeakLocationLimitLow,
                            (int)noiseParams.PeakLocationLimitHigh);
                    }
                }

                double peakWidth =
                    randomInnerLoop.NextDouble(noiseParams.PeakWidthLimitLow, noiseParams.PeakWidthLimitHigh);
                double peakIntensity = randomInnerLoop.NextDouble(noiseParams.PeakIntensityLimitLow,
                    noiseParams.PeakIntensityLimitHigh);

                double[] yArrayToAdd = Enumerable.Range(0, Length)
                    .Select(z => z * _stepSize + _startValue)
                    .ToArray();
                switch (noiseParams.PeakShapeOptions)
                {
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

        private double NextDouble(Random random, int low, int high)
        {
            return random.NextDouble() * ((double)high - (double)low) + (double)low;
        }

        private double NextDouble(Random random, double low, double high)
        {
            return random.NextDouble() * (high - low) + low;
        }
        private double GaussianFunc(double d, double mean, double stddev, double intensity) => intensity * Normal.PDF(mean, stddev, d);
    }
}