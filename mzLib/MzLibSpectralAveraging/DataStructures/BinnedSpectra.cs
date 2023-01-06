using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Collections.Immutable;
using System.Data;
using System.Linq;
using System.Linq.Expressions;
using System.Reflection.Metadata.Ecma335;
using System.Runtime.InteropServices.ComTypes;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using MathNet.Numerics;
using MathNet.Numerics.RootFinding;
using Nett;

namespace MzLibSpectralAveraging
{
    /// <summary>
    /// Object to organize and unify operations performed on a set of spectra to be averaged.
    /// </summary>
    public class BinnedSpectra
    {
        public List<PixelStack> PixelStacks { get; set; }
        // TODO: Concurrent dictionaries should be temporary, and these should be sorted Dictionaries. 
        public SortedDictionary<int, double> NoiseEstimates { get; private set; }
        public SortedDictionary<int, double> ScaleEstimates { get; private set; }
        public SortedDictionary<int, double> Weights { get; private set; }
        public double[] Tics { get; private set; }
        public int NumSpectra { get; set; }
        public List<double[]> RecalculatedSpectra => PixelStackListToSpectra(); 
        public int ReferenceSpectra { get; }

        /// <summary>
        /// Creates a binned spectra object given the number of spectra to be averaged
        /// and the reference spectra. 
        /// </summary>
        /// <param name="numSpectra">The number of spectra to be averaged.</param>
        /// <param name="referenceSpectra">The reference spectra, i.e. the spectra that will be used
        /// to calculate the scales and, subsequently, the weights. Generally, this is the first spectra in the
        /// set to be averaged.</param>
        /// <remarks>Attempting to create a method that will select the "best" reference spectra should require a lot of thought.
        /// If you choose to go by noise values, just know that you cannot differentiate a "good" spectrum vs a "bad" spectrum based on
        /// noise estimates alone. A bad spectrum will have a high noise estimate relative to other spectra in the set. But so will a "good" high
        /// SNR spectra: As signal increases, so too does noise. So pick a criteria that doesn't solely use noise if you want to "intelligently" select
        /// the reference spectra. -AVC </remarks>
        public BinnedSpectra(int numSpectra, int referenceSpectra = 0)
        {
            PixelStacks = new List<PixelStack>();
            NoiseEstimates = new SortedDictionary<int, double>();
            ScaleEstimates = new SortedDictionary<int, double>();
            Weights = new SortedDictionary<int, double>();
            NumSpectra = numSpectra; 
            Tics = new double[numSpectra];
            ReferenceSpectra = referenceSpectra; 
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        private List<double[]> PixelStackListToSpectra()
        {
            List<double[]> results = new(); 
            for (int i = 0; i < NumSpectra; i++)
            {
                results.Add(PopIntensityValuesFromPixelStackList(i)); 
            }

            return results; 
        }

        /// <summary>
        /// Perform rejection based on the SpectralAveragingOptions based as an argument.
        /// </summary>
        /// <remarks>This method is optimized by using Parallel.ForEach. This is the correct choice because we are iterating over many
        /// pixel stacks in the spectra, exceeding my general criteria for surpassing the additional overhead (roughly 1000 operations) that parallelization requires. </remarks>
        /// <param name="options"></param>
        public void RejectOutliers(SpectralAveragingOptions options)
        {
            Parallel.ForEach(PixelStacks, pixelStack =>
            {
                pixelStack.PerformRejection(options);
            }); 
        }

        /// <summary>
        /// Consumer Spectra override for MsDataScan objects
        /// </summary>
        /// <param name="scans">Scans to be added to this bin of spectra</param>
        /// <param name="binSize">size of bins along x array</param>
        public void ConsumeSpectra(IEnumerable<MsDataScan> scans, double binSize)
        {
            ConsumeSpectra(scans.Select(p => p.MassSpectrum), binSize);
        }

        /// <summary>
        /// Consume Spectra override for MzSpectrum objects
        /// </summary>
        /// <param name="spectra">Spectra to be added to this bin of spectra</param>
        /// <param name="binSize">size of bins along x array</param>
        public void ConsumeSpectra(IEnumerable<MzSpectrum> spectra, double binSize)
        {
            ConsumeSpectra(spectra.Select(p => p.XArray).ToArray(),
                spectra.Select(p => p.YArray).ToArray(), spectra.Count(), binSize);
        }

        /// <summary>
        /// Takes jagged arrays of x axis and y axis values and converts them into PixelStack objects. Further ensures that bins
        /// with multiple values in a given spectra are managed appropriately and that empty bins are zero filled in the y array and have the correct
        /// value in the x array. 
        /// </summary>
        /// <remarks>Heavily optimized for speed. Original method took 1.8 minutes in a test case, and this
        /// current iteration takes ~250 ms for the same test. Uses a custom record to store mz value, intensity value,
        /// and bin number so that I could use a BinarySearch override that uses a custom comparer to increase search speed.
        /// Previous version used a .Where() which was the bottleneck of the method.</remarks>
        /// <param name="xArrays">Jagged array of spectra x axes.</param>
        /// <param name="yArrays">Jagged array of intensity values for each spectra.</param>
        /// <param name="numSpectra">Number of spectra to be averaged.</param>
        /// <param name="binSize">Size of the bins, e.g. 1.0, 0.5, 0.0001, etc.</param>
        public void ConsumeSpectra(double[][] xArrays, double[][] yArrays,
            int numSpectra, double binSize)
        {
            double min = 100000;
            double max = 0;
            for (int i = 0; i < numSpectra; i++)
            {
                min = Math.Min(xArrays[i][0], min);
                max = Math.Max(xArrays[i].Max(), max);
            }

            int numberOfBins = (int)Math.Ceiling((max - min) * (1 / binSize));
            // go through each scan and place each (m/z, int) from the spectra into a jagged array

            // 1) find all values of x that fall within a bin.

            List<List<BinValue>> listBinValues = new();

            for (int i = 0; i < numSpectra; i++)
            {
               listBinValues.Add(CreateBinValueList(xArrays[i], yArrays[i], min, binSize));
               listBinValues[i].Sort((m, n) => m.Bin.CompareTo(n.Bin));
            }

            for (int i = 0; i < numberOfBins; i++)
            {
                List<double> xVals = new();
                List<double> yVals = new();
                foreach(var binValList in listBinValues)
                {
                    List<BinValue> binValRecord = new();
                    int index = binValList.BinarySearch(new BinValue() { Bin = i }, new BinValueComparer());
                    if (index < 0)
                    {
                        index = ~index;
                    }
                    int k = 0;
                    while (k + index <= binValList.Count - 1)
                    {
                        // binary search gets just the first index that it finds, so you 
                        // need to check on the left and right side for elements that 
                        // match the bin value. 

                        if (index + k < binValList.Count && binValList[index + k].Bin == i)
                        {
                            binValRecord.Add(binValList[index + k]);
                            if (k == 0)
                            {
                                k++; 
                                continue;
                            } 
                        }

                        if (index - k > 0 && binValList[index - k].Bin == i )
                        {
                            binValRecord.Add(binValList[index - k]);
                        }

                        k++; 
                    }

                    if (binValRecord.Count() > 1)
                    {
                        xVals.Add(binValRecord.Average(m => m.Mz));
                        yVals.Add(binValRecord.Average(m => m.Intensity));
                        continue; 
                    }
                    if (binValRecord.Count == 0)
                    {
                        // TODO: x array values aren't handled correctly. Should be the value of the m/z axis at the particular bin, not a zero. 
                        xVals.Add(0); 
                        yVals.Add(0);
                        continue; 
                    }
                    xVals.Add(binValRecord.First().Mz);
                    yVals.Add(binValRecord.First().Intensity);
                }

                //if (xVals.Count > NumSpectra || yVals.Count > NumSpectra)
                //{

                //}

                PixelStacks.Add(new PixelStack(xVals, yVals));
            }
        }

        /// <summary>
        /// Performs normalization based on the total ion current value for each pixel stack.
        /// </summary>
        /// <remarks>Will require refactoring if this C# library is to be expanded to other applications beyond mass spectrometry.</remarks>
        public void PerformNormalization()
        {
            for (int i = 0; i < PixelStacks.Count; i++)
            {
                for (int j = 0; j < PixelStacks[i].Length; j++)
                {
                    double temp = PixelStacks[i].GetIntensityAtIndex(j) / Tics[j];
                    PixelStacks[i].ModifyPixelIntensity(j, temp);
                }
            }
        }
        /// <summary>
        /// Calculates noise estimates for each spectra using multi-resolution support noise estimation.
        /// </summary>
        /// <param name="waveletType">wavelet to be used in MRS noise estimation.</param>
        /// <param name="epsilon">Noise estimate convergence to be reached before returning the noise estimate.</param>
        /// <param name="maxIterations">Maximum number of iterations to be performed in the MRS noise estimation before returning.</param>
        public void CalculateNoiseEstimates(WaveletType waveletType = WaveletType.Haar, 
            double epsilon = 0.01, int maxIterations = 25)
        {
            ConcurrentDictionary<int, double> tempConcurrentDictionary = new(); 
            RecalculatedSpectra
                .Select((w, i) => new { Index = i, Array = w })
                .AsParallel()
                .ForAll(x =>
                {
                    bool success = MRSNoiseEstimator.MRSNoiseEstimation(x.Array, epsilon, out double noiseEstimate,
                        waveletType: waveletType, maxIterations: maxIterations);
                    if (!success || double.IsNaN(noiseEstimate))
                    {
                        noiseEstimate = BasicStatistics.CalculateStandardDeviation(x.Array);
                    }
                    tempConcurrentDictionary.TryAdd(x.Index, noiseEstimate);
                });
            NoiseEstimates = new SortedDictionary<int, double>(tempConcurrentDictionary); 
        }
        /// <summary>
        /// Calculates the estimates of scale for each spectra in this object. Scale is determined by
        /// taking the square root of the biweight midvariance. 
        /// </summary>
        public void CalculateScaleEstimates()
        {
            ConcurrentDictionary<int, double> tempScaleEstimates = new();
            RecalculatedSpectra
                .Select((w,i) => new {Index = i, Array = w})
                .AsParallel().ForAll(x =>
                {
                    double scale = Math.Sqrt(BiweightMidvariance(x.Array));
                    tempScaleEstimates.TryAdd(x.Index, Math.Sqrt(scale));
                });
            ScaleEstimates = new SortedDictionary<int, double>(tempScaleEstimates); 
        }
        /// <summary>
        /// Calculates the median absolute deviation from the median, which is then used in
        /// calculating the biweight midvariance. Original algorithm found here:
        /// https://pixinsight.com/doc/tools/ImageIntegration/ImageIntegration.html#__equation_26__
        /// </summary>
        /// <param name="array">Array of values to be calculated.</param>
        /// <returns>The median absolute deviation from median.</returns>
        private double MedianAbsoluteDeviationFromMedian(double[] array)
        {
            double arrayMedian = BasicStatistics.CalculateMedian(array);
            double[] results = new double[array.Length];
            for (int j = 0; j < array.Length; j++)
            {
                results[j] = Math.Abs(array[j] - arrayMedian);
            }

            return BasicStatistics.CalculateMedian(results); 
        }
        /// <summary>
        /// Calcultes the biweight midvariance for an array. Algorithm orignally found here:
        /// https://pixinsight.com/doc/tools/ImageIntegration/ImageIntegration.html#__equation_27__
        /// </summary>
        /// <param name="array">Array of doubles.</param>
        /// <returns>The biweight midvariance.</returns>
        private double BiweightMidvariance(double[] array)
        {
            double[] y_i = new double[array.Length];
            double[] a_i = new double[array.Length]; 
            double MAD_X = MedianAbsoluteDeviationFromMedian(array);
            double median = BasicStatistics.CalculateMedian(array); 
            for (int i = 0; i < y_i.Length; i++)
            {
                y_i[i] = (array[i] - median) / (9d * MAD_X);
                if (y_i[i] < 1d)
                {
                    a_i[i] = 1d; 
                }
                else
                {
                    a_i[i] = 0; 
                }
            }

            // biweight midvariance calculation

            double denomSum = 0;
            double numeratorSum = 0; 
            for (int i = 0; i < y_i.Length; i++)
            {
                numeratorSum += a_i[i] * Math.Pow(array[i] - median, 2) * Math.Pow(1 - y_i[i] * y_i[i], 4); 
                denomSum += a_i[i] * (1 - 5 * y_i[i] * y_i[i]) * (1 - y_i[i] * y_i[i]);
            }

            return (double)y_i.Length * numeratorSum / Math.Pow(Math.Abs(denomSum), 2); 
        }
        /// <summary>
        /// Given the noise estimates and the scale estimates, calculates the weight given to
        /// each spectra when averaging using w_i = 1 / (k * noise_estimate)^2,
        /// where k = scaleEstimate_reference / scaleEstimate_i
        /// </summary>
        public void CalculateWeights()
        {
            double referenceScale = ScaleEstimates[ReferenceSpectra]; 
            foreach (var entry in NoiseEstimates)
            {
                var successScale = ScaleEstimates.TryGetValue(entry.Key,
                    out double scale);
                if (!successScale) continue;

                var successNoise = NoiseEstimates.TryGetValue(entry.Key,
                    out double noise);
                if (!successNoise) continue;
                
                double k = referenceScale / scale; 

                double weight = 1d / Math.Pow((k * noise), 2);

                Weights.TryAdd(entry.Key, weight);
            }
        }

        /// <summary>
        /// After spectra are consumed, spectra that had more than one value in a bin have had that value averaged.
        /// Therefore, the original tic value is incorrect and the tic need to be recalculated.
        /// </summary>
        public void RecalculateTics()
        {
            
            for (int i = 0; i < NumSpectra; i++)
            {
                foreach (var pixelStack in PixelStacks)
                {
                    Tics[i] += pixelStack.Intensity[i];
                }
            }
        }
        /// <summary>
        /// Collapses the pixel stack into a single m/z and intensity value. 
        /// </summary>
        /// <remarks>Speed optimized by using a Parallel.Foreach loop.</remarks>
        public void MergeSpectra()
        {
            Parallel.ForEach(PixelStacks, pixelStack =>
            {
                pixelStack.Average(Weights);
            });
        }
        /// <summary>
        /// Utility method to expose the merged x and y array. 
        /// </summary>
        /// <returns>Jagged array where double[0] contains the averaged x array
        /// and double[1] contains the averaged y array. </returns>
        public double[][] GetMergedSpectrum()
        {
            double[] xArray = PixelStacks.Select(i => i.MergedMzValue).ToArray();
            double[] yArray = PixelStacks.Select(i => i.MergedIntensityValue).ToArray();
            return new[] { xArray, yArray };
        }

        /// <summary>
        /// Creates a List of BinValue records given the x and y axis of a spectra. The BinValue record maintains
        /// the x and y values while also recording the bin that each set of values belongs to.  
        /// </summary>
        /// <param name="xArray">m/z values of a spectra.</param>
        /// <param name="yArray">Intensity value of a spectra.</param>
        /// <param name="min">The lowest value in the xArray.</param>
        /// <param name="binSize">The size of each bin.</param>
        /// <returns></returns>
        private static List<BinValue> CreateBinValueList(double[] xArray, double[] yArray,
            double min, double binSize)
        {
            var binIndices = xArray
                .Select((w, i) => 
                    new { Index = i, Bin = (int)Math.Round((w - min) / binSize, MidpointRounding.AwayFromZero) }); 
            List<BinValue> binValues = new List<BinValue>();
            foreach (var bin in binIndices)
            {
                binValues.Add(new BinValue(bin.Bin, xArray[bin.Index], yArray[bin.Index]));
            }
            return binValues;
        }
        /// <summary>
        /// Pixel stack objects make rejection and averaging calculation easy, but makes linear calculations like
        /// noise estimation, scale estimation, and weight calculation difficult. This method transforms the PixelStacks into
        /// linear spectra. 
        /// </summary>
        /// <param name="index">The index of the spectra. Bounded from 0 to less than number of spectra to be averaged.</param>
        /// <returns>A linear spectra from PixelStacks property.</returns>
        private double[] PopIntensityValuesFromPixelStackList(int index)
        {
            double[] results = new double[PixelStacks.Count];
            for (int i = 0; i < PixelStacks.Count; i++)
            {
                results[i] = PixelStacks[i].GetIntensityAtIndex(index); 
            }

            return results; 
        }
    }

    /// <summary>
    /// Record type used to facilitate bin, mz, and intensity matching. 
    /// </summary>
    /// <param name="Bin">Integer bin number.</param>
    /// <param name="Mz">Mz value.</param>
    /// <param name="Intensity">Intensity value.</param>
    internal readonly record struct BinValue(int Bin, double Mz, double Intensity);

    /// <summary>
    /// Custom comparer to use in override for List.BinarySearch() that accepts a custom comparer as an argument. 
    /// </summary>
    internal class BinValueComparer : IComparer<BinValue>
    {
        public int Compare(BinValue x, BinValue y)
        {
            return x.Bin.CompareTo(y.Bin); 
        }
    }
}
