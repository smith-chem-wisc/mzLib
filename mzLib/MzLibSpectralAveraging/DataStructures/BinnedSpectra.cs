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
        public List<PeakBin> PeakBins { get; set; }
        // TODO: Concurrent dictionaries should be temporary, and these should be sorted Dictionaries. 
        public SortedDictionary<int, double> NoiseEstimates { get; private set; }
        public SortedDictionary<int, double> ScaleEstimates { get; private set; }
        public SortedDictionary<int, double> Weights { get;  set; }
        public double[] Tics { get; private set; }
        public int NumSpectra { get; set; }
        public List<double[]> RecalculatedSpectra => PixelStackListToSpectra(); 
        public int ReferenceSpectra { get; }

        /// <summary>
        /// Creates a binned spectra object given the number of spectra to be averaged
        /// and the reference spectra. 
        /// </summary>
        /// <param name="referenceSpectra">The reference spectra, i.e. the spectra that will be used
        /// to calculate the scales and, subsequently, the weights. Generally, this is the first spectra in the
        /// set to be averaged.</param>
        /// <param name="xArrays">Jagged array of spectra x axes.</param>
        /// <param name="yArrays">Jagged array of intensity values for each spectra.</param>
        ///  <param name="binSize">Size of the bins, e.g. 1.0, 0.5, 0.0001, etc.</param>
        /// <remarks>Attempting to create a method that will select the "best" reference spectra should require a lot of thought.
        /// If you choose to go by noise values, just know that you cannot differentiate a "good" spectrum vs a "bad" spectrum based on
        /// noise estimates alone. A bad spectrum will have a high noise estimate relative to other spectra in the set. But so will a "good" high
        /// SNR spectra: As signal increases, so too does noise. So pick a criteria that doesn't solely use noise if you want to "intelligently" select
        /// the reference spectra. -AVC </remarks>
        public BinnedSpectra(double[][] xArrays, double[][] yArrays, double binSize,
            int referenceSpectra = 0)
        {
            PeakBins = new List<PeakBin>();
            NoiseEstimates = new SortedDictionary<int, double>();
            ScaleEstimates = new SortedDictionary<int, double>();
            Weights = new SortedDictionary<int, double>();
            NumSpectra = xArrays.Length;
            Tics = new double[NumSpectra];
            ReferenceSpectra = referenceSpectra;
            ConsumeSpectra(xArrays, yArrays, binSize);
        }

        /// <summary>
        /// Performs normalization based on the total ion current value for each pixel stack
        /// </summary>
        /// <remarks>Will require refactoring if this C# library is to be expanded to other applications beyond mass spectrometry.</remarks>
        public void PerformNormalization()
        {
            for (int i = 0; i < PeakBins.Count; i++)
            {
                for (int j = 0; j < PeakBins[i].Length; j++)
                {
                    var specOfInterestIndex = PeakBins[i].Peaks[j].SpectraId;
                    var tic = Tics[specOfInterestIndex] == 0 ? 1 : Tics[specOfInterestIndex];
                    double temp = PeakBins[i].GetIntensityAtIndex(j) / tic;
                    PeakBins[i].ModifyPixelIntensity(j, temp);
                }
            }
        }

        [Obsolete("Do not want to recalculate tics, need to save this for recapitulating a normalized spectrum and are no longer merging as described in summary comment")]
        /// <summary>
        /// After spectra are consumed, spectra that had more than one value in a bin have had that value averaged.
        /// Therefore, the original tic value is incorrect and the tic need to be recalculated.
        /// </summary>
        public void RecalculateTics()
        {
            for (int i = 0; i < NumSpectra; i++)
            {
                foreach (var pixelStack in PeakBins)
                {
                    Tics[i] += pixelStack.Intensities[i];
                }
            }
        }

        /// <summary>
        /// Collapses the pixel stack into a single m/z and intensity value. 
        /// </summary>
        /// <remarks>Speed optimized by using a Parallel.Foreach loop.</remarks>
        public double[][] MergeSpectra(SpectralAveragingOptions options)
        {
            // average each pixel stack
            double averageTic = Tics.Average();
            Parallel.ForEach(PeakBins, pixelStack =>
            {
                pixelStack.Average(Weights);
            });

            // return averaged values
            double[] xArray = PeakBins.Select(i => i.UnrejectedMzAverage).ToArray();
            double[] yArray = PeakBins.Select(i => i.MergedIntensityValue).ToArray();

            var combined = new List<(double mz, double intensity)>();
            for (int i = 0; i < xArray.Length; i++)
            {
                if (xArray[i] != 0)
                {
                    combined.Add(options.PerformNormalization
                        ? (xArray[i], yArray[i] * averageTic)
                        : (xArray[i], yArray[i]));
                }
            }

            return new[] { combined.Select(p => p.mz).ToArray(),
                combined.Select(p => p.intensity).ToArray() };
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
        ///  <param name="binSize">Size of the bins, e.g. 1.0, 0.5, 0.0001, etc.</param>
        private void ConsumeSpectra(double[][] xArrays, double[][] yArrays,
            double binSize)
        {
            double min = 100000;
            double max = 0;
            for (int i = 0; i < NumSpectra; i++)
            {
                min = Math.Min(xArrays[i][0], min);
                max = Math.Max(xArrays[i].Max(), max);
            }

            int numberOfBins = (int)Math.Ceiling((max - min) * (1 / binSize));
            // go through each scan and place each (m/z, int) from the spectra into a jagged array

            // 1) find all values of x that fall within a bin.

            List<List<BinValue>> listBinValues = new();

            // assign each peak in the spectra to a bin
            for (int i = 0; i < NumSpectra; i++)
            {
                listBinValues.Add(CreateBinValueList(xArrays[i], yArrays[i], min, binSize, i));
                listBinValues[i].Sort((m, n) => m.Bin.CompareTo(n.Bin));
                Tics[i] = yArrays[i].Sum();
            }

            // iterate through each bin
            for (int binIndex = 0; binIndex <= numberOfBins; binIndex++)
            {
                List<double> xVals = new();
                List<double> yVals = new();
                List<Peak> peaks = new();
                // iterate through each spectra
                for (var spectraId = 0; spectraId < listBinValues.Count; spectraId++)
                {
                    var peakBinDataFromOneSpectra = listBinValues[spectraId];
                    List<BinValue> binValRecord = new();
                    int index = peakBinDataFromOneSpectra.BinarySearch(new BinValue() { Bin = binIndex },
                        new BinValueComparer());
                    if (index < 0)
                    {
                        index = ~index;
                    }

                    int k = 0;
                    // collect all mz and intensity values from a spectra that fall into a bin
                    while (k + index <= peakBinDataFromOneSpectra.Count - 1)
                    {
                        // binary search gets just the first index that it finds, so you 
                        // need to check on the left and right side for elements that 
                        // match the bin value. 

                        if (index + k < peakBinDataFromOneSpectra.Count &&
                            peakBinDataFromOneSpectra[index + k].Bin == binIndex)
                        {
                            binValRecord.Add(peakBinDataFromOneSpectra[index + k]);
                            if (k == 0)
                            {
                                k++;
                                continue;
                            }
                        }

                        if (index - k >= 0 && peakBinDataFromOneSpectra[index - k].Bin == binIndex)
                        {
                            binValRecord.Add(peakBinDataFromOneSpectra[index - k]);
                        }


                        k++;
                    }

                    // collapse values within the bin
                    if (binValRecord.Count == 0)
                    {
                        // TODO: x array values aren't handled correctly. Should be the value of the m/z axis at the particular bin, not a zero. 
                        peaks.Add(new Peak(spectraId, 0, 0, false));
                    }
                    else
                    {
                        var averageMz = Math.Round(binValRecord.Average(p => p.Mz), 8);
                        foreach (var peak in binValRecord)
                        {
                            peaks.Add(new Peak(spectraId, averageMz, peak.Intensity, false));
                        }
                    }
                }

                if (peaks.Average(p => p.Mz) != 0)
                    PeakBins.Add(new PeakBin(peaks));
            }

            PeakBins = PeakBins.OrderBy(p => p.UnrejectedMzAverage).ToList();
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
        /// Creates a List of BinValue records given the x and y axis of a spectra. The BinValue record maintains
        /// the x and y values while also recording the bin that each set of values belongs to.  
        /// </summary>
        /// <param name="xArray">m/z values of a spectra.</param>
        /// <param name="yArray">Intensity value of a spectra.</param>
        /// <param name="min">The lowest value in the xArray.</param>
        /// <param name="binSize">The size of each bin.</param>
        /// <returns></returns>
        private static List<BinValue> CreateBinValueList(double[] xArray, double[] yArray,
            double min, double binSize, int spectraId)
        {
            var binIndices = xArray
                .Select((w, i) => 
                    new { Index = i, Bin = (int)Math.Round((w - min) / binSize, MidpointRounding.AwayFromZero) }); 
            List<BinValue> binValues = new List<BinValue>();
            foreach (var bin in binIndices)
            {
                binValues.Add(new BinValue(bin.Bin, xArray[bin.Index], yArray[bin.Index], spectraId));
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
            double[] results = new double[PeakBins.Count];
            for (int i = 0; i < PeakBins.Count; i++)
            {
                results[i] = PeakBins[i].GetIntensityAtIndex(index); 
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
    internal readonly record struct BinValue(int Bin, double Mz, double Intensity, double SpectraId);

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
