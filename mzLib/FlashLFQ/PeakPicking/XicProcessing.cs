using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using Easy.Common.Extensions;
using FlashLFQ.PeakPicking;
using MathNet.Numerics;
using MathNet.Numerics.IntegralTransforms;
using MathNet.Numerics.Interpolation;
using MathNet.Numerics.Statistics;
using System.Runtime.CompilerServices;

[assembly: InternalsVisibleTo("MyTests")]

namespace FlashLFQ
{
    public static class XicProcessing
    {
        /// <summary>
        /// Aligns two XICs and reports the relative shift in time using Fast Fourier transforms.
        /// This function performs better if the XIC contains signal from before and after the
        /// last chromatographic peaks of interest (i.e., longer XICs are better). Alignment will fail if the
        /// magnitude of the RT shift is greater than 1/4 the RT span of either XIC
        /// The XICs are up-sampled to allow for sub-pixel resolution (one Peaks datapoint = one pixel).
        /// </summary>
        /// <param name="refPeaks">List of peaks to be used as a reference, ordered by retention time</param>
        /// <param name="expPeaks">List of peaks to be aligned, ordered by retention time</param>
        /// <param name="resolution">Up-sampling resolution. Higher values allow for more precise shifts</param>
        /// <returns>The retention time correction to align the experimental and reference XICs</returns>
        /// (Suggested reading: Convolution theorem, https://dsp.stackexchange.com/questions/51409/maximum-of-cross-correlation-not-moving)
        public static double AlignPeaks(List<IndexedMassSpectralPeak> refPeaks, List<IndexedMassSpectralPeak> expPeaks, int resolution = 100)
        {
            // First step is to pad the XICs, adding zeros to either end
            // This could be further optimized by padding at the array level w/o creating new peak objects
            refPeaks = PadPeaks(refPeaks);
            expPeaks = PadPeaks(expPeaks);

            // Creates a spline, or a smoothed representation of the Peaks which allows for up-sampling of the data
            // Unsure about whether to use cubic or linear. Can maybe be a delegate argument?
            var referenceSpline = LinearSpline.InterpolateSorted(
                refPeaks.Select(p => p.RetentionTime).ToArray(), refPeaks.Select(p => p.Intensity).ToArray());
            var expSpline = LinearSpline.InterpolateSorted(
                expPeaks.Select(p => p.RetentionTime).ToArray(), expPeaks.Select(p => p.Intensity).ToArray());

            var startTime = Math.Max(refPeaks.First().RetentionTime, expPeaks.First().RetentionTime);
            var endTime = Math.Min(refPeaks.Last().RetentionTime, expPeaks.Last().RetentionTime);
            var timeSpan = endTime - startTime;
            int minArrayLength = Math.Min(expPeaks.Count, refPeaks.Count);
            
            // Create the arrays of interpolated data
            int interpArrayLength = resolution * minArrayLength;
            double[] interpolatedTimes = Generate.LinearSpaced(length: interpArrayLength, startTime, endTime);
            double[] refInterpolatedIntensities = interpolatedTimes.Select(t => referenceSpline.Interpolate(t)).ToArray();
            double[] expInterpolatedIntensities = interpolatedTimes.Select(t => expSpline.Interpolate(t)).ToArray();

            // subtract the mean for both intensity arrays, leaving the zero padded data points intact
            var refIntensityMean = refInterpolatedIntensities.Average();
            var expIntensityMean = expInterpolatedIntensities.Average();
            for (int i = 0; i < interpArrayLength / 2; i++)
            {
                refInterpolatedIntensities[i + interpArrayLength / 4] -= refIntensityMean;
                expInterpolatedIntensities[i + interpArrayLength / 4] -= expIntensityMean;
            }

            // Create arrays of complex numbers from the intensity arrays
            Complex[] refArray = new Complex[refInterpolatedIntensities.Length];
            Complex[] expArray = new Complex[expInterpolatedIntensities.Length];
            for (int i = 0; i < refInterpolatedIntensities.Length; i++)
            {
                refArray[i] = new Complex(refInterpolatedIntensities[i], 0);
                expArray[i] = new Complex(expInterpolatedIntensities[i], 0);
            }

            // The fourier transform of the refArray is multiplied (element-wise) by 
            // the complex conjugate of the fourier transform of the exp array.
            // Then, an inverse Fourier transform is performed on the resulting product array
            Fourier.Forward(refArray);
            Fourier.Forward(expArray);
            Complex[] product = new Complex[refArray.Length];
            for (int i = 0; i < refArray.Length; i++)
            {
                product[i] = Complex.Multiply(refArray[i], expArray[i].Conjugate());
            }
            Fourier.Inverse(product); 

            // The multiplication above results in a rearrangement.
            // This shifts the inverse-product array such that the point corresponding to
            // a shift of zero minutes is in the middle, negative shifts are on the left, and
            // positive shifts are on the right
            double[] shiftedReals = new double[product.Length];
            for (int i = 0; i < shiftedReals.Length / 2; i++)
            {
                shiftedReals[i] = product[i + shiftedReals.Length / 2].Real;
                shiftedReals[i + shiftedReals.Length / 2] = product[i].Real;
            }
            
            var shifts = Generate.LinearSpaced(interpArrayLength, -0.5 * timeSpan, 0.5 * timeSpan);
            return shifts[shiftedReals.IndexOf(shiftedReals.Max())];
        }

        /// <summary>
        /// Returns a new List of indexed mass spectral peaks with added
        /// zero intensity peaks at the beginning and end of
        /// </summary>
        /// <param name="peaks"> Peaks to be padded</param>
        internal static List<IndexedMassSpectralPeak> PadPeaks(List<IndexedMassSpectralPeak> peaks)
        {
            List<IndexedMassSpectralPeak> paddedXICs = new List<IndexedMassSpectralPeak>();
            double mz = peaks.First().Mz;
            int arrayQuarterLength = peaks.Count / 4;
            double startTime = peaks.First().RetentionTime;
            double endTime = peaks.Last().RetentionTime;

            // Calculate the average RT difference between successive scans
            double summedScanSeparations = 0;
            for (int i = 0; i < peaks.Count - 1; i++)
            {
                summedScanSeparations += peaks[i + 1].RetentionTime - peaks[i].RetentionTime;
            }
            double averageScanSeparation = summedScanSeparations / (peaks.Count - 1);

            for (int i = 0; i < arrayQuarterLength; i++)
            {
                paddedXICs.Add(new IndexedMassSpectralPeak(mz: mz, intensity: 0, zeroBasedMs1ScanIndex: 0, 
                    retentionTime: startTime - averageScanSeparation * (arrayQuarterLength - i)));
            }
            paddedXICs.AddRange(peaks);
            for (int i = 0; i < arrayQuarterLength; i++)
            {
                paddedXICs.Add(new IndexedMassSpectralPeak(mz: mz, intensity: 0, zeroBasedMs1ScanIndex: 0,
                    retentionTime: endTime + averageScanSeparation * (i)));
            }

            return paddedXICs;
        }

        // Visual inspection suggests that the true peaks tend to have a min and max relatively close to one another. 
        // A large time delta between minima and maxima suggests a broad peak, or two otherwise poorly behave peaks.
        internal static Extremum[,] ReconcileExtrema(List<Extremum> refExtrema, List<List<Extremum>> expExtremaList)
        {
            // When finding Extrema on a cubic spline, there is no guarantee that it will alternate min + max indices 
            // (e.g., max, min, max, min). So, we need to record the location of each min and max w/in the ref array
            Extremum[] refMaxima = refExtrema.Where(e => e.ExtremumType == ExtremumType.Maximum).ToArray();
            int[] refMaxIndexMap = Enumerable.Range(0, refExtrema.Count)
                .Where(i => refExtrema[i].ExtremumType == ExtremumType.Maximum)
                .ToArray();
            Extremum[] refMinima = refExtrema.Where(e => e.ExtremumType == ExtremumType.Minimum).ToArray();
            int[] refMinIndexMap = Enumerable.Range(0, refExtrema.Count)
                .Where(i => refExtrema[i].ExtremumType == ExtremumType.Minimum)
                .ToArray();

            // 2D array where rows are the extrema from different runs (reference extrema, then exp)
            // and each column maps to a reference extrema
            Extremum[,] extrema2dArray = new Extremum[1 + expExtremaList.Count, refExtrema.Count];
            
            for (int i = 0; i < expExtremaList.Count; i++)
            {
                Extremum[] expMaxima = expExtremaList[i].Where(e => e != null && e.ExtremumType == ExtremumType.Maximum).ToArray();
                if (expMaxima.IsNotNullOrEmpty())
                {
                    var matchedExpMaxima = MatchExtrema(refMaxima, expMaxima);
                    for (int j = 0; j < refMaxIndexMap.Length; j++)
                    {
                        extrema2dArray[i + 1, refMaxIndexMap[j]] = matchedExpMaxima[j];
                    }
                }
                

                Extremum[] expMinima = expExtremaList[i].Where(e => e != null && e.ExtremumType == ExtremumType.Minimum).ToArray();
                if (expMinima.IsNotNullOrEmpty())
                {
                    var matchedExpMinima = MatchExtrema(refMinima, expMinima);
                    for (int j = 0; j < refMinIndexMap.Length; j++)
                    {
                        extrema2dArray[i + 1, refMinIndexMap[j]] = matchedExpMinima[j];
                    }
                }
            }

            // Fill in the first row with all the refExtrema. Some will be removed in the ResolveExtremaArray method
            for (int col = 0; col < extrema2dArray.GetLength(1); col++) 
                extrema2dArray[0, col] = refExtrema[col]; 
            
            ResolveExtremaArray(ref extrema2dArray);

            return extrema2dArray;
        }

        internal static Extremum[] MatchExtrema(Extremum[] refExtrema, Extremum[] expExtrema)
        {
            if (!refExtrema.IsNotNullOrEmpty() || !expExtrema.IsNotNullOrEmpty())
                return null;
            // For every reference Extremum, we search for potential pairs in the exp Extrema. If a 
            // partner is found, it is placed into the paired Exp Extrema array at the same index as 
            // that of the reference Extremum in the refExtrema array
            Extremum[] pairedExpExtrema = new Extremum[refExtrema.Length];

            // For each position in the reference away, stores the index of the closest Extremum in the experimental array
            int[] expIndices = new int[refExtrema.Length];
            for(int j = 0; j < refExtrema.Length; j++)
            {
                expIndices[j] = expExtrema.GetClosestIndex(refExtrema[j]);
            }

            int i = 0;
            while(i < expIndices.Length)
            {
                // It's possible that multiple reference extrema will all be closest
                // to one experimental extremum. In that case, we should report the RT
                // pair where the experimental and reference extrema are closest in time
                if (i < expIndices.Length - 1 && expIndices[i] == expIndices[i + 1])
                {
                    int duplicateIndex = expIndices[i];
                    double diff = Math.Abs(expExtrema[expIndices[i]] - refExtrema[i]);
                    Dictionary<int, double> indexDiffDict = new Dictionary<int, double> { { i, diff } };
                    i++;

                    while (i < expIndices.Length && expIndices[i] == duplicateIndex)
                    {
                        diff = Math.Abs(expExtrema[expIndices[i]] - refExtrema[i]);
                        indexDiffDict.Add(i, diff);
                        i++;
                    }

                    int indexOfClosestPair = indexDiffDict.MinBy(kvp => kvp.Value).Key;
                    pairedExpExtrema[indexOfClosestPair] = expExtrema[expIndices[indexOfClosestPair]];
                }
                else
                {
                    pairedExpExtrema[i] = expExtrema[expIndices[i]];
                    i++;
                }
            }

            return pairedExpExtrema;
        }

        /// <summary>
        /// Creates consensus within the array containing extrema that were paired with the reference extrema array.
        /// Reference array is the first row of the array.
        /// Method iterates through each non-reference column - if half or more ((num rows - 1)/2) have values,
        /// then that reference extremum is retained, all paired non-reference extrema are retained,
        /// and any rows (here, rows are unique XICs from different runs) that don't contain a paired extrema,
        /// an imputed extrema is inserted with retention time = average retention time for the column,
        /// type = reference type, and intensity = -1. If fewer than half have values, all extrema in that
        /// column are discarded, including the reference extrema
        /// </summary>
        public static void ResolveExtremaArray(ref Extremum[,] array)
        {
            for (int col = 0; col < array.GetLength(1); col++)
            {
                double meanTime = array[0, col].RetentionTime;
                int nullCount = 0;
                for(int row = 1; row < array.GetLength(0); row++)
                {
                    if (array[row, col] == null)
                        nullCount++;
                    else
                        meanTime += array[row, col].RetentionTime;
                }

                if (nullCount > array.GetLength(0) / 2)
                    SetColumnToNull(ref array, col);
                else if (nullCount > 0)
                {
                    meanTime /= (double)(array.GetLength(0) - nullCount);
                    Extremum imputedExtremum = new Extremum(meanTime, -1, array[0, col].ExtremumType);
                    for (int row = 1; row < array.GetLength(0); row++)
                    {
                        if (array[row, col] == null)
                            array[row, col] = imputedExtremum;
                    }
                }
            }
        }

        public static void SetColumnToNull(ref Extremum[,] array, int col)
        {
            for(int row = 0; row < array.GetLength(0); row++)
                array[row, col] = null;
        }
        //public static Extremum[] OrderExtrema(List<Extremum> maxima, List<Extremum> minima)
        //{
        //    Extremum[] orderedExtrema = new Extremum[maxima.Count + minima.Count];
        //    int maxPointer = 0;
        //    int minPointer = 0;
        //    while(maxPointer < maxima.Count && minPointer < minima.Count)
        //    {
        //        if (minima[minPointer] < maxima[maxPointer])
        //        {
        //            orderedExtrema[minPointer + maxPointer] = minima[minPointer++];
        //        }
        //        else
        //        {
        //            orderedExtrema[minPointer + maxPointer] = maxima[maxPointer++];
        //        }
        //    }

        //    // Only one of these loops will be hit. After the first while loop, only one list will
        //    // have remaining entries.
        //    while (maxPointer < maxima.Count)
        //    {
        //        orderedExtrema[minPointer + maxPointer] = maxima[maxPointer++];
        //    }
        //    while (minPointer < minima.Count)
        //    {
        //        orderedExtrema[minPointer + maxPointer] = minima[minPointer++];
        //    }

        //    return orderedExtrema;
        //}


        /// <summary>
        /// Modifies lists such that they have an equal length. The longer of the two lists will
        /// have the ends trimmed, yielding a segment from the middle. In cases where an odd number
        /// of entries are removed, one more entry is removed from the end.
        /// </summary>
        /// <exception cref="ArgumentException"> If passed a null or empty list, an exception is thrown</exception>
        public static void EqualizeListLength<T>(ref List<T> list1, ref List<T> list2)
        {
            if (!list1.IsNotNullOrEmpty() || !list2.IsNotNullOrEmpty())
            {
                throw new ArgumentException("Lists can not be null or empty");
            }

            if (list1.Count > list2.Count)
            {
                int trimLength = list1.Count - list2.Count;
                list1 = TrimList(list1, trimLength);
            }
            else if (list2.Count > list1.Count)
            {
                int trimLength = list2.Count - list1.Count;
                list2 = TrimList(list2, trimLength);
            }
        }

        private static List<T> TrimList<T>(List<T> list, int trimLength)
        {
            int startIndex = trimLength / 2;
            int endIndex = list.Count() -  (trimLength - startIndex);
            List<T> trimmedList = new();

            for (int i = startIndex; i <= endIndex; i++)
            {
                trimmedList.Add(list[i]);
            }

            return trimmedList;
        } 
    }

}
