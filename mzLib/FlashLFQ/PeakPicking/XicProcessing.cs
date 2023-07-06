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
        /// <param name="Peaks">Peaks to be padded</param>
        private static List<IndexedMassSpectralPeak> PadPeaks(List<IndexedMassSpectralPeak> Peaks)
        {
            List<IndexedMassSpectralPeak> paddedXICs = new List<IndexedMassSpectralPeak>();
            double mz = Peaks.First().Mz;
            int arrayQuarterLength = Peaks.Count / 4;
            double startTime = Peaks.First().RetentionTime;
            double endTime = Peaks.Last().RetentionTime;

            // Calculate the average RT difference between successive scans
            double summedScanSeparations = 0;
            for (int i = 0; i < Peaks.Count - 1; i++)
            {
                summedScanSeparations += Peaks[i + 1].RetentionTime - Peaks[i].RetentionTime;
            }
            double averageScanSeparation = summedScanSeparations / (Peaks.Count - 1);

            for (int i = 0; i < arrayQuarterLength; i++)
            {
                paddedXICs.Add(new IndexedMassSpectralPeak(mz: mz, intensity: 0, zeroBasedMs1ScanIndex: 0, 
                    retentionTime: startTime - averageScanSeparation * (arrayQuarterLength - i)));
            }
            paddedXICs.AddRange(Peaks);
            for (int i = 0; i < arrayQuarterLength; i++)
            {
                paddedXICs.Add(new IndexedMassSpectralPeak(mz: mz, intensity: 0, zeroBasedMs1ScanIndex: 0,
                    retentionTime: endTime + averageScanSeparation * (i)));
            }

            return paddedXICs;
        }

        // Visual inspection suggests that the true peaks tend to have a min and max relatively close to one another. 
        // A large time delta between minima and maxima suggests a broad peak, or two otherwise poorly behave peaks.
        public static List<Extremum[]> ReconcileExtrema(List<Extremum> refExtrema, List<List<Extremum>> expExtremaList)
        {
            Extremum[] refMaxima = refExtrema.Where(e => e.Type == ExtremumType.Maximum).ToArray();
            Extremum[] refMinima = refExtrema.Where(e => e.Type == ExtremumType.Minimum).ToArray();
            // One array corresponding to the reference extrema paired against each of the extrema
            List<Extremum[]> pairedRefExtremaList = new();
            // One array containing the paired extrema from each experimental array
            List<Extremum[]> pairedExpExtremaList = new(); 
            foreach (List<Extremum> expExtrema in expExtremaList)
            {
                Extremum[] expMaxima = expExtrema.Where(e => e.Type == ExtremumType.Maximum).ToArray();
                Extremum[] expMinima = expExtrema.Where(e => e.Type == ExtremumType.Minimum).ToArray();

                var maxPairs = GetRetentionTimePairs(refMaxima, expMaxima);
                var minPairs = GetRetentionTimePairs(refMinima, expMinima);

                pairedRefExtremaList.Add(maxPairs
                        .Select(t => t.reference)
                        .Union(minPairs.Select(t => t.reference))
                        .OrderBy(t => t.RetentionTime)
                        .ToArray());
                pairedExpExtremaList.Add(maxPairs
                    .Select(t => t.exp)
                    .Union(minPairs.Select(t => t.exp))
                    .OrderBy(t => t.RetentionTime)
                    .ToArray());
            }
            // TODO: Create some consensus list from the different sets of reference extrema

            pairedExpExtremaList.Insert(0,pairedRefExtremaList[0]);
            return pairedExpExtremaList;
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

        public static List<(Extremum reference, Extremum exp)> GetRetentionTimePairs(Extremum[] refExtrema, Extremum[] expExtrema)
        {
            // The reference should be shorter than the experimental (requires fewer iterations)
            // if that isn't the case, we swap the ref and exp arrays.
            bool swapped = false;
            if (refExtrema.Length > expExtrema.Length)
            {
                Extremum[] temp = new Extremum[refExtrema.Length];
                Array.Copy(refExtrema, temp, refExtrema.Length);
                refExtrema = expExtrema;
                expExtrema = temp;
                swapped = true;
            }

            List<(Extremum, Extremum)> rtPairs = new List<(Extremum, Extremum)>();
            List<int> expIndices = new();
            foreach(Extremum refExtremum in refExtrema)
            {
                expIndices.Add(expExtrema.GetClosestIndex(refExtremum));
            }

            for(int i = 0; i < expIndices.Count; i++)
            {
                // It's possible that multiple reference extrema will all be closest
                // to one experimental extremum. In that case, we should report the RT
                // pair where the experimental and reference extrema are closest in time
                if (i < expIndices.Count - 1 && expIndices[i] == expIndices[i+1])
                {
                    int duplicateIndex = expIndices[i];
                    double diff = Math.Abs(expExtrema[expIndices[i]] - refExtrema[i]);
                    Dictionary<int, double> indexDiffDict = new Dictionary<int, double>{ {i, diff} };
                    i++;

                    while(i < expIndices.Count && expIndices[i] == duplicateIndex)
                    {
                        diff = Math.Abs(expExtrema[expIndices[i]] - refExtrema[i]);
                        indexDiffDict.Add(i, diff);
                        i++;
                    }

                    int indexOfClosestPair = indexDiffDict.MinBy(kvp => kvp.Value).Key;
                    rtPairs.Add((refExtrema[indexOfClosestPair], expExtrema[expIndices[indexOfClosestPair]]));
                }
                else
                {
                    rtPairs.Add((refExtrema[i], expExtrema[expIndices[i]]));
                }
            }

            return swapped // If the two were swapped up top, they're swapped back here
                ? rtPairs.Select(pair => (pair.Item2, pair.Item1)).ToList() 
                : rtPairs;
        }

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
