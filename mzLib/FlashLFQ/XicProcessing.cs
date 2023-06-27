using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using Easy.Common.Extensions;
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
        /// The XICs are up-sampled to allow for sub-pixel resolution (one XIC datapoint = one pixel).
        /// </summary>
        /// <param name="refXIC">List of peaks to be used as a reference, ordered by retention time</param>
        /// <param name="expXIC">List of peaks to be aligned, ordered by retention time</param>
        /// <param name="resolution">Up-sampling resolution. Higher values allow for more precise shifts</param>
        /// <returns>The retention time correction to align the experimental and reference XICs</returns>
        /// (Suggested reading: Convolution theorem, https://dsp.stackexchange.com/questions/51409/maximum-of-cross-correlation-not-moving)
        public static double AlignXICs(List<IndexedMassSpectralPeak> refXIC, List<IndexedMassSpectralPeak> expXIC, int resolution = 100)
        {
            // First step is to pad the XICs, adding zeros to either end
            // This could be further optimized by padding at the array level w/o creating new peak objects
            refXIC = PadXICs(refXIC);
            expXIC = PadXICs(expXIC);

            // Creates a spline, or a smoothed representation of the XIC which allows for up-sampling of the data
            // Unsure about whether to use cubic or linear. Can maybe be a delegate argument?
            var referenceSpline = LinearSpline.InterpolateSorted(
                refXIC.Select(p => p.RetentionTime).ToArray(), refXIC.Select(p => p.Intensity).ToArray());
            var expSpline = LinearSpline.InterpolateSorted(
                expXIC.Select(p => p.RetentionTime).ToArray(), expXIC.Select(p => p.Intensity).ToArray());

            var startTime = Math.Max(refXIC.First().RetentionTime, expXIC.First().RetentionTime);
            var endTime = Math.Min(refXIC.Last().RetentionTime, expXIC.Last().RetentionTime);
            var timeSpan = endTime - startTime;
            int minArrayLength = Math.Min(expXIC.Count, refXIC.Count);
            
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
        /// <param name="XIC">XIC to be padded</param>
        private static List<IndexedMassSpectralPeak> PadXICs(List<IndexedMassSpectralPeak> XIC)
        {
            List<IndexedMassSpectralPeak> paddedXICs = new List<IndexedMassSpectralPeak>();
            double mz = XIC.First().Mz;
            int arrayQuarterLength = XIC.Count / 4;
            double startTime = XIC.First().RetentionTime;
            double endTime = XIC.Last().RetentionTime;

            // Calculate the average RT difference between successive scans
            double summedScanSeparations = 0;
            for (int i = 0; i < XIC.Count - 1; i++)
            {
                summedScanSeparations += XIC[i + 1].RetentionTime - XIC[i].RetentionTime;
            }
            double averageScanSeparation = summedScanSeparations / (XIC.Count - 1);

            for (int i = 0; i < arrayQuarterLength; i++)
            {
                paddedXICs.Add(new IndexedMassSpectralPeak(mz: mz, intensity: 0, zeroBasedMs1ScanIndex: 0, 
                    retentionTime: startTime - averageScanSeparation * (arrayQuarterLength - i)));
            }
            paddedXICs.AddRange(XIC);
            for (int i = 0; i < arrayQuarterLength; i++)
            {
                paddedXICs.Add(new IndexedMassSpectralPeak(mz: mz, intensity: 0, zeroBasedMs1ScanIndex: 0,
                    retentionTime: endTime + averageScanSeparation * (i)));
            }

            return paddedXICs;
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
