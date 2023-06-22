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

        // Probably best to pass in more XIC than you need, not less
        public static double AlignXICs(List<IndexedMassSpectralPeak> refXIC, List<IndexedMassSpectralPeak> expXIC, int resolution = 100)
        {
            
            // Creates a spline, or a smoothed representation of the XIC which allows for up-sampling of the data
            // Unsure about whether to use cubic or linear. Can maybe be a delegate argument?
            var referenceSpline = LinearSpline.InterpolateSorted(
                refXIC.Select(p => p.RetentionTime).ToArray(), refXIC.Select(p => p.Intensity).ToArray());
            var expSpline = LinearSpline.InterpolateSorted(
                expXIC.Select(p => p.RetentionTime).ToArray(), expXIC.Select(p => p.Intensity).ToArray());

            var startTime = Math.Max(refXIC.First().RetentionTime, expXIC.First().RetentionTime);
            var endTime = Math.Min(refXIC.Last().RetentionTime, expXIC.Last().RetentionTime);
            var minLength = Math.Min(expXIC.Count, refXIC.Count);
            var timeSpan = endTime - startTime;

            var interpArrayLength = (int)Math.Floor(timeSpan * resolution * minLength);
            double[] refInterpolatedTimes = Generate.LinearSpaced(length: interpArrayLength, startTime, endTime);
            double[] refInterpolatedIntensities =
                refInterpolatedTimes.Select(t => referenceSpline.Interpolate(t)).ToArray();


            EqualizeListLength(ref refXIC, ref expXIC);

            double[] expInterpolatedTimes = Generate.LinearSpaced(length: interpArrayLength, startTime, endTime);
            double[] expInterpolatedIntensities =
                expInterpolatedTimes.Select(t => expSpline.Interpolate(t)).ToArray();


            var refIntensityMean = refInterpolatedIntensities.Average();
            var expIntensityMean = expInterpolatedIntensities.Average();
            for (int i = 0; i < refInterpolatedIntensities.Length; i++)
            {
                refInterpolatedIntensities[i] -= refIntensityMean;
                expInterpolatedIntensities[i] -= expIntensityMean;
            }

            double[][] intensityMatrix = { refInterpolatedIntensities, expInterpolatedIntensities };
            //var test = Correlation.PearsonMatrix(intensityMatrix);

            //double[] complexComponent = new double[refInterpolatedIntensities.Length];
            Complex[] refArray = new Complex[refInterpolatedIntensities.Length];
            Complex[] expArray = new Complex[expInterpolatedIntensities.Length];
            for (int i = 0; i < refInterpolatedIntensities.Length; i++)
            {
                refArray[i] = new Complex(refInterpolatedIntensities[i], 0);
                expArray[i] = new Complex(expInterpolatedIntensities[i], 0);
            }

            Fourier.Forward(refArray);
            Fourier.Forward(expArray);
            Complex[] product = new Complex[refArray.Length];
            for (int i = 0; i < refArray.Length; i++)
            {
                product[i] = Complex.Multiply(refArray[i], expArray[i].Conjugate());
            }

            Fourier.Inverse(product); // Looks like the zero is already in the middle?

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
