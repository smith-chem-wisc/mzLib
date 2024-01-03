using System;
using System.Collections.Generic;
using System.Linq;
using MzLibUtil;

namespace MassSpectrometry
{
    public class ExtensionMethods
    {
        public static void MergeMzAndIntensityValuesForAnyMzsWithinTolerance(ref double[] mzArray, ref double[] intensityArray, PpmTolerance ppmTolerance)
        {
            List<double> newMzArray = new List<double>();
            List<double> newIntensityArray = new List<double>();

            List<int> indiciesToMerge = new List<int>();
            for (int i = 0; i < mzArray.Length; i++)
            {
                if (i < (mzArray.Length - 1) && ppmTolerance.Within(mzArray[i], mzArray[i + 1])) // two mzValues are close enough to merge
                {
                    indiciesToMerge.Add(i);
                    indiciesToMerge.Add(i + 1);
                }
                else // mzValues are NOT close enough to merge so we merge any that we've collected so far.
                {
                    if (indiciesToMerge.Any())
                    {
                        newMzArray.Add(GetWeightedMz(indiciesToMerge.DistinctBy(i => i).ToList(), mzArray.ToList(), intensityArray.ToList()));
                        newIntensityArray.Add(GetIntensityAverage(indiciesToMerge, mzArray, intensityArray));
                        indiciesToMerge.Clear(); //we clear this array because we have merged a contiguous group and can begin again with the next group
                    }
                    else
                    {
                        newMzArray.Add(mzArray[i]);
                        newIntensityArray.Add(intensityArray[i]);
                    }
                }
            }

            mzArray = newMzArray.ToArray();
            intensityArray = newIntensityArray.ToArray();
        }
        private static double GetIntensityAverage(List<int> indiciesToMerge, double[] mzArray, double[] intensityArray)
        {
            double intensitySum = 0;
            foreach (var index in indiciesToMerge)
            {
                intensitySum += intensityArray[index];
            }
            return (intensitySum / (double)indiciesToMerge.Count);
        }

        private static double GetWeightedMz(List<int> indiciesToMerge, List<double> mzArray, List<double> intensityArray)
        {
            double weightedMz = 0;
            double intensitySum = 0;
            foreach (var index in indiciesToMerge)
            {
                weightedMz += mzArray[index] * intensityArray[index];
                intensitySum += intensityArray[index];
            }
            return (weightedMz / intensitySum);
        }

        public static double[] MergeTwoSortedArraysIntoSortedArray(double[] sortedArrayOne, double[] sortedArrayTwo)
        {
            double[] mergedSortedArray = new double[sortedArrayOne.Length + sortedArrayTwo.Length];
            int oneIndex = 0;
            int twoIndex = 0;

            for (int i = 0; i < mergedSortedArray.Length; i++)
            {
                if (oneIndex<sortedArrayOne.Length && sortedArrayOne[oneIndex] < sortedArrayTwo[twoIndex])
                {
                    mergedSortedArray[i] = sortedArrayOne[oneIndex];
                    oneIndex++;
                }
                else
                {
                    if (twoIndex < sortedArrayTwo.Length)
                    {
                        mergedSortedArray[i] = sortedArrayTwo[twoIndex];
                        twoIndex++;
                    }
                }
            }
            return mergedSortedArray;
        }

        public (T[], U[]) MergeSortedArrays<T, U>(T[][] arraysToBeSorted, U[][] valuesToBePreserved) where T : IComparable
        {

        }
    }
}
