using Easy.Common.Extensions;
using MassSpectrometry;
using MzLibUtil;

namespace Readers.Bruker
{
    public static class TofSpectraMerger
    {
        public static readonly double DefaultPpmTolerance = 10;

        internal static MzSpectrum CreateFilteredSpectrum(IList<double> mzs, IList<int> intensities, FilteringParams filteringParams = null)
        {
            // Convert the indices to m/z values
            double[] mzsArray = mzs.ToArray();

            // Convert the intensities to an array
            double[] intensitiesArray = intensities.Select(intensity => (double)intensity).ToArray();

            if (mzsArray.Length != intensitiesArray.Length)
                throw new Exception("Collapsed m/z and intensity arrays are not the same length.");

            if (filteringParams != null
                && mzsArray.Length > 0
                && filteringParams.ApplyTrimmingToMs1)
            {
                WindowModeHelper.Run(ref intensitiesArray,
                    ref mzsArray, filteringParams,
                    mzsArray[0], mzsArray[^0]);
            }
            // TODO: This would be more performant if we kept the intensities as ints
            return new MzSpectrum(mzsArray, intensitiesArray, shouldCopy: false);
        }

        internal static (uint[] Indices, int[] Intensities) MergeArraysToArray(List<uint[]> indexArrays, List<int[]> intensityArrays)
        {
            if (!indexArrays.IsNotNullOrEmpty() || intensityArrays == null || intensityArrays.Count() != indexArrays.Count())
                return (new uint[0], new int[0]);

            // Merge all index arrays and intensity arrays into a single array
            uint[] combinedIndices = indexArrays[0];
            int[] combinedIntensities = intensityArrays[0];
            for (int i = 1; i < indexArrays.Count(); i++)
            {
                var mergeResults = TwoPointerMerge(combinedIndices, indexArrays[i], combinedIntensities, intensityArrays[i]);
                combinedIndices = mergeResults.Indices;
                combinedIntensities = mergeResults.Intensities;
            }

            // Collapse the combined arrays into a single array (centroiding, more or less)
            return CollapseArrays(combinedIndices, combinedIntensities);
        }

        public static (uint[] Indices, int[] Intensities) TwoPointerMerge(uint[] indexArray1, uint[] indexArray2, int[] intensityArray1, int[] intensityArray2)
        {
            int p1 = 0;
            int p2 = 0;

            uint[] mergedIndices = new uint[indexArray1.Length + indexArray2.Length];
            int[] mergedIntensities = new int[intensityArray1.Length + intensityArray2.Length];

            while (p1 < indexArray1.Length || p2 < indexArray2.Length)
            {
                if (p1 == indexArray1.Length)
                {
                    while (p2 < indexArray2.Length)
                    {
                        mergedIndices[p1 + p2] = indexArray2[p2];
                        mergedIntensities[p1 + p2] = intensityArray2[p2];
                        p2++;
                    }
                }
                else if (p2 == indexArray2.Length)
                {
                    while (p1 < indexArray1.Length)
                    {
                        mergedIndices[p1 + p2] = indexArray1[p1];
                        mergedIntensities[p1 + p2] = intensityArray1[p1];
                        p1++;
                    }
                }
                else if (indexArray1[p1] < indexArray2[p2])
                {
                    mergedIndices[p1 + p2] = indexArray1[p1];
                    mergedIntensities[p1 + p2] = intensityArray1[p1];
                    p1++;
                }
                else
                {
                    mergedIndices[p1 + p2] = indexArray2[p2];
                    mergedIntensities[p1 + p2] = intensityArray2[p2];
                    p2++;
                }
            }

            return (mergedIndices, mergedIntensities);
        }

        public static (uint[] Indices, int[] Intensities) CollapseArrays(uint[] indexArray, int[] intensityArray)
        {
            // This is a quick and dirty implementation of the collapse function. There are some edge cases that are not handled here.

            // Define lists to store the collapsed indices and intensities
            List<uint> collapsedIndices = new();
            List<int> collapsedIntensities = new();

            // Initialize pointers to the first two elements in the index array
            int p1 = 0;
            int p2 = 1;
            while (p1 < indexArray.Length)
            {
                uint currentIdx = indexArray[p1];

                // Find clusters of indices that are close together
                // increment pointer 2 until the cluster ends and we're further than 3 indices away
                while (p2 < indexArray.Length && (2 + currentIdx) >= indexArray[p2])
                {
                    p2++;
                }
                p2--; // Move the pointer back by one
                int medianPointer = (p1 + p2) / 2;
                // Use the median index in each cluster as the collapsed index
                collapsedIndices.Add(indexArray[medianPointer]);

                // Sum the intensities in each cluster to get the collapsed intensity
                int summedIntensity = 0;
                for (int i = p1; i <= p2; i++)
                {
                    summedIntensity += intensityArray[i];
                }
                collapsedIntensities.Add(summedIntensity);

                // Move the pointers forward
                p1 = p2 + 1;
                p2 = p1 + 1;
            }

            return (collapsedIndices.ToArray(), collapsedIntensities.ToArray());
        }

        #region MzLevelOperations

        internal static MzSpectrum MergeArraysToSpectrum(
            List<uint[]> indexArrays, 
            List<int[]> intensityArrays, 
            FrameProxyFactory proxyFactory,
            FilteringParams filteringParams = null)
        {
            if (!indexArrays.IsNotNullOrEmpty() || intensityArrays == null || intensityArrays.Count() != indexArrays.Count())
                return null;

            // Merge all index arrays and intensity arrays into a single array
            uint[] combinedIndices = indexArrays[0];
            int[] combinedIntensities = intensityArrays[0];
            for (int i = 1; i < indexArrays.Count(); i++)
            {
                var mergeResults = TwoPointerMerge(combinedIndices, indexArrays[i], combinedIntensities, intensityArrays[i]);
                combinedIndices = mergeResults.Indices;
                combinedIntensities = mergeResults.Intensities;
            }

            // Collapse the combined arrays into a single array (centroiding, more or less)
            var centroidedResults = CollapseArrays(proxyFactory.ConvertIndicesToMz(combinedIndices), combinedIntensities);

            return CreateFilteredSpectrum(
                centroidedResults.Mzs,
                centroidedResults.Intensities,
                filteringParams);
        }

        internal static (double[] Mzs, int[] Intensities) CollapseArrays(double[] mzArray, int[] intensityArray, double ppmTolerance = -1)
        {
            // Define lists to store the collapsed indices and intensities
            List<double> collapsedMzs = new();
            List<int> collapsedIntensities = new();

            PpmTolerance tol = new(ppmTolerance < 1 ? DefaultPpmTolerance : ppmTolerance);

            // Initialize pointers to the first two elements in the index array
            int p1 = 0;
            int p2 = 1;
            while (p1 < mzArray.Length)
            {
                double currentMz = mzArray[p1];
                double upperBoundMz = tol.GetMaximumValue(currentMz);

                // Find clusters of indices that are close together
                // increment pointer 2 until the cluster ends and we're further than 3 indices away
                while (p2 < mzArray.Length && upperBoundMz >= mzArray[p2])
                {
                    upperBoundMz = tol.GetMaximumValue(mzArray[p2]);
                    p2++;
                }
                p2--; // Move the pointer back by one

                if(p1 == p2)
                {
                    collapsedIntensities.Add(intensityArray[p1]);
                    collapsedMzs.Add(mzArray[p1]);
                }
                else
                {
                    // Calculate the summed intensity in the cluster
                    int summedIntensity = 0;
                    for (int i = p1; i <= p2; i++)
                    {
                        summedIntensity += intensityArray[i];
                    }
                    collapsedIntensities.Add(summedIntensity);

                    // weighted averaging to determine the collapsed m/z of the cluster
                    double collapsedMz = 0;
                    for (int i = p1; i <= p2; i++)
                    {
                        double weight = (double)intensityArray[i] / (double)summedIntensity;
                        collapsedMz += weight * mzArray[i];
                    }
                    collapsedMzs.Add(collapsedMz);
                }

                // Move the pointers forward
                p1 = p2 + 1;
                p2 = p1 + 1;
            }

            return (collapsedMzs.ToArray(), collapsedIntensities.ToArray());
        }

        internal static (double[] Mzs, int[] Intensities) MergeArraysToMzArray(List<uint[]> indexArrays, List<int[]> intensityArrays, FrameProxyFactory proxyFactory)
        {
            if (!indexArrays.IsNotNullOrEmpty() || intensityArrays == null || intensityArrays.Count() != indexArrays.Count())
                return (new double[0], new int[0]);

            // Merge all index arrays and intensity arrays into a single array
            uint[] combinedIndices = indexArrays[0];
            int[] combinedIntensities = intensityArrays[0];
            for (int i = 1; i < indexArrays.Count(); i++)
            {
                var mergeResults = TwoPointerMerge(combinedIndices, indexArrays[i], combinedIntensities, intensityArrays[i]);
                combinedIndices = mergeResults.Indices;
                combinedIntensities = mergeResults.Intensities;
            }
            double[] mzsArray = proxyFactory.ConvertIndicesToMz(combinedIndices);

            // Collapse the combined arrays into a single array (centroiding, more or less)
            return CollapseArrays(mzsArray, combinedIntensities);
        }

        internal static MzSpectrum MergeArraysToSpectrum(
            List<double[]> mzArrays, 
            List<int[]> intensityArrays, 
            FilteringParams filteringParams = null, 
            double ppmTolerance = -1)
        {
            if (!mzArrays.IsNotNullOrEmpty() || intensityArrays == null || intensityArrays.Count() != mzArrays.Count())
                return null;

            // Merge all index arrays and intensity arrays into a single array
            double[] combinedMzs = mzArrays[0];
            int[] combinedIntensities = intensityArrays[0];
            for (int i = 1; i < mzArrays.Count(); i++)
            {
                var mergeResults = TwoPointerMerge(combinedMzs, mzArrays[i], combinedIntensities, intensityArrays[i]);
                combinedMzs = mergeResults.Mzs;
                combinedIntensities = mergeResults.Intensities;
            }

            // Collapse the combined arrays into a single array (centroiding, more or less)
            var centroidedResults = CollapseArrays(combinedMzs, combinedIntensities, ppmTolerance);

            return CreateFilteredSpectrum(centroidedResults.Mzs, centroidedResults.Intensities, filteringParams);
        }

        public static (double[] Mzs, int[] Intensities) TwoPointerMerge(double[] mzArray1, double[] mzArray2, int[] intensityArray1, int[] intensityArray2)
        {
            int p1 = 0;
            int p2 = 0;

            double[] mergedMzs = new double[mzArray1.Length + mzArray2.Length];
            int[] mergedIntensities = new int[intensityArray1.Length + intensityArray2.Length];

            while (p1 < mzArray1.Length || p2 < mzArray2.Length)
            {
                if (p1 == mzArray1.Length)
                {
                    while (p2 < mzArray2.Length)
                    {
                        mergedMzs[p1 + p2] = mzArray2[p2];
                        mergedIntensities[p1 + p2] = intensityArray2[p2];
                        p2++;
                    }
                }
                else if (p2 == mzArray2.Length)
                {
                    while (p1 < mzArray1.Length)
                    {
                        mergedMzs[p1 + p2] = mzArray1[p1];
                        mergedIntensities[p1 + p2] = intensityArray1[p1];
                        p1++;
                    }
                }
                else if (mzArray1[p1] < mzArray2[p2])
                {
                    mergedMzs[p1 + p2] = mzArray1[p1];
                    mergedIntensities[p1 + p2] = intensityArray1[p1];
                    p1++;
                }
                else
                {
                    mergedMzs[p1 + p2] = mzArray2[p2];
                    mergedIntensities[p1 + p2] = intensityArray2[p2];
                    p2++;
                }
            }

            return (mergedMzs, mergedIntensities);
        }


        #endregion

    }
}