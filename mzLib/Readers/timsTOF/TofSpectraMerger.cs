using Easy.Common.Extensions;
using MassSpectrometry;
using MzLibUtil;

namespace Readers
{
    public static class TofSpectraMerger
    {
        public static readonly double DefaultPpmTolerance = 10;

        #region NewCentroidingApproach
        // The following methods are used to merge and collapse index arrays and intensity arrays
        // The timsTOF data format doesn't store m/z values directly, but rather indices in a lookup table where the mz values are stored 
        // Keeping these indices as ints allows for more efficient storage and processing of the data

        public const int NoiseFloor = 75;

        // The general workflow is to read in individual scans as index and intensity arrays
        // Then, merge these arrays into a single array sorted by tofIndex, ascending
        // Then, collapse these arrays by combining entries with the same tofIndex and removing entries that fall below a certain intensity threshold (NoiseFloor)
        // Finally, perform centroiding by grouping adjacent tofIndices, summing their intensities, and then calculating a weighted average of the m/z values in the cluster


        /// <summary>
        /// Merges two index and intensity arrays using a two-pointer technique.
        /// The merged arrays are sorted by index, ascending
        /// </summary>
        /// <param name="indexArray1">First index array.</param>
        /// <param name="indexArray2">Second index array.</param>
        /// <param name="intensityArray1">First intensity array.</param>
        /// <param name="intensityArray2">Second intensity array.</param>
        /// <returns>A tuple containing the merged indices and intensities.</returns>
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

        /// <summary>
        /// This combines the index and intensity arrays into a TimsSpectrum object
        /// If the same index is present multiple times, the intensities are summed
        /// </summary>
        /// <param name="indexArray"></param>
        /// <param name="intensityArray"></param>
        /// <param name="zeroIndexedTimsScanNumber"></param>
        /// <returns></returns>
        internal static (uint[] Indices, int[] Intensities) CollapseArrays(uint[] indexArray, int[] intensityArray, bool removeLowIntensityPeaks = true)
        {
            // Define lists to store the collapsed indices and intensities
            List<uint> collapsedIndices = new List<uint>(indexArray.Length);
            List<int> collapsedIntensities = new List<int>(intensityArray.Length);

            // Initialize pointers to the first two elements in the index array
            int p1 = 0;
            int p2 = 1;
            while (p1 < indexArray.Length)
            {
                uint currentIdx = indexArray[p1];

                // Combine intensities if array contains same idx multiple times
                while (p2 < indexArray.Length && currentIdx == indexArray[p2])
                {
                    p2++;
                }
                p2--; // Move the pointer back by one

                // Sum the intensities in each cluster to get the collapsed intensity
                int summedIntensity = 0;
                for (int i = p1; i <= p2; i++)
                {
                    summedIntensity += intensityArray[i];
                }

                if (!removeLowIntensityPeaks || summedIntensity > NoiseFloor)
                {
                    collapsedIndices.Add(currentIdx);
                    collapsedIntensities.Add(summedIntensity);
                }

                // Move the pointers forward
                p1 = p2 + 1;
                p2 = p1 + 1;
            }

            collapsedIndices.TrimExcess();
            collapsedIntensities.TrimExcess();

            return (collapsedIndices.ToArray(), collapsedIntensities.ToArray());
        }

        /// <summary>
        /// This centroids a TOF spectrum by grouping peaks with adjacent indices, summing their intensities
        /// and reporting the weighted average of the m/z values in the cluster
        /// </summary>
        /// <param name="indexArray"></param>
        /// <param name="intensityArray"></param>
        /// <param name="proxyFactory"></param>
        /// <returns></returns>
        internal static (double[] Mzs, int[] Intensities) Centroid(uint[] indexArray, int[] intensityArray, FrameProxyFactory proxyFactory)
        {
            // Define lists to store the collapsed indices and intensities
            List<double> collapsedMzs = new();
            List<int> collapsedIntensities = new();

            // Initialize pointers to the first two elements in the index array
            int p1 = 0;
            int p2 = 1;
            while (p1 < indexArray.Length)
            {
                // We could do this based on tolerance
                // Now, i'm testing what happens if we say grouped points must have adjacent tof indices
                uint upperBoundTofIndex = indexArray[p1] + 2;

                // Find clusters of indices that are close together
                // increment pointer 2 until the cluster ends and we're further than 3 indices away
                while (p2 < indexArray.Length && indexArray[p2] <= upperBoundTofIndex)
                {
                    upperBoundTofIndex = indexArray[p2] + 2;
                    p2++;
                }
                p2--; 

                if (p1 == p2)
                {
                    collapsedIntensities.Add(intensityArray[p1]);
                    collapsedMzs.Add(proxyFactory.ConvertIndexToMz(indexArray[p1]));
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
                        collapsedMz += weight * proxyFactory.ConvertIndexToMz(indexArray[i]);
                    }
                    collapsedMzs.Add(collapsedMz);
                }

                // Move the pointers forward
                p1 = p2 + 1;
                p2 = p1 + 1;
            }

            return (collapsedMzs.ToArray(), collapsedIntensities.ToArray());
        }

        /// <summary>
        /// Combines multiple scans into a single TimsSpectrum object without performing centroiding
        /// This is called when analyzing MS2 scans, where the same precursor is selected for fragmentation over multiple frames
        /// Each framge gets one TimsSpectrum
        /// </summary>
        /// <param name="indexArrays"></param>
        /// <param name="intensityArrays"></param>
        /// <param name="zeroIndexedTimsScanNumber"></param>
        /// <returns></returns>
        internal static TimsSpectrum CreateTimsSpectrum(List<uint[]> indexArrays, List<int[]> intensityArrays)
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

            var collapsedResults = CollapseArrays(combinedIndices, combinedIntensities, removeLowIntensityPeaks: false);
            return new TimsSpectrum(collapsedResults.Indices, collapsedResults.Intensities);
        }

        /// <summary>
        /// Merges multiple index and intensity arrays into an mzSpectrum.
        /// This operation is somewhere between averaging and centroiding
        /// In the TimsTofFileReader, MS1 scans are kept as index arrays and intensity arrays.
        /// </summary>
        /// <param name="indexArrays">List of index arrays.</param>
        /// <param name="intensityArrays">List of intensity arrays.</param>
        /// <param name="proxyFactory">Frame proxy factory.</param>
        /// <param name="filteringParams">Filtering parameters (optional).</param>
        /// <returns>A merged MS1 spectrum.</returns>
        internal static MzSpectrum CreateMzSpectrum(
            List<TimsSpectrum> timsSpectra,
            FrameProxyFactory proxyFactory,
            int msnLevel = 2,
            FilteringParams filteringParams = null)
        {
            List<uint[]> indices = new List<uint[]>(timsSpectra.Count);
            List<int[]> intensities = new List<int[]>(timsSpectra.Count);
            foreach(var timsSpectrum in timsSpectra)
            {
                indices.Add(timsSpectrum.XArray);
                intensities.Add(timsSpectrum.YArray);
            }

            return CreateMzSpectrum(indices, intensities, proxyFactory, msnLevel, filteringParams);
        }

        /// <summary>
        /// Merges multiple index and intensity arrays into an MzSpectrum
        /// This operation centroids and averages multiple component scans
        /// </summary>
        /// <param name="indexArrays">List of index arrays.</param>
        /// <param name="intensityArrays">List of intensity arrays.</param>
        /// <param name="proxyFactory">Frame proxy factory.</param>
        /// <param name="msnLevel">MSn level of the spectrum (1 or 2).</param>
        /// <param name="filteringParams">Filtering parameters (optional).</param>
        /// <returns>A merged MS1 spectrum.</returns>
        internal static MzSpectrum CreateMzSpectrum(
        List<uint[]> indexArrays,
        List<int[]> intensityArrays,
        FrameProxyFactory proxyFactory,
        int msnLevel,
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
            var collapsedResults = CollapseArrays(combinedIndices, combinedIntensities);
            var centroidedResults = Centroid(collapsedResults.Indices, collapsedResults.Intensities, proxyFactory);

            return CreateFilteredSpectrum(
                centroidedResults.Mzs,
                centroidedResults.Intensities,
                filteringParams,
                msnLevel: msnLevel);
        }

        
        internal static MzSpectrum CreateFilteredSpectrum(IList<double> mzs, IList<int> intensities,
            FilteringParams filteringParams = null, int msnLevel = 1)
        {
            double[] mzsArray;
            if (mzs is double[])
                mzsArray = (double[])mzs;
            else
                mzsArray = mzs.ToArray();

            // Convert the intensities to an array
            double[] intensitiesArray = intensities.Select(intensity => (double)intensity).ToArray();

            if (mzsArray.Length != intensitiesArray.Length)
                throw new Exception("Collapsed m/z and intensity arrays are not the same length.");

            if (filteringParams != null
                && mzsArray.Length > 0
                && ((filteringParams.ApplyTrimmingToMs1 && msnLevel == 1)
                    || (filteringParams.ApplyTrimmingToMsMs && msnLevel > 1)))
            {
                WindowModeHelper.Run(ref intensitiesArray,
                    ref mzsArray, filteringParams,
                    mzsArray[0], mzsArray[^1]);
            }
            // TODO: This would be more performant if we kept the intensities as ints
            return new MzSpectrum(mzsArray, intensitiesArray, shouldCopy: false);
        }

        #endregion

        /// <summary>
        /// Merges two m/z and intensity arrays using a two-pointer technique.
        /// Used when merging component spectra into one MS2 spectrum
        /// </summary>
        /// <param name="mzArray1">First m/z array.</param>
        /// <param name="mzArray2">Second m/z array.</param>
        /// <param name="intensityArray1">First intensity array.</param>
        /// <param name="intensityArray2">Second intensity array.</param>
        /// <returns>A tuple containing the merged m/z values and intensities.</returns>
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

        /// <summary>
        /// Collapses the given mz and intensity arrays. 
        /// mz values within ppmTolerance (and their corresponding intensity values) are merged. 
        /// The idea here is to centroid a spectrum
        /// </summary>
        /// <param name="mzArray">The mz array to collapse.</param>
        /// <param name="intensityArray">The intensity array to collapse.</param>
        /// /// <param name="ppmTolerance">PPM tolerance value (default is 10).</param>
        /// <returns>A tuple containing the collapsed mz and intensities.</returns>
        internal static (double[] Mzs, int[] Intensities) CollapseArrays(double[] mzArray, int[] intensityArray, double ppmTolerance = 10)
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

                if (p1 == p2)
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

        #region MzLevelOperations



        /// <summary>
        /// Merges multiple index and intensity arrays into an MS1 spectrum.
        /// This operation is somewhere between averaging and centroiding
        /// In the TimsTofFileReader, MS1 scans are kept as index arrays and intensity arrays.
        /// </summary>
        /// <param name="indexArrays">List of index arrays.</param>
        /// <param name="intensityArrays">List of intensity arrays.</param>
        /// <param name="proxyFactory">Frame proxy factory.</param>
        /// <param name="filteringParams">Filtering parameters (optional).</param>
        /// <returns>A merged MS1 spectrum.</returns>
        internal static MzSpectrum MergeArraysToMs1Spectrum(
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
                filteringParams,
                msnLevel: 1);
        }

        /// <summary>
        /// Merges multiple m/z and intensity arrays into an MS2 spectrum.
        /// This operation is somewhere between averaging and centroiding.
        /// In the TimsTofFileReader, MS2 component spectrum are stored as 
        /// double[] m/z arrays and int[] intensity arrays.
        /// Each component scan 
        /// </summary>
        /// <param name="mzArrays">List of m/z arrays.</param>
        /// <param name="intensityArrays">List of intensity arrays.</param>
        /// <param name="filteringParams">Filtering parameters (optional).</param>
        /// <param name="ppmTolerance">PPM tolerance value (default is -1).</param>
        /// <returns>A merged MS2 spectrum.</returns>
        internal static MzSpectrum MergeArraysToMs2Spectrum(
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

            return CreateFilteredSpectrum(
                centroidedResults.Mzs,
                centroidedResults.Intensities,
                filteringParams,
                msnLevel: 2);
        }

        /// <summary>
        /// Merges multiple index and intensity arrays into an m/z array.
        /// Used when building the component spectra for an MS2 scan
        /// </summary>
        /// <param name="indexArrays">List of index arrays.</param>
        /// <param name="intensityArrays">List of intensity arrays.</param>
        /// <param name="proxyFactory">Frame proxy factory.</param>
        /// <returns>A tuple containing the merged m/z values and intensities.</returns>
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


        #endregion

    }
}