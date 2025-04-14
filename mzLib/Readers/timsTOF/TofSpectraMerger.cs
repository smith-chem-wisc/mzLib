using Easy.Common.Extensions;
using MassSpectrometry;
using MzLibUtil;
using System.Collections.Immutable;

namespace Readers
{
    public class TofSpectraMerger
    {

        #region NewCentroidingApproach
        // The following methods are used to merge and collapse index arrays and intensity arrays
        // The timsTOF data format doesn't store m/z values directly, but rather indices in a lookup table where the mz values are stored 
        // Keeping these indices as ints allows for more efficient storage and processing of the data

        // The general workflow is to read in individual scans as index and intensity arrays
        // Then, combine these arrays into a single array sorted by tofIndex, ascending
        // Then, collapse the combined array by merging entries with the same tofIndex and removing entries that fall below a certain intensity threshold (NoiseFloor)
        // Finally, perform centroiding by grouping adjacent tofIndices, summing their intensities, and then calculating a weighted average of the m/z values in the cluster

        public int? Ms1NoiseFloor { get; private set; }
        public int? Ms2NoiseFloor { get; private set; }

        public void SetNoiseFloor(List<TimsSpectrum> spectra, int msnLevel)
        {
            // Merge all index arrays and intensity arrays into a single array
            uint[] combinedIndices = spectra[0].XArray;
            int[] combinedIntensities = spectra[0].YArray;
            for (int i = 1; i < spectra.Count(); i++)
            {
                var mergeResults = TwoPointerMerge(combinedIndices, spectra[i].XArray, combinedIntensities, spectra[i].YArray);
                combinedIndices = mergeResults.Indices;
                combinedIntensities = mergeResults.Intensities;
            }

            Array.Sort(combinedIntensities);

            int median = combinedIntensities[combinedIntensities.Length / 2];
            int firstQuartile = combinedIntensities[combinedIntensities.Length / 4];
            int firstQuintile = combinedIntensities[combinedIntensities.Length / 5];
            int firstDecile = combinedIntensities[combinedIntensities.Length / 10];
            if (msnLevel == 1)
                Ms1NoiseFloor = median;
            else
                Ms2NoiseFloor = median;
        }

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
        internal (uint[] Indices, int[] Intensities) CollapseArrays(uint[] indexArray, int[] intensityArray, int msnLevel = 1, bool removeLowIntensityPeaks = true)
        {
            // Define lists to store the collapsed indices and intensities
            List<uint> collapsedIndices = new List<uint>(indexArray.Length);
            List<int> collapsedIntensities = new List<int>(intensityArray.Length);
            int noiseFloor = 0;
            if(removeLowIntensityPeaks && (Ms1NoiseFloor == null || Ms2NoiseFloor == null))
                throw new Exception("Cannot remove low intensity peaks without a noise floor. Set the noise floor before calling this method.");
            if(removeLowIntensityPeaks)
                noiseFloor = msnLevel == 1 ? (int)Ms1NoiseFloor : (int)Ms2NoiseFloor;

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

                if (!removeLowIntensityPeaks || summedIntensity > noiseFloor)
                { 
                    collapsedIndices.Add(currentIdx);
                    collapsedIntensities.Add(summedIntensity);
                }

                // Move the pointers forward
                p1 = p2 + 1;
                p2 = p1 + 1;
            }

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
        /// Each frame gets one TimsSpectrum
        /// </summary>
        /// <param name="indexArrays"></param>
        /// <param name="intensityArrays"></param>
        /// <param name="zeroIndexedTimsScanNumber"></param>
        /// <returns></returns>
        internal TimsSpectrum CreateTimsSpectrum(List<uint[]> indexArrays, List<int[]> intensityArrays, int msnLevel = 2, int timsScanIndex = -1, bool? removeLowIntensityPeaks = null)
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

            var collapsedResults = CollapseArrays(combinedIndices, combinedIntensities, msnLevel: msnLevel, removeLowIntensityPeaks: removeLowIntensityPeaks ?? msnLevel == 1);
            return new TimsSpectrum(collapsedResults.Indices, collapsedResults.Intensities, timsScanIndex);
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
        internal MzSpectrum CreateMzSpectrum(
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
        /// <returns>A merged  spectrum.</returns>
        internal MzSpectrum CreateMzSpectrum(
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
            var collapsedResults = CollapseArrays(combinedIndices, combinedIntensities, msnLevel);
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

    }
}