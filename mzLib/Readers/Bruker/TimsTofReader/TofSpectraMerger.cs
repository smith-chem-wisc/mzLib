using Easy.Common.Extensions;
using MassSpectrometry;
using MzLibUtil;

namespace Readers.Bruker
{
    // The TofSpectraMerger implements an algorithm that merges an arbitray number of ordered lists.
    // The implementation is based off the code found here: https://leetcode.com/problems/merge-k-sorted-lists/solutions/3286058/image-explanation-5-methods-divide-conquer-priority-queue-complete-intuition/
    public static class TofSpectraMerger
    {
        public static readonly double DefaultPpmTolerance = 5;

        internal static MzSpectrum MergeArraysToSpectrum(List<uint[]> indexArrays, List<int[]> intensityArrays, FrameProxyFactory proxyFactory,
            FilteringParams filteringParams = null)
        {
            if (!indexArrays.IsNotNullOrEmpty() || intensityArrays == null || intensityArrays.Count() != indexArrays.Count())
                return null;

            // Merge all index arrays and intensity arrays into a single array
            uint[] combinedIndices = indexArrays[0];
            int[] combinedIntensities = intensityArrays[0];
            for(int i = 1; i < indexArrays.Count(); i++)
            {
                var mergeResults = TwoPointerMerge(combinedIndices, indexArrays[i], combinedIntensities, intensityArrays[i]);
                combinedIndices = mergeResults.Indices;
                combinedIntensities = mergeResults.Intensities;
            }

            // Collapse the combined arrays into a single array (centroiding, more or less)
            var centroidedResults = CollapseArrays(combinedIndices, combinedIntensities);

            return ConvertAndFilterSpectrum(centroidedResults.Indices, centroidedResults.Intensities, proxyFactory, filteringParams);
        }

        internal static MzSpectrum ConvertAndFilterSpectrum(IList<uint> indices, IList<int> intensities,
            FrameProxyFactory proxyFactory, FilteringParams filteringParams = null)
        {
            // Convert the indices to m/z values
            double[] mzsArray = proxyFactory.ConvertIndicesToMz(indices);

            // Convert the intensities to an array
            double[] intensitiesArray = intensities.Select(intensity => (double)intensity).ToArray();

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

        public static MzSpectrum MergesMs1Spectra(List<double[]> mzArrays, List<int[]> intensityArrays,
            FilteringParams filteringParams = null, Tolerance tolerance = null) 
        {
            tolerance ??= new PpmTolerance(DefaultPpmTolerance);
            if (!mzArrays.IsNotNullOrEmpty() || intensityArrays == null || intensityArrays.Count() != mzArrays.Count())
                return null;

            int preMergePeaks = 0;
            List<ListNode<TofPeak>> peakNodes = new();
            for (int i = 0; i < mzArrays.Count(); i++)
            {
                if (mzArrays[i].Length < 1) continue;

                ListNode<TofPeak> peakNode = new ListNode<TofPeak>(
                    new TofPeak(mzArrays[i][0], intensityArrays[i][0]));
                peakNodes.Add(peakNode);
                preMergePeaks++;
                for (int j = 1; j < mzArrays[i].Length; j++)
                {
                    peakNode.Next = new ListNode<TofPeak>(
                        new TofPeak(mzArrays[i][j], intensityArrays[i][j]));
                    peakNode = peakNode.Next;
                    preMergePeaks++;
                }
            }
            var mergedListHead = MergeSpectraHelper(peakNodes, 0, peakNodes.Count-1, tolerance, out int mergedPeaksCount);
            int finalPeaksCount = preMergePeaks - mergedPeaksCount;

            double[] mzs = new double[finalPeaksCount];
            double[] intensities = new double[finalPeaksCount];
            for(int i = 0; i < finalPeaksCount; i++)
            {
                if (mergedListHead == null) throw new Exception("Disagreement in peak count after recursive merge");
                mzs[i] = mergedListHead.Value.Mz;
                intensities[i] = mergedListHead.Value.Intensity;
                mergedListHead = mergedListHead.Next;
            }

            if (filteringParams != null
                && mzs.Length > 0
                && filteringParams.ApplyTrimmingToMs1)
            {
                WindowModeHelper.Run(ref intensities,
                    ref mzs, filteringParams,
                    mzs[0], mzs[^0]);
            }
            return new MzSpectrum(mzs, intensities, shouldCopy: false);
        }


        public static (uint[] Indices, int[] Intensities) TwoPointerMerge(uint[] indexArray1, uint[] indexArray2, int[] intensityArray1, int[] intensityArray2)
        {
            int p1 = 0;
            int p2 = 0;

            uint[] mergedIndices = new uint[indexArray1.Length + indexArray2.Length];
            int[] mergedIntensities = new int[intensityArray1.Length + intensityArray2.Length];

            while(p1 < indexArray1.Length || p2 < indexArray2.Length)
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
            while(p1 < indexArray.Length)
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

        internal static MzSpectrum MergeMsmsSpectra(List<ListNode<TofPeak>> spectrumHeadNodes,
            int allSpectraPeakCount, FilteringParams filteringParams = null, Tolerance tolerance = null)
        {
            tolerance ??= new PpmTolerance(DefaultPpmTolerance);

            var mergedListHead = MergeSpectraHelper(spectrumHeadNodes, 0, spectrumHeadNodes.Count - 1, tolerance, out int mergedPeaksCount);
            int finalPeaksCount = allSpectraPeakCount - mergedPeaksCount;

            double[] mzs = new double[finalPeaksCount];
            double[] intensities = new double[finalPeaksCount];
            for (int i = 0; i < finalPeaksCount; i++)
            {
                if (mergedListHead == null) throw new Exception("Disagreement in peak count after recursive merge");
                mzs[i] = mergedListHead.Value.Mz;
                intensities[i] = mergedListHead.Value.Intensity;
                mergedListHead = mergedListHead.Next;
            }

            if (filteringParams != null
                && mzs.Length > 0
                && filteringParams.ApplyTrimmingToMsMs)
            {
                WindowModeHelper.Run(ref intensities,
                    ref mzs, filteringParams,
                    mzs[0], mzs[finalPeaksCount-1]); // The unary operator was throwing "index out of range" errors here
            }
            return new MzSpectrum(mzs, intensities, shouldCopy: false);
        }

        internal static ListNode<TofPeak> MergeSpectraToLinkedList(List<double[]> mzArrays, List<int[]> intensityArrays, out int listLength,
           Tolerance tolerance = null)
        {
            tolerance ??= new PpmTolerance(DefaultPpmTolerance);
            listLength = 0;
            if (!mzArrays.IsNotNullOrEmpty() || intensityArrays == null || intensityArrays.Count() != mzArrays.Count())
                return null;

            int preMergePeaks = 0;
            List<ListNode<TofPeak>> peakNodes = new();
            for (int i = 0; i < mzArrays.Count(); i++)
            {
                if (mzArrays[i].Length < 1) continue;

                ListNode<TofPeak> peakNode = new ListNode<TofPeak>(
                    new TofPeak(mzArrays[i][0], intensityArrays[i][0]));
                peakNodes.Add(peakNode);
                preMergePeaks++;
                for (int j = 1; j < mzArrays[i].Length; j++)
                {
                    peakNode.Next = new ListNode<TofPeak>(
                        new TofPeak(mzArrays[i][j], intensityArrays[i][j]));
                    peakNode = peakNode.Next;
                    preMergePeaks++;
                }
            }
            var mergedListHead = MergeSpectraHelper(peakNodes, 0, peakNodes.Count - 1, tolerance, out int mergedPeaksCount);
            listLength = preMergePeaks - mergedPeaksCount;

            return mergedListHead;
        }

        /// <summary>
        /// This is gonna coerce all doubles in the intensity arrays to ints and back again. Sorry babes
        /// </summary>
        /// <param name="mzArrays"></param>
        /// <param name="intensityArrays"></param>
        /// <param name="tolerance"></param>
        /// <returns></returns>
        /// <exception cref="Exception"></exception>
        public static MzSpectrum MergeSpectra(List<double[]> mzArrays, List<double[]> intensityArrays,
           Tolerance tolerance = null)
        {
            tolerance ??= new PpmTolerance(DefaultPpmTolerance);
            if (!mzArrays.IsNotNullOrEmpty() || intensityArrays == null || intensityArrays.Count() != mzArrays.Count())
                return null;

            int preMergePeaks = 0;
            List<ListNode<TofPeak>> peakNodes = new();
            for (int i = 0; i < mzArrays.Count(); i++)
            {
                if (mzArrays[i].Length < 1) continue;

                ListNode<TofPeak> peakNode = new ListNode<TofPeak>(
                    new TofPeak(mzArrays[i][0], (int)intensityArrays[i][0]));
                peakNodes.Add(peakNode);
                preMergePeaks++;
                for (int j = 1; j < mzArrays[i].Length; j++)
                {
                    peakNode.Next = new ListNode<TofPeak>(
                        new TofPeak(mzArrays[i][j], (int)intensityArrays[i][j]));
                    peakNode = peakNode.Next;
                    preMergePeaks++;
                }
            }
            var mergedListHead = MergeSpectraHelper(peakNodes, 0, peakNodes.Count - 1, tolerance, out int mergedPeaksCount);
            int finalPeaksCount = preMergePeaks - mergedPeaksCount;

            double[] mzs = new double[finalPeaksCount];
            double[] intensities = new double[finalPeaksCount];
            for (int i = 0; i < finalPeaksCount; i++)
            {
                if (mergedListHead == null) throw new Exception("Disagreement in peak count after recursive merge");
                mzs[i] = mergedListHead.Value.Mz;
                intensities[i] = mergedListHead.Value.Intensity;
                mergedListHead = mergedListHead.Next;
            }

            return new MzSpectrum(mzs, intensities, shouldCopy: false);
        }

        internal static ListNode<TofPeak> MergeSpectraHelper(
            List<ListNode<TofPeak>> peakNodes,
            int nodeListStart,
            int nodeListEnd,
            Tolerance tolerance,
            out int mergedPeaksCount)
        {
            mergedPeaksCount = 0;
            if (nodeListStart > nodeListEnd) return null;
            if (nodeListStart == nodeListEnd) return peakNodes[nodeListStart];

            int nodeListMid = nodeListStart + (nodeListEnd - nodeListStart)/2;
            ListNode<TofPeak> left = MergeSpectraHelper(peakNodes, nodeListStart, nodeListMid, tolerance, out int partialMergedPeakCount);
            mergedPeaksCount += partialMergedPeakCount;
            ListNode<TofPeak> right = MergeSpectraHelper(peakNodes, nodeListMid+1, nodeListEnd, tolerance, out partialMergedPeakCount);
            mergedPeaksCount += partialMergedPeakCount;
            ListNode<TofPeak> mergedHeadNode = Merge(left, right, tolerance, out partialMergedPeakCount);
            mergedPeaksCount += partialMergedPeakCount;
            return mergedHeadNode; // add tolerance here!
        }

        internal static ListNode<TofPeak> Merge(
            ListNode<TofPeak> peak1,
            ListNode<TofPeak> peak2,
            Tolerance tolerance,
            out int mergedPeakCount)
        {
            // Set temporary head/tail with null values
            ListNode<TofPeak> dummyHead = new ListNode<TofPeak>(initElement: null);
            ListNode<TofPeak> tail = dummyHead;
            mergedPeakCount = 0;

            while (peak1 != null && peak2 != null)
            {
                int comparison = peak1.Value.CompareTo(peak2.Value, tolerance);
                if (comparison == 0)
                {
                    tail.Next = new ListNode<TofPeak>(peak1.Value + peak2.Value);
                    peak1 = peak1.Next;
                    peak2 = peak2.Next;
                    mergedPeakCount++;
                }
                else if (comparison < 0)
                {
                    tail.Next = peak1;
                    peak1 = peak1.Next;
                }
                else
                {
                    tail.Next = peak2;
                    peak2 = peak2.Next;
                }
                tail = tail.Next;
            }
            // If either list has remaining elements, append them to the end of the linked list
            tail.Next = peak1 ?? peak2;
            // Return the head node of the list
            return dummyHead.Next;
        }
    }

    public class TofPeak : IComparable
    {
        public double Mz { get; }
        public int Intensity { get; }

        public TofPeak(double mz, int intensity)
        {
            Mz = mz;
            Intensity = intensity;
        }

        public static TofPeak operator +(TofPeak peak1, TofPeak peak2)
        { 
            int summedIntensity = peak1.Intensity + peak2.Intensity;
            return new TofPeak( (peak1.Mz * peak1.Intensity + peak2.Mz* peak2.Intensity) / (summedIntensity), summedIntensity);
        }

        public int CompareTo(Object obj)
        {
            if (obj == null) return 1;

            TofPeak otherPeak = obj as TofPeak;
            if (otherPeak == null) throw new ArgumentException("Object is not a TofPeak");
            return this.Mz.CompareTo(otherPeak.Mz);
        }

        public int CompareTo(TofPeak otherPeak, Tolerance tol)
        {
            if (tol.Within(this.Mz, otherPeak.Mz)) return 0; // peaks are equal within tolerance
            return this.Mz.CompareTo(otherPeak.Mz);
        }

        public override string ToString()
        {
            return string.Format("({0:G7},{1:G7})", Mz, Intensity);
        }
    }

    // This class is a custom implementation of a Linked List. The LinkedList class built into 
    // C# is not compatible with the SpectraMerger algorithm implemented above. 
    // this class is based on the code found here: https://stackoverflow.com/questions/40221059/how-do-i-set-the-next-property-of-a-node-in-a-linkedlist-to-another-node-in-anot
    class ListNode<T> : IEnumerable<T>
    {
        public T Value;
        public ListNode<T> Next;
        private ListNode() { }
        public ListNode(IEnumerable<T> init) 
        {
            ListNode<T> current = null;
            foreach (T elem in init)
            {
                if (current == null) current = this; else current = current.Next = new ListNode<T>();
                current.Value = elem;
            }
        }

        public ListNode(T initElement)
        {
            ListNode<T> current = this;
            current.Value = initElement;
        }

        class ListEnum : IEnumerator<T>
        {
            private ListNode<T> first;
            private ListNode<T> current;
            bool more;
            public ListEnum(ListNode<T> first) { this.first = first; more = true; }
            public T Current { get { return current.Value; } }
            public void Dispose() { }
            object System.Collections.IEnumerator.Current { get { return current.Value; } }
            public void Reset() { current = null; more = true; }
            public bool MoveNext()
            {
                if (!more)
                    return false;
                else if (current == null)
                    return more = ((current = first) != null);
                else
                    return more = ((current = current.Next) != null);
            }
        }

        public IEnumerator<T> GetEnumerator()
        {
            return new ListEnum(this);
        }

        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }
    }
}
