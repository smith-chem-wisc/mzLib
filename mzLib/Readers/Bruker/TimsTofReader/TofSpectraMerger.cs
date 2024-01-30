using Easy.Common.Extensions;
using MassSpectrometry;
using MzLibUtil;

namespace Readers.Bruker.TimsTofReader
{
    // The TofSpectraMerger implements an algorithm that merges an arbitray number of ordered lists.
    // The implementation is based off the code found here: https://leetcode.com/problems/merge-k-sorted-lists/solutions/3286058/image-explanation-5-methods-divide-conquer-priority-queue-complete-intuition/
    public static class TofSpectraMerger
    {
        public static MzSpectrum MergeSpectra(List<MzSpectrum> spectra, FilteringParams filteringParams = null, Tolerance tolerance = null)
        {
            List<double[]> mzArrays = new();
            List<double[]> intensityArrays = new();
            foreach (var spectrum in spectra)
            {
                mzArrays.Add(spectrum.XArray);
                intensityArrays.Add(spectrum.YArray);
            }

            return MergeSpectra(mzArrays, intensityArrays, tolerance);
        }

        // This isfor ms1
        public static MzSpectrum MergeSpectra(List<double[]> mzArrays, List<int[]> intensityArrays,
            FilteringParams filteringParams = null, Tolerance tolerance = null) 
        {
            tolerance ??= new PpmTolerance(5);
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

        // This is for msms
        internal static MzSpectrum MergeSpectra(List<ListNode<TofPeak>> spectrumHeadNodes,
            int allSpectraPeakCount, FilteringParams filteringParams = null, Tolerance tolerance = null)
        {
            tolerance ??= new PpmTolerance(5);

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
            tolerance ??= new PpmTolerance(5);
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
            tolerance ??= new PpmTolerance(5);
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
            => new TofPeak((peak1.Mz + peak2.Mz) / 2, peak1.Intensity + peak2.Intensity);

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
