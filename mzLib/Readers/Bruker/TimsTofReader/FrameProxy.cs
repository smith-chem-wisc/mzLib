using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers.Bruker.TimsTofReader
{

    /// <summary>
    /// From Bruker documentation:
    /// Raw data layout: (N = scan_end - scan_begin = number of requested scans)
    ///   N x uint32_t: number of peaks in each of the N requested scans
    ///   N x (two uint32_t arrays: first indices, then intensities)
    /// </summary>
    internal class FrameProxy
    {
        private int[] _scanOffsets;
        private uint[] _rawData;
        public UInt64 FileHandle { get; }
        internal int NumberOfScans { get; }
        internal int TotalNumberOfPeaks => _scanOffsets[NumberOfScans - 1];

        internal FrameProxy(uint[] rawData, int numScans, UInt64 fileHandle)
        {
            NumberOfScans = numScans;
            FileHandle = fileHandle;
            _rawData = rawData;
            _scanOffsets = PartialSum(rawData, 0, numScans);
        }

        internal int GetNumberOfPeaks(int zeroIndexedScanNumber)
        {
            ThrowIfInvalidScanNumber(zeroIndexedScanNumber);
            return (int)_rawData[zeroIndexedScanNumber];
        }

        /// <summary>
        /// Returns a range containing the start(inclusive) and end (exclusive) indices
        /// for the segment of the _rawData array corresponding to the m/z lookup values for 
        /// a given scan
        /// </summary>
        internal Range GetXRange(int zeroIndexedScanNumber)
        {
            ThrowIfInvalidScanNumber(zeroIndexedScanNumber);
            return GetReadScanRange(zeroIndexedScanNumber, offset: 0);
        }

        /// <summary>
        /// Returns a range containing the start(inclusive) and end (exclusive) indices
        /// for the segment of the _rawData array corresponding to raw intensity values for a given scan
        /// </summary>
        internal Range GetYRange(int zeroIndexedScanNumber)
        {
            ThrowIfInvalidScanNumber(zeroIndexedScanNumber);
            return GetReadScanRange(zeroIndexedScanNumber, offset: (int)_rawData[zeroIndexedScanNumber]);
        }

        private void ThrowIfInvalidScanNumber(int zeroIndexedScanNumber)
        {
            if (zeroIndexedScanNumber < 0 || zeroIndexedScanNumber >= NumberOfScans)
                throw new ArgumentException("Scan number out of range.");
        }

        private Range GetReadScanRange(int zeroIndexedScanNumber, int offset)
        {
            int start = NumberOfScans + 2*_scanOffsets[zeroIndexedScanNumber] + offset;
            return new Range(start, start + (int)_rawData[zeroIndexedScanNumber]);
        }

        /// <summary>
        /// Calculates the running total of an array, beginning with 
        /// the start index (inclusive) and ending with the end index (exclusive).
        /// Used for determining scan offsets.
        /// </summary>
        /// <param name="array">Array to be summed</param>
        /// <param name="start">Where to begin summing</param>
        /// <param name="end">Where summing ends (exclusive) </param>
        /// <returns> An array of length (end - start) containing the 
        /// partial sums at each index of the input array</returns>
        public static int[] PartialSum(uint[] array, int start, int end)
        {
            int runningTotal = 0;
            int[] sums = new int[end - start];
            for(int i = 0; i < end; i++)
            {
                runningTotal += (int)array[i];
                sums[i] = runningTotal;
            }
            return sums;
        }
    }
}
