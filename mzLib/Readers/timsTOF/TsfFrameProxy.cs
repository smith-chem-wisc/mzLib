using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace Readers
{
    //internal class TsfFrameProxyFactory: FrameProxyFactory
    //{
    //    internal TsfFrameProxyFactory(FrameTable framesTable, UInt64 fileHandle, Object fileLock, int maxIndex) 
    //        : base(framesTable, fileHandle, fileLock, maxIndex)
    //    {
    //        // Initialize the factory with the provided parameters
    //        // Additional initialization logic can be added here if needed
    //    }

    //    internal override TsfFrameProxy GetFrameProxy(long frameId) 
    //    {
    //        // Create a new instance of TsfFrameProxy with the provided parameters
    //        return new TsfFrameProxy(FileHandle, frameId, FramesTable.NumScans[frameId - 1], FileLock, Converter);
    //    }
    //}

    internal class TsfFrameProxy : FrameProxy
    {
        public TsfFrameProxy(UInt64 fileHandle, long frameId, int numScans, Object fileLock, TimsConversion converter) 
            : base(fileHandle, frameId, numScans, fileLock, converter)
        {
            // Initialize the frame proxy with the provided file path and native ID format
            // Additional initialization logic can be added here if needed
        }

        //public override string GetFileType()
        //{
        //    return "TsfFrame"; // Return the type of frame this proxy represents
        //}

        /// <summary>
        /// Read a range of scans from a single frame.
        ///
        /// Output layout: (N = scan_end - scan_begin = number of requested scans)
        ///   N x uint32_t: number of peaks in each of the N requested scans
        ///   N x (two uint32_t arrays: first indices, then intensities)
        ///
        /// Note: different threads must not read scans from the same storage handle
        /// concurrently.
        /// </summary> 
        internal static uint[] GetLineSpectrumData(UInt64 fileHandle, long frameId, UInt32 numScans, Object fileLock)
        {
            int bufferSize = _defaultBufferSize;
            // buffer expansion loop
            while (true)
            {
                IntPtr indices = Marshal.AllocHGlobal(bufferSize * Marshal.SizeOf<double>());
                IntPtr intensities = Marshal.AllocHGlobal(bufferSize * Marshal.SizeOf<float>());
                try
                {
                    uint outputLength;

                    lock (fileLock)
                    {
                        outputLength = tims_read_line_spectrum_v2(
                            fileHandle,
                            frameId,
                            indices,
                            intensities,
                            length: (uint)(bufferSize * 4));
                    }

                    if (4 * bufferSize > outputLength)
                    {
                        var dataArray = new double[bufferSize];
                        CopyToManaged(indices, dataArray, 0, bufferSize);

                        //return dataArray;
                    }

                    if (outputLength > 16777216) // Arbitrary 16 mb frame limit
                    {
                        throw new Exception("Maximum frame size exceeded");
                    }

                    // Increase buffer size if necessary
                    bufferSize = ((int)outputLength / 4) + 1;
                }
                finally { 
                    Marshal.FreeHGlobal(indices);
                    Marshal.FreeHGlobal(intensities);
                }
            }
        }

        /// <summary>
        /// Read a range of scans from a single frame.
        ///
        /// Output layout: (N = scan_end - scan_begin = number of requested scans)
        ///   N x uint32_t: number of peaks in each of the N requested scans
        ///   N x (two uint32_t arrays: first indices, then intensities)
        ///
        /// Note: different threads must not read scans from the same storage handle
        /// concurrently.
        /// </summary> 
        /// <param name="handle"> Unique Handle of .d file ( returned on tims_open() )</param>
        /// <param name="frame_id"> From .tdf SQLite: Frames.Id </param>
        /// <param name="scan_begin"> first scan number to read (inclusive) </param>
        /// <param name="scan_end"> Last scan number (exclusive) </param>
        /// <param name="buffer"> Destination buffer allocated by user </param>
        /// <param name="length"> Length of the buffer (in bytes, i.e. 4 * buffer.length) </param>
        /// <returns> 0 on error, otherwise the number of buffer bytes necessary for the output
        /// of this call (if this is larger than the provided buffer length, the result is not
        /// complete). </returns>
        [DllImport("timsdata.dll", CallingConvention = CallingConvention.Cdecl)]
        unsafe static extern UInt32 tims_read_line_spectrum_v2
            (UInt64 handle, Int64 spectrum_id, IntPtr index_array, IntPtr intensity_array, UInt32 length);

        /// <summary>
        /// Read a range of scans from a single frame.
        ///
        /// Output layout: (N = scan_end - scan_begin = number of requested scans)
        ///   N x uint32_t: number of peaks in each of the N requested scans
        ///   N x (two uint32_t arrays: first indices, then intensities)
        ///
        /// Note: different threads must not read scans from the same storage handle
        /// concurrently.
        /// </summary> 
        /// <param name="handle"> Unique Handle of .d file ( returned on tims_open() )</param>
        /// <param name="spectrum_id"> From .tsf SQLite: Frames.Id </param>
        /// <param name="buffer"> Destination buffer allocated by user </param>
        /// <param name="length"> Length of the buffer (in bytes, i.e. 4 * buffer.length) </param>
        /// <returns> 0 on error, otherwise the number of buffer bytes necessary for the output
        /// of this call (if this is larger than the provided buffer length, the result is not
        /// complete). </returns>
        [DllImport("timsdata.dll", CallingConvention = CallingConvention.Cdecl)]
        unsafe static extern UInt32 tims_read_profile_spectrum_v2
            (UInt64 handle, Int64 spectrum_id, IntPtr buffer, UInt32 length);

    }
}
