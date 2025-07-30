using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace Readers
{
    internal class TsfFrameProxy : FrameProxy
    {
        public TsfFrameProxy(UInt64 fileHandle, long frameId, int numScans, Object fileLock, TimsConversion converter) 
            : base()
        {
            NumberOfScans = numScans;
            FileHandle = fileHandle;
            FrameId = frameId;
            Converter = converter;
            FileLock = fileLock;

            var lineSpectrumTuple = GetLineSpectrumData(FileHandle, FrameId, (uint)NumberOfScans, FileLock);
            IndexArray = lineSpectrumTuple.Item1;
            IntensityArray = lineSpectrumTuple.Item2;
        }

        public double[] IndexArray { get; }
        public float[] IntensityArray { get; private set; }


        /// <summary>
        /// Read a centroided spectrum from a single tsf frame.
        /// </summary> 
        internal static (double[], float[]) GetLineSpectrumData(UInt64 fileHandle, long frameId, UInt32 numScans, Object fileLock)
        {
            int bufferSize = 1024;
            // buffer expansion loop
            while (true)
            {
                IntPtr indices = Marshal.AllocHGlobal(bufferSize * Marshal.SizeOf<double>());
                IntPtr intensities = Marshal.AllocHGlobal(bufferSize * Marshal.SizeOf<float>());
                try
                {
                    int outputLength;

                    lock (fileLock)
                    {
                        outputLength = tsf_read_line_spectrum_v2(
                            fileHandle,
                            frameId,
                            indices,
                            intensities,
                            length: bufferSize);
                    }

                    if (bufferSize >= outputLength)
                    {
                        var indexArray = new double[outputLength];
                        CopyToManaged(indices, indexArray, 0, outputLength);

                        var intensityArray = new float[outputLength];
                        CopyToManaged(intensities, intensityArray, 0, outputLength);

                        return (indexArray, intensityArray);
                    }

                    if (outputLength > 16777216) // Arbitrary 16 mb frame limit
                    {
                        throw new Exception("Maximum frame size exceeded");
                    }

                    // Increase buffer size if necessary
                    bufferSize = outputLength + 1;
                }
                finally { 
                    Marshal.FreeHGlobal(indices);
                    Marshal.FreeHGlobal(intensities);
                }
            }
        }

        /// <summary>
        /// Read a centroided spectrum from a single tsf frame.
        ///
        /// Note: different threads must not read scans from the same storage handle
        /// concurrently.
        /// </summary> 
        /// <param name="handle"> Unique Handle of .d file ( returned on tims_open() )</param>
        /// <param name="spectrum_id"> From .tsf SQLite: Frames.Id </param>
        /// <param name="index_array"> pointer to the address of a double array to store the indices of the peaks </param>
        /// <param name="intensity_array"> pointer to the address of a float array to store the intensities of the peaks </param>
        /// <param name="length"> number of peaks/ entries in the arrays </param>
        /// <returns> 0 on error, otherwise the number of entries necessary for the output arrays
        /// of this call (if this is larger than the provided output array length, the result is not
        /// complete). </returns>
        [DllImport("timsdata.dll", CallingConvention = CallingConvention.Cdecl)]
        unsafe static extern Int32 tsf_read_line_spectrum_v2
            (UInt64 handle, Int64 spectrum_id, IntPtr index_array, IntPtr intensity_array, Int32 length);

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
        unsafe static extern UInt32 tsf_read_profile_spectrum_v2
            (UInt64 handle, Int64 spectrum_id, IntPtr buffer, UInt32 length);

    }
}
