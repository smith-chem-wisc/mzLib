using Easy.Common.Extensions;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection.Metadata;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace Readers
{
    internal enum ConversionFunctions
    {
        IndexToMz,
        MzToIndex,
        ScanToOneOverK0,
        OneOverK0ToScan,
        ScanToVoltage,
        VoltageToScan,
        IndexToMzTsf,
        MzToIndexTsf
    }

    internal unsafe class TimsConversion
    {

        private UInt64 _fileHandle;
        private Object _fileLock;

        internal TimsConversion(UInt64 fileHandle, Object fileLock)
        {
            _fileHandle = fileHandle;
            _fileLock = fileLock;
        }

        /// <summary>
        /// Takes an array of raw values and converts them according to the specified conversion function, 
        /// returning an equal length array containing the transformed values
        /// </summary>
        /// <param name="fileHandle"> Unique identifier associated with the open timsTof .d data file </param>
        /// <param name="frameId"> Frame identified</param>
        /// <returns> Double array containing the transformed input values </returns>
        internal unsafe double[] DoTransformation(UInt64 fileHandle, long frameId, double[] input, ConversionFunctions function)
        {
            if(!input.IsNotNullOrEmpty())
            {
                return Array.Empty<double>();
            }
            double[] transformedValues = new double[input.Length];
            fixed (double* inputPtr = &input[0])
            {
                IntPtr outPtr = Marshal.AllocHGlobal(input.Length * Marshal.SizeOf<double>());
                try
                {
                    lock (_fileLock)
                    {
                        switch (function)
                        {
                            case ConversionFunctions.IndexToMz:
                                tims_index_to_mz(fileHandle, frameId, inputPtr, (double*)outPtr, (UInt32)input.Length);
                                break;
                            case ConversionFunctions.MzToIndex:
                                tims_mz_to_index(fileHandle, frameId, inputPtr, (double*)outPtr, (UInt32)input.Length);
                                break;
                            case ConversionFunctions.ScanToOneOverK0:
                                tims_scannum_to_oneoverk0(fileHandle, frameId, inputPtr, (double*)outPtr, (UInt32)input.Length);
                                break;
                            case ConversionFunctions.OneOverK0ToScan:
                                tims_oneoverk0_to_scannum(fileHandle, frameId, inputPtr, (double*)outPtr, (UInt32)input.Length);
                                break;
                            case ConversionFunctions.ScanToVoltage:
                                tims_scannum_to_voltage(fileHandle, frameId, inputPtr, (double*)outPtr, (UInt32)input.Length);
                                break;
                            case ConversionFunctions.VoltageToScan:
                                tims_voltage_to_scannum(fileHandle, frameId, inputPtr, (double*)outPtr, (UInt32)input.Length);
                                break;
                            case ConversionFunctions.IndexToMzTsf:
                                tsf_index_to_mz(fileHandle, frameId, inputPtr, (double*)outPtr, (UInt32)input.Length);
                                break;
                            case ConversionFunctions.MzToIndexTsf:
                                tsf_mz_to_index(fileHandle, frameId, inputPtr, (double*)outPtr, (UInt32)input.Length);
                                break;
                        }
                    }

                    Marshal.Copy(outPtr, transformedValues, 0, input.Length);
                }
                finally { Marshal.FreeHGlobal(outPtr); }
            }
            return transformedValues;
        }

        [DllImport("timsdata.dll", CharSet = CharSet.Ansi, CallingConvention = CallingConvention.Cdecl)]
        unsafe static extern void tims_index_to_mz
              (UInt64 fileHandle, Int64 frame_id, double* inputPtr, double* outPtr, UInt32 count);
        [DllImport("timsdata.dll", CharSet = CharSet.Ansi, CallingConvention = CallingConvention.Cdecl)]
        unsafe static extern void tims_mz_to_index
              (UInt64 fileHandle, Int64 frame_id, double* inputPtr, double* outPtr, UInt32 count);
        [DllImport("timsdata.dll", CharSet = CharSet.Ansi, CallingConvention = CallingConvention.Cdecl)]
        unsafe static extern void tims_scannum_to_oneoverk0
              (UInt64 fileHandle, Int64 frame_id, double* inputPtr, double* outPtr, UInt32 count);
        [DllImport("timsdata.dll", CharSet = CharSet.Ansi, CallingConvention = CallingConvention.Cdecl)]
        unsafe static extern void tims_oneoverk0_to_scannum
              (UInt64 fileHandle, Int64 frame_id, double* inputPtr, double* outPtr, UInt32 count);
        [DllImport("timsdata.dll", CharSet = CharSet.Ansi, CallingConvention = CallingConvention.Cdecl)]
        unsafe static extern void tims_scannum_to_voltage
              (UInt64 fileHandle, Int64 frame_id, double* inputPtr, double* outPtr, UInt32 count);
        [DllImport("timsdata.dll", CharSet = CharSet.Ansi, CallingConvention = CallingConvention.Cdecl)]
        unsafe static extern void tims_voltage_to_scannum
              (UInt64 fileHandle, Int64 frame_id, double* inputPtr, double* outPtr, UInt32 count);
        [DllImport("timsdata.dll", CharSet = CharSet.Ansi, CallingConvention = CallingConvention.Cdecl)]
        unsafe static extern void tsf_index_to_mz
            (UInt64 fileHandle, Int64 frame_id, double* inputPtr, double* outPtr, UInt32 count);
        [DllImport("timsdata.dll", CharSet = CharSet.Ansi, CallingConvention = CallingConvention.Cdecl)]
        unsafe static extern void tsf_mz_to_index
            (UInt64 fileHandle, Int64 frame_id, double* inputPtr, double* outPtr, UInt32 count);
    }
}
