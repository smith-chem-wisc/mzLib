using System;
using System.Runtime.InteropServices;
using System.Text;
using MassSpectrometry;
using System.Data.SQLite;
using Easy.Common.Extensions;
using MzLibUtil;
using UsefulProteomicsDatabases;
using System.Data.Common;

namespace Readers.Bruker
{ 
    public static class TimsTofFileReader
    {

        public const int InitialFrameBufferSize = 128;

        public unsafe static void ReadData(string pathToDFile)
        {
            byte[] binBytePath = BrukerFileReader.ConvertStringToUTF8ByteArray(pathToDFile);
            UInt64 fileHandle = tims_open(binBytePath, 0);

            var connection = new SQLiteConnection();
            connection.ConnectionString = "DataSource=" + Path.Combine(pathToDFile, "analysis.tdf");

            connection.Open();

            using var command = new SQLiteCommand(connection);
            command.CommandText = @"SELECT COUNT(*) FROM Frames;";
            using var sqliteReader = command.ExecuteReader();
            int count = 0;
            while (sqliteReader.Read())
            {
                count = sqliteReader.GetInt32(0);
                break;
            }


            sqliteReader.Close();

            int placeholder = 0;

            command.CommandText = @"SELECT * FROM Frames;";
            using var sqliteReader2 = command.ExecuteReader();

            var columns = Enumerable.Range(0, sqliteReader2.FieldCount)
                .Select(sqliteReader2.GetName).ToList();

            int idIndex = sqliteReader2.GetOrdinal("Id");
            int scanCountIndex = sqliteReader2.GetOrdinal("NumScans");
            int timeIndex = sqliteReader2.GetOrdinal("Time");
            int scanModeIndex = sqliteReader2.GetOrdinal("ScanMode");
            int msMsTypeIndex = sqliteReader2.GetOrdinal("MsMsType");
            int numPeaksIndex = sqliteReader2.GetOrdinal("NumPeaks");

            long[] ids = new long[count];
            int[] numScans = new int[count];
            float[] times = new float[count];
            int[] scanMode = new int[count];
            int[] msMsTypes = new int[count];
            int[] numPeaks = new int[count];
            int i = 0;
            while(sqliteReader2.Read())
            {
                ids[i] = sqliteReader2.GetInt64(idIndex);
                numScans[i] = sqliteReader2.GetInt32(scanCountIndex);
                times[i] = sqliteReader2.GetFloat(timeIndex);
                scanMode[i] = sqliteReader2.GetInt32(scanModeIndex);
                msMsTypes[i] = sqliteReader2.GetInt32(msMsTypeIndex);
                numPeaks[i++] = sqliteReader2.GetInt32(numPeaksIndex);

                if(i == 100)
                {
                    placeholder++;
                }
            }
            int placeholded2 = 1;


            var scancountMean = numScans.Select(u => (int)u).Average();
            var totalScans = scancountMean * count;

            placeholded2++;

            float finalTime = times[count - 1];
            var scansPerHour = totalScans / (finalTime / 60);

            placeholded2++;

            //var shortPtrs = new IntPtr[95578];
            
            // Note, the data is actually stored as UInt32, will need 
            // to do conversion via type punning after marshalling
            //var dataArray = new Int32[95578];
            //GCHandle hData = GCHandle.Alloc(dataArray, GCHandleType.Pinned);
            //IntPtr pData = hData.AddrOfPinnedObject();

            int dataLength = 2^12; // Default allocation ~ 16 kB
            IntPtr pData = Marshal.AllocHGlobal(dataLength * Marshal.SizeOf<Int32>());

            var outputLength = tims_read_scans_v2(
                fileHandle,
                ids[1000],
                scan_begin: 0,
                scan_end: (UInt32)numScans[1000],
                buffer: pData,
                length: (UInt32)(dataLength * 4));

            var dataArray = new Int32[95578];
            Marshal.Copy(pData, dataArray, 0, dataLength);

            var max = dataArray.Max();

            placeholded2++;

        }

        /// <summary>
        /// This is reimplementation of the Marshal.Copy method that allows for arbitrary types
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="source"></param>
        /// <param name="destination"></param>
        /// <param name="startIndex"></param>
        /// <param name="length"></param>
        /// <exception cref="ArgumentNullException"></exception>
        /// <exception cref="ArgumentOutOfRangeException"></exception>
        private static unsafe void CopyToManaged<T>(IntPtr source, T[] destination, int startIndex, int length)
        {
            if (source == IntPtr.Zero) throw new ArgumentNullException(nameof(source));
            if (destination is null) throw new ArgumentNullException(nameof(destination));
            if (startIndex < 0) throw new ArgumentOutOfRangeException(nameof(startIndex));
            if (length < 0) throw new ArgumentOutOfRangeException(nameof(length));

            void* sourcePtr = (void*)source;
            Span<T> srcSpan = new Span<T>(sourcePtr, length);
            Span<T> destSpan = new Span<T>(destination, startIndex, length);

            srcSpan.CopyTo(destSpan);
        }

        internal unsafe static uint[] GetScanRawData(UInt64 fileHandle, long frameId, UInt32 numScans)
        {
            int dataLength = 4096; // Default allocation ~ 16 kB
            
            // buffer expansion loop
            while(true)
            {
                IntPtr pData = Marshal.AllocHGlobal(dataLength * Marshal.SizeOf<Int32>());

                var outputLength = tims_read_scans_v2(
                    fileHandle,
                    frameId,
                    scan_begin: 0,
                    scan_end: numScans,
                    buffer: pData,
                    length: (uint)(dataLength * 4));

                if (4 * dataLength > outputLength)
                {
                    var dataArray = new uint[dataLength];
                    CopyToManaged(pData, dataArray, 0, dataLength);
                    Marshal.FreeHGlobal(pData);

                    return dataArray;
                }

                if(outputLength > 16777216 ) // Arbitrary 16 mb frame limit
                {
                    throw new Exception("Maximum frame size exceeded");
                }

                // Increase buffer size if necessary
                dataLength = ((int)outputLength / 4) + 1;

                Marshal.FreeHGlobal(pData);
            }
        }

        internal unsafe static void ParseScan(uint[] rawData)
        {

        }


        #region Bruker Dll Functions 

        /// <summary>
        /// Returns a unique handle that references an open timsTOF data file
        /// </summary>
        /// <param name="analysis_directory_name_utf8"></param>
        /// <param name="use_recalibrated_state"></param>
        /// <returns></returns>
        [DllImport("Bruker/TimsTofReader/timsdata.dll")]
        public static extern UInt64 tims_open
              (byte[] analysis_directory_name_utf8, UInt32 use_recalibrated_state);

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
        [DllImport("Bruker/TimsTofReader/timsdata.dll", CallingConvention = CallingConvention.Cdecl)]
        unsafe static extern UInt32 tims_read_scans_v2
              (UInt64 handle, Int64 frame_id, UInt32 scan_begin, UInt32 scan_end, IntPtr buffer, UInt32 length);

        #endregion Bruker Dll Functions

    }
}
