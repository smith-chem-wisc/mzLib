using System;
using System.Runtime.InteropServices;
using System.Text;
using MassSpectrometry;
using System.Data.SQLite;
using Easy.Common.Extensions;
using MzLibUtil;
using UsefulProteomicsDatabases;
using System.Data.Common;
using Readers.Bruker.TimsTofReader;

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

            //int dataLength = 2^12; // Default allocation ~ 16 kB
            //IntPtr pData = Marshal.AllocHGlobal(dataLength * Marshal.SizeOf<Int32>());

            //var outputLength = tims_read_scans_v2(
            //    fileHandle,
            //    ids[1000],
            //    scan_begin: 0,
            //    scan_end: (UInt32)numScans[1000],
            //    buffer: pData,
            //    length: (UInt32)(dataLength * 4));

            //var dataArray = new Int32[95578];
            //Marshal.Copy(pData, dataArray, 0, dataLength);

            int frameIdIndex = 65365;
            long frameId = ids[frameIdIndex];

            var dataArray = FrameProxy.GetScanRawData(fileHandle, ids[frameIdIndex], (UInt32)numScans[frameIdIndex]);
            FrameProxy test = new FrameProxy(dataArray, numScans[frameIdIndex], fileHandle);

            var max = dataArray.Max();

            int maxScan = dataArray[0..900].IndexOf(dataArray[0..900].Max());
            //int maxScan = 237;
            Range scanMzRange = test.GetXRange(maxScan);
            
            var mzSlice = dataArray[scanMzRange];
            double[] mzDubSlice = new double[mzSlice.Length];
            for(int j = 0; j < mzSlice.Length; j++)
            {
                mzDubSlice[j] = (double)mzSlice[j];
            }

            placeholded2++;

            var conversionTest = TimsConversion.DoTransformation(fileHandle,
                ids[frameIdIndex],
                mzDubSlice,
                TimsConversion.ConversionFunctions.IndexToMz);

            placeholded2++;

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

        

        #endregion Bruker Dll Functions

    }
}
