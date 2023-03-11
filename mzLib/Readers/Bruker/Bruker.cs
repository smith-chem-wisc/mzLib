using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using System.Data.SQLite;
using Easy.Common.Extensions;

namespace Readers.Bruker
{
	public class Bruker : MsDataFile
	{
		private SQLiteConnection _connection;
		private ulong _handle; 
		public override MsDataFile LoadAllStaticData(FilteringParams filteringParams = null, int maxThreads = 1)
		{
			throw new NotImplementedException();
		}

		public override SourceFile GetSourceFile()
		{
			throw new NotImplementedException();
		}

		public override MsDataScan GetOneBasedScanFromDynamicConnection(int oneBasedScanNumber, IFilteringParams filterParams = null)
		{
			throw new NotImplementedException();
		}

		public override void CloseDynamicConnection()
		{
			throw new NotImplementedException();
		}

		public override void InitiateDynamicConnection()
		{
			throw new NotImplementedException();
		}

		private const string GetFullSpectraTableString = "SELECT * FROM Spectra ORDER BY Rt";
		private const int numberColumnsSpectraTable = 18; 
		private List<SpectraTableEntry> GetFullSpectraTable()
		{
			List<SpectraTableEntry> spectraList = new();
			using var command = new SQLiteCommand(_connection); 
			command.CommandText = GetFullSpectraTableString;

			using var sqliteReader = command.ExecuteReader();
			while (sqliteReader.Read())
			{
				SpectraTableEntry entry = new SpectraTableEntry();

				for (int i = 0; i < numberColumnsSpectraTable; i++)
				{
					if (!sqliteReader.IsDBNull(i))
					{
						continue; 
					}
					switch (i)
					{
						case 0:
							entry.Id = sqliteReader.GetInt32(i);
							break;
						case 1: 
							entry.Rt = sqliteReader.GetDouble(i);
							break;
						case 2:
							entry.Segment = sqliteReader.GetInt32(i);
							break;
						case 3:
							entry.AcquisitionKey = sqliteReader.GetInt32(i);
							break;
						case 4:
							entry.ParentSpectrum = sqliteReader.GetInt32(i);
							break;
						case 5:
							entry.MzAcqRangeLower = sqliteReader.GetInt32(i);
							break;
						case 6: 
							entry.MzAcqRangeUpper = sqliteReader.GetInt32(i);
							break;
						case 7: 
							entry.SumIntensity = sqliteReader.GetDouble(i);
							break;
						case 8:
							entry.MaxIntensity = sqliteReader.GetDouble(i);
							break;
						case 9: 
							entry.TransformatorId = sqliteReader.GetInt32(i);
							break;
						case 10: 
							entry.ProfileMzId = sqliteReader.GetInt32(i);
							break;
						case 11: 
							entry.ProfileIntensityId = sqliteReader.GetInt32(i);
							break;
						case 12: 
							entry.LineIndexId = sqliteReader.GetInt32(i);
							break;
						case 13: 
							entry.LineMzId = sqliteReader.GetInt32(i);
							break;
						case 14: 
							entry.LineIntensityId = sqliteReader.GetInt32(i);
							break;
						case 15: 
							entry.LineIndexWidthId = sqliteReader.GetInt32(i);
							break;
						case 16: 
							entry.LinePeakAreaId = sqliteReader.GetInt32(i);
							break;
						case 17: 
							entry.LineSnrId = sqliteReader.GetInt32(i);
							break; 
						default:
							continue; 
					}
					spectraList.Add(entry);
				}
			}

			return spectraList;
		}

		private void GetFullMetaDataTable()
		{

		}

		private void OpenFileConnection(string path)
		{
			string sqlite_fn = GetSQLiteCacheFilename(path);
			_handle = baf2sql_array_open_storage(1, ConvertStringToUTF8ByteArray(sqlite_fn));
			if (_handle == 0)
			{
				ThrowLastBaf2SqlError();
			}

			_connection = new SQLiteConnection();
			_connection.ConnectionString = "DataSource=" + sqlite_fn; 
			_connection.Open();
		}

		private void CloseFileConnection()
		{
			baf2sql_array_close_storage(_handle);
		}

		#region Bruker Dll Functions 
		[DllImport("baf2sql_c", CallingConvention = CallingConvention.Cdecl)]
		public static extern UInt32 baf2sql_get_sqlite_cache_filename
			  (byte[] sql_filename_buf_utf8, UInt32 sql_filename_buflen, byte[] baf_filename_utf8);

		[DllImport("baf2sql_c", CallingConvention = CallingConvention.Cdecl)]
		public static extern UInt64 baf2sql_array_open_storage
			   (int ignore_calibrator_ami, byte[] filename_utf8);

		[DllImport("baf2sql_c", CallingConvention = CallingConvention.Cdecl)]
		public static extern void baf2sql_array_close_storage
			   (UInt64 handle);

		[DllImport("baf2sql_c", CallingConvention = CallingConvention.Cdecl)]
		public static extern void baf2sql_array_get_num_elements
			   (UInt64 handle, UInt64 id, ref UInt64 num_elements);

		[DllImport("baf2sql_c", CallingConvention = CallingConvention.Cdecl)]
		public static extern int baf2sql_array_read_double
			   (UInt64 handle, UInt64 id, double[] buf);

		[DllImport("baf2sql_c", CallingConvention = CallingConvention.Cdecl)]
		public static extern int baf2sql_array_read_float
			   (UInt64 handle, UInt64 id, float[] buf);

		[DllImport("baf2sql_c", CallingConvention = CallingConvention.Cdecl)]
		public static extern int baf2sql_array_read_uint32
			   (UInt64 handle, UInt64 id, UInt32[] buf);

		[DllImport("baf2sql_c", CallingConvention = CallingConvention.Cdecl)]
		public static extern UInt32 baf2sql_get_last_error_string(StringBuilder buf, UInt32 len);

		[DllImport("baf2sql_c", CallingConvention = CallingConvention.Cdecl)]
		public static extern void baf2sql_set_num_threads(UInt32 n);

		/* ----------------------------------------------------------------------------------------------- */

		[Serializable()]
		public class Baf2SqlException : System.Exception
		{
			public Baf2SqlException() : base() { }
			public Baf2SqlException(string message) : base(message) { }
			public Baf2SqlException(string message, System.Exception inner) : base(message, inner) { }
			protected Baf2SqlException(System.Runtime.Serialization.SerializationInfo info,
				System.Runtime.Serialization.StreamingContext context)
			{ }
		}

		/* Throw last error string as an exception. */
		public static void ThrowLastBaf2SqlError()
		{
			StringBuilder buf = new StringBuilder("");
			UInt32 len = baf2sql_get_last_error_string(buf, 0);
			buf.EnsureCapacity((int)(len + 1));
			baf2sql_get_last_error_string(buf, len);
			throw new Baf2SqlException(buf.ToString());
		}

		/* ----------------------------------------------------------------------------------------------- */
		private static byte[] ConvertStringToUTF8ByteArray(String input)
		{
			byte[] utf8 = Encoding.UTF8.GetBytes(input);
			var result = new byte[utf8.Length + 1];
			Array.Copy(utf8, result, utf8.Length);

			return result;
		}
		/* Find out the file name of the SQL cache corresponding to the specified BAF file.
         * (If the SQL cache doesn't exist yet, it will be created.) */
		public static String GetSQLiteCacheFilename(String baf_filename)
		{
			byte[] buf = new byte[1];
			byte[] baf_filename_utf8 = ConvertStringToUTF8ByteArray(baf_filename);

			UInt32 len = baf2sql_get_sqlite_cache_filename(buf, 0, baf_filename_utf8);
			if (len == 0) ThrowLastBaf2SqlError();

			buf = new byte[len];
			len = baf2sql_get_sqlite_cache_filename(buf, len, baf_filename_utf8);
			if (len == 0) ThrowLastBaf2SqlError();

			return Encoding.UTF8.GetString(buf, 0, buf.Length - 1);
		}

		/* ----------------------------------------------------------------------------------------------- */

		/* Given the Id of one spectral component (e.g., a 'ProfileMzId' from the SQL cache),
         * load the binary data from the BAF (returning a double array). */
		public static double[] GetBafDoubleArray(UInt64 handle, UInt64 id)
		{
			UInt64 n = 0;
			baf2sql_array_get_num_elements(handle, id, ref n);

			double[] myArray = new double[n];
			int rc = baf2sql_array_read_double(handle, id, myArray);
			if (rc == 0) ThrowLastBaf2SqlError();

			return myArray;
		}

		/* Return array 'id', converting to float format */
		public static float[] GetBafFloatArray(UInt64 handle, UInt64 id)
		{
			UInt64 n = 0;
			baf2sql_array_get_num_elements(handle, id, ref n);

			float[] myArray = new float[n];
			int rc = baf2sql_array_read_float(handle, id, myArray);
			if (rc == 0) ThrowLastBaf2SqlError();

			return myArray;
		}

		/* Return array 'id', converting to UInt32 format */
		public static UInt32[] GetBafUInt32Array(UInt64 handle, UInt64 id)
		{
			UInt64 n = 0;
			baf2sql_array_get_num_elements(handle, id, ref n);

			UInt32[] myArray = new UInt32[n];
			int rc = baf2sql_array_read_uint32(handle, id, myArray);
			if (rc == 0) ThrowLastBaf2SqlError();

			return myArray;
		}
		#endregion
	}

	internal class SpectraTableEntry
	{
		public int Id { get; set; }
		public double Rt { get; set; }
		public int Segment { get; set; }
		public int AcquisitionKey { get; set; }
		public int? ParentSpectrum { get; set; }
		public int MzAcqRangeLower { get; set; }
		public int MzAcqRangeUpper { get; set; }
		public double SumIntensity { get; set; }
		public double MaxIntensity { get; set; }
		public int TransformatorId { get; set; }
		public int ProfileMzId { get; set; }

		public int ProfileIntensityId { get; set; }
		public int LineIndexId { get; set; }
		public int LineMzId { get; set; }
		public int LineIntensityId { get; set; }
		public int LineIndexWidthId { get; set; }
		public int LinePeakAreaId { get; set; }
		public int LineSnrId { get; set; }
	}
}
