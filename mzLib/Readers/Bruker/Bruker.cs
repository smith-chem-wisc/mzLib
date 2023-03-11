using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using System.Data.SQLite;
using Easy.Common.Extensions;
using MzLibUtil;
using UsefulProteomicsDatabases;

namespace Readers.Bruker
{
	public class Bruker : MsDataFile
	{
		private SQLiteConnection _connection;
		private ulong _handle; 
		public override MsDataFile LoadAllStaticData(FilteringParams filteringParams = null, int maxThreads = 1)
		{
			if (!File.Exists(FilePath))
			{
				throw new FileNotFoundException(); 
			}
			Loaders.LoadElements();

			List<MsDataScan> scans = new(); 
			OpenFileConnection(FilePath);
			int totalSpectra = GetTotalSpectraCount();
			LoadTablesInMemory();
			for (int i = 1; i <= totalSpectra; i++)
			{
				var tempScan = GetMsDataScan(i, filteringParams);
				scans.Add(tempScan); 
			}
			CloseFileConnection();
			Scans = scans.OrderBy(x => x.OneBasedPrecursorScanNumber).ToArray();
			SourceFile = GetSourceFile();
			return this; 
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

		// ids are the one based scan numbers
		private List<AcqKeyTable> _acqKeyTable = new List<AcqKeyTable>();
		private List<SpectraTableEntry> _spectraTable = new List<SpectraTableEntry>();
		private List<StepsRow> _stepsTable = new List<StepsRow>();

		private const string GetTotalSpectraCountString =
			@"SELECT COUNT(*) FROM Spectra"; 
		private int GetTotalSpectraCount()
		{
			using var command = new SQLiteCommand(_connection); 
			command.CommandText = GetTotalSpectraCountString;
			using var sqliteReader = command.ExecuteReader();
			int count = 0; 
			while (sqliteReader.Read())
			{
				count = sqliteReader.GetInt32(0);
				break; 
			}

			return count; 
		}

		private void LoadTablesInMemory()
		{
			_acqKeyTable = GetAcqKeyTable();
			_spectraTable = GetFullSpectraTable(); 
			_stepsTable = GetFullStepsTable();
		}

		private MsDataScan GetMsDataScan(int id, IFilteringParams? filteringParams)
		{
			// if the id fields are null, then it doesn't exist.
			MzSpectrum spectrumData = GetSpectraData(_spectraTable[id], filteringParams);
			int? oneBasedPrecursorScanNumber =
				_spectraTable[id].ParentSpectrum is null ? null : _spectraTable[id].ParentSpectrum;
			// need the plus one because Bruker codes MS scans as 0, but we code them as 1. 
			int msOrder = _acqKeyTable[id].MsLevel + 1; 
			bool isCentroid = _spectraTable[id].LineIndexId is not null;  
			var polarity = _acqKeyTable[id].Polarity == 0 ? Polarity.Positive : Polarity.Negative;
			double? selectedIonMz = null;

			MZAnalyzerType analyzer = _acqKeyTable[id].AcquisitionMode == 33 ? MZAnalyzerType.FTICR : MZAnalyzerType.TOF;  

			DissociationType? dissociationType = null; 
			if (msOrder > 1)
			{
				selectedIonMz = _stepsTable.First(i => i.TargetSpectrum == id).ParentMass;
				switch (_acqKeyTable[id].ScanMode)
				{
					case 0:
						dissociationType = null;
						break;
					case 2:
						dissociationType = DissociationType.CID;
						break;
					case 4:
						dissociationType = DissociationType.ISCID;
						break;
					case 255:
						dissociationType = DissociationType.Unknown;
						break; 
					default:
						dissociationType = DissociationType.Unknown;
						break; 
				}
			}

			// need to get the isolation width from the variables table, but the variables table doesn't 
			// have an easily accessible name. I guess I don't need to but it is kind of nice. 
			MzRange scanWindowRange = new MzRange((double)_spectraTable[id].MzAcqRangeLower,
				(double)_spectraTable[id].MzAcqRangeUpper);
			string scanFilter = "f";
			string nativeId = "scan=" + id; 

			return new MsDataScan(spectrumData,
				oneBasedScanNumber: id, msnOrder: msOrder, isCentroid: isCentroid,
				polarity: polarity, retentionTime: _spectraTable[id].Rt,
				scanWindowRange: scanWindowRange, scanFilter: scanFilter,
				mzAnalyzer: analyzer, totalIonCurrent: _spectraTable[id].SumIntensity,
				injectionTime: null, noiseData: new double[,] { { 0 }, { 0 } },
				nativeId: nativeId, selectedIonMz: selectedIonMz, dissociationType:dissociationType,
				oneBasedPrecursorScanNumber:oneBasedPrecursorScanNumber, 
				isolationMZ:selectedIonMz);
		}

		private MzSpectrum GetSpectraData(SpectraTableEntry spectraInfo, IFilteringParams filteringParams)
		{
			// get centroided data if available. Otherwise, default to profile.
			// if neither exists, return a spectra consisting of two zeroes. This probably shouldn't 
			// happen, but can't rule it out. 
			if (spectraInfo.LineMzId == null && spectraInfo.ProfileMzId == null)
			{
				return new MzSpectrum(new double[,] {{0}, {0}}); 
			}
			if (spectraInfo.LineMzId == null)
			{
				double[] profileMzs = GetBafDoubleArray(_handle, (ulong)spectraInfo.ProfileMzId);
				double[] profileInts = GetBafDoubleArray(_handle, (ulong)spectraInfo.ProfileIntensityId);
				if (filteringParams != null
				    && profileMzs.Length > 0
				    && filteringParams.ApplyTrimmingToMsMs)
				{
					WindowModeHelper.Run(ref profileInts,
						ref profileMzs, filteringParams, 
						profileMzs[0], profileMzs[^0]);
				}

				return new MzSpectrum(profileMzs, profileInts, true); 
			}

			double[] lineMzs = GetBafDoubleArray(_handle, (ulong)spectraInfo.LineMzId);
			double[] lineInt = GetBafDoubleArray(_handle, (ulong)spectraInfo.LineIntensityId);
			if (filteringParams != null
			    && lineMzs.Length > 0
			    && filteringParams.ApplyTrimmingToMsMs)
			{
				WindowModeHelper.Run(ref lineInt,
					ref lineMzs, filteringParams,
					lineMzs[0], lineMzs[^0]);
			}
			return new MzSpectrum(lineMzs, lineInt, true); 
		}

		private const string GetAcqKeyTableString = @"SELECT * FROM AcquisitionKeys";
		private const int AcqKeyTableColumns = 5; 
		private List<AcqKeyTable> GetAcqKeyTable()
		{
			List<AcqKeyTable> spectraList = new();
			using var command = new SQLiteCommand(_connection);
			command.CommandText = GetAcqKeyTableString;

			using var sqliteReader = command.ExecuteReader();
			while (sqliteReader.Read())
			{
				AcqKeyTable table = new AcqKeyTable();
				for (int i = 0; i < AcqKeyTableColumns; i++)
				{
					if (sqliteReader.IsDBNull(i))
					{
						continue; 
					}

					switch (i)
					{
						case 0: 
							table.Id = sqliteReader.GetInt32(i);
							break;
						case 1: 
							table.Polarity = sqliteReader.GetInt32(i);
							break;
						case 2: 
							table.ScanMode = sqliteReader.GetInt32(i);
							break;
						case 3: 
							table.AcquisitionMode = sqliteReader.GetInt32(i);
							break;
						case 4: 
							table.MsLevel = sqliteReader.GetInt32(i);
							break; 
					}
				}
				spectraList.Add(table);
			}
			return spectraList;
		}

		private const string GetStepsTableString = @"SELECT * FROM Steps";
		private const int StepsTableColumns = 6;

		private List<StepsRow> GetFullStepsTable()
		{
			List<StepsRow> stepsTableList = new();
			using var command = new SQLiteCommand(_connection);
			command.CommandText = GetStepsTableString; 
			using var sqliteReader = command.ExecuteReader();
			while (sqliteReader.Read())
			{
				StepsRow stepsRow = new StepsRow();
				var propertiesArray = stepsRow.GetPropertyNames();
				stepsTableList.Add(SqlColumnReader<StepsRow>(sqliteReader));
			}
			return stepsTableList;
		}

		private T SqlColumnReader<T>(SQLiteDataReader reader) where T: new()
		{
			// get all the property names, then iterate over that. 
			// The objects should be exact 1:1 column corresponding so as 
			// long as that's guaranteed, then it should work no problem. 
			var propertiesArray = typeof(T).GetPropertyNames();
			if (propertiesArray == null)
			{
				throw new ArgumentNullException(); 
			}

			T newSqlMappedCsObject = new(); 
			for (int i = 0; i < propertiesArray.Length; i++)
			{
				// check type of property about to be returned
				// then use the correct method from SqliteDataReader to parse 
				// the value. 
				Type propertyType = typeof(T).GetProperty(propertiesArray[i])!.PropertyType;
				if (propertyType == typeof(int))
				{
					typeof(T).GetProperty(propertiesArray[i])!
						.SetValue(newSqlMappedCsObject, reader.GetInt32(i));
				}else if (propertyType == typeof(double))
				{
					typeof(T).GetProperty(propertiesArray[i])!
						.SetValue(newSqlMappedCsObject, reader.GetDouble(i));
				}
				else if (propertyType == typeof(string))
				{
					typeof(T).GetProperty(propertiesArray[i])!
						.SetValue(newSqlMappedCsObject, reader.GetString(i));
				}
			}

			return newSqlMappedCsObject;
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
		public static extern void baf2sql_array_close_storage(UInt64 handle);

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

	internal class StepsRow
	{
		public int TargetSpectrum { get; set; }
		public int Number { get; set; }
		public int IsolationType { get; set; }

		public int ReactionType { get; set; }
		// 0 = MS; 1 = MS/MS
		public int MsLevel { get; set; }
		public double ParentMass { get; set; }
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
		public int? TransformatorId { get; set; }
		public int? ProfileMzId { get; set; }

		public int? ProfileIntensityId { get; set; }
		public int? LineIndexId { get; set; }
		public int? LineMzId { get; set; }
		public int? LineIntensityId { get; set; }
		public int? LineIndexWidthId { get; set; }
		public int? LinePeakAreaId { get; set; }
		public int? LineSnrId { get; set; }
	}

	public class AcqKeyTable
	{
		public int Id { get; set; }
		// 0 = positive, 1 = negative 
		public int Polarity { get; set; }
		// 0 = MS
		// 2 = MS/MS
		// 4 = in-source CID
		// 5 = broadband CID
		// 255 = unknown
		public int ScanMode { get; set; }
		// 1 = (axial or orthogonal) TOF, linear detection mode
		// 2 = (axial or orthogonal) TOF, reflector detection mode
		// 33 = FTMS
		// 255 = unknown
		public int AcquisitionMode { get; set; }
		// 0 = MS (no fragmentation),
		// 1 = MS/MS
		public int MsLevel { get; set; }
	}
}
