using System.Globalization;
using System.Runtime.InteropServices;
using System.Text;
using MassSpectrometry;
using System.Data.SQLite;
using Easy.Common.Extensions;
using MzLibUtil;

namespace Readers
{
	/// <summary>
	/// Reads Bruker data files (but not TimsTOF data). The major quirk of the
	/// Burker file format is that it relies on an SQL database in order to retrieve data.
	/// Therefore, I had to write a bunch of fairly simple SQL statements in order to access the data.
	/// These are stored in const strings right above the methods they are used in. Please don't move
	/// or change them without understanding SQL. -AVC
	/// </summary>
	public class BrukerFileReader : MsDataFile
	{
		public BrukerFileReader(string filePath) : base(filePath) { }

		/* There are three tables that contain all of the information required to create an MsDataScan object: Steps (which has isolation information),
		 Spectra (contains the spectrum metadata), and AcqKeys (contains more metadata). Unfortunately, the metadata is spread over the three tables, 
		making it rather complicated to cleanly return an MsDataScan. LoadAllStatic data loads all the scans into memory using a single SQL query. 
		LoadDynamic uses an sqlite statement to get the requested scan each time it is requested. 
		 */ 
		/// <summary>
		/// Connection to sqlite database that is initialized when this object is created.
		/// </summary>
		private SQLiteConnection? _connection;
		/// <summary>
		/// This is an identification for the requested sqlite database. This is very important. Don't touch it. 
		/// </summary>
		private ulong? _handle;
		
		// tables used for LoadAllStaticData.
		// These tables reflect the sqlite tables found in the Bruker file support.
		private List<AcqKeyRow> _acqKeyTable = new List<AcqKeyRow>();
		private List<SpectraTableRow> _spectraTable = new List<SpectraTableRow>();
		private List<StepsRow> _stepsTable = new List<StepsRow>();
		// Per-spectrum precursor metadata (isolation width, charge, collision energy) keyed by spectrum Id.
		// Bruker stores these in the (Per)SpectrumVariables table, NOT in Steps/Spectra/AcquisitionKeys.
		private Dictionary<int, PrecursorVariables> _precursorVariables = new();

		// PermanentName values in the SupportedVariables table for the precursor metadata we need.
		// The isolation width is REQUIRED downstream: MsDataScan.IsolationRange (and therefore
		// GetIsolatedMassesAndCharges) returns null/empty unless both IsolationMz and IsolationWidth are set,
		// which means a consumer such as MetaMorpheus cannot deconvolute a precursor and drops every MS2 scan.
		private const string IsolationWidthVariableName = "Quadrupole_IsolationResolution_Act";
		private const string PrecursorChargeVariableName = "MSMS_PreCursorChargeState";
		private const string CollisionEnergyVariableName = "Collision_Energy_Act";

		public override MsDataFile LoadAllStaticData(FilteringParams? filteringParams = null, int maxThreads = 1)
		{
			// special notes regarding the Bruker file import: 
			/*
			 * 1) Database file connection is opened.
			 * 2) The tables from the bruker sqlite databases are loaded into C# in memory.
			 * 3) The tables are parsed. 
			 * 4) File connection is closed.
			 * 5) Scans are ordered and returned. 
			 */ 

			if (!Directory.Exists(FilePath))
			{
				throw new FileNotFoundException(); 
			}

			List<MsDataScan> scans = new(); 
			OpenFileConnection(FilePath+@"\analysis.baf");
			int totalSpectra = GetTotalSpectraCount();
			LoadTablesInMemory();
			for (int i = 0; i < totalSpectra; i++)
			{
				// check to see if there is a matching element in the steps tables 
                _precursorVariables.TryGetValue(_spectraTable[i].Id, out var precursorVariables);
                if (!_spectraTable[i].Parent.HasValue)
                {
                    // get the index to the row in steps table that corresponds to the id in spectraTable
                    var tempScan = GetMsDataScan(_spectraTable[i], _acqKeyTable, null, filteringParams, precursorVariables);
                    scans.Add(tempScan);
                    continue;
                }

                var targetStepsRow = _stepsTable.Find(z => z.TargetSpectrum == _spectraTable[i].Id);
                var scan = GetMsDataScan(_spectraTable[i], _acqKeyTable, targetStepsRow, filteringParams, precursorVariables);
                scans.Add(scan);

            }
			// close the file connection. At this point, you don't need to be connected to the sqlite database anymore. You have all the data 
			// you need. 
			CloseFileConnection();
            Scans = scans.OrderBy(x => x.OneBasedScanNumber).ToArray();
			SourceFile = GetSourceFile();
            return this; 
		}

		private const string nativeIdFormat = "scan number only nativeID format";
		private const string massSpecFileFormat = "mzML format"; 
		public override SourceFile GetSourceFile()
        {
			// append the analysis.baf because the constructor for SourceFile will look for the 
			// parent directory. 
            string fileName = FilePath + @"\analysis.baf"; 
            return new SourceFile(nativeIdFormat, massSpecFileFormat,
				null, null, id: null, filePath: fileName);
        }

        public override MsDataScan GetOneBasedScanFromDynamicConnection(int oneBasedScanNumber, IFilteringParams? filterParams = null)
		{
            if (CheckIfScansLoaded() && oneBasedScanNumber <= Scans.Length)
                return GetOneBasedScan(oneBasedScanNumber);

            lock (DynamicReadingLock)
            {
                return GetMsDataScanDynamic(oneBasedScanNumber, filterParams);
            }
        }

		public override void CloseDynamicConnection()
		{
			CloseFileConnection();
		}

		public override void InitiateDynamicConnection()
		{
            if (!File.Exists(FilePath + @"\analysis.baf"))
            {
                throw new FileNotFoundException();
            }
            OpenFileConnection(FilePath + @"\analysis.baf");
		}

		private const string GetFullSpectraTableString = "SELECT * FROM Spectra ORDER BY Rt";
		private List<SpectraTableRow> GetFullSpectraTable()
		{
			List<SpectraTableRow> spectraList = new();
			using var command = new SQLiteCommand(_connection); 
			command.CommandText = GetFullSpectraTableString;

			using var sqliteReader = command.ExecuteReader();
			while (sqliteReader.Read())
			{
				SpectraTableRow row = SqlColumnReader<SpectraTableRow>(sqliteReader);
				spectraList.Add(row);
			}
            return spectraList;
		}
		private const string GetTotalSpectraCountString =
			@"SELECT COUNT(*) FROM Spectra"; 
		// Total spectra count required for determining iteration values in LoadAllStaticData. 
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
			_precursorVariables = GetAllPrecursorVariables();
		}

		// Sqlite commands to fetch tables and load them into C#. 
		private const string GetSingleSpectrumString = @"SELECT * FROM Spectra " +
		                                               "WHERE Id = ";

		private const string GetSingleAcqKeys = @"SELECT * FROM AcquisitionKeys";

		private const string GetSingleStepsKey = @"SELECT * FROM Steps " +
		                                         "WHERE TargetSpectrum = "; 
		// Executes as single sqlite queries. 
		private MsDataScan GetMsDataScanDynamic(int id, IFilteringParams? filteringParams)
		{
			// use a SQL statement with a variable to filter the tables
			string spectrumQuery = new StringBuilder().Append(GetSingleSpectrumString).Append(id).ToString();
			string acqKeyQuery = new StringBuilder().Append(GetSingleAcqKeys).Append(id).ToString();
			string stepsQuery = new StringBuilder().Append(GetSingleStepsKey).Append(id).ToString();

			string[] queryStringsArray = new[]
			{
				spectrumQuery,
                stepsQuery,
				acqKeyQuery
            }; 

			SpectraTableRow? spectraTableRow = new SpectraTableRow();
			List<AcqKeyRow> acqKey = new List<AcqKeyRow>();
			StepsRow? stepsRow = new StepsRow();

			// only need to get the acquisition key list a single time. 
            var acqKeyList = GetAcqKeyTable(); 

			// load in the rows that correspond to the one based spectrum count that you want.
			for (int i = 0; i < 2; i++)
			{
				using var sqliteCommand = new SQLiteCommand(_connection);
				sqliteCommand.CommandText = queryStringsArray[i]; 

				using var sqliteReader = sqliteCommand.ExecuteReader();
				// there should only be one scan with a given id. 
				 
				while (sqliteReader.Read())
				{
					if (sqliteReader == null)
					{
						break; 
					}
					switch (i)
					{
						case 0:
							spectraTableRow = SqlColumnReader<SpectraTableRow>(sqliteReader); 
							break;
                        case 1:
							stepsRow = SqlColumnReader<StepsRow>(sqliteReader);
							break;
					}
				}
			}

			var precursorVariables = GetPrecursorVariablesForSpectrum(id);
			return GetMsDataScan(spectraTableRow, acqKeyList, stepsRow, filteringParams, precursorVariables);
		}

		/// <summary>
		/// Converts relevant information from the three main sqlite datatables into
		/// an MsDataScan object. 
		/// </summary>
		/// <param name="spectraRow"></param>
		/// <param name="acqKeyRow"></param>
		/// <param name="stepsRow"></param>
		/// <param name="filterParams"></param>
		/// <returns></returns>
		private MsDataScan GetMsDataScan(SpectraTableRow spectraRow,
			List<AcqKeyRow> acqKeyRow, StepsRow? stepsRow, IFilteringParams? filterParams,
			PrecursorVariables? precursorVariables = null)
        {

            var acquisitionKeyMsLevelLink = acqKeyRow
                .ToDictionary(i => i.Id, i => (i.MsLevel, i.Polarity));

			int? oneBasedPrecursorScanNumber = spectraRow.Parent is null ? null : spectraRow.Parent;
			// need the plus one because Bruker codes MS scans as 0, but we code them as 1.
			int msOrder = acquisitionKeyMsLevelLink[spectraRow.AcquisitionKey].MsLevel + 1;

            // Bruker data can record profile, centroid or both. So GetSpectraData return an integer
			// value indicating the scenario. msOrder is needed so peak trimming honors the per-MS-level
			// FilteringParams flags (trimming MS1 would strip precursor isotopes and break deconvolution).
			int isCentroidIntSwitch = GetSpectraData(spectraRow, filterParams, msOrder, out MzSpectrum spectrumData);
            bool isCentroid = false;
            switch (isCentroidIntSwitch)
            {
				case -1:
                    isCentroid = false;
                    break;
				case 0:
                    isCentroid = false;
                    break;
				case 1:
                    isCentroid = true;
                    break;
            }

			var polarity = acquisitionKeyMsLevelLink[spectraRow.AcquisitionKey].Polarity == 0 ? Polarity.Positive : Polarity.Negative;
			double? selectedIonMz = null;

			MZAnalyzerType analyzer = acqKeyRow[spectraRow.AcquisitionKey - 1].AcquisitionMode == 33 ? MZAnalyzerType.FTICR : MZAnalyzerType.TOF;


			// Isolation width, precursor charge and collision energy live in the (Per)SpectrumVariables table.
			// The isolation width is essential: without it MsDataScan.IsolationRange is null, so precursor
			// deconvolution yields nothing and downstream search engines drop every MS2 scan.
			double? isolationWidth = null;
			int? selectedIonChargeStateGuess = null;
			string? hcdEnergy = null;
			if (msOrder > 1)
			{
                if (spectraRow.Parent != null && stepsRow != null)
                {
                    selectedIonMz = stepsRow.ParentMass;
                }

                if (precursorVariables != null)
                {
                    // A zero/non-positive width is Bruker's "not set" sentinel (seen on MS1 rows); ignore it.
                    if (precursorVariables.IsolationWidth is > 0)
                        isolationWidth = precursorVariables.IsolationWidth;
                    // Charge state is 0 when the instrument could not determine it; treat that as unknown.
                    if (precursorVariables.PrecursorChargeState is > 0)
                        selectedIonChargeStateGuess = precursorVariables.PrecursorChargeState;
                    if (precursorVariables.CollisionEnergy.HasValue)
                        hcdEnergy = precursorVariables.CollisionEnergy.Value.ToString(CultureInfo.InvariantCulture);
                }
            }
            // ScanMode (0 = MS, 2 = MS/MS, 4 = in-source CID, 255 = unknown) is what the dissociation map is keyed on;
            // AcquisitionMode (1/2 = TOF, 33 = FTMS) is the analyzer and was previously used here by mistake.
            int scanMode = acqKeyRow[spectraRow.AcquisitionKey - 1].ScanMode;
            DissociationType dissociationType = _dissociationTypes.TryGetValue(scanMode, out var dt)
                ? dt
                : DissociationType.Unknown;

            MzRange scanWindowRange = new MzRange((double)spectraRow.MzAcqRangeLower,
				(double)spectraRow.MzAcqRangeUpper);
			string scanFilter = "f";
			string nativeId = "scan=" + spectraRow.Id;

			// Bruker records retention time in seconds; mzLib's convention (mzML/Thermo) is minutes.
			return new MsDataScan(spectrumData,
				oneBasedScanNumber: spectraRow.Id, msnOrder: msOrder, isCentroid: isCentroid,
				polarity: polarity, retentionTime: spectraRow.Rt / 60.0,
				scanWindowRange: scanWindowRange, scanFilter: scanFilter,
				mzAnalyzer: analyzer, totalIonCurrent: spectraRow.SumIntensity,
				injectionTime: null, noiseData: new double[,] { { 0 }, { 0 } },
				nativeId: nativeId, selectedIonMz: selectedIonMz,
				selectedIonChargeStateGuess: selectedIonChargeStateGuess,
				dissociationType: dissociationType,
				oneBasedPrecursorScanNumber: oneBasedPrecursorScanNumber,
				isolationMZ: selectedIonMz, isolationWidth: isolationWidth,
				hcdEnergy: hcdEnergy);
		}

		// maps the integer codes in the Bruker tables to the enum in mzlib. 
		private static Dictionary<int, DissociationType> _dissociationTypes = new()
        {
            {0, DissociationType.Unknown},
            {2, DissociationType.CID},
            {4, DissociationType.ISCID},
            {5, DissociationType.CID},   // broadband CID
            {255, DissociationType.Unknown}
        };
		// Converst SQLite to C# parseable objects. Returns 0 upon success.
		private int GetSpectraData(SpectraTableRow spectraInfo, IFilteringParams? filteringParams, int msOrder, out MzSpectrum spectrum)
		{
			// get centroided data if available. Otherwise, default to profile.
			// if neither exists, return a spectra consisting of two zeroes. This probably shouldn't
			// happen, but can't rule it out.
			if (spectraInfo.LineMzId == null && spectraInfo.ProfileMzId == null)
			{
				spectrum = new MzSpectrum(new double[,] {{0}, {0}});
                return -1;
            }
			if (spectraInfo.LineMzId == null && spectraInfo.ProfileMzId != null)
			{
				double[] profileMzs = GetBafDoubleArray(_handle!.Value, (ulong)spectraInfo.ProfileMzId);
				double[] profileInts = GetBafDoubleArray(_handle!.Value, (ulong)spectraInfo.ProfileIntensityId);
				if (filteringParams != null
				    && profileMzs.Length > 1
				    && ShouldTrim(filteringParams, msOrder))
				{
					WindowModeHelper.Run(ref profileInts,
						ref profileMzs, filteringParams,
						profileMzs[0], profileMzs[^1]);
				}

				spectrum = new MzSpectrum(profileMzs, profileInts, true);
                return 0;
            }

			double[] lineMzs = GetBafDoubleArray(_handle!.Value, (ulong)spectraInfo.LineMzId!);
			double[] lineInt = GetBafDoubleArray(_handle!.Value, (ulong)spectraInfo.LineIntensityId!);
			if (filteringParams != null
			    && lineMzs.Length > 1
			    && ShouldTrim(filteringParams, msOrder))
			{
				WindowModeHelper.Run(ref lineInt,
					ref lineMzs, filteringParams,
					lineMzs[0], lineMzs[^1]);
			}
			spectrum = new MzSpectrum(lineMzs, lineInt, true);
            return 1;
        }

		// Peak trimming must honor the per-MS-level FilteringParams flags, exactly as the mzML/Thermo/timsTOF
		// readers do. In particular MS1 scans must only be trimmed when ApplyTrimmingToMs1 is set; trimming an
		// MS1 precursor scan strips its isotope envelopes and precursor deconvolution then finds nothing.
		private static bool ShouldTrim(IFilteringParams filteringParams, int msOrder)
		{
			return (filteringParams.ApplyTrimmingToMs1 && msOrder == 1)
			       || (filteringParams.ApplyTrimmingToMsMs && msOrder == 2)
			       || (filteringParams.ApplyTrimmingToMsN && msOrder > 2);
		}

		private const string GetAcqKeyTableString = @"SELECT * FROM AcquisitionKeys";
        private List<AcqKeyRow> GetAcqKeyTable()
		{
			List<AcqKeyRow> spectraList = new();
			using var command = new SQLiteCommand(_connection);
			command.CommandText = GetAcqKeyTableString;

			using var sqliteReader = command.ExecuteReader();
			
            while (sqliteReader.Read())
            {
                AcqKeyRow row = SqlColumnReader<AcqKeyRow>(sqliteReader);
                spectraList.Add(row);
            }
            
            return spectraList;
		}

		private const string GetStepsTableString = @"SELECT * FROM Steps";
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

		/// <summary>
		/// Holds the per-spectrum precursor metadata that Bruker stores in the (Per)SpectrumVariables table
		/// (rather than in Spectra/Steps/AcquisitionKeys). The isolation width in particular is required for
		/// precursor deconvolution: <see cref="MsDataScan.IsolationRange"/> is null unless it is set.
		/// </summary>
		private class PrecursorVariables
		{
			public double? IsolationWidth { get; set; }
			public int? PrecursorChargeState { get; set; }
			public double? CollisionEnergy { get; set; }
		}

		/// <summary>
		/// Resolves the integer Variable ids for the precursor metadata we care about by their stable
		/// PermanentName, since the numeric ids are not guaranteed to be identical across files/instruments.
		/// </summary>
		private Dictionary<string, int> GetVariableIdsByName()
		{
			var map = new Dictionary<string, int>();
			try
			{
				using var command = new SQLiteCommand(_connection);
				command.CommandText = "SELECT Variable, PermanentName FROM SupportedVariables";
				using var reader = command.ExecuteReader();
				while (reader.Read())
				{
					if (reader.IsDBNull(1)) continue;
					map[reader.GetString(1)] = reader.GetInt32(0);
				}
			}
			catch (SQLiteException)
			{
				// Table absent (older/unusual acquisition) - callers fall back to no precursor variables.
			}
			return map;
		}

		/// <summary>
		/// Loads isolation width, precursor charge and collision energy for every spectrum from the
		/// PerSpectrumVariables table, keyed by spectrum Id. Returns an empty dictionary if the table
		/// (or the required variables) are not present.
		/// </summary>
		private Dictionary<int, PrecursorVariables> GetAllPrecursorVariables()
		{
			var result = new Dictionary<int, PrecursorVariables>();
			var variableIds = GetVariableIdsByName();
			if (!variableIds.TryGetValue(IsolationWidthVariableName, out int widthId))
				return result; // without an isolation width there is nothing worth loading

			variableIds.TryGetValue(PrecursorChargeVariableName, out int chargeId);
			variableIds.TryGetValue(CollisionEnergyVariableName, out int collisionEnergyId);

			try
			{
				using var command = new SQLiteCommand(_connection);
				command.CommandText =
					$"SELECT Spectrum, Variable, Value FROM PerSpectrumVariables " +
					$"WHERE Variable IN ({widthId}, {chargeId}, {collisionEnergyId})";
				using var reader = command.ExecuteReader();
				while (reader.Read())
				{
					if (reader.IsDBNull(2)) continue;
					int spectrumId = reader.GetInt32(0);
					int variableId = reader.GetInt32(1);
					if (!result.TryGetValue(spectrumId, out var pv))
					{
						pv = new PrecursorVariables();
						result[spectrumId] = pv;
					}
					AssignVariable(pv, variableId, reader.GetValue(2), widthId, chargeId, collisionEnergyId);
				}
			}
			catch (SQLiteException)
			{
				// PerSpectrumVariables absent - leave the dictionary empty.
			}
			return result;
		}

		/// <summary>
		/// Single-spectrum equivalent of <see cref="GetAllPrecursorVariables"/> for the dynamic-connection path.
		/// </summary>
		private PrecursorVariables? GetPrecursorVariablesForSpectrum(int spectrumId)
		{
			var variableIds = GetVariableIdsByName();
			if (!variableIds.TryGetValue(IsolationWidthVariableName, out int widthId))
				return null;

			variableIds.TryGetValue(PrecursorChargeVariableName, out int chargeId);
			variableIds.TryGetValue(CollisionEnergyVariableName, out int collisionEnergyId);

			PrecursorVariables? pv = null;
			try
			{
				using var command = new SQLiteCommand(_connection);
				command.CommandText =
					$"SELECT Variable, Value FROM PerSpectrumVariables " +
					$"WHERE Spectrum = {spectrumId} AND Variable IN ({widthId}, {chargeId}, {collisionEnergyId})";
				using var reader = command.ExecuteReader();
				while (reader.Read())
				{
					if (reader.IsDBNull(1)) continue;
					pv ??= new PrecursorVariables();
					AssignVariable(pv, reader.GetInt32(0), reader.GetValue(1), widthId, chargeId, collisionEnergyId);
				}
			}
			catch (SQLiteException)
			{
				// PerSpectrumVariables absent.
			}
			return pv;
		}

		// The Value column has no declared affinity and may come back as a boxed long or double; convert defensively.
		private static void AssignVariable(PrecursorVariables pv, int variableId, object rawValue,
			int widthId, int chargeId, int collisionEnergyId)
		{
			double value;
			try
			{
				value = Convert.ToDouble(rawValue, CultureInfo.InvariantCulture);
			}
			catch (Exception e) when (e is FormatException or InvalidCastException or OverflowException)
			{
				return; // e.g. the byte[] exclusion-mass blob; nothing useful for us
			}

			if (variableId == widthId) pv.IsolationWidth = value;
			else if (variableId == chargeId) pv.PrecursorChargeState = (int)value;
			else if (variableId == collisionEnergyId) pv.CollisionEnergy = value;
		}
		/// <summary>
		/// Generic function to read an sql table. Currently supports the types that exist in
		/// the sql tables, which are int32, string, and double. 
		/// </summary>
		/// <typeparam name="T">Object ("receiving object") of a type that has a 1:1 col to property mapping. The order in the receiving object
		/// must match the order of the columns in the sql database. 
		/// </typeparam>
		/// <param name="reader">SQLiteReader object, initialized after the execution of a command.</param>
		/// <returns>Return null exception if there is an error in the data format of the baf file.</returns>
		/// <exception cref="ArgumentNullException"></exception>
		public static T SqlColumnReader<T>(SQLiteDataReader reader) where T: new()
		{
			// get all the property names, then iterate over that. 
			// The objects should be exact 1:1 column corresponding so as 
			// long as that's guaranteed, then it should work no problem. 
			var propertiesArray = typeof(T).GetProperties()
                .Select(i => i.Name)
                .ToArray();
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
				}else if (propertyType == typeof(long?))
                {
                    long? value = null;
                    if (reader.IsDBNull(i))
                    {
                        typeof(T).GetProperty(propertiesArray[i])!
                            .SetValue(newSqlMappedCsObject, value);
                        continue; 
                    }
                    typeof(T).GetProperty(propertiesArray[i])!
                        .SetValue(newSqlMappedCsObject, reader.GetInt64(i));
                }else if (propertyType == typeof(int?))
                {
                    int? value = null;
                    if (reader.IsDBNull(i))
                    {
						typeof(T).GetProperty(propertiesArray[i])!
                            .SetValue(newSqlMappedCsObject, value);
						continue;
                    }
					typeof(T).GetProperty(propertiesArray[i])!
                        .SetValue(newSqlMappedCsObject, reader.GetInt32(i));
                }
			}

			return newSqlMappedCsObject;
		}
		private void OpenFileConnection(string path)
		{
			string sqlite_fn = GetSQLiteCacheFilename(path);
			_handle = baf2sql_array_open_storage(1, ConvertStringToUTF8ByteArray(path));
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
			baf2sql_array_close_storage(_handle!.Value);
		}

		#region Bruker Dll Functions 
		[DllImport("baf2sql_c.dll", CharSet = CharSet.Unicode,
            CallingConvention = CallingConvention.Cdecl)]
		private static extern UInt32 baf2sql_get_sqlite_cache_filename
			  (byte[] sql_filename_buf_utf8, UInt32 sql_filename_buflen, byte[] baf_filename_utf8);

		[DllImport("baf2sql_c.dll", CallingConvention = CallingConvention.Cdecl)]
		private static extern UInt64 baf2sql_array_open_storage
			   (int ignore_calibrator_ami, byte[] filename_utf8);

		[DllImport("baf2sql_c.dll", CallingConvention = CallingConvention.Cdecl)]
		// bruker doesn't actually provide a way to determine if the sqlite database is closed correctly. 
		private static extern void baf2sql_array_close_storage(UInt64 handle);

		[DllImport("baf2sql_c", CallingConvention = CallingConvention.Cdecl)]
		private static extern void baf2sql_array_get_num_elements
			(UInt64 handle, UInt64 id, ref UInt64 num_elements);

		[DllImport("baf2sql_c", CallingConvention = CallingConvention.Cdecl)]
		private static extern int baf2sql_array_read_double
			   (UInt64 handle, UInt64 id, double[] buf);

		[DllImport("baf2sql_c", CallingConvention = CallingConvention.Cdecl)]
		private static extern int baf2sql_array_read_float
			   (UInt64 handle, UInt64 id, float[] buf);

		[DllImport("baf2sql_c", CallingConvention = CallingConvention.Cdecl)]
		private static extern int baf2sql_array_read_uint32
			   (UInt64 handle, UInt64 id, UInt32[] buf);

		[DllImport("baf2sql_c", CallingConvention = CallingConvention.Cdecl)]
		private static extern UInt32 baf2sql_get_last_error_string(StringBuilder buf, UInt32 len);

		[DllImport("baf2sql_c", CallingConvention = CallingConvention.Cdecl)]
		private static extern void baf2sql_set_num_threads(UInt32 n);

		/* ----------------------------------------------------------------------------------------------- */

		[Serializable()]
		private class Baf2SqlException : System.Exception
		{
			public Baf2SqlException() : base() { }
			public Baf2SqlException(string message) : base(message) { }
			public Baf2SqlException(string message, System.Exception inner) : base(message, inner) { }
			protected Baf2SqlException(System.Runtime.Serialization.SerializationInfo info,
				System.Runtime.Serialization.StreamingContext context)
			{ }
		}

		/* Throw last error string as an exception. */
		private static void ThrowLastBaf2SqlError()
		{
			StringBuilder buf = new StringBuilder("");
			UInt32 len = baf2sql_get_last_error_string(buf, 0);
			buf.EnsureCapacity((int)(len + 1));
			baf2sql_get_last_error_string(buf, len);
			throw new Baf2SqlException(buf.ToString());
		}

		/* ----------------------------------------------------------------------------------------------- */
		public static byte[] ConvertStringToUTF8ByteArray(String input)
		{
			byte[] utf8 = Encoding.UTF8.GetBytes(input);
			var result = new byte[utf8.Length + 1];
			Array.Copy(utf8, result, utf8.Length);

			return result;
		}
		/* Find out the file name of the SQL cache corresponding to the specified BAF file.
         * (If the SQL cache doesn't exist yet, it will be created.) */
		private static String GetSQLiteCacheFilename(String baf_filename)
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
		private static double[] GetBafDoubleArray(UInt64 handle, UInt64 id)
		{
			UInt64 n = 0;
			baf2sql_array_get_num_elements(handle, id, ref n);

			double[] myArray = new double[n];
			int rc = baf2sql_array_read_double(handle, id, myArray);
			if (rc == 0) ThrowLastBaf2SqlError();

			return myArray;
		}

		/* Return array 'id', converting to float format */
		private static float[] GetBafFloatArray(UInt64 handle, UInt64 id)
		{
			UInt64 n = 0;
			baf2sql_array_get_num_elements(handle, id, ref n);

			float[] myArray = new float[n];
			int rc = baf2sql_array_read_float(handle, id, myArray);
			if (rc == 0) ThrowLastBaf2SqlError();

			return myArray;
		}

		/* Return array 'id', converting to UInt32 format */
		private static UInt32[] GetBafUInt32Array(UInt64 handle, UInt64 id)
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
}
