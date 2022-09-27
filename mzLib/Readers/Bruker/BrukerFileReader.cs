using System.Data.SQLite; 
using System.Data; 
using MassSpectrometry;
using MzLibUtil;

namespace Readers
{
    public class BrukerFileConnection : IDisposable
    {
        public readonly string PathToBaf;
        public int TotalScans { get; private set; }
        private UInt64 Handle { get; set; }
        private String SqliteFn { get; set; }
        private byte[] FnUtf8 { get; set; }
        public SQLiteConnection SqlConnection { get; private set; }
        private DataTable MetaDataTable { get; set; } = new();
        private DataTable LineSpectraTable { get; set; } = new();
        private DataTable SupportedVarsTable { get; set; } = new();


        public BrukerFileConnection(string pathToBaf)
        {
            PathToBaf = pathToBaf;
        }

        public void InitializeConnection()
        {
            SqliteFn = BafReaderExtern.GetSQLiteCacheFilename(PathToBaf);
            FnUtf8 = BafReaderExtern.ConvertStringToUTF8ByteArray(PathToBaf);
            Handle = BafReaderExtern.baf2sql_array_open_storage(1, FnUtf8);
            SqlConnection = new SQLiteConnection();
            SetConnectionString();
            SqlConnection.Open();
        }

        private void SetConnectionString()
        {
            SqlConnection.ConnectionString = "DataSource=" + SqliteFn;
        }

        public int GetTotalSpectraTableRows()
        {
            return LineSpectraTable.Rows.Count;
        }
        public MsDataScan CreateMsDataScan(int rowNumber)
        {
            int sqliteRowNumber = rowNumber;
            double[] xArray = GetLineMzId(sqliteRowNumber);
            double[] yArray = GetLineIntensityId(sqliteRowNumber);

            MzSpectrum mzSpectrum = new MzSpectrum(xArray, yArray, false);

            GetSpectraAssociatedData(sqliteRowNumber,
                out int oneBasedScanNumber,
                out double retentionTime,
                out int msnOrder,
                out int? parentScan,
                out int lowMz,
                out int upperMz,
                out double totalIonCurrent);
            var metaDataDict = GetScanMetaData(sqliteRowNumber);
            string scanFilter = string.Empty;
            string nativeId = string.Empty;
            return new MsDataScan(
                mzSpectrum,
                oneBasedScanNumber: oneBasedScanNumber,
                msnOrder: msnOrder,
                isCentroid: true,
                polarity: Polarity.Positive,
                retentionTime: retentionTime,
                scanWindowRange: new MzRange(lowMz, upperMz),
                scanFilter: scanFilter,
                mzAnalyzer: MZAnalyzerType.TOF,
                totalIonCurrent: totalIonCurrent,
                injectionTime: null,
                noiseData: new double[,] { { 0, 0 }, { 0, 0 } },
                nativeId: nativeId,
                selectedIonMz: metaDataDict["MSMS_IsolationMass_Act"],
                selectedIonChargeStateGuess: (int)metaDataDict["MSMS_PreCursorChargeState"],
                selectedIonIntensity: null,
                isolationMZ: metaDataDict["MSMS_IsolationMass_Act"],
                isolationWidth: metaDataDict["Quadrupole_IsolationResolution_Act"],
                DissociationType.CID,
                oneBasedPrecursorScanNumber: parentScan,
                selectedIonMonoisotopicGuessMz: null,
                hcdEnergy: metaDataDict["Collision_Energy_Act"].ToString()
            );
        }

        #region Data
        public void GetDataTable()
        {
            try
            {
                using (SQLiteCommand command = new SQLiteCommand(SqlConnection))
                {
                    command.CommandText =
                        "SELECT * FROM Spectra";
                    using (SQLiteDataReader reader = command.ExecuteReader())
                    {
                        LineSpectraTable.Load(reader);
                    }
                }
            }
            catch (Exception e)
            {
                Console.WriteLine(e);
                throw;
            }
        }

        public void GetSupportedVariablesTable()
        {
            try
            {
                using (SQLiteCommand command = new SQLiteCommand(SqlConnection))
                {
                    command.CommandText =
                        "SELECT * FROM SupportedVariables";
                    using (SQLiteDataReader reader = command.ExecuteReader())
                    {
                        SupportedVarsTable.Load(reader);
                    }
                }
            }
            catch (Exception e)
            {
                Console.WriteLine(e);
                throw;
            }
        }

        public double[] GetLineMzId(int scanId)
        {
            var row = LineSpectraTable.Rows[scanId];
            var id = (UInt64)row.Field<Int64>("LineMzId");
            return BafReaderExtern.GetBafDoubleArray(Handle, id);
        }

        public double[] GetLineIndexId(int scanId)
        {
            var row = LineSpectraTable.Rows[scanId];
            var id = (UInt64)row.Field<Int64>("LineIndexId");
            return BafReaderExtern.GetBafDoubleArray(Handle, id);
        }

        public double[] GetLineIntensityId(int scanId)
        {
            var row = LineSpectraTable.Rows[scanId];
            var id = (UInt64)row.Field<Int64>("LineIntensityId");
            return BafReaderExtern.GetBafDoubleArray(Handle, id);
        }

        public double[] GetLineIndexWidthId(int scanId)
        {
            var row = LineSpectraTable.Rows[scanId];
            var id = (UInt64)row.Field<Int64>("LineIndexWidthId");
            return BafReaderExtern.GetBafDoubleArray(Handle, id);
        }

        public double[] GetLineSnrId(int scanId)
        {
            var row = LineSpectraTable.Rows[scanId];
            var id = (UInt64)row.Field<Int64>("LineSnrId");
            return BafReaderExtern.GetBafDoubleArray(Handle, id);
        }

        public double[] GetLinePeakAreaId(int scanId)
        {
            var row = LineSpectraTable.Rows[scanId];
            var id = (UInt64)row.Field<Int64>("LinePeakAreaId");
            return BafReaderExtern.GetBafDoubleArray(Handle, id);
        }
        // Needs to get Id, Rt, AcquisitionKey(MsnOrder), Parent (parentScan if Ms2), 
        // MzAcqRangeLower, MzAcqRangeUpper, SumIntensity
        public void GetSpectraAssociatedData(int scanNumber, out int oneBasedScanNumber,
            out double retentionTime, out int msnOrder, out int? parentScan, out int lowMz,
            out int upperMz, out double totalIonCurrent)
        {
            var row = LineSpectraTable.Rows[scanNumber];

            oneBasedScanNumber = (int)row.Field<long>("Id");
            retentionTime = row.Field<double>("Rt");
            msnOrder = (int)row.Field<long>("AcquisitionKey");
            parentScan = (int?)row.Field<long?>("Parent");
            lowMz = (int)row.Field<long>("MzAcqRangeLower");
            upperMz = (int)row.Field<long>("MzAcqRangeUpper");
            totalIonCurrent = row.Field<double>("SumIntensity");
        }
        #endregion

        // TODO: Refactor so each method has a datarow parameter, so you only get one datarow ever
        #region METADATA
        public void GetMetadataTable()
        {
            // copy the PerSpectrumVariables 
            // combine PerSpectrumVariables with the SupportedVariables table 
            // load into DataTable 
            MetaDataTable = new DataTable();
            MetaDataTable.Columns.Add("Spectrum", typeof(Int64));
            MetaDataTable.Columns.Add("Variable", typeof(Int64));
            MetaDataTable.Columns.Add("Value", typeof(double));

            // column 1 is Spectrum, 2 is Variable, 3 is Value. 
            // 1 is an integer, 2 is an integer, 3 is an integer column, but it 
            // has double values inside it. 
            // CREATE NEW TABLE AND COPY INTO NEW TABLE 
            // DO NOT MODIFY ORIGINAL TABLES. 
            using (SQLiteCommand command = new SQLiteCommand(SqlConnection))
            {
                command.CommandText =
                    command.CommandText =
                        "CREATE TABLE IF NOT EXISTS copied(" +
                        "Spectrum INTEGER," +
                        "Variable INTEGER, " +
                        "Value REAL, " +
                        "PRIMARY KEY (Spectrum, Variable)" +
                        ");" +
                        "INSERT INTO copied(Spectrum, Variable, Value) " +
                        "SELECT Spectrum, Variable, Value " +
                        "FROM PerSpectrumVariables;" +
                        "SELECT * FROM copied";
                try
                {
                    using (SQLiteDataReader reader = command.ExecuteReader())
                    {
                        while (reader.Read())
                        {
                            DataRow dr = MetaDataTable.NewRow();
                            dr[0] = reader.GetInt64(0);
                            dr[1] = reader.GetInt64(1);
                            dr[2] = reader.GetDouble(2);
                            MetaDataTable.Rows.Add(dr);
                        }
                    }
                }
                catch (Exception e)
                {
                    Console.WriteLine(e.ToString());
                    throw;
                }
                finally // cleans up the copied table to prevent errors on repeated runs. 
                {
                    command.CommandText = "DROP TABLE IF EXISTS copied";
                    using (SQLiteDataReader reader = command.ExecuteReader()) ;
                }
            }
        }
        /// <summary>
        /// Outputs the scan metadata for a given scan ID as a Dictionary
        /// </summary>
        /// <param name="scanId"></param> One based scan number. 
        /// <returns>Dictionary with string of the variable name and the associated value.
        /// Note that any doubles are returned as integers cast to doubles. 
        /// </returns>
        public Dictionary<string, double> GetScanMetaData(int scanId)
        {
            int sqliteRowIndex = scanId + 1;
            // filter out all but scan of interest. 
            var rows = MetaDataTable.AsEnumerable()
                .Where(i => (int)i.Field<Int64>("Spectrum") == sqliteRowIndex)
                .Select(i => i);
            var metaDataDict = new Dictionary<string, double>();
            foreach (var row in rows)
            {
                int keyInt = (int)row.Field<Int64>("Variable");
                string keyString = GetMetaDataVariableString(keyInt);
                double value = (double)row.Field<double>("Value");
                metaDataDict.Add(keyString, value);
            }

            return metaDataDict;
        }

        public string GetMetaDataVariableString(int variable)
        {
            return ((BrukerBafVariableNames)variable).ToString();
        }
        #endregion

        public void Dispose()
        {
            BafReaderExtern.baf2sql_array_close_storage(Handle);
        }
    }
}