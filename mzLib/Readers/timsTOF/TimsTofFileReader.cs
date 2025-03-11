﻿using System;
using System.Runtime.InteropServices;
using System.Text;
using MassSpectrometry;
using System.Data.SQLite;
using Easy.Common.Extensions;
using MzLibUtil;
using UsefulProteomicsDatabases;
using System.Data.Common;
using Readers;
using System.Data.SqlClient;
using System.Data;
using ThermoFisher.CommonCore.Data.Business;
using Polarity = MassSpectrometry.Polarity;
using System.Security.AccessControl;
using System.Collections.Concurrent;
using System.Diagnostics;
using System.Security.Permissions;
using System.ComponentModel;

namespace Readers
{ 
    public class TimsTofFileReader : MsDataFile, IDisposable
    {
        // timsTOF instruments collect frames, packets of ions collected by the tims, then analyzed 
        // over multiple scans with each scan corresponding to the same retention time but different
        // ion mobility valuess. When reading the file, multiple scans from the same frame are collapsed into 
        // a single spectrum

        public TimsTofFileReader(string filePath) : base (filePath) { }

        private UInt64? _fileHandle;
        private Object _fileLock;
        private SQLiteConnection? _sqlConnection;
        private int _maxThreads;
        public int NumberOfFrames { get; private set; }
        public List<long> Ms1FrameIds { get; private set; }
        internal FrameProxyFactory FrameProxyFactory { get; private set; }
        
        // I don't know what the default scan range is, and at this point I'm too afraid to ask...
        private MzRange? _scanWindow;
        public MzRange ScanWindow => _scanWindow ??= new MzRange(20, 2000);
        public const string ScanFilter = "f";

        public override void InitiateDynamicConnection()
        {
            if (!File.Exists(FilePath + @"\analysis.tdf") | !File.Exists(FilePath + @"\analysis.tdf_bin"))
            {
                throw new FileNotFoundException("Data file is missing .tdf and/or .tdf_bin file");
            }

            OpenSqlConnection();

            if(_fileHandle != null) tims_close((UInt64)_fileHandle);
            OpenBinaryFileConnection();
            _fileLock = new();

            CountFrames();
            BuildProxyFactory();
        }

        internal void OpenSqlConnection()
        {
            if (_sqlConnection?.State == ConnectionState.Open)
                return;

            _sqlConnection = new SQLiteConnection("Data Source=" +
                Path.Combine(FilePath, "analysis.tdf") +
                "; Version=3");
            try
            {
                _sqlConnection.Open();
            }
            catch (Exception e)
            {
                throw new MzLibException("Error opening the .tdf file: " + e.Message);
            }
        }

        internal void OpenBinaryFileConnection()
        {
            byte[] binaryFileBytePath = BrukerFileReader.ConvertStringToUTF8ByteArray(FilePath);
            _fileHandle = tims_open(binaryFileBytePath, 0);
            if (_fileHandle == null || _fileHandle == 0)
                throw new MzLibException("Could not open the analysis.tdf_bin file");
        }

        public override void CloseDynamicConnection()
        {
            if (_sqlConnection?.State == ConnectionState.Open) _sqlConnection.Close();
            _sqlConnection?.Dispose();
            if (_fileHandle != null)
            {
                tims_close((UInt64)_fileHandle);
                _fileHandle = null;
            }   
        }

        public void Dispose()
        {
            CloseDynamicConnection();
        }

        /// <summary>
        /// WARNING! This method reads in the entire data file before
        /// returning the requested scan! It is recommended to call the 
        /// GetScanFromPrecursorAndFrameIdFromDynamicConnection()
        /// </summary>
        public override MsDataScan GetOneBasedScanFromDynamicConnection(int oneBasedScanNumber, IFilteringParams filterParams = null)
        {
            if(oneBasedScanNumber <= 0)
                throw new IndexOutOfRangeException("Invalid one-based index given when accessing data scans. Index: " + oneBasedScanNumber);    
            if (Scans != null && Scans.Length >= oneBasedScanNumber && Scans[oneBasedScanNumber - 1] != null)
                return Scans[oneBasedScanNumber - 1];

            LoadAllStaticData(filteringParams: (FilteringParams)filterParams);
            if (oneBasedScanNumber > Scans.Length)
                throw new IndexOutOfRangeException("Invalid one-based index given when accessing data scans. Index: " + oneBasedScanNumber);
            return Scans[oneBasedScanNumber - 1];
        }

        /// <summary>
        /// Returns a TimsDataScan with the specified frame and precursor id
        /// WARNING! The returned data scan will have a OneBasedScanNumber of -1
        /// </summary>
        /// <param name="precursorId"></param>
        /// <param name="frameId"></param>
        /// <param name="filteringParams"></param>
        /// <returns></returns>
        /// <exception cref="MzLibException"></exception>
        /// <exception cref="ArgumentException"></exception>
        /// <exception cref="NotImplementedException"></exception>
        public TimsDataScan GetScanFromPrecursorAndFrameIdFromDynamicConnection(int precursorId, int frameId, IFilteringParams filteringParams = null)
        {
            if(_fileHandle == null || _fileHandle == 0 || _sqlConnection.IsCanceled() || FrameProxyFactory == null)
            {
                throw new MzLibException("The dynamic connection has not been created yet!");
            }

            var frameType = FrameProxyFactory.FramesTable.GetAnalysisType(frameId);
            switch(frameType)
            {
                case TimsTofMsMsType.MS:
                    var records = GetMs1Records(frameId);
                    var recordForPrecursor = records.FirstOrDefault(x => x.PrecursorId == precursorId);
                    return GetMs1Scan(recordForPrecursor, FrameProxyFactory.GetFrameProxy(frameId), (FilteringParams)filteringParams);
                case TimsTofMsMsType.PASEF:
                    return BuildPasefScanFromPrecursor([precursorId], (FilteringParams)filteringParams).FirstOrDefault();
                default:
                    throw new NotImplementedException();
            }
        }

        internal void CountFrames()
        {
            if (_sqlConnection == null) return;
            using var command = new SQLiteCommand(_sqlConnection);
            command.CommandText = @"SELECT COUNT(*) FROM Frames;";
            using var sqliteReader = command.ExecuteReader();
            int count = 0;
            while (sqliteReader.Read())
            {
                count = sqliteReader.GetInt32(0);
                break;
            }
            NumberOfFrames = count;
        }

        internal void CountMS1Frames()
        {
            if (_sqlConnection == null) return;
            using var command = new SQLiteCommand(_sqlConnection);
            command.CommandText = @"SELECT f.Id FROM Frames f WHERE f.MsMsType = 0;";
            using var sqliteReader = command.ExecuteReader();
            Ms1FrameIds = new();

            while (sqliteReader.Read())
            {
                Ms1FrameIds.Add(sqliteReader.GetInt64(0));
            }
            
        }

        /// <summary>
        /// Builds a new FrameProxyFactory to pull frames from the timsTOF data file
        /// and sets the FrameProxyFactory property 
        /// </summary>
        /// <exception cref="MzLibException"></exception>
        internal void BuildProxyFactory()
        {
            if (_sqlConnection == null || _fileHandle == null) return;
            var framesTable = new FrameTable(_sqlConnection, NumberOfFrames);
            if (framesTable == null)
                throw new MzLibException("Something went wrong while loading the Frames table from the analysis.tdf database.");

            int numberOfIndexedMzs = GetNumberOfDigitizerSamples();
            FrameProxyFactory = new FrameProxyFactory(framesTable, (ulong)_fileHandle, _fileLock, numberOfIndexedMzs);
        }

        internal void CountPrecursors()
        {
            if (_sqlConnection == null) return;
            using var command = new SQLiteCommand(_sqlConnection);
            command.CommandText = @"SELECT MAX(Id) FROM Precursors;";
            using var sqliteReader = command.ExecuteReader();
            var columns = Enumerable.Range(0, sqliteReader.FieldCount)
                .Select(sqliteReader.GetName).ToList();
            long maxPrecursorId = 0;
            while (sqliteReader.Read())
            {
                maxPrecursorId = sqliteReader.GetInt64(0);
            }
            Ms1ScanArray = new TimsDataScan[maxPrecursorId];
            PasefScanArray = new TimsDataScan[maxPrecursorId];
        }

        public ConcurrentBag<TimsDataScan> Ms1ScansNoPrecursorsBag { internal get; set; }
        public TimsDataScan[] Ms1ScanArray { internal get; set; }
        public TimsDataScan[] PasefScanArray { internal get; set; }

        internal int GetNumberOfDigitizerSamples()
        {
            using var command = new SQLiteCommand(_sqlConnection);
            command.CommandText = @"SELECT value FROM GlobalMetadata" +
                " WHERE GlobalMetadata.Key = 'DigitizerNumSamples'";
            using var reader = command.ExecuteReader();
            reader.Read();
            return Int32.Parse(reader.GetString(0));
        }

        public override MsDataFile LoadAllStaticData(FilteringParams filteringParams = null, int maxThreads = 1)
        {
            InitiateDynamicConnection();

            CountMS1Frames();
            CountPrecursors();
            
            _maxThreads = maxThreads;
            Ms1ScansNoPrecursorsBag = new();
            Parallel.ForEach(
                Partitioner.Create(0, Ms1FrameIds.Count),
                new ParallelOptions() { MaxDegreeOfParallelism = _maxThreads },
                (range) =>
            {
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    BuildAllScans(Ms1FrameIds[i], filteringParams);
                }
            });

            CloseDynamicConnection();
            AssignOneBasedPrecursorsToPasefScans();
            SourceFile = GetSourceFile();
            return this;
        }

        internal void AssignOneBasedPrecursorsToPasefScans()
        {
            var localMs1Scans = this.Ms1ScanArray.Where(scan => scan != null).OrderBy(scan => scan.FrameId).ThenBy(scan => scan.PrecursorId).ToList();
            var localPasefScans = this.PasefScanArray.Where(scan => scan != null).OrderBy(scan => scan.PrecursorId).ToList();
            var localMs1ScansNoPrecursor = Ms1ScansNoPrecursorsBag.OrderBy(scan => scan.FrameId).ToList();
            TimsDataScan[] scanArray = new TimsDataScan[localMs1Scans.Count*2 + localMs1ScansNoPrecursor.Count];

            int oneBasedScanIndex = 1;
            int pasefScanIndex = 0;
            int ms1NoPrecursorIndex = 0;
            TimsDataScan? ms1ScanNoPrecursor = localMs1ScansNoPrecursor.IsNotNullOrEmpty() ? localMs1ScansNoPrecursor[ms1NoPrecursorIndex] : null;
            //Write the scans to the scanArray and assign scan indices
            for (int i = 0; i < localMs1Scans.Count; i++)
            {
                var ms1Scan = localMs1Scans[i];
                while (ms1ScanNoPrecursor != null && ms1ScanNoPrecursor.FrameId < ms1Scan.FrameId)
                {
                    ms1ScanNoPrecursor.SetOneBasedScanNumber(oneBasedScanIndex);
                    scanArray[oneBasedScanIndex - 1] = ms1ScanNoPrecursor;
                    ms1NoPrecursorIndex++;
                    oneBasedScanIndex++;
                    ms1ScanNoPrecursor = ms1NoPrecursorIndex < localMs1ScansNoPrecursor.Count ? localMs1ScansNoPrecursor[ms1NoPrecursorIndex] : null;
                }
                ms1Scan.SetOneBasedScanNumber(oneBasedScanIndex);
                scanArray[oneBasedScanIndex - 1] = ms1Scan;
                oneBasedScanIndex++;
                if (ms1Scan.PrecursorId == -1) continue; // Continue if the scan didn't have any precursors (as there will be no MS2 scans)
                
                // This assumes that there is a one to one correspondence between the MS1 scans and the PASEF scans
                var pasefScan = localPasefScans[pasefScanIndex];
                while(pasefScan.PrecursorId < ms1Scan.PrecursorId)
                {
                    pasefScanIndex++;
                    pasefScan = localPasefScans[pasefScanIndex];
                }
                if(pasefScan.PrecursorId == ms1Scan.PrecursorId)
                {
                    pasefScan.SetOneBasedPrecursorScanNumber(ms1Scan.OneBasedScanNumber);
                    pasefScan.SetOneBasedScanNumber(oneBasedScanIndex);
                    scanArray[oneBasedScanIndex - 1] = pasefScan;
                    pasefScanIndex++;
                    oneBasedScanIndex++;
                }
            }

            // Some MS1 scans contain no peaks where the precursor was identified, so they are not included in the scanArray
            Scans = scanArray.Where(scan => scan != null).ToArray();
        }

        /// <summary>
        /// This function will create multiple MS1 scans from each MS1 frame in the timsTOF data file
        /// One Ms1 Scan per precursor
        /// It will then create an Ms2 scan for each precursor by averaging MS2 spectra for the precursor
        /// collected over multiple frames
        /// Created scans are then added to their respective scan arrays
        /// </summary>
        /// <param name="frameId"></param>
        /// <param name="filteringParams"></param>
        internal void BuildAllScans(long frameId, FilteringParams filteringParams)
        {
            FrameProxy frame = FrameProxyFactory.GetFrameProxy(frameId);
            var records = GetMs1Records(frameId);
            foreach(Ms1Record record in records)
            {
                TimsDataScan? dataScan = GetMs1Scan(record, frame, filteringParams);
                if (dataScan != null)
                {
                    if (dataScan.PrecursorId > 0)
                        Ms1ScanArray[(int)dataScan.PrecursorId - 1] = dataScan;
                    else
                        Ms1ScansNoPrecursorsBag.Add(dataScan);
                }
            }

            // Then, build ONE MS2 scan by averaging every PASEF frame that sampled that precursor
            var pasefScans = BuildPasefScanFromPrecursor(precursorIds: records.Select(r => r.PrecursorId), filteringParams);
            foreach (var scan in pasefScans)
            {
                if (scan?.PrecursorId != null)
                    PasefScanArray[(int)scan.PrecursorId - 1] = scan;
            }
        }

        internal List<Ms1Record> GetMs1Records(long frameId)
        {
            List<Ms1Record> records = new List<Ms1Record>();
            // Only do this if we have valid precursors (which we don't for like SRM/inclusion list type stuff) 
            using (var command = new SQLiteCommand(_sqlConnection))
            {
                // This command finds all the precursors identified and fragmented in each MS/MS Pasef scan
                // It is used to take an MS1 frame and create multiple "MsDataScans" by averaging the 
                // spectra from each scan within a given Ion Mobility (i.e. ScanNum) range
                command.CommandText =
                    @"SELECT MIN(m.ScanNumBegin), MAX(m.ScanNumEnd), p.ScanNumber, p.Id" +
                    " FROM Precursors p" +
                    " INNER JOIN PasefFrameMsMsInfo m on m.Precursor = p.Id" +
                    " WHERE p.Parent = " + frameId.ToString() +
                    " GROUP BY p.Id;";
                using var sqliteReader = command.ExecuteReader();

                while (sqliteReader.Read())
                {
                    var scanStart = sqliteReader.GetInt32(0);
                    var scanEnd = sqliteReader.GetInt32(1);
                    var scanMedian = sqliteReader.GetFloat(2);
                    int precursorId = sqliteReader.GetInt32(3);
                    records.Add(new Ms1Record(precursorId, scanStart, scanEnd, (double)scanMedian));
                }
            }
            // If no records were returned, then no precursors were observed in the frame
            // In that case, create a record that contains every scan and a precursorID of -1
            if (records.Count == 0)
                records.Add(new Ms1Record(-1, 1, FrameProxyFactory.FramesTable.NumScans[frameId - 1], FrameProxyFactory.FramesTable.NumScans[frameId - 1]));
            return records;
        }

        internal TimsDataScan? GetMs1Scan(Ms1Record record, FrameProxy frame, FilteringParams filteringParams)
        {
            List<uint[]> indexArrays = new();
            List<int[]> intensityArrays = new();
            for (int scan = record.ScanStart; scan < record.ScanEnd; scan++)
            {
                indexArrays.Add(frame.GetScanIndices(scan-1));
                intensityArrays.Add(frame.GetScanIntensities(scan-1));
            }
            // Step 2: Average those suckers
            MzSpectrum averagedSpectrum = TofSpectraMerger.MergeArraysToMs1Spectrum(indexArrays, intensityArrays, FrameProxyFactory, filteringParams: filteringParams);
            if (averagedSpectrum.Size < 1)
            {
                return null;
            }
            // Step 3: Make an MsDataScan bby
            var dataScan = new TimsDataScan(
                massSpectrum: averagedSpectrum,
                oneBasedScanNumber: -1, // This gets adjusted once all data has been read
                msnOrder: 1,
                isCentroid: true,
                polarity: FrameProxyFactory.GetPolarity(frame.FrameId),
                retentionTime: FrameProxyFactory.GetRetentionTime(frame.FrameId),
                scanWindowRange: ScanWindow,
                scanFilter: ScanFilter,
                mzAnalyzer: MZAnalyzerType.TOF,
                totalIonCurrent: intensityArrays.Sum(array => array.Sum()),
                injectionTime: FrameProxyFactory.GetInjectionTime(frame.FrameId),
                noiseData: null,
                nativeId: "frame=" + frame.FrameId.ToString() +
                    ";scans=" + record.ScanStart.ToString() + "-" + record.ScanEnd.ToString() +
                    ";precursor=" + record.PrecursorId.ToString(),
                frameId: frame.FrameId,
                scanNumberStart: record.ScanStart,
                scanNumberEnd: record.ScanEnd,
                medianOneOverK0: FrameProxyFactory.GetOneOverK0(record.ScanMedian),
                precursorId: record.PrecursorId);

            return dataScan;
        }

        internal List<TimsDataScan> BuildPasefScanFromPrecursor(IEnumerable<int> precursorIds, FilteringParams filteringParams)
        {
            HashSet<long> allFrames = new();
            List<TimsDataScan> pasefScans = new();

            // Create  TimsDataScans with all relevant metadata, but without an mzSpectrum
            foreach (PasefRecord record in GetPasefRecords(precursorIds))
            {
                allFrames.UnionWith(record.FrameList);
                var dataScan = new TimsDataScan(
                    massSpectrum: null,
                    oneBasedScanNumber: -1, // This will be adjusted once all scans have been read
                    msnOrder: 2,
                    isCentroid: true,
                    polarity: FrameProxyFactory.GetPolarity(record.FrameList.First()),
                    retentionTime: FrameProxyFactory.GetRetentionTime(record.FrameList.First()),
                    scanWindowRange: ScanWindow,
                    scanFilter: ScanFilter,
                    mzAnalyzer: MZAnalyzerType.TOF,
                    totalIonCurrent: -1, // Will be set later
                    injectionTime: FrameProxyFactory.GetInjectionTimeSum(record.FrameList.First(), record.FrameList.Last()),
                    noiseData: null,
                    nativeId: "frames=" + record.FrameList.First().ToString() + "-" + record.FrameList.Last().ToString() +
                              ";scans=" + record.ScanStart.ToString() + "-" + record.ScanEnd.ToString(),
                    frameId: record.FrameList.First(),
                    scanNumberStart: record.ScanStart,
                    scanNumberEnd: record.ScanEnd,
                    medianOneOverK0: FrameProxyFactory.GetOneOverK0(record.ScanMedian), // Needs to be set later
                    precursorId: record.PrecursorId,
                    selectedIonMz: record.MostAbundantPrecursorMz,
                    selectedIonChargeStateGuess: record.Charge,
                    selectedIonIntensity: record.PrecursorIntensity,
                    isolationMZ: record.IsolationMz,
                    isolationWidth: record.IsolationWidth,
                    dissociationType: DissociationType.CID,
                    oneBasedPrecursorScanNumber: -1, // This will be set later
                    selectedIonMonoisotopicGuessMz: record.PrecursorMonoisotopicMz,
                    hcdEnergy: record.CollisionEnergy.ToString(),
                    frames: record.FrameList.ToList());
                pasefScans.Add(dataScan);
            }

            // Grab all fragmentation spectra for each precursor
            // Each TimsDataScan in pasefScans corresponds to one precursor.
            // A precursor can be isolated and fragmented in multiple pasef frames
            // Here, we iterate through each frame, averaging the scans that correspond to each precursor
            foreach (long frameId in allFrames)
            {
                FrameProxy frame = FrameProxyFactory.GetFrameProxy(frameId);
                //Iterate through all the datascans created above with this frame
                foreach (var scan in pasefScans)
                {
                    if (scan.FrameIds.Contains(frameId))
                    {
                        List<uint[]> indexArrays = new();
                        List<int[]> intensityArrays = new();
                        for (int mobilityScanIdx = scan.ScanNumberStart; mobilityScanIdx < scan.ScanNumberEnd; mobilityScanIdx++)
                        {
                            indexArrays.Add(frame.GetScanIndices(mobilityScanIdx-1));
                            intensityArrays.Add(frame.GetScanIntensities(mobilityScanIdx-1));
                        }
                        // Perform frame level averaging, where all scans from one frame associated with a given precursor are merged and centroided
                        // Need to convert indexArrays to one uint[] and intensityArrays to one int[]
                        (double[] Mzs, int[] Intensities) summedArrays = TofSpectraMerger.MergeArraysToMzArray(indexArrays, intensityArrays, FrameProxyFactory);
                        scan.AddComponentArrays(summedArrays.Mzs, summedArrays.Intensities);
                    }
                }
            }

            // Now, we average the fragmentation spectra (each spectra originating in a different frame)
            // to yield one spectrum per precursor
            foreach (TimsDataScan scan in pasefScans)
            {
                scan.AverageComponentSpectra(FrameProxyFactory, filteringParams);
            }

            return pasefScans;
        }

        internal IEnumerable<PasefRecord> GetPasefRecords(IEnumerable<int> precursorIds)
        {
            using (var command = new SQLiteCommand(_sqlConnection))
            {
                string multiplePrecursorString = "(" +
                    String.Join(',', precursorIds.Select(id => "\'" + id.ToString() + "\'")) +
                    ")";
                // SQL Command for getting some info from both PasefFrameMsMsInfo table and
                // Precursors table
                command.CommandText =
                    @"SELECT GROUP_CONCAT(m.Frame), m.ScanNumBegin, m.ScanNumEnd, m.IsolationMz, m.IsolationWidth," +
                    " m.CollisionEnergy, p.LargestPeakMz, p.MonoisotopicMz, p.Charge, p.Intensity, p.ScanNumber, p.Id" +
                    " FROM PasefFrameMsMsInfo m" +
                    " INNER JOIN Precursors p on m.Precursor = p.Id" +
                    " WHERE m.Precursor IN " + multiplePrecursorString +
                    " GROUP BY m.Precursor;";

                using var sqliteReader = command.ExecuteReader();

                // Each call to read returns the information associated with a given precursor
                while (sqliteReader.Read())
                {
                    var frameList = sqliteReader.GetString(0).Split(',').Select(id => Int64.Parse(id));
                    var scanStart = sqliteReader.GetInt32(1);
                    var scanEnd = sqliteReader.GetInt32(2);
                    var isolationMz = sqliteReader.GetFloat(3);
                    var isolationWidth = sqliteReader.GetFloat(4);
                    var collisionEnergy = sqliteReader.GetFloat(5);
                    var mostAbundantPrecursorPeak = sqliteReader.GetFloat(6);
                    float precursorMonoisotopicMz = sqliteReader.IsDBNull(7) ? isolationMz : sqliteReader.GetFloat(7);
                    int charge = sqliteReader.IsDBNull(8) ? 1 : sqliteReader.GetInt32(8);
                    var precursorIntensity = sqliteReader.GetFloat(9);
                    var scanMedian = sqliteReader.GetFloat(10);
                    var precursorId = sqliteReader.GetInt32(11);

                    yield return new PasefRecord(frameList, precursorId, scanStart, scanEnd, scanMedian, isolationMz, isolationWidth, collisionEnergy, mostAbundantPrecursorPeak, precursorMonoisotopicMz, charge, precursorIntensity);
                }
            }
        }

        private const string nativeIdFormat = "Frame ID + scan number range format";
        private const string massSpecFileFormat = ".D format";
        public override SourceFile GetSourceFile()
        {
            // append the analysis.baf because the constructor for SourceFile will look for the 
            // parent directory. 
            string fileName = FilePath + @"\analysis.tdf";
            return new SourceFile(nativeIdFormat, massSpecFileFormat,
                null, null, id: null, filePath: fileName);
        }

        #region Bruker Dll Functions 

        /// <summary>
        /// Returns a unique handle that references an open timsTOF data file
        /// </summary>
        /// <param name="analysis_directory_name_utf8"></param>
        /// <param name="use_recalibrated_state"></param>
        /// <returns></returns>
        [DllImport("timsdata.dll", CharSet = CharSet.Ansi, CallingConvention = CallingConvention.Cdecl)]
        public static extern UInt64 tims_open
              (byte[] analysis_directory_name_utf8, UInt32 use_recalibrated_state);

        /// <summary>
        /// Closes a file connection to a .tdf binary file
        /// </summary>
        [DllImport("timsdata.dll", CharSet = CharSet.Ansi, CallingConvention = CallingConvention.Cdecl)]
        public static extern void tims_close
              (UInt64 fileHandle);

        #endregion Bruker Dll Functions

    }
}
