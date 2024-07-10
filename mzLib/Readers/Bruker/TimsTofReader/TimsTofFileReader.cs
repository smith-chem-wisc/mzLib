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
using System.Data.SqlClient;
using System.Data;
using ThermoFisher.CommonCore.Data.Business;
using Polarity = MassSpectrometry.Polarity;
using System.Security.AccessControl;
using System.Collections.Concurrent;
using System.Diagnostics;
using System.Security.Permissions;

namespace Readers.Bruker
{ 
    public class TimsTofFileReader : MsDataFile
    {
        // timsTOF instruments collect frames, packets of ions collected by the tims, then analyzed 
        // over multiple scans with each scan corresponding to the same retention time but different
        // ion mobility valuess. When reading the file, multiple scans from the same frame are collapsed into 

        public TimsTofFileReader(string filePath) : base (filePath) 
        {
            PrecursorToOneBasedParentScanIndex = new();
        }

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

        /// <summary>
        /// Each precursor is uniquely linked to one MS1 scan
        /// </summary>
        public Dictionary<int, int> PrecursorToOneBasedParentScanIndex { get; private set; }

        public override void InitiateDynamicConnection()
        {
            if (!File.Exists(FilePath + @"\analysis.tdf") | !File.Exists(FilePath + @"\analysis.tdf_bin"))
            {
                throw new FileNotFoundException();
            }
            OpenSqlConnection();
            OpenBinaryFileConnection();
            _fileLock = new();
        }

        public long TotalScans { get; private set; }

        internal void OpenSqlConnection()
        {
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
        }

        public override void CloseDynamicConnection()
        {
            _sqlConnection?.Close();
            if (_fileHandle != null)
                tims_close((UInt64)_fileHandle);
        }

        public override MsDataScan GetOneBasedScanFromDynamicConnection(int oneBasedScanNumber, IFilteringParams filterParams = null)
        {
            throw new NotImplementedException();
        }

        internal void CountFrames()
        {
            if (_sqlConnection == null) return;
            using var command = new SQLiteCommand(_sqlConnection);
            command.CommandText = @"SELECT COUNT(*) FROM Frames;";
            using var sqliteReader = command.ExecuteReader();
            int count = 0;
            var columns = Enumerable.Range(0, sqliteReader.FieldCount)
                .Select(sqliteReader.GetName).ToList();
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
            var columns = Enumerable.Range(0, sqliteReader.FieldCount)
                .Select(sqliteReader.GetName).ToList();
            while (sqliteReader.Read())
            {
                Ms1FrameIds.Add(sqliteReader.GetInt64(0));
            }
        }

        //public List<long> PasefFrameIds { get; private set; }

        //internal void CountPasefFrames()
        //{
        //    if (_sqlConnection == null) return;
        //    using var command = new SQLiteCommand(_sqlConnection);
        //    command.CommandText = @"SELECT f.Id FROM Frames f WHERE f.MsMsType = 8;";
        //    using var sqliteReader = command.ExecuteReader();
        //    PasefFrameIds = new();
        //    var columns = Enumerable.Range(0, sqliteReader.FieldCount)
        //        .Select(sqliteReader.GetName).ToList();
        //    while (sqliteReader.Read())
        //    {
        //        PasefFrameIds.Add(sqliteReader.GetInt64(0));
        //    }
        //}

        //internal void BuildPrecursorToOneBasedParentScanDict()
        //{
        //    PrecursorToOneBasedParentScanIndex = new();
        //    if (_sqlConnection == null) return;
        //    using var command = new SQLiteCommand(_sqlConnection);
        //    command.CommandText = @"SELECT p.Id FROM Precursors p";
        //    using var sqliteReader = command.ExecuteReader();
        //    var columns = Enumerable.Range(0, sqliteReader.FieldCount)
        //        .Select(sqliteReader.GetName).ToList();
        //    while (sqliteReader.Read())
        //    {
        //        PrecursorToOneBasedParentScanIndex.Add(sqliteReader.GetInt32(0), 0);
        //    }
        //}

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
            FrameProxyFactory = new FrameProxyFactory(framesTable, (ulong)_fileHandle, _fileLock);
            if (FrameProxyFactory == null)
                throw new MzLibException("Something went wrong constructing the FrameProxyFactory.");
        }

        public override MsDataFile LoadAllStaticData(FilteringParams filteringParams = null, int maxThreads = 1)
        {
            InitiateDynamicConnection();
            if (_fileHandle == null || _fileHandle == 0)
                throw new MzLibException("Could not open the analysis.tdf_bin file");
            CountFrames();
            CountMS1Frames();
            //CountPasefFrames();
            //BuildPrecursorToOneBasedParentScanDict(); // This only needs to be done for DDA data
            BuildProxyFactory();
            
            _maxThreads = maxThreads;

            if (false)
            {
                ConcurrentDictionary<int, List<MsDataScan>> scanListDictionary = new();
                //Parallel.ForEach(Partitioner.Create(0, Ms1FrameIds.Count),
                Parallel.ForEach(Partitioner.Create(0, Ms1FrameIds.Count),
                        new ParallelOptions { MaxDegreeOfParallelism = maxThreads }, (range, loopState) =>
                        {
                            //Debugger.Break();
                            List<MsDataScan> scansInRange = new();
                            //var newConnection = new SQLiteConnection("Data Source=" +
                            //    Path.Combine(FilePath, "analysis.tdf") +
                            //    "; Version=3");
                            //newConnection.Open();
                            for (int i = range.Item1; i < range.Item2; i++)
                            {
                                long frameId = Ms1FrameIds[i];
                                //Debugger.Break();

                                //BuildMS1Scans(scansInRange, frameId, newConnection);
                                //BuildMS1Scans(scansInRange, frameId);
                            }
                            //Debugger.Break();
                            scanListDictionary[range.Item1] = scansInRange;
                            //newConnection.Close();
                        });
                List<MsDataScan> combinedScanList = new();
                foreach (var kvp in scanListDictionary.OrderBy(kvp => kvp.Key))
                {
                    combinedScanList.AddRange(kvp.Value);
                }

                CloseDynamicConnection();

                Scans = combinedScanList.ToArray();
                SourceFile = GetSourceFile();
                return this;
            }
            else
            {
                List<TimsDataScan> scanList = new();
                for (int i = 0; i < Ms1FrameIds.Count; i++)
                {
                    long frameId = Ms1FrameIds[i];
                    long nextFrameId = i < (Ms1FrameIds.Count -1) ? Ms1FrameIds[i+1] : long.MaxValue;
                    BuildMS1Scans(scanList, frameId, nextFrameId, filteringParams);

                    //if (FramesTable.MsMsType[i].ToEnum<TimsTofMsMsType>(out var msMsType))
                    //{
                    //    switch (msMsType)
                    //    {
                    //        case TimsTofMsMsType.MS:
                    //            BuildMS1Scans(scanList, currentFrame, i);
                    //            break;
                    //        case TimsTofMsMsType.PASEF:
                    //            BuildPasefScans(scanList, currentFrame, i);
                    //            break;
                    //        default:
                    //            throw new NotImplementedException("Only PASEF data is currently supported");
                    //    }
                    //}

                }
                CloseDynamicConnection();
                AssignOneBasedPrecursorsToPasefScans(scanList);
                Scans = scanList.ToArray();
                SourceFile = GetSourceFile();
                return this;
            }
        }

        internal void AssignOneBasedPrecursorsToPasefScans(List<TimsDataScan> scanList)
        {
            var pasefScans = scanList.Where(scan => scan.MsnOrder > 1);
            foreach (var scan in pasefScans)
            {
                if (scan.PrecursorId != null 
                    && PrecursorToOneBasedParentScanIndex.TryGetValue((int)scan.PrecursorId, out int oneBasedParentScanIndex))
                {
                    scan.SetOneBasedPrecursorScanNumber(oneBasedParentScanIndex);
                }
            }
        }

        /// <summary>
        /// This function will create multiple MS1 scans from each MS1 frame in the timsTOF data file
        /// The spectra of every scan in a given Ion Mobility range will be averaged to create a single spectrum
        /// </summary>
        /// <param name="scanList"></param>
        /// <param name="frameId"></param>
        /// <param name="filteringParams"></param>
        /// <param name="sqLiteConnection"></param>
        internal void BuildMS1Scans(List<TimsDataScan> scanList, long frameId, long nextFrameId, FilteringParams filteringParams)
        {
            FrameProxy frame = FrameProxyFactory.GetFrameProxy(frameId);
            List<Ms1Record> records = new();

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
                    " WHERE p.Parent = " + frame.FrameId.ToString() +
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
            if(records.Count == 0)
            {
                BuildMs1ScanNoPrecursor(scanList, frame, filteringParams);
                return;
            }

            ConcurrentBag<TimsDataScan> scanBag = new();
            Parallel.ForEach(records, record =>
            {
                TimsDataScan dataScan = GetMs1Scan(record, frame, scanList.Count + 1, filteringParams);
                if(dataScan != null)
                {
                    scanBag.Add(dataScan);
                }
            });

            foreach (TimsDataScan scan in scanBag.OrderBy(scan => scan.PrecursorId))
            {
                scan.SetOneBasedScanNumber(scanList.Count + 1);
                if(scan.PrecursorId != null)
                {
                    PrecursorToOneBasedParentScanIndex[(int)scan.PrecursorId] = scan.OneBasedScanNumber;
                }
                scanList.Add(scan);
            }

            // Then, build ONE MS2 scan from every PASEF frame that sampled that precursor
            BuildPasefScanFromPrecursor(
                scanList,
                precursorIds: scanBag.Where(scan => scan.PrecursorId != null).Select(scan => (int)scan.PrecursorId).Distinct(),
                filteringParams);
        }

        internal TimsDataScan GetMs1Scan(Ms1Record record, FrameProxy frame, int oneBasedScanNumber, FilteringParams filteringParams)
        {
            List<double[]> mzArrays = new();
            List<int[]> intensityArrays = new();
            for (int scan = record.ScanStart; scan < record.ScanEnd; scan++)
            {
                mzArrays.Add(frame.GetScanMzs(scan));
                intensityArrays.Add(frame.GetScanIntensities(scan));
            }
            // Step 2: Average those suckers
            MzSpectrum averagedSpectrum = SumScans(mzArrays, intensityArrays, filteringParams);
            if(averagedSpectrum.Size < 1)
            {
                return null;
            }
            // Step 3: Make an MsDataScan bby
            var dataScan = new TimsDataScan(
                massSpectrum: averagedSpectrum,
                //massSpectrum: null,
                oneBasedScanNumber: oneBasedScanNumber,
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
                medianOneOverK0: frame.GetOneOverK0((double)record.ScanMedian),
                precursorId: record.PrecursorId);

            return dataScan;
        }

        internal void BuildPasefScanFromPrecursor(List<TimsDataScan> scanList, IEnumerable<int> precursorIds, FilteringParams filteringParams)
        {
            HashSet<long> allFrames = new();
            List<TimsDataScan> pasefScans = new();

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
                int runningTotal = 0;
                
                // Each call to read returns the information associated with a given precursor
                while (sqliteReader.Read())
                {
                    var frameList = sqliteReader.GetString(0).Split(',').Select(id => Int64.Parse(id));
                    allFrames.UnionWith(frameList);
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
                    runningTotal++;

                    // Now, we build data scans with null MzSpectra. The MassSpectrum will be constructed later.
                    var dataScan = new TimsDataScan(
                            massSpectrum: null,
                            oneBasedScanNumber: scanList.Count + pasefScans.Count + 1,
                            msnOrder: 2,
                            isCentroid: true,
                            polarity: FrameProxyFactory.GetPolarity(frameList.First()),
                            retentionTime: FrameProxyFactory.GetRetentionTime(frameList.First()),
                            scanWindowRange: ScanWindow,
                            scanFilter: ScanFilter,
                            mzAnalyzer: MZAnalyzerType.TOF,
                            totalIonCurrent: -1, // Will be set later
                            injectionTime: FrameProxyFactory.GetInjectionTimeSum(frameList.First(), frameList.Last()),
                            noiseData: null,
                            nativeId: "frame=" + frameList.First().ToString() + "-" + frameList.Last().ToString() +
                                ";scans=" + scanStart.ToString() + "-" + scanEnd.ToString(),
                            frameId: frameList.First(),
                            scanNumberStart: scanStart,
                            scanNumberEnd: scanEnd,
                            medianOneOverK0: -1, // Needs to be set later
                            precursorId: precursorId,
                            selectedIonMz: mostAbundantPrecursorPeak,
                            selectedIonChargeStateGuess: charge,
                            selectedIonIntensity: precursorIntensity,
                            isolationMZ: isolationMz,
                            isolationWidth: isolationWidth,
                            dissociationType: DissociationType.CID, // I think? Not totally sure about this
                                                                    // We don't really have one based precursors in timsTOF data, so it's not immediately clear how to handle this field
                                                                    //oneBasedPrecursorScanNumber: PrecursorIdToZeroBasedScanIndex.TryGetValue(precursorId, out int precursorScanIdx)
                                                                    //    ? precursorScanIdx
                                                                    //    : null,
                            oneBasedPrecursorScanNumber: null,
                            selectedIonMonoisotopicGuessMz: precursorMonoisotopicMz,
                            hcdEnergy: collisionEnergy.ToString(),
                            frames: frameList.ToList());
                    pasefScans.Add(dataScan);
                }
            }


            // For scan 1, we have 8 unique precursors, each of which is sampled 10-11 times
            // We need a way of iteratively building an mzSpectrum, 
            foreach (long frameId in allFrames)
            {
                FrameProxy frame = FrameProxyFactory.GetFrameProxy(frameId);
                //Iterate through all the datascans created above with this frame
                //foreach (TimsDataScan scan in pasefScans)
                Parallel.ForEach(pasefScans, scan =>
                {
                    if (scan.FrameIds.Contains(frameId))
                    {
                        List<double[]> mzArrays = new();
                        List<int[]> intensityArrays = new();
                        for (int mobilityScanIdx = scan.ScanNumberStart; mobilityScanIdx < scan.ScanNumberEnd; mobilityScanIdx++)
                        {
                            mzArrays.Add(frame.GetScanMzs(mobilityScanIdx));
                            intensityArrays.Add(frame.GetScanIntensities(mobilityScanIdx));
                        }
                        // Step 2: Average those suckers
                        ListNode<TofPeak> spectrumHeadNode = SumScansToLinkedList(mzArrays, intensityArrays, out int listLength);
                        scan.AddComponentSpectrum(spectrumHeadNode, listLength);
                    }
                }); 
            }

            //foreach (TimsDataScan scan in pasefScans)
            Parallel.ForEach(pasefScans, scan =>
            {
                scan.AverageComponentSpectra(filteringParams);
            });

            scanList.AddRange(pasefScans);
        }

        internal void BuildMs1ScanNoPrecursor(List<TimsDataScan> scanList, FrameProxy frame, FilteringParams filteringParams)
        {
            List<double[]> mzArrays = new();
            List<int[]> intensityArrays = new();
            // I don't know if scans are zero-indexed or one-based and at this point I'm too afraid to ask
            for (int scan = 1; scan < frame.NumberOfScans; scan++)
            {
                mzArrays.Add(frame.GetScanMzs(scan));
                intensityArrays.Add(frame.GetScanIntensities(scan));
            }
            // Step 2: Average those suckers
            MzSpectrum averagedSpectrum = SumScans(mzArrays, intensityArrays, filteringParams);
            // Step 3: Make an MsDataScan bby
            var dataScan = new TimsDataScan(
                massSpectrum: averagedSpectrum,
                oneBasedScanNumber: scanList.Count + 1,
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
                    ";scans=1-" + (frame.NumberOfScans - 1).ToString() +
                    ";precursor=NULL",
                frameId: frame.FrameId,
                scanNumberStart: 1,
                scanNumberEnd: frame.NumberOfScans,
                medianOneOverK0: frame.GetOneOverK0(frame.NumberOfScans/2),
                precursorId: null);

            scanList.Add(dataScan);
        }

        internal MzSpectrum SumScans(List<double[]> mzArrays, List<int[]> intensityArrays, FilteringParams filteringParams)
        {
            return TofSpectraMerger.MergesMs1Spectra(mzArrays, intensityArrays, filteringParams: filteringParams);
        }

        internal ListNode<TofPeak> SumScansToLinkedList(List<double[]> mzArrays, List<int[]> intensityArrays, out int listLength)
        {
            return TofSpectraMerger.MergeSpectraToLinkedList(mzArrays, intensityArrays, out listLength);
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
        [DllImport("Bruker/TimsTofReader/timsdata.dll")]
        public static extern UInt64 tims_open
              (byte[] analysis_directory_name_utf8, UInt32 use_recalibrated_state);

        /// <summary>
        /// Closes a file connection to a .tdf binary file
        /// </summary>
        [DllImport("Bruker/TimsTofReader/timsdata.dll")]
        public static extern void tims_close
              (UInt64 fileHandle);

        #endregion Bruker Dll Functions

    }
}
