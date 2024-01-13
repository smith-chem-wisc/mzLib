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

namespace Readers.Bruker
{ 
    public class TimsTofFileReader : MsDataFile
    {
        // timsTOF instruments collect frames, packets of ions collected by the tims, then analyzed 
        // over multiple scans with each scan corresponding to the same retention time but different
        // ion mobility valuess. When reading the file, multiple scans from the same frame are collapsed into 

        public TimsTofFileReader(string filePath) : base (filePath) 
        {
            PrecursorIdToZeroBasedScanIndex = new();
        }

        private UInt64? _fileHandle;
        private Object _fileLock;
        private SQLiteConnection? _sqlConnection;
        public int NumberOfFrames { get; private set; }
        public List<long> Ms1FrameIds { get; private set; }
        internal FrameProxyFactory? FrameProxyFactory { get; private set; }
        
        // I don't know what the default scan range is, and at this point I'm too afraid to ask...
        private MzRange? _scanWindow;
        public MzRange ScanWindow => _scanWindow ??= new MzRange(20, 2000);
        public const string ScanFilter = "f";

        public ConcurrentDictionary<int, int> PrecursorIdToZeroBasedScanIndex { get; }

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

        internal void OpenSqlConnection()
        {
            // You really only need to do this if you're doing a static load
            // If dynamic, this is a worse way to do it
            //_sqlConnection = new SQLiteConnection("Data Source=:memory:;Mode=Memory");
            //_sqlConnection.Open();

            //SQLiteConnection onDiskDBConnection = new SQLiteConnection("Data Source=" + 
            //    Path.Combine(FilePath, "analysis.tdf") +
            //    "; Version=3");
            //onDiskDBConnection.Open();
            //onDiskDBConnection.BackupDatabase(_sqlConnection, "main", "main", -1, null, -1);
            //onDiskDBConnection.Close();

            _sqlConnection = new SQLiteConnection("Data Source=" +
            Path.Combine(FilePath, "analysis.tdf") +
            "; Version=3");
            _sqlConnection.Open();
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
            BuildProxyFactory();

            ConcurrentDictionary<int, List<MsDataScan>> scanListDictionary = new();


            //Parallel.ForEach(Partitioner.Create(0, Ms1FrameIds.Count),
            Parallel.ForEach(Partitioner.Create(0, 100),
                    new ParallelOptions { MaxDegreeOfParallelism = maxThreads }, (range, loopState) =>
                    {
                        //Debugger.Break();
                        List<MsDataScan> scansInRange = new();
                        for (int i = range.Item1; i < range.Item2; i++)
                        {
                            long frameId = Ms1FrameIds[i];
                            //Debugger.Break();
                            
                            BuildMS1Scans(scansInRange, frameId);   
                        }
                        //Debugger.Break();
                        scanListDictionary[range.Item1] = scansInRange;
                    });

            //for (int i = 0; i < NumberOfFrames; i++)
            //{
            //    long oneBasedFrameNumber = FramesTable.OneBasedFrameIndex[i];
            //    FrameProxy currentFrame = new FrameProxy((ulong)_fileHandle, oneBasedFrameNumber, FramesTable.NumScans[i]);

            //    if (FramesTable.MsMsType[i].ToEnum<TimsTofMsMsType>(out var msMsType))
            //    {
            //        switch(msMsType)
            //        {
            //            case TimsTofMsMsType.MS:
            //                BuildMS1Scans(scanList, currentFrame, i);
            //                break;
            //            case TimsTofMsMsType.PASEF:
            //                //BuildPasefScans(scanList, currentFrame, i);
            //                break;
            //            default:
            //                throw new NotImplementedException("Only PASEF data is currently supported");
            //        }

            //    }
            //}
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

        internal void BuildMS1Scans(List<MsDataScan> scanList, long frameId)
        {
            FrameProxy frame = FrameProxyFactory.GetFrameProxy(frameId);
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
                List<int> precursorIdsInFrame = new();
                while (sqliteReader.Read())
                {
                    var scanStart = sqliteReader.GetInt32(0);
                    var scanEnd = sqliteReader.GetInt32(1);
                    var scanMedian = sqliteReader.GetFloat(2);
                    int precursorId = sqliteReader.GetInt32(3);
                    precursorIdsInFrame.Add(precursorId);
                    // Now we're gonna try and build some scans!
                    // Step 1: Pull the m/z + intensity arrays for each scan in range
                    List<double[]> mzArrays = new();
                    List<int[]> intensityArrays = new();
                    for (int scan = scanStart; scan < scanEnd; scan++)
                    {
                        mzArrays.Add(frame.GetScanMzs(scan));
                        intensityArrays.Add(frame.GetScanIntensities(scan));
                    }
                    // Step 2: Average those suckers
                    MzSpectrum averagedSpectrum = SumScans(mzArrays, intensityArrays);
                    // Step 3: Make an MsDataScan bby
                    var dataScan = new TimsDataScan(
                        massSpectrum: averagedSpectrum,
                        oneBasedScanNumber: scanList.Count + 1,
                        msnOrder: 1,
                        isCentroid: true,
                        polarity: FrameProxyFactory.GetPolarity(frameId),
                        retentionTime: FrameProxyFactory.GetRetentionTime(frameId),
                        scanWindowRange: ScanWindow,
                        scanFilter: ScanFilter,
                        mzAnalyzer: MZAnalyzerType.TOF,
                        totalIonCurrent: intensityArrays.Sum(array => array.Sum()),
                        injectionTime: FrameProxyFactory.GetInjectionTime(frameId),
                        noiseData: null,
                        nativeId: "frame=" + frame.FrameId.ToString() + ";scans=" + scanStart.ToString() + "-" + scanEnd.ToString(),
                        frameId: frame.FrameId,
                        scanNumberStart: scanStart,
                        scanNumberEnd: scanEnd,
                        medianOneOverK0: frame.GetOneOverK0((double)scanMedian),
                        precursorId: precursorId);

                    //PrecursorIdToZeroBasedScanIndex.TryAdd(precursorId, scanList.Count());
                    scanList.Add(dataScan);

                }
                // Then, build ONE MS2 scan from every PASEF frame that sampled that precursor
                BuildPasefScanFromPrecursor(scanList, precursorIdsInFrame.Distinct(), frameId);
            }
        }

        internal void BuildPasefScanFromPrecursor(List<MsDataScan> scanList, IEnumerable<int> precursorIds, long parentFrameId)
        {
            using var command = new SQLiteCommand(_sqlConnection);
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
            HashSet<long> allFrames = new();
            List<MsDataScan> pasefScans = new();
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
                var precursorMonoisotopicMz = sqliteReader.GetFloat(7);
                var charge = sqliteReader.GetInt32(8);
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

            // For scan 1, we have 8 unique precursors, each of which is sampled 10-11 times
            // We need a way of iteratively building an mzSpectrum, 
            foreach(long frameId in allFrames)
            {
                FrameProxy frame = FrameProxyFactory.GetFrameProxy(frameId);
                //Iterate through all the datascans created above with this frame
                foreach(TimsDataScan scan in pasefScans)
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
                        MzSpectrum summedSpectrum = SumScans(mzArrays, intensityArrays);
                        scan.AddComponentSpectrum(summedSpectrum);
                    }
                }
            }

            foreach(TimsDataScan scan in pasefScans)
            {
                scan.AverageComponentSpectra();
            }

            scanList.AddRange(pasefScans);
        }

        //internal void BuildPasefScans(List<MsDataScan> scanList, FrameProxy frame, int frameIndex)
        //{
        //    using var command = new SQLiteCommand(_sqlConnection);
        //    // SQL Command for getting some info from both PasefFrameMsMsInfo table and
        //    // Precursors table
        //    command.CommandText =
        //        @"SELECT m.ScanNumBegin, m.ScanNumEnd, m.IsolationMz, m.IsolationWidth," +
        //        " m.CollisionEnergy, p.LargestPeakMz, p.MonoisotopicMz, p.Charge, p.Intensity, p.ScanNumber, p.Id" +
        //        " FROM PasefFrameMsMsInfo m" +
        //        " INNER JOIN Precursors p on p.Id = m.Precursor" +
        //        " WHERE m.Frame = " + frame.FrameId.ToString() +
        //        ";";

        //    using var sqliteReader = command.ExecuteReader();
        //    while (sqliteReader.Read())
        //    {
        //        var scanStart = sqliteReader.GetInt32(0);
        //        var scanEnd = sqliteReader.GetInt32(1);
        //        var isolationMz = sqliteReader.GetFloat(2);
        //        var isolationWidth = sqliteReader.GetFloat(3);
        //        var collisionEnergy = sqliteReader.GetFloat(4);
        //        var mostAbundantPrecursorPeak = sqliteReader.GetFloat(5);
        //        var precursorMonoisotopicMz = sqliteReader.GetFloat(6);
        //        var charge = sqliteReader.GetInt32(7);
        //        var precursorIntensity = sqliteReader.GetFloat(8);
        //        var scanMedian = sqliteReader.GetFloat(9);
        //        var precursorId = sqliteReader.GetInt32(10);

        //        // Now we're gonna try and build some scans!
        //        // Step 1: Pull the m/z + intensity arrays for each scan in range
        //        List<double[]> mzArrays = new();
        //        List<int[]> intensityArrays = new();
        //        for (int scan = scanStart; scan < scanEnd; scan++)
        //        {
        //            mzArrays.Add(frame.GetScanMzs(scan));
        //            intensityArrays.Add(frame.GetScanIntensities(scan));
        //        }
        //        // Step 2: Average those suckers
        //        MzSpectrum averagedSpectrum = SumScans(mzArrays, intensityArrays);
        //        // Step 3: Make an MsDataScan bby
        //        var dataScan = new TimsDataScan(
        //                massSpectrum: averagedSpectrum,
        //                oneBasedScanNumber: scanList.Count + 1,
        //                msnOrder: 2,
        //                isCentroid: true,
        //                polarity: FramesTable.Polarity[frameIndex] == '+' ? Polarity.Positive : Polarity.Negative,
        //                retentionTime: (double)FramesTable.RetentionTime[frameIndex],
        //                scanWindowRange: ScanWindow,
        //                scanFilter: ScanFilter,
        //                mzAnalyzer: MZAnalyzerType.TOF,
        //                totalIonCurrent: intensityArrays.Sum(array => array.Sum()),
        //                injectionTime: FramesTable.FillTime[frameIndex],
        //                noiseData: null,
        //                nativeId: "frame=" + frame.FrameId.ToString() + ";scans=" + scanStart.ToString() + "-" + scanEnd.ToString(),
        //                frameId: frame.FrameId,
        //                scanNumberStart: scanStart,
        //                scanNumberEnd: scanEnd,
        //                medianOneOverK0: frame.GetOneOverK0((double)scanMedian),
        //                precursorId: precursorId,
        //                selectedIonMz: mostAbundantPrecursorPeak,
        //                selectedIonChargeStateGuess: charge,
        //                selectedIonIntensity: precursorIntensity,
        //                isolationMZ: isolationMz,
        //                isolationWidth: isolationWidth,
        //                dissociationType: DissociationType.CID, // I think? Not totally sure about this
        //                // We don't really have one based precursors in timsTOF data, so it's not immediately clear how to handle this field
        //                oneBasedPrecursorScanNumber: PrecursorIdToZeroBasedScanIndex.TryGetValue(precursorId, out int precursorScanIdx) 
        //                    ? precursorScanIdx 
        //                    : null,
        //                selectedIonMonoisotopicGuessMz: precursorMonoisotopicMz,
        //                hcdEnergy: collisionEnergy.ToString());

        //        scanList.Add(dataScan);
        //    }
        //}

        internal MzSpectrum SumScans(List<double[]> mzArrays, List<int[]> intensityArrays)
        {
            return TofSpectraMerger.MergeSpectra(mzArrays, intensityArrays);

            //int mostIntenseScanIndex = 0;
            //int maxIntensity = 0;
            //for (int i = 0; i < intensityArrays.Count; i++)
            //{
            //    int[] intensityArray = intensityArrays[i];
            //    int intensitySum = intensityArray.Sum();
            //    if(intensitySum > maxIntensity)
            //    {
            //        mostIntenseScanIndex = i;
            //        maxIntensity = intensitySum;
            //    }
            //}

            //return new MzSpectrum(
            //    mz: mzArrays[mostIntenseScanIndex], 
            //    intensities: Array.ConvertAll(intensityArrays[mostIntenseScanIndex], entry => (double)entry),
            //    shouldCopy: false);
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
