﻿using MassSpectrometry;

namespace Readers
{
    /// <summary>
    /// Adapter class to allow for MsDataFile loading structure to integrate into all file reading
    /// </summary>
    public sealed class MsDataFileToResultFileAdapter : MsDataFile, IResultFile
    {
        private MsDataFile _dataFile;


        #region MsDataFile Members

        public override MsDataFile LoadAllStaticData(FilteringParams filteringParams = null, int maxThreads = 1) =>
            _dataFile.LoadAllStaticData(filteringParams, maxThreads);

        public override SourceFile GetSourceFile() => _dataFile.GetSourceFile();

        public override MsDataScan GetOneBasedScanFromDynamicConnection(int oneBasedScanNumber, IFilteringParams filterParams = null) =>
            _dataFile.GetOneBasedScanFromDynamicConnection(oneBasedScanNumber, filterParams);

        public override void CloseDynamicConnection() => _dataFile.CloseDynamicConnection();

        public override void InitiateDynamicConnection() => _dataFile.InitiateDynamicConnection();

        #endregion

        #region IResultFile Members

        public new string FilePath { get; set; }
        public SupportedFileType FileType => FilePath.ParseFileType();
        public Software Software { get; set; } = Software.MassSpecFile;
        public List<MsDataScan> Results { get; set; }
        public void LoadResults()
        {
            _dataFile = MsDataFileReader.GetDataFile(FilePath).LoadAllStaticData();
            Results = _dataFile.GetAllScansList();
            this.Scans = _dataFile.Scans;
            this.SourceFile = _dataFile.GetSourceFile();
        }

        public void WriteResults(string outputPath)
        {
            if (!outputPath.EndsWith(".mzml", StringComparison.InvariantCultureIgnoreCase))
                outputPath += SupportedFileType.MzML.GetFileExtension();
            _dataFile.ExportAsMzML(outputPath, true);
        }

        #endregion

        public MsDataFileToResultFileAdapter(string filePath) : base(filePath)
        {
            FilePath = filePath;
        }

        public MsDataFileToResultFileAdapter() : base("")
        {
        }
    }
}
