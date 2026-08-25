using MassSpectrometry;

namespace Readers
{
    /// <summary>
    /// Adapter class to allow for MsDataFile loading structure to integrate into all file reading
    /// </summary>
    public sealed class MsDataFileToResultFileAdapter : MsDataFile, IResultFile
    {
        private MsDataFile _dataFile;

        private static readonly string MzMlExtension = SupportedFileType.MzML.GetFileExtension();

        /// <summary>
        /// Extensions of the spectra formats this adapter reads but cannot write. Named by enum member
        /// so the list stays greppable against the entries that resolve to this type in
        /// <see cref="SupportedFileTypeExtensions.GetResultFileType(SupportedFileType)"/>, but the
        /// extension strings themselves are pulled from
        /// <see cref="SupportedFileTypeExtensions.GetFileExtension"/> so no literal here can go stale.
        /// Path.GetExtension is what reduces the compound suffixes ("_ms1.msalign") to the extension a
        /// caller's path would actually end in (".msalign").
        /// </summary>
        private static readonly HashSet<string> UnwritableSpectraExtensions =
            new[]
            {
                SupportedFileType.ThermoRaw,
                SupportedFileType.Mgf,
                SupportedFileType.BrukerD,
                SupportedFileType.BrukerTimsTof,
                SupportedFileType.Ms1Align,
                SupportedFileType.Ms2Align,
            }
            .Select(type => Path.GetExtension(type.GetFileExtension()))
            .ToHashSet(StringComparer.InvariantCultureIgnoreCase);


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

        /// <summary>
        /// Writes the loaded scans to <paramref name="outputPath"/> as mzML, whatever format was read.
        /// </summary>
        /// <remarks>
        /// This adapter is a converter -- it reads every spectra format that maps to it in
        /// <see cref="SupportedFileTypeExtensions.GetResultFileType(SupportedFileType)"/> and writes
        /// only mzML -- so the extension it is handed and the format it emits can disagree. A path
        /// with no extension, or one that is not a spectra format, still gets ".mzML" appended; that
        /// is the long-standing behaviour and the reason a dotted sample name such as
        /// "gradient_1.5uL" is not mistaken for a format request. A path that names a *different*
        /// spectra format is refused instead: "output.mgf" used to become "output.mgf.mzML", an mzML
        /// file wearing two extensions, neither of which the caller chose. Refusing states what the
        /// silent rename only implied, and a caller who genuinely wants MGF can now call
        /// <see cref="MsDataFileExtensions.ExportAsMgf"/> directly.
        /// </remarks>
        /// <param name="outputPath">
        /// Destination path. Must end in ".mzML", or carry no spectra-format extension, in which case
        /// ".mzML" is appended.
        /// </param>
        /// <exception cref="ArgumentException">
        /// <paramref name="outputPath"/> names a spectra format this adapter cannot write.
        /// </exception>
        public void WriteResults(string outputPath)
        {
            if (!Path.GetExtension(outputPath).Equals(MzMlExtension, StringComparison.InvariantCultureIgnoreCase))
            {
                if (UnwritableSpectraExtensions.Contains(Path.GetExtension(outputPath)))
                    throw new ArgumentException(
                        $"This is an mzML writer, so it cannot write '{Path.GetFileName(outputPath)}'. Pass a path " +
                        $"ending in {MzMlExtension}, or one with no extension and {MzMlExtension} will be appended.",
                        nameof(outputPath));

                outputPath += MzMlExtension;
            }

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
