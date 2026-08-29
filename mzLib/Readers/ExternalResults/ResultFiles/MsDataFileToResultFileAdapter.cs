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
        private static readonly string MgfExtension = SupportedFileType.Mgf.GetFileExtension();

        /// <summary>
        /// Extensions of the spectra formats this adapter reads but cannot write. Mgf is deliberately
        /// absent: #1165 added <c>ExportAsMgf</c>, so a caller asking for one now gets one.
        ///
        /// Named by enum member
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
        /// <summary>
        /// Writes the spectra in the format the path's extension names.
        /// </summary>
        /// <remarks>
        /// The extension chooses the writer rather than being overruled by it, which is the whole point
        /// of #1181: the old guard asked "is this already .mzML" instead of "what did the caller ask
        /// for", so <c>WriteResults("output.mgf")</c> produced <c>output.mgf.mzML</c> -- a file wearing
        /// two extensions, in a format nobody chose.
        ///
        /// Three outcomes, and the middle one is the reason this is not simply
        /// <c>Path.ChangeExtension</c>:
        ///
        ///   * a format this adapter can write (mzML, and MGF since #1165) is written as asked;
        ///   * a spectra format it can only READ is refused, because silently handing back mzML under
        ///     the caller's chosen name is worse than saying no; and
        ///   * anything else has .mzML appended, which keeps a bare path working and also keeps a
        ///     dotted sample name working -- Path.GetExtension("gradient_1.5uL") is ".5uL", and that is
        ///     a naming habit, not a format request.
        /// </remarks>
        public void WriteResults(string outputPath)
        {
            string extension = Path.GetExtension(outputPath);

            if (extension.Equals(MzMlExtension, StringComparison.InvariantCultureIgnoreCase))
            {
                _dataFile.ExportAsMzML(outputPath, true);
                return;
            }

            if (extension.Equals(MgfExtension, StringComparison.InvariantCultureIgnoreCase))
            {
                _dataFile.ExportAsMgf(outputPath);
                return;
            }

            if (UnwritableSpectraExtensions.Contains(extension))
                throw new ArgumentException(
                    $"This writer handles {MzMlExtension} and {MgfExtension}, so it cannot write " +
                    $"'{Path.GetFileName(outputPath)}'. Pass a path ending in one of those, or one with " +
                    $"no extension and {MzMlExtension} will be appended.",
                    nameof(outputPath));

            _dataFile.ExportAsMzML(outputPath + MzMlExtension, true);
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
