namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// One file found by walking a PRIDE project's FTP directory tree — the COMPLETE listing. PRIDE's
    /// REST manifest (<see cref="PrideArchiveClient.GetProjectFilesAsync"/>) is knowingly incomplete:
    /// for PXD000001 it returns 8 files while the FTP tree holds 13 — it omits five, including the two
    /// largest. When a caller must know everything a project actually contains — or its true size —
    /// this is the authoritative list.
    /// </summary>
    public sealed class PrideFtpFile
    {
        /// <param name="relativePath">The path relative to the project's FTP root (see <see cref="RelativePath"/>).</param>
        /// <param name="approximateSizeBytes">PRIDE's rounded index size in bytes (see <see cref="ApproximateSizeBytes"/>).</param>
        /// <param name="url">The file's HTTPS download URL.</param>
        public PrideFtpFile(string relativePath, long approximateSizeBytes, string url)
        {
            RelativePath = relativePath;
            ApproximateSizeBytes = approximateSizeBytes;
            Url = url;
        }

        /// <summary>
        /// The path relative to the project's FTP root, e.g. "run1.raw" or, for a file in a
        /// subdirectory, "generated/summary.mztab".
        /// </summary>
        public string RelativePath { get; }

        /// <summary>The bare file name — the last segment of <see cref="RelativePath"/>.</summary>
        public string FileName => RelativePath.Substring(RelativePath.LastIndexOf('/') + 1);

        /// <summary>
        /// The file's size as PRIDE's FTP directory index reports it, ROUNDED to about three
        /// significant figures (its "429M" becomes 449 839 104). It is good for a project-size
        /// estimate but is NOT exact — deliberately named to say so. For the precise transfer size of a
        /// single file, issue an HTTP HEAD against <see cref="Url"/> and read its Content-Length.
        /// </summary>
        public long ApproximateSizeBytes { get; }

        /// <summary>The HTTPS URL the file can be downloaded from.</summary>
        public string Url { get; }
    }
}
