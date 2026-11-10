using System;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// One entry in a MassIVE dataset's file listing, as <see cref="MassiveClient.GetDatasetFilesAsync"/> reads it from
    /// MassIVE's FTP server: a file, or a symbolic link the server will not follow.
    /// </summary>
    /// <remarks>
    /// The listing is MassIVE's FTP server's own (decision D3, 2026-10-18), not the GNPS2 dataset cache, whose
    /// <c>collection</c> column is unreliable (MSV000080763's <c>ccms_peak/RAW/</c> copies are filed under <c>RAW</c>).
    /// <para>
    /// A dataset update is a folder <c>updates/&lt;YYYY-MM-DD&gt;_&lt;user&gt;_&lt;id&gt;/</c> holding collections of
    /// its own (measured 2026-11-10 on MSV000090383 and MSV000090777), so <see cref="UpdateName"/> and
    /// <see cref="Collection"/> are read from the path alone.
    /// </para>
    /// </remarks>
    public sealed class MassiveFile
    {
        /// <summary>The folder a dataset's updates are kept in, at the dataset's root.</summary>
        internal const string UpdatesFolder = "updates";

        /// <param name="datasetAccession">The dataset's MassIVE accession, e.g. "MSV000096226".</param>
        /// <param name="datasetRoot">The dataset's folder on MassIVE's FTP server, e.g. "/v07/MSV000096226".</param>
        /// <param name="path">The entry's path below the dataset's folder, '/'-separated, e.g. "raw/Heart/x.raw".</param>
        /// <param name="sizeBytes">The size the listing gives; null for a link, which has none worth reporting.</param>
        /// <param name="modified">The date the listing gives, or null when it gives none.</param>
        /// <param name="isLink">Whether the entry is a symbolic link rather than a file.</param>
        /// <exception cref="ArgumentException">An argument is blank, the path or root is not '/'-relative as described, or a file has no size.</exception>
        public MassiveFile(string datasetAccession, string datasetRoot, string path, long? sizeBytes,
            DateTime? modified = null, bool isLink = false)
        {
            if (string.IsNullOrWhiteSpace(datasetAccession))
                throw new ArgumentException("A MassIVE dataset accession is required.", nameof(datasetAccession));
            if (string.IsNullOrWhiteSpace(datasetRoot) || !datasetRoot.StartsWith("/", StringComparison.Ordinal) || datasetRoot.EndsWith("/", StringComparison.Ordinal))
                throw new ArgumentException("The dataset's folder must be an absolute FTP path with no trailing '/', e.g. \"/v07/MSV000096226\".", nameof(datasetRoot));
            if (string.IsNullOrWhiteSpace(path) || path.StartsWith("/", StringComparison.Ordinal) || path.EndsWith("/", StringComparison.Ordinal))
                throw new ArgumentException("The path must be relative to the dataset's folder, e.g. \"raw/x.raw\".", nameof(path));
            if (!isLink && sizeBytes is not >= 0)
                throw new ArgumentException("A file needs its size from the listing.", nameof(sizeBytes));

            DatasetAccession = datasetAccession;
            DatasetRoot = datasetRoot;
            Path = path;
            SizeBytes = isLink ? null : sizeBytes;
            Modified = modified;
            IsLink = isLink;

            // The folders the entry sits in. A link counts as a folder itself: ccms_peak, a collection, is one.
            string[] segments = path.Split('/');
            int folders = isLink ? segments.Length : segments.Length - 1;
            bool inUpdate = folders >= 2 && segments[0] == UpdatesFolder;
            UpdateName = inUpdate ? segments[1] : null;
            int collectionAt = inUpdate ? 2 : 0;
            Collection = folders > collectionAt ? segments[collectionAt] : null;
        }

        /// <summary>The dataset's MassIVE accession, e.g. "MSV000096226".</summary>
        public string DatasetAccession { get; }

        /// <summary>The dataset's folder on MassIVE's FTP server, e.g. "/v07/MSV000096226".</summary>
        public string DatasetRoot { get; }

        /// <summary>
        /// The path below the dataset's folder, '/'-separated, e.g. "raw/Heart/x.raw". A download keeps it below the
        /// destination folder, so "raw/a/x.raw" and "raw/b/x.raw" never overwrite each other.
        /// </summary>
        public string Path { get; }

        /// <summary>The path on MassIVE's FTP server: <see cref="DatasetRoot"/>, '/', <see cref="Path"/>.</summary>
        public string ServerPath => DatasetRoot + "/" + Path;

        /// <summary>The leaf name, the last segment of <see cref="Path"/>.</summary>
        public string FileName => Path.Substring(Path.LastIndexOf('/') + 1);

        /// <summary>
        /// The MassIVE collection the entry is in: its first folder ("raw", "peak", "search", ...), or, inside an
        /// update, the first folder below the update's own. A link is its own collection ("ccms_peak"). Null for a
        /// file with no folder above it (a file at the dataset's root, or directly in an update's folder).
        /// </summary>
        public string Collection { get; }

        /// <summary>
        /// The update the entry was added in (the folder name below "updates/", e.g. "2022-09-27_RGHuffman_25d0cd1a"),
        /// or null when it is not in an update. Whether an update replaces an earlier file is not visible from the
        /// listing.
        /// </summary>
        public string UpdateName { get; }

        /// <summary>
        /// The exact size in bytes the FTP listing gives, known before any download, and the size every download is
        /// checked against. Null only for a link.
        /// </summary>
        public long? SizeBytes { get; }

        /// <summary>
        /// The date the FTP listing shows, in the server's time: to the minute for an entry changed in the last six
        /// months, otherwise to the day. Null when the listing gives none. Downloads do not rely on it; they compare
        /// the server's own MDTM.
        /// </summary>
        public DateTime? Modified { get; }

        /// <summary>
        /// Whether the entry is a symbolic link. MassIVE's FTP server will not open one: its <c>ccms_peak</c> folder
        /// is a link on every dataset checked (2026-11-10), and listing or fetching through it answers
        /// "permission denied". A link is listed so a caller can tell "MassIVE holds files that cannot be listed" from
        /// "MassIVE holds nothing there"; <see cref="MassiveClient.DownloadFileAsync"/> refuses it.
        /// </summary>
        public bool IsLink { get; }

        /// <inheritdoc />
        public override string ToString() => $"{DatasetAccession}:{Path}";
    }
}
