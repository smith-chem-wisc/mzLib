using System;
using System.Collections.Generic;
using MzLibUtil;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// One file belonging to a PRIDE Archive project, as returned by the PRIDE Archive REST API
    /// (v3 <c>projects/{accession}/files</c>). A plain data object populated by JSON deserialization;
    /// controlled-vocabulary fields are modeled with <see cref="CvParam"/>.
    /// </summary>
    /// <remarks>
    /// To obtain a download location, match <see cref="PublicFileLocations"/> on the protocol's stable
    /// CV accession — e.g. "PRIDE:0000469" for FTP, "PRIDE:0000468" for Aspera — rather than on the
    /// display name. Not every file exposes an HTTPS location.
    /// </remarks>
    public class PrideArchiveFile
    {
        /// <summary>The file name (e.g. "run1.raw").</summary>
        public string FileName { get; set; } = string.Empty;

        /// <summary>The file size in bytes.</summary>
        public long FileSizeBytes { get; set; }

        /// <summary>The file checksum, if the repository provides one; otherwise empty.</summary>
        public string Checksum { get; set; } = string.Empty;

        /// <summary>The file category as a controlled-vocabulary term (e.g. value "RAW", "SEARCH", "PEAK").</summary>
        public CvParam FileCategory { get; set; } = new();

        /// <summary>
        /// The public download locations as controlled-vocabulary terms; each term's
        /// <see cref="CvParam.Value"/> is a URL and its <see cref="CvParam.Accession"/> identifies the
        /// protocol (FTP "PRIDE:0000469", Aspera "PRIDE:0000468").
        /// </summary>
        public List<CvParam> PublicFileLocations { get; set; } = new();

        /// <summary>The date the file was originally submitted.</summary>
        public DateTimeOffset SubmissionDate { get; set; }

        /// <summary>The date the file was made public.</summary>
        public DateTimeOffset PublicationDate { get; set; }

        /// <summary>The date the file was last updated.</summary>
        public DateTimeOffset UpdatedDate { get; set; }
    }
}
