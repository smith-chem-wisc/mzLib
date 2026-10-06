namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// One row of a PRIDE project's checksum list (<see cref="PrideArchiveClient.GetFileChecksumsAsync"/>):
    /// the MD5 and exact size PRIDE recorded for a file the submitter deposited.
    /// </summary>
    /// <remarks>
    /// <see cref="SizeBytes"/> is the size of the bytes a download delivers, which is what makes it useful for
    /// checking one. <see cref="PrideArchiveFile.FileSizeBytes"/> is not: for a <c>.gz</c> file it can be more
    /// than ten times larger (PXD015239's <c>190712-VL-BR-LFQ-1.mzid.gz</c> is 21 027 335 bytes here and
    /// 294 548 296 there; measured 2026-09-24).
    /// </remarks>
    public sealed class PrideFileChecksum
    {
        /// <param name="fileName">The bare file name (see <see cref="FileName"/>).</param>
        /// <param name="md5">The MD5 as 32 hexadecimal characters (see <see cref="Md5"/>).</param>
        /// <param name="sizeBytes">The exact size in bytes (see <see cref="SizeBytes"/>).</param>
        public PrideFileChecksum(string fileName, string md5, long sizeBytes)
        {
            FileName = fileName;
            Md5 = md5;
            SizeBytes = sizeBytes;
        }

        /// <summary>The bare file name, as <see cref="PrideArchiveFile.FileName"/> spells it. Never a path.</summary>
        public string FileName { get; }

        /// <summary>The file's MD5 as 32 hexadecimal characters, as PRIDE wrote it (lower case in every list seen).</summary>
        public string Md5 { get; }

        /// <summary>The file's exact size in bytes: the length a download of it comes to.</summary>
        public long SizeBytes { get; }
    }
}
