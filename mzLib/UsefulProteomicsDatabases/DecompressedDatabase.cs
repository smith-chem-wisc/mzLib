using System;
using System.IO;
using System.IO.Compression;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// A gzipped database decompressed to a sibling of the input and deleted when disposed. The name
    /// is unique per instance: a fixed one collides when two processes load the same database at once,
    /// because File.Create takes FileShare.None while the reader holds FileShare.Read.
    /// </summary>
    internal sealed class DecompressedDatabase : IDisposable
    {
        private readonly string decompressed;

        private DecompressedDatabase(string location, string decompressed)
        {
            Location = location;
            this.decompressed = decompressed;
        }

        /// <summary>The path to read from: the input itself unless it was gzipped.</summary>
        internal string Location { get; }

        /// <summary>
        /// A sibling rather than the system temp directory: the caller already needs the input's
        /// directory writable and large enough for the whole database, so this adds no requirement.
        /// The extension is cosmetic, and only keeps the copy recognisable while it exists.
        /// </summary>
        internal static DecompressedDatabase For(string dbLocation, string extension)
        {
            if (!dbLocation.EndsWith(".gz"))
            {
                return new DecompressedDatabase(dbLocation, null);
            }

            //we had trouble decompressing and streaming on the fly so we decompress completely first, then stream the file, then delete the decompressed file
            string copy = Path.Combine(Path.GetDirectoryName(dbLocation), $"temp_{Guid.NewGuid():N}{extension}");
            var database = new DecompressedDatabase(copy, copy);
            try
            {
                using var stream = new FileStream(dbLocation, FileMode.Open, FileAccess.Read, FileShare.Read);
                using FileStream outputFileStream = File.Create(copy);
                using var decompressor = new GZipStream(stream, CompressionMode.Decompress);
                decompressor.CopyTo(outputFileStream);
            }
            catch
            {
                //a corrupt gz throws part way through the copy, and a unique name would otherwise
                //orphan one partial file per attempt rather than reusing a single one
                database.Dispose();
                throw;
            }

            return database;
        }

        public void Dispose()
        {
            if (decompressed != null)
            {
                File.Delete(decompressed);
            }
        }
    }
}
