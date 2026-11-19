using System;
using System.IO;
using System.IO.Compression;
using System.Security.Cryptography;

namespace UsefulProteomicsDatabases.Ensembl
{
    /// <summary>
    /// The two steps every Ensembl reader here takes before parsing: hash the file as published, then
    /// open its text. The sha256 is of the bytes as read -- the compressed bytes for a .gz, which is what
    /// Ensembl publishes and checksums -- and is taken in its own pass, before any row is trusted.
    /// </summary>
    internal static class EnsemblFile
    {
        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        internal static string Sha256(string path, string description)
        {
            if (!File.Exists(path))
            {
                throw new FileNotFoundException($"{description} not found.", path);
            }

            using var stream = File.OpenRead(path);
            return Convert.ToHexString(SHA256.HashData(stream)).ToLowerInvariant();
        }

        /// <summary>Opens the file's text, decompressing when the name ends in .gz.</summary>
        internal static StreamReader OpenText(string path)
        {
            Stream file = File.OpenRead(path);
            Stream content = path.EndsWith(".gz", StringComparison.OrdinalIgnoreCase)
                ? new GZipStream(file, CompressionMode.Decompress)
                : file;
            return new StreamReader(content);
        }
    }
}
