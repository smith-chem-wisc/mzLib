using System;
using System.Collections.Generic;
using System.Linq;
using MzLibUtil;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// Pure, offline helpers over the PRIDE Archive manifest: resolving a <see cref="PrideArchiveFile"/> to an
    /// HTTPS download URL, and filtering / summarizing a manifest (a <see cref="IEnumerable{T}"/> of
    /// <see cref="PrideArchiveFile"/>). No network access — the download itself lives on
    /// <see cref="PrideArchiveClient"/>.
    /// </summary>
    public static class PrideArchiveExtensions
    {
        /// <summary>The stable CV accession of a PRIDE file's FTP location ("PRIDE:0000469").</summary>
        public const string FtpLocationAccession = "PRIDE:0000469";

        /// <summary>The stable CV accession of a PRIDE file's Aspera location ("PRIDE:0000468").</summary>
        public const string AsperaLocationAccession = "PRIDE:0000468";

        /// <summary>The PRIDE FTP host ("ftp.pride.ebi.ac.uk") known to also serve the identical path over HTTPS.</summary>
        public const string PrideFtpHost = "ftp.pride.ebi.ac.uk";

        /// <summary>
        /// Resolves the best HTTPS URL from which <paramref name="file"/> can be downloaded, if one exists.
        /// PRIDE files advertise FTP and Aspera locations (rarely an explicit HTTPS one); its FTP host
        /// (<c>ftp.pride.ebi.ac.uk</c>) serves the identical path over HTTPS, so an FTP location on that host
        /// is scheme-upgraded from <c>ftp://</c> to <c>https://</c>. The upgrade is restricted to that host,
        /// since it is the only one known to serve the same path over HTTPS. Preference order: an explicit
        /// HTTPS location, then the FTP location upgraded to HTTPS. An Aspera-only file cannot be reached over HTTPS.
        /// </summary>
        /// <param name="file">The file whose <see cref="PrideArchiveFile.PublicFileLocations"/> are examined.</param>
        /// <param name="httpsUrl">On success, the resolved HTTPS URL; otherwise null.</param>
        /// <returns>True if an HTTPS URL was resolved; otherwise false.</returns>
        /// <exception cref="ArgumentNullException">The file is null.</exception>
        public static bool TryGetHttpsDownloadUrl(this PrideArchiveFile file, out string httpsUrl)
        {
            if (file == null)
                throw new ArgumentNullException(nameof(file));

            httpsUrl = null;
            IReadOnlyList<CvParam> locations = file.PublicFileLocations;
            if (locations == null || locations.Count == 0)
                return false;

            // 1) An explicit HTTPS location, should the repository ever provide one, is used as-is.
            foreach (CvParam location in locations)
            {
                if (location?.Value != null && location.Value.StartsWith("https://", StringComparison.OrdinalIgnoreCase))
                {
                    httpsUrl = location.Value;
                    return true;
                }
            }

            // 2) Otherwise upgrade a PRIDE FTP location to HTTPS. Only ftp.pride.ebi.ac.uk is known to serve
            //    the identical path over HTTPS, so the scheme-upgrade is restricted to that host. Selecting by
            //    a usable ftp:// value (not by accession alone) means a mislabeled or empty FTP-accession entry
            //    cannot shadow a well-formed ftp:// URL elsewhere in the list.
            foreach (CvParam location in locations)
            {
                string value = location?.Value;
                if (value == null || !value.StartsWith("ftp://", StringComparison.OrdinalIgnoreCase))
                    continue;
                if (!Uri.TryCreate(value, UriKind.Absolute, out Uri ftpUri)
                    || !string.Equals(ftpUri.Host, PrideFtpHost, StringComparison.OrdinalIgnoreCase))
                    continue;

                httpsUrl = "https://" + value.Substring("ftp://".Length);
                return true;
            }

            return false;
        }

        /// <summary>
        /// Returns the HTTPS URL from which <paramref name="file"/> can be downloaded, upgrading the FTP
        /// location if necessary (see <see cref="TryGetHttpsDownloadUrl"/>).
        /// </summary>
        /// <exception cref="ArgumentNullException">The file is null.</exception>
        /// <exception cref="NotSupportedException">The file exposes no HTTPS-reachable location (e.g. Aspera-only).</exception>
        public static string GetHttpsDownloadUrl(this PrideArchiveFile file)
        {
            if (!file.TryGetHttpsDownloadUrl(out string httpsUrl))
                throw new NotSupportedException(
                    $"PRIDE file '{file.FileName}' exposes no HTTPS-reachable download location (only Aspera, or none).");
            return httpsUrl;
        }

        /// <summary>
        /// Filters a manifest to files whose <see cref="PrideArchiveFile.FileCategory"/> value equals
        /// <paramref name="category"/> (case-insensitive), e.g. "RAW", "SEARCH", "PEAK", "RESULT".
        /// </summary>
        /// <exception cref="ArgumentNullException">The sequence is null.</exception>
        /// <exception cref="ArgumentException">The category is blank.</exception>
        public static IEnumerable<PrideArchiveFile> WhereCategory(this IEnumerable<PrideArchiveFile> files, string category)
        {
            if (files == null)
                throw new ArgumentNullException(nameof(files));
            if (string.IsNullOrWhiteSpace(category))
                throw new ArgumentException("A file category is required.", nameof(category));

            return files.Where(f => f?.FileCategory != null
                && string.Equals(f.FileCategory.Value, category, StringComparison.OrdinalIgnoreCase));
        }

        /// <summary>
        /// Filters a manifest to files whose name ends with any of <paramref name="extensions"/> (case-insensitive).
        /// A leading dot is optional: "raw" and ".raw" are equivalent.
        /// </summary>
        /// <exception cref="ArgumentNullException">The sequence is null.</exception>
        /// <exception cref="ArgumentException">No extension is supplied.</exception>
        public static IEnumerable<PrideArchiveFile> WhereExtension(this IEnumerable<PrideArchiveFile> files, params string[] extensions)
        {
            if (files == null)
                throw new ArgumentNullException(nameof(files));
            if (extensions == null || extensions.Length == 0)
                throw new ArgumentException("At least one file extension is required.", nameof(extensions));

            string[] normalized = extensions
                .Where(e => !string.IsNullOrWhiteSpace(e))
                .Select(e => e.StartsWith(".") ? e : "." + e)
                .ToArray();
            if (normalized.Length == 0)
                throw new ArgumentException("At least one non-blank file extension is required.", nameof(extensions));

            return files.Where(f => f?.FileName != null
                && normalized.Any(ext => f.FileName.EndsWith(ext, StringComparison.OrdinalIgnoreCase)));
        }

        /// <summary>The total size in bytes of the files in a manifest (nulls ignored).</summary>
        /// <exception cref="ArgumentNullException">The sequence is null.</exception>
        public static long TotalSizeBytes(this IEnumerable<PrideArchiveFile> files)
        {
            if (files == null)
                throw new ArgumentNullException(nameof(files));

            return files.Where(f => f != null).Sum(f => f.FileSizeBytes);
        }
    }
}
