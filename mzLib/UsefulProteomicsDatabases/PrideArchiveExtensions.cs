using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// Pure, offline helpers over PRIDE data: resolving a <see cref="PrideArchiveFile"/> to an HTTPS download
    /// URL, filtering / summarizing a manifest (a <see cref="IEnumerable{T}"/> of <see cref="PrideArchiveFile"/>),
    /// and reading a <see cref="PrideProxiSpectrum"/>'s controlled-vocabulary metadata into an mzLib
    /// <see cref="MsDataScan"/>. No network access — the fetch/download themselves live on
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

        // ---- PROXI spectrum metadata -------------------------------------------------------------------
        //
        // The PSI-MS accessions PROXI uses to describe a spectrum. Matched on accession, never on name (see
        // CvParam's own doc). Verified against a live PRIDE PROXI response for
        // mzspec:PXD000561:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555:VLHPLEGAVVIIFK/2.

        /// <summary>"MS:1000511" — ms level (1 for a survey scan, 2 for a fragment spectrum, ...).</summary>
        private const string MsLevelAccession = "MS:1000511";

        /// <summary>"MS:1003057" — the scan's number within its originating run.</summary>
        private const string ScanNumberAccession = "MS:1003057";

        /// <summary>"MS:1000041" — the precursor's charge state.</summary>
        private const string ChargeStateAccession = "MS:1000041";

        /// <summary>"MS:1000744" — the selected (precursor) ion's m/z.</summary>
        private const string SelectedIonMzAccession = "MS:1000744";

        /// <summary>"MS:1000042" — the selected ion's peak intensity.</summary>
        private const string SelectedIonIntensityAccession = "MS:1000042";

        /// <summary>"MS:1000827" — the isolation window's target m/z.</summary>
        private const string IsolationWindowTargetMzAccession = "MS:1000827";

        /// <summary>"MS:1000828" — the isolation window's lower offset from its target.</summary>
        private const string IsolationWindowLowerOffsetAccession = "MS:1000828";

        /// <summary>"MS:1000829" — the isolation window's upper offset from its target.</summary>
        private const string IsolationWindowUpperOffsetAccession = "MS:1000829";

        /// <summary>"MS:1000016" — the scan's start (retention) time.</summary>
        private const string ScanStartTimeAccession = "MS:1000016";

        /// <summary>"MS:1000927" — the ion injection time, in milliseconds.</summary>
        private const string IonInjectionTimeAccession = "MS:1000927";

        /// <summary>"MS:1000512" — the vendor filter string describing the scan.</summary>
        private const string FilterStringAccession = "MS:1000512";

        /// <summary>"MS:1000465" — scan polarity, whose value names the polarity ("positive scan"/"negative scan").</summary>
        private const string ScanPolarityAccession = "MS:1000465";

        /// <summary>"MS:1000130" — positive scan, when polarity is given as a bare term rather than a value.</summary>
        private const string PositiveScanAccession = "MS:1000130";

        /// <summary>"MS:1000129" — negative scan, when polarity is given as a bare term rather than a value.</summary>
        private const string NegativeScanAccession = "MS:1000129";

        /// <summary>"MS:1000525" — spectrum representation, whose value is "centroid spectrum" or "profile spectrum".</summary>
        private const string SpectrumRepresentationAccession = "MS:1000525";

        /// <summary>"MS:1000127" — centroid spectrum, when representation is given as a bare term.</summary>
        private const string CentroidSpectrumAccession = "MS:1000127";

        /// <summary>"UO:0000031" — the unit term for minutes.</summary>
        private const string MinuteUnitAccession = "UO:0000031";

        /// <summary>
        /// Reads a <see cref="PrideProxiSpectrum"/>'s controlled-vocabulary <see cref="PrideProxiSpectrum.Attributes"/>
        /// into a fully-populated mzLib <see cref="MsDataScan"/> — scan number, ms level, polarity, retention time,
        /// centroid-ness, precursor m/z and charge, isolation window, injection time and filter string — carrying
        /// the spectrum itself through as the scan's <see cref="MsDataScan.MassSpectrum"/> (a
        /// <see cref="PrideProxiSpectrum"/> already is an <see cref="MzSpectrum"/>). Pure and offline; the fetch
        /// lives on <see cref="PrideArchiveClient.GetProxiSpectrumAsync"/>.
        /// </summary>
        /// <remarks>
        /// PROXI describes far less than a vendor file does, so several <see cref="MsDataScan"/> fields are filled
        /// with the honest "not stated" value rather than a guess:
        /// <list type="bullet">
        /// <item><description><see cref="MZAnalyzerType.Unknown"/> — PROXI names the <i>instrument</i>
        /// ("MS:1000463"), not the mass analyzer, so there is no analyzer term to map. (An mzML-sourced analyzer
        /// accession would map through <c>Readers.Mzml.AnalyzerDictionary</c>, which is not reachable from this
        /// project: Readers references UsefulProteomicsDatabases, not the reverse.)</description></item>
        /// <item><description><see cref="DissociationType.Unknown"/> — PROXI carries no dissociation term; the
        /// method is recorded only inside the vendor filter string, which is too vendor-specific to parse
        /// reliably. The filter string is preserved verbatim as <see cref="MsDataScan.ScanFilter"/>.</description></item>
        /// <item><description>No noise data, and no precursor scan number — PROXI reports neither.</description></item>
        /// </list>
        /// Retention time is normalized to <b>minutes</b>, mzLib's convention everywhere. "MS:1000016" is
        /// interpreted per its unit term when one is present, and as seconds when it is absent — which is what
        /// PRIDE serves (a live scan reports 4662.32, i.e. 77.7 min, not 4662 min).
        /// Ms level defaults to 2 and the scan number to 1 when the respective term is absent; a USI addresses a
        /// fragment spectrum, so 2 is the meaningful default (the same one <c>Mgf</c> hardcodes).
        /// </remarks>
        /// <param name="spectrum">The PROXI spectrum to convert.</param>
        /// <returns>An <see cref="MsDataScan"/> wrapping <paramref name="spectrum"/> and its CV metadata.</returns>
        /// <exception cref="ArgumentNullException">The spectrum is null.</exception>
        public static MsDataScan ToMsDataScan(this PrideProxiSpectrum spectrum)
        {
            if (spectrum == null)
                throw new ArgumentNullException(nameof(spectrum));

            int scanNumber = GetInt(spectrum, ScanNumberAccession) ?? 1;
            double? selectedIonMz = GetDouble(spectrum, SelectedIonMzAccession);

            // PROXI states the isolation window as a target m/z plus offsets; MsDataScan wants a centre and a
            // total width. Fall back to the precursor m/z as the centre, as the readers do when only one is known.
            double? lowerOffset = GetDouble(spectrum, IsolationWindowLowerOffsetAccession);
            double? upperOffset = GetDouble(spectrum, IsolationWindowUpperOffsetAccession);
            double? isolationWidth = lowerOffset.HasValue || upperOffset.HasValue
                ? (lowerOffset ?? 0) + (upperOffset ?? 0)
                : null;

            // PROXI states no acquisition scan window, so bound it by the peaks actually returned — the same
            // synthesis Mgf performs for a format that omits it.
            MzRange scanWindowRange = spectrum.Size > 0
                ? new MzRange(spectrum.XArray[0], spectrum.XArray[spectrum.Size - 1])
                : null;

            return new MsDataScan(
                massSpectrum: spectrum,
                oneBasedScanNumber: scanNumber,
                msnOrder: GetInt(spectrum, MsLevelAccession) ?? 2,
                isCentroid: IsCentroid(spectrum),
                polarity: GetPolarity(spectrum),
                retentionTime: GetRetentionTimeInMinutes(spectrum),
                scanWindowRange: scanWindowRange,
                scanFilter: GetValue(spectrum, FilterStringAccession),
                mzAnalyzer: MZAnalyzerType.Unknown,
                totalIonCurrent: spectrum.SumOfAllY,
                injectionTime: GetDouble(spectrum, IonInjectionTimeAccession),
                noiseData: null,
                nativeId: $"scan={scanNumber}",
                selectedIonMz: selectedIonMz,
                selectedIonChargeStateGuess: GetInt(spectrum, ChargeStateAccession),
                selectedIonIntensity: GetDouble(spectrum, SelectedIonIntensityAccession),
                isolationMZ: GetDouble(spectrum, IsolationWindowTargetMzAccession) ?? selectedIonMz,
                isolationWidth: isolationWidth,
                dissociationType: DissociationType.Unknown,
                oneBasedPrecursorScanNumber: null,
                selectedIonMonoisotopicGuessMz: selectedIonMz);
        }

        /// <summary>
        /// The value of the first attribute with <paramref name="accession"/>, or null if absent or blank. First,
        /// not single: PRIDE legitimately repeats a term (one "instrument" per instrument in the project).
        /// </summary>
        private static string GetValue(PrideProxiSpectrum spectrum, string accession)
        {
            foreach (CvParam attribute in spectrum.Attributes)
            {
                if (attribute != null && string.Equals(attribute.Accession, accession, StringComparison.Ordinal))
                    return string.IsNullOrEmpty(attribute.Value) ? null : attribute.Value;
            }
            return null;
        }

        /// <summary>True if an attribute with <paramref name="accession"/> is present, whatever its value.</summary>
        private static bool HasAccession(PrideProxiSpectrum spectrum, string accession) =>
            spectrum.Attributes.Any(a => a != null && string.Equals(a.Accession, accession, StringComparison.Ordinal));

        /// <summary>The first attribute with <paramref name="accession"/> parsed as a double, or null.</summary>
        private static double? GetDouble(PrideProxiSpectrum spectrum, string accession) =>
            double.TryParse(GetValue(spectrum, accession), NumberStyles.Float, CultureInfo.InvariantCulture, out double value)
                ? value
                : null;

        /// <summary>The first attribute with <paramref name="accession"/> parsed as an int, or null.</summary>
        private static int? GetInt(PrideProxiSpectrum spectrum, string accession) =>
            int.TryParse(GetValue(spectrum, accession), NumberStyles.Integer, CultureInfo.InvariantCulture, out int value)
                ? value
                : null;

        /// <summary>
        /// The scan's polarity. PRIDE states it as "MS:1000465" carrying the value "positive scan"/"negative
        /// scan"; other PROXI servers may instead send the bare PSI terms "MS:1000130"/"MS:1000129". Both forms
        /// are read; an unstated polarity is <see cref="Polarity.Unknown"/>.
        /// </summary>
        private static Polarity GetPolarity(PrideProxiSpectrum spectrum)
        {
            if (HasAccession(spectrum, PositiveScanAccession))
                return Polarity.Positive;
            if (HasAccession(spectrum, NegativeScanAccession))
                return Polarity.Negative;

            string polarity = GetValue(spectrum, ScanPolarityAccession);
            if (polarity == null)
                return Polarity.Unknown;
            if (polarity.Contains("positive", StringComparison.OrdinalIgnoreCase))
                return Polarity.Positive;
            if (polarity.Contains("negative", StringComparison.OrdinalIgnoreCase))
                return Polarity.Negative;
            return Polarity.Unknown;
        }

        /// <summary>
        /// Whether the peaks are centroided. PRIDE states this as "MS:1000525" carrying "centroid spectrum";
        /// other servers may send the bare term "MS:1000127". Absent either, false (profile) is assumed, matching
        /// the PSI default.
        /// </summary>
        private static bool IsCentroid(PrideProxiSpectrum spectrum)
        {
            if (HasAccession(spectrum, CentroidSpectrumAccession))
                return true;

            string representation = GetValue(spectrum, SpectrumRepresentationAccession);
            return representation != null && representation.Contains("centroid", StringComparison.OrdinalIgnoreCase);
        }

        /// <summary>
        /// The scan's retention time in minutes, mzLib's convention. "MS:1000016" is converted per its unit term
        /// when one is present; with no unit term the value is read as seconds, which is what PRIDE serves.
        /// Zero when the term is absent.
        /// </summary>
        private static double GetRetentionTimeInMinutes(PrideProxiSpectrum spectrum)
        {
            double? scanStartTime = GetDouble(spectrum, ScanStartTimeAccession);
            if (!scanStartTime.HasValue)
                return 0;

            CvParam term = spectrum.Attributes.FirstOrDefault(
                a => a != null && string.Equals(a.Accession, ScanStartTimeAccession, StringComparison.Ordinal));
            bool isMinutes = string.Equals(term?.UnitAccession, MinuteUnitAccession, StringComparison.Ordinal)
                || string.Equals(term?.UnitName, "minute", StringComparison.OrdinalIgnoreCase);

            return isMinutes ? scanStartTime.Value : scanStartTime.Value / 60.0;
        }
    }
}
