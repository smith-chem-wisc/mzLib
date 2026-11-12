using System.Globalization;
using System.Text.RegularExpressions;

namespace Readers
{
    /// <summary>
    /// Reads the scan number out of a spectrum's nativeID, the PSI-MS identifier that mzIdentML, mzTab and mzML
    /// all use to name a spectrum within its file.
    /// </summary>
    internal static class SpectrumNativeId
    {
        private static readonly Regex Index = new(@"(^|\s)index=(\d+)($|\s)", RegexOptions.Compiled);

        /// <summary>
        /// The nativeID's scan number: "scan=N" gives N; "index=N", a zero-based position in a peak list, gives
        /// N + 1, as CasanovoMzTabFile does; anything else gives -1.
        /// </summary>
        internal static int OneBasedScanNumberOf(string? nativeId)
        {
            if (string.IsNullOrEmpty(nativeId))
            {
                return -1;
            }

            var scan = Mzml.nativeIdScanNumberParser.Match(nativeId);
            if (scan.Success && int.TryParse(scan.Groups[2].Value, NumberStyles.None, CultureInfo.InvariantCulture, out int scanNumber))
            {
                return scanNumber;
            }

            var index = Index.Match(nativeId);
            if (index.Success && int.TryParse(index.Groups[2].Value, NumberStyles.None, CultureInfo.InvariantCulture, out int position))
            {
                return position + 1;
            }

            return -1;
        }
    }
}
