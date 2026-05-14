namespace MassSpectrometry;

/// <summary>
/// Lightweight, immutable snapshot of scan and precursor metadata extracted from an MS2 scan.
/// Designed to be shared across spectral matches (PSMs) from the same scan/precursor,
/// avoiding duplication of scalar metadata while allowing the heavyweight scan objects
/// (MsDataScan, MzSpectrum, IsotopicEnvelope[]) to be released from memory after scoring.
///
/// Scan-level properties (OneBasedScanNumber through NativeId) are identical for all
/// precursors deconvoluted from the same raw scan. Precursor-level properties
/// (PrecursorCharge through OneOverK0) are specific to a single deconvoluted precursor
/// and may differ across chimeric identifications from the same scan.
/// </summary>
/// <param name="OneBasedScanNumber">One-based scan number from the raw file.</param>
/// <param name="OneBasedPrecursorScanNumber">One-based scan number of the precursor (MS1) scan, if available.</param>
/// <param name="RetentionTime">Retention time in minutes.</param>
/// <param name="NumPeaks">Number of peaks in the MS2 spectrum at the time of extraction.</param>
/// <param name="TotalIonCurrent">Total ion current of the MS2 scan.</param>
/// <param name="NativeId">Vendor-native scan identifier string.</param>
/// <param name="FullFilePath">Absolute or relative path to the originating spectra file.</param>
/// <param name="PrecursorCharge">Charge state assigned to the deconvoluted precursor.</param>
/// <param name="PrecursorMonoisotopicPeakMz">Monoisotopic m/z of the deconvoluted precursor.</param>
/// <param name="PrecursorMass">Neutral monoisotopic mass of the precursor, derived from m/z and charge.</param>
/// <param name="PrecursorIntensity">MS1 intensity of the precursor ion.</param>
/// <param name="PrecursorEnvelopePeakCount">Number of peaks in the precursor isotopic envelope.</param>
/// <param name="PrecursorFractionalIntensity">Fraction of precursor intensity relative to envelope total. -1 if unavailable.</param>
/// <param name="OneOverK0">Inverse reduced ion mobility (1/K0) for TIMS data; null for non-IMS instruments.</param>
public record ScanMetadata(
    // Scan-level properties
    int OneBasedScanNumber,
    int? OneBasedPrecursorScanNumber,
    double RetentionTime,
    int NumPeaks,
    double TotalIonCurrent,
    string NativeId,
    string FullFilePath,

    // Precursor-level properties
    int PrecursorCharge,
    double PrecursorMonoisotopicPeakMz,
    double PrecursorMass,
    double PrecursorIntensity,
    int PrecursorEnvelopePeakCount,
    double PrecursorFractionalIntensity,
    double? OneOverK0 = null)
{
    /// <summary>
    /// Convenience property deriving the file name without extension from <see cref="FullFilePath"/>.
    /// </summary>
    public string FilenameWithoutExtension =>
        System.IO.Path.GetFileNameWithoutExtension(FullFilePath);
}
