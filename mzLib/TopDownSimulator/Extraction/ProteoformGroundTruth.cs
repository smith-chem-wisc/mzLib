namespace TopDownSimulator.Extraction;

/// <summary>
/// The raw MS1 signal around a single proteoform, laid out in the shape the
/// per-factor fitter consumes: one entry per (charge, isotopologue, scan).
///
/// Zero intensities represent "peak not found within tolerance" and are kept
/// (rather than dropped) so the fitter sees the correct scan count per charge.
/// </summary>
public sealed class ProteoformGroundTruth
{
    public required double MonoisotopicMass { get; init; }
    public required double RetentionTimeCenter { get; init; }
    public required int MinCharge { get; init; }
    public required int MaxCharge { get; init; }

    /// <summary>Zero-based indices of the MS1 scans covered by this ground truth.</summary>
    public required int[] ZeroBasedScanIndices { get; init; }

    /// <summary>Retention time (minutes) of each covered scan, parallel to <see cref="ZeroBasedScanIndices"/>.</summary>
    public required double[] ScanTimes { get; init; }

    /// <summary>
    /// Theoretical centroid m/z per isotopologue for each charge.
    /// Indexed as [chargeOffset][isotopologue], where chargeOffset = z - MinCharge.
    /// </summary>
    public required double[][] CentroidMzs { get; init; }

    /// <summary>
    /// Observed peak intensity per (chargeOffset, isotopologue, scanIndex).
    /// This is the intensity of the single peak closest to the theoretical centroid
    /// (within the m/z window). Missing peaks are stored as 0.
    /// </summary>
    public required double[][][] IsotopologueIntensities { get; init; }

    /// <summary>
    /// All peaks observed within ±<see cref="MzWindowHalfWidth"/> m/z of each theoretical
    /// centroid, per (chargeOffset, isotopologue, scanIndex). Innermost array may be empty
    /// if no peaks fell in the window. Supplies the raw m/z profile the envelope-width fitter
    /// consumes.
    /// </summary>
    public required PeakSample[][][][] IsotopologuePeakWindows { get; init; }

    /// <summary>m/z half-width used when collecting <see cref="IsotopologuePeakWindows"/>.</summary>
    public required double MzWindowHalfWidth { get; init; }

    /// <summary>
    /// Summed-isotopologue intensity per (chargeOffset, scanIndex). Shape [nCharges][nScans].
    /// This is the per-charge XIC the RT-profile and charge-distribution fitters read.
    /// </summary>
    public required double[][] ChargeXics { get; init; }

    public int ChargeCount => MaxCharge - MinCharge + 1;
    public int ScanCount => ZeroBasedScanIndices.Length;
}
