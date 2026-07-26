using System;

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
    /// <remarks>
    /// When produced by <c>GroundTruthExtractor</c>, only each charge's own apex scan
    /// (<see cref="ApexScanIndexByCharge"/>) is populated; every other cell is an empty array.
    /// Materializing the full tensor costs hundreds of megabytes across a run, while one column per
    /// charge gives the width fitter a width measurement at several m/z spanning the charge
    /// envelope, which is what fitting an m/z-dependent width law needs. Objects built by hand may
    /// of course populate whatever they like.
    /// </remarks>
    public required PeakSample[][][][] IsotopologuePeakWindows { get; init; }

    /// <summary>
    /// Index of the experimental peak that supplied each entry of
    /// <see cref="IsotopologueIntensities"/>, or -1 where no peak matched. The index is into that
    /// scan's <c>MzSpectrum</c>, so (<see cref="ZeroBasedScanIndices"/>[s], this) identifies a
    /// physical peak uniquely across the whole run.
    /// </summary>
    /// <remarks>
    /// Two proteoforms that overlap in m/z and RT are matched to the <i>same</i> experimental peak,
    /// once each. Without this key any energy accounting over sample sets counts that peak twice.
    /// Empty on hand-built truths; see <see cref="HasSourcePeakIdentity"/>.
    /// </remarks>
    public int[][][] SourcePeakIndices { get; init; } = Array.Empty<int[][]>();

    /// <summary>
    /// |observed m/z − theoretical centroid| for the peak named by <see cref="SourcePeakIndices"/>,
    /// or NaN where none matched. Lets a consumer tell which of several claimants sits closest to
    /// the real peak.
    /// </summary>
    /// <remarks>
    /// Single precision, and a delta rather than the m/z itself. This is a dense tensor in a class
    /// that otherwise works hard to avoid them, and every truth is retained for the length of a run;
    /// storing the delta directly is both what the tie-break needs and half the width of the m/z it
    /// would otherwise be derived from. Distances here are ≤ the extraction window, so a float's
    /// ~1e-7 relative precision is several orders finer than any tie it has to break.
    /// </remarks>
    public float[][][] SourcePeakDeltas { get; init; } = Array.Empty<float[][]>();

    /// <summary>
    /// Whether <see cref="SourcePeakIndices"/> and <see cref="SourcePeakDeltas"/> were populated.
    /// </summary>
    public bool HasSourcePeakIdentity =>
        SourcePeakIndices.Length == ChargeCount && SourcePeakDeltas.Length == ChargeCount;

    /// <summary>
    /// Scan index at which each charge's XIC peaks, or -1 for a charge with no positive signal.
    /// These are the columns of <see cref="IsotopologuePeakWindows"/> the extractor populates.
    /// </summary>
    public int[] ApexScanIndexByCharge { get; init; } = Array.Empty<int>();

    /// <summary>m/z half-width used when collecting <see cref="IsotopologuePeakWindows"/>.</summary>
    public required double MzWindowHalfWidth { get; init; }

    /// <summary>
    /// Charge offset of the highest point of <see cref="ChargeXics"/>, or -1 when no positive
    /// signal was found. This is the column of <see cref="IsotopologuePeakWindows"/> the extractor
    /// populates.
    /// </summary>
    public int ApexChargeIndex { get; init; } = -1;

    /// <summary>
    /// Scan index of the highest point of <see cref="ChargeXics"/>, or -1 when no positive signal
    /// was found. Parallel to <see cref="ApexChargeIndex"/>.
    /// </summary>
    public int ApexScanIndex { get; init; } = -1;

    /// <summary>
    /// Summed-isotopologue intensity per (chargeOffset, scanIndex). Shape [nCharges][nScans].
    /// This is the per-charge XIC the RT-profile and charge-distribution fitters read.
    /// </summary>
    public required double[][] ChargeXics { get; init; }

    public int ChargeCount => MaxCharge - MinCharge + 1;
    public int ScanCount => ZeroBasedScanIndices.Length;
}
