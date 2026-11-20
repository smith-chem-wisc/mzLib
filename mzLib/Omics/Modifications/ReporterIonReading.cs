namespace Omics.Modifications
{
    /// <summary>
    /// What one scan recorded for one reporter ion. Only <see cref="Value"/> is a measured intensity; the other two are
    /// different kinds of "no value" and must not be read as 0.
    /// </summary>
    public enum ReporterIonReadingKind
    {
        /// <summary>A peak lies within tolerance of the reporter ion; the reading's intensity is the most intense one.</summary>
        Value,

        /// <summary>
        /// The scan recorded the reporter ion's m/z and no peak lies within tolerance: the intensity is below the scan's
        /// detection floor there. It is censored, not zero and not missing.
        /// </summary>
        BelowFloor,

        /// <summary>
        /// The scan could not have recorded the reporter ion: its m/z lies outside what the scan recorded (an ion-trap scan
        /// with a low-mass cutoff, say), or the spectrum is empty. Nothing is known about it.
        /// </summary>
        NotRecorded
    }

    /// <summary>Where a reading's detection floor came from.</summary>
    public enum DetectionFloorSource
    {
        /// <summary>No floor: the reporter ion was not recorded, or the scan gives nothing to estimate one from.</summary>
        None,

        /// <summary>The instrument's noise level at the reporter ion's m/z, from the scan's noise table (Thermo Orbitrap scans have one).</summary>
        NoiseTable,

        /// <summary>
        /// No noise table: the weakest peak the scan recorded between the lowest and highest reporter ions. A reporter ion
        /// with no peak there was weaker than every peak the scan kept in that narrow region.
        /// </summary>
        ReporterRegionMinimum,

        /// <summary>
        /// No noise table and no peak in the reporter region: the weakest peak anywhere in the scan. The weakest guide of
        /// the three, because the floor can differ across the scan's m/z range.
        /// </summary>
        ScanMinimum
    }

    /// <summary>
    /// One reporter ion as one scan recorded it: a value, a value below the scan's detection floor, or not recorded.
    /// </summary>
    /// <param name="Intensity">The reporter ion's intensity when <paramref name="Kind"/> is
    /// <see cref="ReporterIonReadingKind.Value"/>; NaN otherwise, never 0.</param>
    /// <param name="Kind">What the scan recorded for this reporter ion.</param>
    /// <param name="Floor">The scan's detection floor at this reporter ion's m/z: an estimate of the intensity below which
    /// the scan shows no peak there. Null when <see cref="FloorSource"/> is <see cref="DetectionFloorSource.None"/>.</param>
    /// <param name="FloorSource">Where <paramref name="Floor"/> came from.</param>
    public readonly record struct ReporterIonReading(double Intensity, ReporterIonReadingKind Kind, double? Floor,
        DetectionFloorSource FloorSource);
}
