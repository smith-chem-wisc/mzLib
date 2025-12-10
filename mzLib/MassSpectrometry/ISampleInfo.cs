namespace MassSpectrometry
{
    /// <summary>
    /// Common interface for sample information used in quantification.
    /// Provides properties needed for grouping, display, and identification of samples.
    /// </summary>
    public interface ISampleInfo
    {
        /// <summary>
        /// Full path or identifier for the source file. May be empty for non-file-based samples.
        /// </summary>
        string FullFilePathWithExtension { get; }
        /// <summary>
        /// Display name for the sample (used in headers and reports).
        /// </summary>
        string SampleIdentifier { get; }
        /// <summary>
        /// The condition or experimental group this sample belongs to.
        /// </summary>
        string Condition { get; }

        /// <summary>
        /// Biological replicate identifier within a condition.
        /// </summary>
        int BiologicalReplicate { get; }

        /// <summary>
        /// Technical replicate identifier for repeated measurements.
        /// </summary>
        int TechnicalReplicate { get; }

        /// <summary>
        /// Fraction identifier for fractionated workflows. Returns 0 if not applicable.
        /// </summary>
        int Fraction { get; }
    }
}
