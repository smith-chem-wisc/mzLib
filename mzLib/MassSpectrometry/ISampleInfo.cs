using System;

namespace MassSpectrometry
{
    /// <summary>
    /// Common interface for sample information used in quantification.
    /// Provides properties needed for grouping, display, and identification of samples.
    /// </summary>
    public interface ISampleInfo : IComparable<ISampleInfo>, IEquatable<ISampleInfo>
    {
        /// <summary>
        /// Full path or identifier for the source file. May be empty for non-file-based samples.
        /// </summary>
        string FullFilePathWithExtension { get; }
        /// <summary>
        /// "The condition of the sample (e.g., 'Control' or 'Treatment')") and IsobaricQuantSampleInfo (as "The condition or experimental group this sample belongs to"). 
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
