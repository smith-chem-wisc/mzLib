namespace FlashLFQ
{
    /// <summary>
    /// Defines the interface for spectra file information
    /// </summary>
    public interface IQuantifiableSpectraFile
    {
        /// <summary>
        /// The filename without extension
        /// </summary>
        string FilenameWithoutExtension { get; }

        /// <summary>
        /// The experimental condition
        /// </summary>
        string Condition { get; }

        /// <summary>
        /// The biological replicate number
        /// </summary>
        int BiologicalReplicate { get; }

        /// <summary>
        /// The technical replicate number
        /// </summary>
        int TechnicalReplicate { get; }

        /// <summary>
        /// The fraction number
        /// </summary>
        int Fraction { get; }
    }
}
