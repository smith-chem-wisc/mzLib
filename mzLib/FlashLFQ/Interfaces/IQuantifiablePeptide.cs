using System.Collections.Generic;

namespace FlashLFQ
{
    /// <summary>
    /// Defines the interface for peptides that can be quantified
    /// </summary>
    public interface IQuantifiablePeptide
    {
        /// <summary>
        /// The peptide sequence
        /// </summary>
        string Sequence { get; }

        /// <summary>
        /// Whether this peptide should be used for protein quantification
        /// </summary>
        bool UseForProteinQuant { get; }

        /// <summary>
        /// The protein groups this peptide belongs to
        /// </summary>
        IEnumerable<IQuantifiableProteinGroup> ProteinGroups { get; }

        /// <summary>
        /// Gets the intensity for a specific file
        /// </summary>
        double GetIntensity(IQuantifiableSpectraFile fileInfo);

        /// <summary>
        /// Gets the detection type for a specific file
        /// </summary>
        DetectionType GetDetectionType(IQuantifiableSpectraFile fileInfo);

        /// <summary>
        /// Determines if this peptide has unambiguous quantification
        /// </summary>
        bool UnambiguousPeptideQuant();
    }
}
