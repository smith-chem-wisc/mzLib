using System.Collections.Generic;

namespace FlashLFQ
{
    /// <summary>
    /// Defines the interface for protein groups that can be quantified
    /// </summary>
    public interface IQuantifiableProteinGroup
    {
        /// <summary>
        /// The protein group name/identifier
        /// </summary>
        string ProteinGroupName { get; }

        /// <summary>
        /// Gets the intensity for a specific file
        /// </summary>
        double GetIntensity(IQuantifiableSpectraFile fileInfo);

        /// <summary>
        /// Sets the intensity for a specific file
        /// </summary>
        void SetIntensity(IQuantifiableSpectraFile fileInfo, double intensity);
    }
}
