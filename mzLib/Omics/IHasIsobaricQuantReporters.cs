using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Omics
{
    /// <summary>
    /// Capability interface for objects that can expose isobaric-quantification reporter intensities
    /// (for example TMT or iTRAQ reporter channels) and be ordered for comparisons.
    /// Implementers provide a fixed-order array of reporter intensities and a boolean indicating
    /// whether meaningful reporter data is present.
    /// </summary>
    internal interface IHasIsobaricQuantReporters : IComparable<IHasIsobaricQuantReporters>
    {
        /// <summary>
        /// True when the object contains valid isobaric reporter intensities that can be used
        /// for quantification. When false, <see cref="IsobaricQuantReporterIntensities"/> should
        /// be ignored by quantification logic.
        /// </summary>
        bool HasIsobaricQuantReporters { get; }

        /// <summary>
        /// Per-channel reporter intensities for isobaric quantification.
        /// The array is ordered by reporter channel index (implementations must document the
        /// channel ordering they use). Values are raw or pre-processed intensity values
        /// (implementers should document whether values are normalized).
        /// If a channel has no measured intensity, implementations may return 0.0 for that index.
        /// </summary>
        double[] IsobaricQuantReporterIntensities { get; }
    }
}