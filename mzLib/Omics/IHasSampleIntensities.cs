using MassSpectrometry;

namespace Omics
{
    /// <summary>
    /// Implemented by entities that can carry per-sample quantification values written back by a
    /// quantification engine — for example an <see cref="Omics.BioPolymerGroup.IBioPolymerGroup"/>.
    ///
    /// The two members are deliberately kept on their own interface rather than on the entity
    /// interfaces themselves, so that an engine can populate any entity that opts in without every
    /// implementer of a broad interface (such as <see cref="IBioPolymerWithSetMods"/>) being forced
    /// to carry quantification state it does not use.
    ///
    /// <see cref="SamplesForQuantification"/> defines the column order; <see cref="IntensitiesBySample"/>
    /// holds the values. Both are null until a quantification engine populates them.
    /// </summary>
    public interface IHasSampleIntensities
    {
        /// <summary>
        /// Samples that contribute quantification data for this entity, in output column order.
        /// Supports <see cref="SpectraFileInfo"/> (label-free) and <see cref="IsobaricQuantSampleInfo"/>
        /// (TMT/iTRAQ). May be null when no experimental design is available.
        /// </summary>
        List<ISampleInfo>? SamplesForQuantification { get; set; }

        /// <summary>
        /// Measured intensity values for this entity, keyed by sample.
        /// Supports both <see cref="SpectraFileInfo"/> and <see cref="IsobaricQuantSampleInfo"/> as keys.
        /// May be null when no intensity data is available.
        /// </summary>
        Dictionary<ISampleInfo, double>? IntensitiesBySample { get; set; }
    }
}
