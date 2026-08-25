using System.Collections.Generic;
using MassSpectrometry;
using Omics;
using Omics.BioPolymerGroup;

namespace Quantification
{
    /// <summary>
    /// The result of a <see cref="QuantificationEngine"/> run.
    ///
    /// Protein-level values are additionally written back onto the <see cref="IBioPolymerGroup"/>
    /// objects themselves (via <see cref="IHasSampleIntensities"/>) so that existing writers —
    /// <see cref="BioPolymerGroup.PopulateSampleGroupResults"/>, <see cref="BioPolymerGroup.ToString"/>,
    /// and <see cref="BioPolymerGroup.GetTabSeparatedHeader"/> — can render them without any
    /// knowledge of this class.
    ///
    /// Peptide-level values are carried here rather than written onto the peptide objects.
    /// <see cref="Proteomics.ProteolyticDigestion.PeptideWithSetModifications"/> is
    /// <c>[Serializable]</c> and is serialized in bulk during indexing, so adding stored
    /// quantification state to it would enlarge every indexed peptide. A side table costs nothing
    /// when quantification is not run.
    /// </summary>
    public class QuantificationResults
    {
        public QuantificationResults() { }

        /// <summary>
        /// Human-readable outcome. On failure this describes what was wrong with the engine inputs
        /// and the remaining properties are empty.
        /// </summary>
        public string Summary { get; internal set; }

        /// <summary>
        /// True when the engine ran to completion. False when validation rejected the inputs.
        /// </summary>
        public bool Success { get; internal set; }

        /// <summary>
        /// The samples that form the columns of <see cref="PeptideIntensities"/> and
        /// <see cref="ProteinIntensities"/>, in output order. Empty on failure.
        /// </summary>
        public IReadOnlyList<ISampleInfo> Samples { get; internal set; } = new List<ISampleInfo>();

        /// <summary>
        /// Final peptide-level intensities, keyed by peptide then by sample. Samples with no
        /// measured value for a peptide are absent from the inner dictionary rather than present
        /// with a zero, so callers can distinguish "not observed" from "observed at zero".
        /// </summary>
        public IReadOnlyDictionary<IBioPolymerWithSetMods, Dictionary<ISampleInfo, double>> PeptideIntensities { get; internal set; }
            = new Dictionary<IBioPolymerWithSetMods, Dictionary<ISampleInfo, double>>();

        /// <summary>
        /// Final protein-group-level intensities, keyed by group then by sample. The same values are
        /// written onto each group's <see cref="IHasSampleIntensities.IntensitiesBySample"/>.
        /// </summary>
        public IReadOnlyDictionary<IBioPolymerGroup, Dictionary<ISampleInfo, double>> ProteinIntensities { get; internal set; }
            = new Dictionary<IBioPolymerGroup, Dictionary<ISampleInfo, double>>();

        /// <summary>
        /// Paths of the files the engine wrote, in the order they were written. Empty when all of
        /// the <c>Write*Information</c> parameters were left off.
        /// </summary>
        public IReadOnlyList<string> WrittenFiles { get; internal set; } = new List<string>();

        /// <summary>
        /// The directory the engine wrote to. This is <see cref="QuantificationParameters.OutputDirectory"/>
        /// when the caller set one, and otherwise the directory derived from the source data files, so a
        /// caller who relied on the default can find out where output went. Null when the run wrote nothing.
        /// </summary>
        public string OutputDirectory { get; internal set; }

        internal static QuantificationResults Failure(string summary) =>
            new QuantificationResults { Summary = summary, Success = false };
    }
}
