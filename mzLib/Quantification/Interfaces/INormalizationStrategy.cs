using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using Omics;
using Omics.BioPolymerGroup;

namespace Quantification.Interfaces
{
    public enum NormalizationStrategyType
    {
        None,
        Median,
        Quantile,
        VarianceStabilization
    }

    /// <summary>
    /// Normalization is responsible for normalizing data, and then modifying the IntensitiesBySample dictionary in each IBioPolymerWithSetMods or IBioPolymerGroup in place.
    /// </summary>
    public interface INormalizationStrategy
    {
        string Name { get; }

        /// <summary>
        /// Normalize peptide intensities in place by modifying the IntensitiesBySample dictionary in each IBioPolymerWithSetMods.
        /// </summary>
        /// <param name="peptides"></param>
        PeptideMatrix NormalizePeptideIntensities(PeptideMatrix peptideMatrix, List<IBioPolymerWithSetMods> peptides);

        /// <summary>
        /// Normalize protein intensities in place by modifying the IntensitiesBySample dictionary in each IBioPolymerGroup.
        /// </summary>
        /// <param name="proteins"></param>
        ProteinMatrix NormalizeProteinIntensities(ProteinMatrix proteinMatrix, List<IBioPolymerGroup> proteins);
    }

}
