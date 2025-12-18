using MassSpectrometry;
using Omics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Omics.BioPolymerGroup;
using System.Numerics;

namespace Quantification.Interfaces
{
    
    public enum RollUpStrategyType
    {
        /// <summary>
        /// Simple summation of intensities.
        /// </summary>
        Sum,
        /// <summary>
        /// Median of intensities.
        /// </summary>
        Median,
        /// <summary>
        /// Weighted average based on spectral quality.
        /// </summary>
        WeightedAverage
    }

    public interface IRollUpStrategy
    {
        string Name { get; }

        /// <summary>
        /// Takes in a list of spectral matches, then rolls up the quantification values to the peptide level.
        /// The rolled-up values are written to the IntensitiesBySample dictionary in IBioPolymerWithSetMods.
        /// </summary>
        /// <param name="experimentalDesign"></param>
        /// <param name="spectralMatches"></param>
        /// <param name="peptides"></param>
        public PeptideMatrix  RollUpSpectralMatches(IExperimentalDesign experimentalDesign, List<ISpectralMatch> spectralMatches, List<IBioPolymerWithSetMods> peptides);

        /// <summary>
        /// Takes in a list of peptides, then rolls up the quantification values to the protein level.
        /// Rolled-up values are written to the IntensitiesBySample dictionary in IBioPolymerGroup.
        /// </summary>
        /// <param name="peptides"></param>
        /// <param name="proteins"></param>
        public ProteinMatrix RollUpPeptides(PeptideMatrix peptides, List<IBioPolymerGroup> proteins);
    }
}
