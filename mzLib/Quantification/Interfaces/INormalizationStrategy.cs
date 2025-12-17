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
    internal interface INormalizationStrategy
    {
        string Name { get; }

        /// <summary>
        /// Normalize peptide intensities in place by modifying the IntensitiesBySample dictionary in each IBioPolymerWithSetMods.
        /// </summary>
        /// <param name="peptides"></param>
        public void NormalizePeptideIntensities(List<IBioPolymerWithSetMods> peptides);

        /// <summary>
        /// Normalize protein intensities in place by modifying the IntensitiesBySample dictionary in each IBioPolymerGroup.
        /// </summary>
        /// <param name="proteins"></param>
        public void NormalizeProteinIntensities(List<IBioPolymerGroup> proteins);
    }

}
