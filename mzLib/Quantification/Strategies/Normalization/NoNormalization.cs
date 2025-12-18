using Quantification.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Omics;
using Omics.BioPolymerGroup;

namespace Quantification.Strategies
{
    public class NoNormalization : INormalizationStrategy
    {
        public NoNormalization()
        {
        }

        public string Name => "No Normalization";

        public PeptideMatrix NormalizePeptideIntensities(PeptideMatrix peptideMatrix, List<IBioPolymerWithSetMods> peptides)
        {
            //TODO: Populate the IntensitiesBySample dictionaries in each peptide in place without changing the values.
            return peptideMatrix;
        }

        public ProteinMatrix NormalizeProteinIntensities(ProteinMatrix proteinMatrix, List<IBioPolymerGroup> proteins)
        {
            //TODO: Populate the IntensitiesBySample dictionaries in each protein in place without changing the values.
            return proteinMatrix;
        }
    }
}
