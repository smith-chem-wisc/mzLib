using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using Omics;
using Omics.BioPolymerGroup;
using Quantification.Interfaces;

namespace Quantification.Strategies
{
    public class SumRollUp : IRollUpStrategy
    {
        public string Name => "Sum Roll-Up";
        // Implement roll-up methods here

        public SumRollUp()
        {
        }

        public PeptideMatrix RollUpSpectralMatches(IExperimentalDesign experimentalDesign, List<ISpectralMatch> spectralMatches, List<IBioPolymerWithSetMods> peptides)
        {
            throw new NotImplementedException();
            //PeptideMatrix result = new PeptideMatrix(peptides, experimentalDesign.Samples);
            //foreach (var peptide in peptides)
            //{ }
        }

        public ProteinMatrix RollUpPeptides(PeptideMatrix peptides, List<IBioPolymerGroup> proteins)
        {
            throw new NotImplementedException();
            //ProteinMatrix result = new ProteinMatrix();
            //foreach (var peptide in peptides)
            //{


            //}
            //return result;
        }
    }
}
