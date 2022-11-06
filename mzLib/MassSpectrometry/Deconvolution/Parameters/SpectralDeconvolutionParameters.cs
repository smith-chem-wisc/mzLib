using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry.Proteomics;
using MassSpectrometry.Proteomics.ProteolyticDigestion;

namespace MassSpectrometry.Deconvolution.Parameters
{

    public class SpectralDeconvolutionParameters : DeconvolutionParameters
    {
        public List<Protein> Proteins;
        public List<Modification> FixedModifications;
        public List<Modification> VariableModifications;
        public Protease Protease;
        private bool FindNonDatabasePeaks; // This should be linked to a method that generates Averagine envelopes


        public SpectralDeconvolutionParameters(int minAssumedChargeState, int maxAssumedChargeState,
            double deconvolutionTolerancePpm, List<Protein> proteins, List<Modification> fixedModifications,
            List<Modification> variableModifications, Protease protease, bool findNonDatabasePeaks = false) :
            base(minAssumedChargeState, maxAssumedChargeState, deconvolutionTolerancePpm)
        {
            Proteins = proteins;
            FixedModifications = fixedModifications;
            VariableModifications = variableModifications;
            Protease = protease;
            FindNonDatabasePeaks = findNonDatabasePeaks;
        }
    }
}
