using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Proteomics;

namespace MassSpectrometry.Deconvolution.Parameters
{
    public class SpectralDeconvolutionParameters : DeconvolutionParameters
    {
        private List<Protein> Proteins;
        private bool FindNonDatabasePeaks;

        public SpectralDeconvolutionParameters() : base()
        {

        }
    }
}
