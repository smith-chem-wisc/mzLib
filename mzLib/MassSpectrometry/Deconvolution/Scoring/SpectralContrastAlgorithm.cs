using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry.Deconvolution.Scoring
{
    public class SpectralContrastAlgorithm : ScoringAlgorithm
    {
        public SpectralContrastAlgorithm(IScoreArgs arguments) : base(arguments)
        {

        }

        public override double Score()
        {
            return 0;
        }
    }
}
