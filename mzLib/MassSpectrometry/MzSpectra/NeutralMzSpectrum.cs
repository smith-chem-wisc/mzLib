using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry.MzSpectra
{
    public class NeutralMzSpectrum : MzSpectrum
    {
        public NeutralMzSpectrum(double[,] mzintensities) : base(mzintensities)
        {
        }

        public NeutralMzSpectrum(double[] monoisotopicMasses, double[] intensities, int[] charges, bool shouldCopy)
            : base(monoisotopicMasses, intensities, shouldCopy)
        {
            if (monoisotopicMasses.Length != intensities.Length || monoisotopicMasses.Length != charges.Length)
                throw new ArgumentException("The lengths of monoisotopicMasses, intensities, and charges must be the same.");
            
        }
    }
}
