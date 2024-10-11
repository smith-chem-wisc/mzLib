using System;
using Chemistry;

namespace MassSpectrometry
{
    public class NeutralMassSpectrum : MzSpectrum
    {
        public int[] Charges { get; init; }
        public NeutralMassSpectrum(double[,] monoisotopicMassesIntensities, int[] charges) : base(monoisotopicMassesIntensities)
        {
            if (monoisotopicMassesIntensities.Length != charges.Length || monoisotopicMassesIntensities.Length != monoisotopicMassesIntensities.Rank)
                throw new ArgumentException("The lengths of monoisotopicMasses, intensities, and charges must be the same.");

            Charges = charges;
        }

        public NeutralMassSpectrum(double[] monoisotopicMasses, double[] intensities, int[] charges, bool shouldCopy)
            : base(monoisotopicMasses, intensities, shouldCopy)
        {
            if (monoisotopicMasses.Length != intensities.Length || monoisotopicMasses.Length != charges.Length)
                throw new ArgumentException("The lengths of monoisotopicMasses, intensities, and charges must be the same.");

            Charges = charges;
        }

        /// <summary>
        /// Converts to a charged spectrum
        /// </summary>
        protected override MzPeak GeneratePeak(int index)
        {
            return new MzPeak(XArray[index].ToMz(Charges[index]), YArray[index]);
        }
    }
}
