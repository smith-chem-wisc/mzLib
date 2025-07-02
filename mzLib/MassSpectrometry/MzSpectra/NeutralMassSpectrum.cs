using System;
using Chemistry;

namespace MassSpectrometry
{
    public class NeutralMassSpectrum : MzSpectrum
    {
        public int[] Charges { get; init; }
        public NeutralMassSpectrum(double[,] monoisotopicMassesIntensities, int[] charges) : base(monoisotopicMassesIntensities)
        {
            if (monoisotopicMassesIntensities.GetLength(0) != charges.Length)
                throw new ArgumentException("The lengths of monoisotopicMasses, intensities, and charges must be the same.");

            Charges = charges;

            double minMz = double.MaxValue;
            double maxMz = double.MinValue;
            for (int i = 0; i < monoisotopicMassesIntensities.GetLength(0); i++)
            {
                var mz = monoisotopicMassesIntensities[i,0].ToMz(charges[i]);
                if (mz < minMz)
                    minMz = mz;
                if (mz > maxMz)
                    maxMz = mz;
            }

            FirstX = minMz;
            LastX = maxMz;
        }

        public NeutralMassSpectrum(double[] monoisotopicMasses, double[] intensities, int[] charges, bool shouldCopy)
            : base(monoisotopicMasses, intensities, shouldCopy)
        {
            if (monoisotopicMasses.GetLength(0) != intensities.Length || monoisotopicMasses.Length != charges.Length)
                throw new ArgumentException("The lengths of monoisotopicMasses, intensities, and charges must be the same.");

            Charges = charges;

            double minMz = double.MaxValue;
            double maxMz = double.MinValue;
            for (int i = 0; i < monoisotopicMasses.Length; i++)
            {
                var mz = monoisotopicMasses[i].ToMz(charges[i]);
                if (mz < minMz)
                    minMz = mz;
                if (mz > maxMz)
                    maxMz = mz;
            }

            FirstX = minMz;
            LastX = maxMz;
        }

        public override double? FirstX { get; } // in m/z
        public override double? LastX { get; } // in m/z

        /// <summary>
        /// Converts to a charged spectrum
        /// </summary>
        protected override MzPeak GeneratePeak(int index)
        {
            return new MzPeak(XArray[index].ToMz(Charges[index]), YArray[index]);
        }
    }
}
