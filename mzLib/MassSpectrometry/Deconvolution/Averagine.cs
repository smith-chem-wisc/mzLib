using Chemistry;
using System;
using System.Linq;

namespace MassSpectrometry.Deconvolution
{
    public abstract class AverageResidue
    {
        protected const double fineRes = 0.125;
        protected const double minRes = 1e-8;
        protected static readonly int NumAveraginesToGenerate = 1500;
        public abstract int GetMostIntenseMassIndex(double testMass);
        public abstract double[] GetAllTheoreticalMasses(int index);
        public abstract double[] GetAllTheoreticalIntensities(int index);
        public abstract double GetDiffToMonoisotopic(int index);
    }
    public sealed class Averagine : AverageResidue
    {
        static readonly double[][] AllMasses = new double[NumAveraginesToGenerate][];
        static readonly double[][] AllIntensities = new double[NumAveraginesToGenerate][];
        static readonly double[] MostIntenseMasses = new double[NumAveraginesToGenerate];
        static readonly double[] DiffToMonoisotopic = new double[NumAveraginesToGenerate];
        static Averagine() 
        {
            double averageC = 4.9384;
            double averageH = 7.7583;
            double averageO = 1.4773;
            double averageN = 1.3577;
            double averageS = 0.0417;

            for (int i = 0; i < NumAveraginesToGenerate; i++)
            {
                double averagineMultiplier = (i + 1) / 2.0;
                //Console.Write("numAveragines = " + numAveragines);
                ChemicalFormula chemicalFormula = new ChemicalFormula();
                chemicalFormula.Add("C", Convert.ToInt32(averageC * averagineMultiplier));
                chemicalFormula.Add("H", Convert.ToInt32(averageH * averagineMultiplier));
                chemicalFormula.Add("O", Convert.ToInt32(averageO * averagineMultiplier));
                chemicalFormula.Add("N", Convert.ToInt32(averageN * averagineMultiplier));
                chemicalFormula.Add("S", Convert.ToInt32(averageS * averagineMultiplier));

                var chemicalFormulaReg = chemicalFormula;
                IsotopicDistribution ye = IsotopicDistribution.GetDistribution(chemicalFormulaReg, fineRes, minRes);
                var masses = ye.Masses.ToArray();
                var intensities = ye.Intensities.ToArray();
                Array.Sort(intensities, masses);
                Array.Reverse(intensities);
                Array.Reverse(masses);

                MostIntenseMasses[i] = masses[0];
                DiffToMonoisotopic[i] = masses[0] - chemicalFormulaReg.MonoisotopicMass;
                AllMasses[i] = masses;
                AllIntensities[i] = intensities;
            }
        }
        public override int GetMostIntenseMassIndex(double testMass)
        {
            return MostIntenseMasses.GetClosestIndex(testMass);
        }
        public override double[] GetAllTheoreticalMasses(int index)
        {
            return AllMasses[index];
        }
        public override double[] GetAllTheoreticalIntensities(int index)
        {
            return AllIntensities[index];
        }

        public override double GetDiffToMonoisotopic(int index)
        {
            return DiffToMonoisotopic[index];
        }
    }

    public sealed class Averatide : AverageResidue
    {
        static readonly double[][] AllMasses = new double[NumAveraginesToGenerate][];
        static readonly double[][] AllIntensities = new double[NumAveraginesToGenerate][];
        static readonly double[] MostIntenseMasses = new double[NumAveraginesToGenerate];
        static readonly double[] DiffToMonoisotopic = new double[NumAveraginesToGenerate];
        static Averatide()
        {
            double averageC = 9.5;
            double averageH = 12.75;
            double averageO = 3.75;
            double averageN = 5;

            for (int i = 0; i < NumAveraginesToGenerate; i++)
            {
                double averagineMultiplier = (i + 1) / 2.0;
                //Console.Write("numAveragines = " + numAveragines);
                ChemicalFormula chemicalFormula = new ChemicalFormula();
                chemicalFormula.Add("C", Convert.ToInt32(averageC * averagineMultiplier));
                chemicalFormula.Add("H", Convert.ToInt32(averageH * averagineMultiplier));
                chemicalFormula.Add("O", Convert.ToInt32(averageO * averagineMultiplier));
                chemicalFormula.Add("N", Convert.ToInt32(averageN * averagineMultiplier));

                {
                    var chemicalFormulaReg = chemicalFormula;
                    IsotopicDistribution ye = IsotopicDistribution.GetDistribution(chemicalFormulaReg, fineRes, minRes);
                    var masses = ye.Masses.ToArray();
                    var intensities = ye.Intensities.ToArray();
                    Array.Sort(intensities, masses);
                    Array.Reverse(intensities);
                    Array.Reverse(masses);

                    MostIntenseMasses[i] = masses[0];
                    DiffToMonoisotopic[i] = masses[0] - chemicalFormulaReg.MonoisotopicMass;
                    AllMasses[i] = masses;
                    AllIntensities[i] = intensities;
                }
            }
        }
        public override int GetMostIntenseMassIndex(double testMass)
        {
            return MostIntenseMasses.GetClosestIndex(testMass);
        }
        public override double[] GetAllTheoreticalMasses(int index)
        {
            return AllMasses[index];
        }
        public override double[] GetAllTheoreticalIntensities(int index)
        {
            return AllIntensities[index];
        }

        public override double GetDiffToMonoisotopic(int index)
        {
            return DiffToMonoisotopic[index];
        }
    }
}
