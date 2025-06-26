using Chemistry;
using System;
using System.Linq;

namespace MassSpectrometry;

/// <summary>
/// Represents the average Amino Acid and is used for most abundant isotopic peak to monoisotopic peak difference.
/// </summary>
/// <remarks>All instance methods return a reference to its static precalculated values</remarks>
public sealed class Averagine : AverageResidue
{
    public static readonly double[][] AllMasses = new double[NumAveraginesToGenerate][];
    public static readonly double[][] AllIntensities = new double[NumAveraginesToGenerate][];
    public static readonly double[] MostIntenseMasses = new double[NumAveraginesToGenerate];
    public static readonly double[] DiffToMonoisotopic = new double[NumAveraginesToGenerate];
    public override int GetMostIntenseMassIndex(double testMass) => MostIntenseMasses.GetClosestIndex(testMass);
    public override double[] GetAllTheoreticalMasses(int index) => AllMasses[index];
    public override double[] GetAllTheoreticalIntensities(int index) => AllIntensities[index];
    public override double GetDiffToMonoisotopic(int index) => DiffToMonoisotopic[index];
    static Averagine() 
    {
        // Magic numbers determined by https://pmc.ncbi.nlm.nih.gov/articles/PMC6166224/
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
            IsotopicDistribution ye = IsotopicDistribution.GetDistribution(chemicalFormulaReg, FineRes, MinRes);
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