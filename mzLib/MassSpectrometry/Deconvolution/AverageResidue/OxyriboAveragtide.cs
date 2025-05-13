using System;
using System.Linq;
using Chemistry;

namespace MassSpectrometry;

/// <summary>
/// Represents the average RNA nucleotide and is used for most abundant isotopic peak to monoisotopic peak difference
/// </summary>
/// <remarks>All instance methods return a reference to its static precalculated values</remarks>
public sealed class OxyriboAveragtide : AverageResidue
{
    public static readonly double[][] AllMasses = new double[NumAveraginesToGenerate][];
    public static readonly double[][] AllIntensities = new double[NumAveraginesToGenerate][];
    public static readonly double[] MostIntenseMasses = new double[NumAveraginesToGenerate];
    public static readonly double[] DiffToMonoisotopic = new double[NumAveraginesToGenerate];
    public override int GetMostIntenseMassIndex(double testMass) => MostIntenseMasses.GetClosestIndex(testMass);

    public override double[] GetAllTheoreticalMasses(int index) => AllMasses[index];

    public override double[] GetAllTheoreticalIntensities(int index) => AllIntensities[index];

    public override double GetDiffToMonoisotopic(int index) => DiffToMonoisotopic[index];
    static OxyriboAveragtide()
    {
        // Magic numbers determined by counting atoms in the main 4 canonical RNA bases. 
        double averageC = 9.5;
        double averageH = 12.75;
        double averageO = 3.75;
        double averageN = 5;

        for (int i = 0; i < NumAveraginesToGenerate; i++)
        {
            double averagineMultiplier = (i + 1) / 4.0;
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
}