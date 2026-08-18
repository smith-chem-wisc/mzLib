using System;
using System.Linq;
using Chemistry;

namespace MassSpectrometry;

/// <summary>
/// Represents the average RNA nucleotide and is used for most abundant isotopic peak to monoisotopic peak difference
/// </summary>
/// <remarks>All instance methods return a reference to its static precalculated values</remarks>
public sealed class OxyriboAveragine : AverageResidue
{
    public static readonly double[][] AllMasses = new double[NumAveraginesToGenerate][];
    public static readonly double[][] AllIntensities = new double[NumAveraginesToGenerate][];
    public static readonly double[] MostIntenseMasses = new double[NumAveraginesToGenerate];
    public static readonly double[] DiffToMonoisotopic = new double[NumAveraginesToGenerate];
    public override int GetMostIntenseMassIndex(double testMass) => MostIntenseMasses.GetClosestIndex(testMass);

    public override double[] GetAllTheoreticalMasses(int index) => AllMasses[index];

    public override double[] GetAllTheoreticalIntensities(int index) => AllIntensities[index];

    public override double GetDiffToMonoisotopic(int index) => DiffToMonoisotopic[index];
    static OxyriboAveragine()
    {
        // Magic numbers determined by counting atoms in the main 4 canonical RNA bases. 
        // This is not the best approach and future work should refine these numbers. 
        // One possible approach is to also incorporate the residue frequency in RNA sequences. 

        var water = ChemicalFormula.ParseFormula("H2O");
        var phosphate = ChemicalFormula.ParseFormula("H3PO4");
        var ribose = ChemicalFormula.ParseFormula("C5H10O5");
        var a = ChemicalFormula.ParseFormula("C5H5N5");
        var c = ChemicalFormula.ParseFormula("C4H5N3O");
        var g = ChemicalFormula.ParseFormula("C5H5N5O");
        var u = ChemicalFormula.ParseFormula("C4H4N2O2");

        var amp = a + ribose - water + phosphate - water;
        var cmp = c + ribose - water + phosphate - water;
        var gmp = g + ribose - water + phosphate - water;
        var ump = u + ribose - water + phosphate - water;

        var combined = amp + cmp + gmp + ump;

        var averageC = combined.CountWithIsotopes(PeriodicTable.GetElement("C")) / 4.0;
        var averageH = combined.CountWithIsotopes(PeriodicTable.GetElement("H")) / 4.0;
        var averageO = combined.CountWithIsotopes(PeriodicTable.GetElement("O")) / 4.0;
        var averageN = combined.CountWithIsotopes(PeriodicTable.GetElement("N")) / 4.0;
        var averageP = combined.CountWithIsotopes(PeriodicTable.GetElement("P")) / 4.0;

        for (int i = 0; i < NumAveraginesToGenerate; i++)
        {
            double averagineMultiplier = (i + 1) / 4.0;
            ChemicalFormula chemicalFormula = new ChemicalFormula();
            chemicalFormula.Add("C", Convert.ToInt32(averageC * averagineMultiplier));
            chemicalFormula.Add("H", Convert.ToInt32(averageH * averagineMultiplier));
            chemicalFormula.Add("O", Convert.ToInt32(averageO * averagineMultiplier));
            chemicalFormula.Add("N", Convert.ToInt32(averageN * averagineMultiplier));
            chemicalFormula.Add("P", Convert.ToInt32(averageP * averagineMultiplier));

            {
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
}