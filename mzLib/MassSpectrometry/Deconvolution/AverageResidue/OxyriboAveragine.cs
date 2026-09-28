using System;
using System.Collections.Generic;
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

    /// <summary>
    /// The average number of each element in a single ribonucleotide averagine unit
    /// (e.g. {'C', 9.5}), determined by averaging the four canonical ribonucleotide monophosphates.
    /// </summary>
    private static readonly Dictionary<char, double> AverageComposition = BuildAverageComposition();

    public override int GetMostIntenseMassIndex(double testMass) => MostIntenseMasses.GetClosestIndex(testMass);

    public override double[] GetAllTheoreticalMasses(int index) => AllMasses[index];

    public override double[] GetAllTheoreticalIntensities(int index) => AllIntensities[index];

    public override double GetDiffToMonoisotopic(int index) => DiffToMonoisotopic[index];

    /// <summary>
    /// Returns the average elemental composition of a single ribonucleotide averagine unit as a map
    /// from element symbol to average atom count (e.g. {'C', 9.5}). Callers can scale this by a mass
    /// to derive an approximate chemical formula for an arbitrary species.
    /// </summary>
    public Dictionary<char, double> GetAverageChemicalFormula() => new(AverageComposition);

    /// <summary>
    /// Builds the average ribonucleotide composition by counting atoms in the four canonical RNA
    /// bases. This is not the best approach and future work should refine these numbers. One possible
    /// approach is to also incorporate the residue frequency in RNA sequences.
    /// </summary>
    private static Dictionary<char, double> BuildAverageComposition()
    {
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

        return new Dictionary<char, double>
        {
            { 'C', combined.CountWithIsotopes(PeriodicTable.GetElement("C")) / 4.0 },
            { 'H', combined.CountWithIsotopes(PeriodicTable.GetElement("H")) / 4.0 },
            { 'O', combined.CountWithIsotopes(PeriodicTable.GetElement("O")) / 4.0 },
            { 'N', combined.CountWithIsotopes(PeriodicTable.GetElement("N")) / 4.0 },
            { 'P', combined.CountWithIsotopes(PeriodicTable.GetElement("P")) / 4.0 },
        };
    }

    static OxyriboAveragine()
    {
        for (int i = 0; i < NumAveraginesToGenerate; i++)
        {
            double averagineMultiplier = (i + 1) / 4.0;
            ChemicalFormula chemicalFormula = new ChemicalFormula();
            foreach (var (element, count) in AverageComposition)
            {
                chemicalFormula.Add(element.ToString(), Convert.ToInt32(count * averagineMultiplier));
            }

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