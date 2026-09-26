using Chemistry;
using System;
using System.Collections.Generic;
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

    /// <summary>
    /// The average number of each element in a single amino-acid averagine unit
    /// (e.g. {'C', 4.9384}). Magic numbers determined by https://pmc.ncbi.nlm.nih.gov/articles/PMC6166224/.
    /// </summary>
    private static readonly Dictionary<char, double> AverageComposition = new()
    {
        { 'C', 4.9384 },
        { 'H', 7.7583 },
        { 'O', 1.4773 },
        { 'N', 1.3577 },
        { 'S', 0.0417 },
    };

    public override int GetMostIntenseMassIndex(double testMass) => MostIntenseMasses.GetClosestIndex(testMass);
    public override double[] GetAllTheoreticalMasses(int index) => AllMasses[index];
    public override double[] GetAllTheoreticalIntensities(int index) => AllIntensities[index];
    public override double GetDiffToMonoisotopic(int index) => DiffToMonoisotopic[index];

    /// <summary>
    /// Returns the average elemental composition of a single averagine unit as a map from element
    /// symbol to average atom count (e.g. {'C', 4.9384}). Callers can scale this by a mass to derive
    /// an approximate chemical formula for an arbitrary species.
    /// </summary>
    public Dictionary<char, double> GetAverageChemicalFormula() => new(AverageComposition);

    static Averagine()
    {
        for (int i = 0; i < NumAveraginesToGenerate; i++)
        {
            double averagineMultiplier = (i + 1) / 2.0;
            //Console.Write("numAveragines = " + numAveragines);
            ChemicalFormula chemicalFormula = new ChemicalFormula();
            foreach (var (element, count) in AverageComposition)
            {
                chemicalFormula.Add(element.ToString(), Convert.ToInt32(count * averagineMultiplier));
            }

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