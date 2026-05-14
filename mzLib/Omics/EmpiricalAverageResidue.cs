using Chemistry;
using MassSpectrometry;
using Omics.Digestion;
using Omics.Modifications;

namespace Omics;

public record struct AverageResidueComposition(double C, double H, double O, double N, double P, double S, double Se);

public class EmpiricalAverageResidue : AverageResidue
{
    public readonly double[][] AllMasses = new double[NumAveraginesToGenerate][];
    public readonly double[][] AllIntensities = new double[NumAveraginesToGenerate][];
    public readonly double[] MostIntenseMasses = new double[NumAveraginesToGenerate];
    public readonly double[] DiffToMonoisotopic = new double[NumAveraginesToGenerate];
    public override int GetMostIntenseMassIndex(double testMass) => MostIntenseMasses.GetClosestIndex(testMass);

    public override double[] GetAllTheoreticalMasses(int index) => AllMasses[index];

    public override double[] GetAllTheoreticalIntensities(int index) => AllIntensities[index];

    public override double GetDiffToMonoisotopic(int index) => DiffToMonoisotopic[index];

    public EmpiricalAverageResidue(IEnumerable<IBioPolymerWithSetMods> withSetMods) : this(GetAverageChemicalFormula(withSetMods))
    {
    }

    public EmpiricalAverageResidue(IEnumerable<IBioPolymer> bioPolymers, IDigestionParams digParams, List<Modification> fixedMods, List<Modification> variableMods) : this(GetAverageChemicalFormula(bioPolymers, digParams, fixedMods, variableMods))
    {
    }

    public EmpiricalAverageResidue(AverageResidueComposition composition)
    {
        for (int i = 0; i < NumAveraginesToGenerate; i++)
        {
            double averagineMultiplier = (i + 1) / 4.0;

            var chemicalFormulaReg = new ChemicalFormula();
            chemicalFormulaReg.Add("C", Convert.ToInt32(composition.C * averagineMultiplier));
            chemicalFormulaReg.Add("H", Convert.ToInt32(composition.H * averagineMultiplier));
            chemicalFormulaReg.Add("O", Convert.ToInt32(composition.O * averagineMultiplier));
            chemicalFormulaReg.Add("N", Convert.ToInt32(composition.N * averagineMultiplier));
            chemicalFormulaReg.Add("P", Convert.ToInt32(composition.P * averagineMultiplier));
            chemicalFormulaReg.Add("S", Convert.ToInt32(composition.S * averagineMultiplier));
            chemicalFormulaReg.Add("Se", Convert.ToInt32(composition.Se * averagineMultiplier));

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

    private static AverageResidueComposition GetAverageChemicalFormula(IEnumerable<IBioPolymer> bioPolymers, IDigestionParams digParams, List<Modification> fixedMods, List<Modification> variableMods)
    {
        var withSetMods = bioPolymers.SelectMany(b => b.Digest(digParams, fixedMods, variableMods));
        return GetAverageChemicalFormula(withSetMods);
    }

    private static AverageResidueComposition GetAverageChemicalFormula(IEnumerable<IBioPolymerWithSetMods> withSetMods)
    {
        var totalChemicalFormula = new ChemicalFormula();
        long totalResidues = 0;
        foreach (var bioPolymer in withSetMods)
        {
            totalChemicalFormula += bioPolymer.ThisChemicalFormula;
            totalResidues += bioPolymer.Length;
        }

        double averageC = totalChemicalFormula.CountWithIsotopes(PeriodicTable.GetElement("C")) / (double)totalResidues;
        double averageH = totalChemicalFormula.CountWithIsotopes(PeriodicTable.GetElement("H")) / (double)totalResidues;
        double averageO = totalChemicalFormula.CountWithIsotopes(PeriodicTable.GetElement("O")) / (double)totalResidues;
        double averageN = totalChemicalFormula.CountWithIsotopes(PeriodicTable.GetElement("N")) / (double)totalResidues;
        double averageP = totalChemicalFormula.CountWithIsotopes(PeriodicTable.GetElement("P")) / (double)totalResidues;
        double averageS = totalChemicalFormula.CountWithIsotopes(PeriodicTable.GetElement("S")) / (double)totalResidues;
        double averageSe = totalChemicalFormula.CountWithIsotopes(PeriodicTable.GetElement("Se")) / (double)totalResidues;

        return new AverageResidueComposition(averageC, averageH, averageO, averageN, averageP, averageS, averageSe);
    }
}
