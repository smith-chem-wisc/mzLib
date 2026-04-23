using Chemistry;
using NUnit.Framework;
using Omics;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Test.Omics;

[TestFixture]
[ExcludeFromCodeCoverage]
public class EmpiricalAverageResidueTests
{
    [Test]
    public void Constructor_Composition_PopulatesExpectedCollectionsAndAccessors()
    {
        var composition = new AverageResidueComposition(4.8, 7.7, 1.4, 1.3, 0.2);
        var model = new EmpiricalAverageResidue(composition);

        Assert.That(model.AllMasses.Length, Is.EqualTo(1500));
        Assert.That(model.AllIntensities.Length, Is.EqualTo(1500));
        Assert.That(model.MostIntenseMasses.Length, Is.EqualTo(1500));
        Assert.That(model.DiffToMonoisotopic.Length, Is.EqualTo(1500));

        int index = 300;
        double[] masses = model.GetAllTheoreticalMasses(index);
        double[] intensities = model.GetAllTheoreticalIntensities(index);

        Assert.That(masses, Is.Not.Null.And.Not.Empty);
        Assert.That(intensities, Is.Not.Null.And.Not.Empty);
        Assert.That(masses.Length, Is.EqualTo(intensities.Length));
        Assert.That(intensities.First(), Is.GreaterThanOrEqualTo(intensities.Last()));
        Assert.That(masses[0], Is.EqualTo(model.MostIntenseMasses[index]).Within(1e-9));
        Assert.That(model.GetDiffToMonoisotopic(index), Is.GreaterThanOrEqualTo(0));
        Assert.That(model.GetMostIntenseMassIndex(model.MostIntenseMasses[index]), Is.EqualTo(index));
    }

    [Test]
    public void Constructor_WithSetMods_MatchesModelFromManuallyComputedAverageComposition()
    {
        var precursors = new List<IBioPolymerWithSetMods>
        {
            new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>()),
            new PeptideWithSetModifications("MKWVTFISLL", new Dictionary<string, Modification>()),
            new PeptideWithSetModifications("ACDEFGHIK", new Dictionary<string, Modification>())
        };

        AverageResidueComposition expectedComposition = ComputeAverageComposition(precursors);
        var fromSetMods = new EmpiricalAverageResidue(precursors);
        var fromComposition = new EmpiricalAverageResidue(expectedComposition);

        AssertModelsMatchAtIndices(fromSetMods, fromComposition, new[] { 0, 10, 100, 500, 1499 });
    }

    [Test]
    public void Constructor_BioPolymersAndDigestionParams_MatchesConstructorFromDigestedSetMods()
    {
        var proteins = new List<IBioPolymer>
        {
            new Protein("MPEPTIDER", "P1"),
            new Protein("ACDEFGHIKLMNPQRST", "P2")
        };

        var digestionParams = new DigestionParams("top-down");
        var fixedMods = new List<Modification>();
        var variableMods = new List<Modification>();

        var fromBioPolymers = new EmpiricalAverageResidue(proteins, digestionParams, fixedMods, variableMods);
        IEnumerable<IBioPolymerWithSetMods> digested = proteins.SelectMany(p => p.Digest(digestionParams, fixedMods, variableMods));
        var fromSetMods = new EmpiricalAverageResidue(digested);

        AssertModelsMatchAtIndices(fromBioPolymers, fromSetMods, new[] { 1, 25, 250, 1000, 1499 });
    }

    private static AverageResidueComposition ComputeAverageComposition(IEnumerable<IBioPolymerWithSetMods> withSetMods)
    {
        var totalChemicalFormula = new ChemicalFormula();
        long totalResidues = 0;

        foreach (IBioPolymerWithSetMods bioPolymer in withSetMods)
        {
            totalChemicalFormula += bioPolymer.ThisChemicalFormula;
            totalResidues += bioPolymer.Length;
        }

        return new AverageResidueComposition(
            totalChemicalFormula.CountWithIsotopes(PeriodicTable.GetElement("C")) / (double)totalResidues,
            totalChemicalFormula.CountWithIsotopes(PeriodicTable.GetElement("H")) / (double)totalResidues,
            totalChemicalFormula.CountWithIsotopes(PeriodicTable.GetElement("O")) / (double)totalResidues,
            totalChemicalFormula.CountWithIsotopes(PeriodicTable.GetElement("N")) / (double)totalResidues,
            totalChemicalFormula.CountWithIsotopes(PeriodicTable.GetElement("P")) / (double)totalResidues);
    }

    private static void AssertModelsMatchAtIndices(
        EmpiricalAverageResidue expected,
        EmpiricalAverageResidue observed,
        IEnumerable<int> indices)
    {
        foreach (int index in indices)
        {
            Assert.That(observed.MostIntenseMasses[index], Is.EqualTo(expected.MostIntenseMasses[index]).Within(1e-9));
            Assert.That(observed.GetDiffToMonoisotopic(index), Is.EqualTo(expected.GetDiffToMonoisotopic(index)).Within(1e-9));

            double[] expectedIntensities = expected.GetAllTheoreticalIntensities(index);
            double[] observedIntensities = observed.GetAllTheoreticalIntensities(index);
            Assert.That(observedIntensities.Length, Is.EqualTo(expectedIntensities.Length));
            Assert.That(observedIntensities[0], Is.EqualTo(expectedIntensities[0]).Within(1e-9));
            Assert.That(observedIntensities[^1], Is.EqualTo(expectedIntensities[^1]).Within(1e-9));
        }
    }
}
