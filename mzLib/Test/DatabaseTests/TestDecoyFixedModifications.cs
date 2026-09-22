using System.Collections.Generic;
using System.Linq;
using Chemistry;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using Transcriptomics;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests;

[TestFixture]
public class TestDecoyFixedModifications
{
    [Test]
    public void ReverseRnaDecoyRemapsFixedModifications()
    {
        var fixedModification = CreateModification('U');
        var rna = new RNA(
            "GUACUG",
            "RNA1",
            oneBasedFixedModifications: new Dictionary<int, Modification>
            {
                [2] = fixedModification,
                [5] = fixedModification,
            });

        var decoy = RnaDecoyGenerator.GenerateDecoys([rna], DecoyType.Reverse, 1).Single();

        Assert.That(decoy.BaseSequence, Is.EqualTo("GUCAUG"));
        Assert.That(decoy.OneBasedFixedModifications.Keys, Is.EquivalentTo(new[] { 2, 5 }));
        Assert.That(decoy.OneBasedFixedModifications.Values, Is.All.EqualTo(fixedModification));
    }

    [Test]
    public void ReverseProteinDecoyRemapsFixedModifications()
    {
        var firstModification = CreateModification('P');
        var secondModification = CreateModification('T');
        var protein = new Protein(
            "MPEPTIDE",
            "P1",
            oneBasedFixedModifications: new Dictionary<int, Modification>
            {
                [2] = firstModification,
                [5] = secondModification,
            });

        var decoy = DecoyProteinGenerator.GenerateDecoys([protein], DecoyType.Reverse, 1).Single();

        Assert.That(decoy.BaseSequence, Is.EqualTo("MEDITPEP"));
        Assert.That(decoy.OneBasedFixedModifications.Keys, Is.EquivalentTo(new[] { 5, 8 }));
        Assert.That(decoy.OneBasedFixedModifications[8], Is.EqualTo(firstModification));
        Assert.That(decoy.OneBasedFixedModifications[5], Is.EqualTo(secondModification));
    }

    [Test]
    public void SlideProteinDecoyRemapsFixedModificationsWithLocalizedModifications()
    {
        var modification = CreateModification('P');
        var localizedModifications = new Dictionary<int, List<Modification>>
        {
            [2] = [modification],
        };
        var protein = new Protein(
            "MPEPTIDE",
            "P1",
            oneBasedModifications: localizedModifications,
            oneBasedFixedModifications: localizedModifications.ToDictionary(pair => pair.Key, pair => pair.Value[0]));

        var decoy = DecoyProteinGenerator.GenerateDecoys([protein], DecoyType.Slide, 1).Single();

        Assert.That(decoy.OneBasedFixedModifications.Keys,
            Is.EquivalentTo(decoy.OneBasedPossibleLocalizedModifications.Keys));
        Assert.That(decoy.OneBasedFixedModifications.Values, Is.All.EqualTo(modification));
    }

    private static Modification CreateModification(char target)
    {
        Assert.That(ModificationMotif.TryGetMotif(target.ToString(), out var motif), Is.True);
        return new Modification(
            _originalId: $"Fixed {target}",
            _modificationType: "Test",
            _target: motif,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("CH2"));
    }
}
