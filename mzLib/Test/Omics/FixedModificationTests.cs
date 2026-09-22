using System.Collections.Generic;
using Chemistry;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using Transcriptomics;

namespace Test.Omics;

[TestFixture]
public class FixedModificationTests
{
    [Test]
    public void RnaStoresAndClonesOneBasedFixedModifications()
    {
        var modification = CreateModification('U');
        var fixedModifications = new Dictionary<int, Modification>
        {
            [3] = modification,
        };

        var rna = new RNA("GUACUG", oneBasedFixedModifications: fixedModifications);
        var clone = (RNA)rna.CloneWithNewSequenceAndMods("GUACUG");

        Assert.That(rna.OneBasedFixedModifications[3], Is.EqualTo(modification));
        Assert.That(clone.OneBasedFixedModifications[3], Is.EqualTo(modification));
    }

    [Test]
    public void ProteinStoresAndClonesOneBasedFixedModifications()
    {
        var modification = CreateModification('M');
        var fixedModifications = new Dictionary<int, Modification>
        {
            [1] = modification,
        };

        var protein = new Protein("MPEPTIDE", "P1", oneBasedFixedModifications: fixedModifications);
        var clone = (Protein)protein.CloneWithNewSequenceAndMods("MPEPTIDE");

        Assert.That(protein.OneBasedFixedModifications[1], Is.EqualTo(modification));
        Assert.That(clone.OneBasedFixedModifications[1], Is.EqualTo(modification));
    }

    private static Modification CreateModification(char target)
    {
        Assert.That(ModificationMotif.TryGetMotif(target.ToString(), out var motif), Is.True);
        return new Modification(
            _originalId: "Fixed test modification",
            _modificationType: "Test",
            _target: motif,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("CH2"));
    }
}
