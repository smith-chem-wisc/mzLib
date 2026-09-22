using System.Collections.Generic;
using System.Linq;
using Chemistry;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Transcriptomics;
using Transcriptomics.Digestion;

namespace Test.Omics;

[TestFixture]
public class FixedModificationDigestionTests
{
    [Test]
    public void RnaDigestAppliesParentFixedBeforeVariableAndGlobalFixedMods()
    {
        var anchored = CreateModification("Anchored", 'U');
        var variable = CreateModification("Variable", 'U');
        var globalFixed = CreateModification("Global fixed", 'U');
        var rna = new RNA(
            "GUACUG",
            oneBasedFixedModifications: new Dictionary<int, Modification>
            {
                [2] = anchored,
            });

        var products = rna.Digest(
                new RnaDigestionParams { MaxMods = 2 },
                [globalFixed],
                [variable])
            .ToList();

        Assert.That(products, Is.Not.Empty);
        var productsContainingFixedResidue = products
            .Where(product => product.OneBasedStartResidue <= 2 && product.OneBasedEndResidue >= 2)
            .ToList();
        Assert.That(productsContainingFixedResidue, Is.Not.Empty);
        Assert.That(productsContainingFixedResidue.Any(product =>
            product.AllModsOneIsNterminus[2 - product.OneBasedStartResidue + 2] == anchored), Is.True);
        Assert.That(products.Any(product => product.NumFixedMods >= 1), Is.True);
    }

    [Test]
    public void ProteinDigestAppliesParentFixedBeforeVariableAndGlobalFixedMods()
    {
        var anchored = CreateModification("Anchored", 'P');
        var variable = CreateModification("Variable", 'P');
        var globalFixed = CreateModification("Global fixed", 'P');
        var protein = new Protein(
            "MPEPTIDE",
            "P1",
            oneBasedFixedModifications: new Dictionary<int, Modification>
            {
                [2] = anchored,
            });

        var products = protein.Digest(
                new DigestionParams("top-down") { MaxMods = 2 },
                [globalFixed],
                [variable])
            .ToList();

        Assert.That(products, Is.Not.Empty);
        var productsContainingFixedResidue = products
            .Where(product => product.OneBasedStartResidueInProtein <= 2 && product.OneBasedEndResidueInProtein >= 2)
            .ToList();
        Assert.That(productsContainingFixedResidue, Is.Not.Empty);
        Assert.That(productsContainingFixedResidue.Any(product =>
            product.AllModsOneIsNterminus[2 - product.OneBasedStartResidueInProtein + 2] == anchored), Is.True);
        Assert.That(products.Any(product => product.NumFixedMods >= 1), Is.True);
    }

    private static Modification CreateModification(string id, char target)
    {
        Assert.That(ModificationMotif.TryGetMotif(target.ToString(), out var motif), Is.True);
        return new Modification(
            _originalId: id,
            _modificationType: "Test",
            _target: motif,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("CH2"));
    }
}
