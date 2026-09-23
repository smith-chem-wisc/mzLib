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

    [Test]
    public void ProteinDigestAppliesParentFixedModsOnlyWhenResidueIsInsideProduct()
    {
        var modAtResidueFour = CreateModification("Anchored four", 'R');
        var modAtResidueSix = CreateModification("Anchored six", 'K');
        var protein = new Protein(
            "AKTRTKTR",
            "P1",
            oneBasedFixedModifications: new Dictionary<int, Modification>
            {
                [4] = modAtResidueFour,
                [6] = modAtResidueSix,
            });

        var products = protein.Digest(
                new DigestionParams("trypsin", maxMissedCleavages: 0, minPeptideLength: 1,
                    initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain),
                new List<Modification>(),
                new List<Modification>())
            .ToList();

        Assert.That(products, Is.Not.Empty);

        var ak = products.Single(product => product.OneBasedStartResidueInProtein == 1);
        Assert.That(ak.AllModsOneIsNterminus.Values, Does.Not.Contain(modAtResidueFour));
        Assert.That(ak.AllModsOneIsNterminus.Values, Does.Not.Contain(modAtResidueSix));

        var tr = products.Single(product => product.OneBasedStartResidueInProtein == 3);
        Assert.That(tr.AllModsOneIsNterminus[4 - 3 + 2], Is.EqualTo(modAtResidueFour));
        Assert.That(tr.AllModsOneIsNterminus.Values, Does.Not.Contain(modAtResidueSix));

        var tk = products.Single(product => product.OneBasedStartResidueInProtein == 5);
        Assert.That(tk.AllModsOneIsNterminus[6 - 5 + 2], Is.EqualTo(modAtResidueSix));
        Assert.That(tk.AllModsOneIsNterminus.Values, Does.Not.Contain(modAtResidueFour));

        var secondTr = products.Single(product => product.OneBasedStartResidueInProtein == 7);
        Assert.That(secondTr.AllModsOneIsNterminus.Values, Does.Not.Contain(modAtResidueFour));
        Assert.That(secondTr.AllModsOneIsNterminus.Values, Does.Not.Contain(modAtResidueSix));
    }

    [Test]
    public void RnaDigestAppliesParentFixedModsOnlyWhenResidueIsInsideProduct()
    {
        var modAtResidueFour = CreateModification("Anchored four", 'A');
        var modAtResidueFive = CreateModification("Anchored five", 'C');
        var rna = new RNA(
            "GAUACG",
            oneBasedFixedModifications: new Dictionary<int, Modification>
            {
                [4] = modAtResidueFour,
                [5] = modAtResidueFive,
            });

        var products = rna.Digest(
                new RnaDigestionParams("RNase U2", minLength: 1),
                new List<Modification>(),
                new List<Modification>())
            .Cast<OligoWithSetMods>()
            .ToList();

        Assert.That(products, Is.Not.Empty);

        var g = products.Single(product => product.OneBasedStartResidue == 1);
        Assert.That(g.AllModsOneIsNterminus.Values, Does.Not.Contain(modAtResidueFour));
        Assert.That(g.AllModsOneIsNterminus.Values, Does.Not.Contain(modAtResidueFive));

        var a = products.Single(product => product.OneBasedStartResidue == 2);
        Assert.That(a.AllModsOneIsNterminus.Values, Does.Not.Contain(modAtResidueFour));
        Assert.That(a.AllModsOneIsNterminus.Values, Does.Not.Contain(modAtResidueFive));

        var ua = products.Single(product => product.OneBasedStartResidue == 3);
        Assert.That(ua.AllModsOneIsNterminus[4 - 3 + 2], Is.EqualTo(modAtResidueFour));
        Assert.That(ua.AllModsOneIsNterminus.Values, Does.Not.Contain(modAtResidueFive));

        var cg = products.Single(product => product.OneBasedStartResidue == 5);
        Assert.That(cg.AllModsOneIsNterminus[5 - 5 + 2], Is.EqualTo(modAtResidueFive));
        Assert.That(cg.AllModsOneIsNterminus.Values, Does.Not.Contain(modAtResidueFour));
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
