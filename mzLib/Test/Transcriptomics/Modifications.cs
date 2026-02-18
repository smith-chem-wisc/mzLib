using NUnit.Framework;
using Omics.Modifications;
using System.Linq;
using Chemistry;
using Omics.Fragmentation;

namespace Test.Transcriptomics;

[TestFixture]
public static class Modifications
{
    [Test]
    public static void BackboneModification_LoadsCorrectly()
    {
        // First mod should be phosphorothiolate
        // ID   Phosphorothioate
        // TG   X
        // FT   Backbone
        // PP   Anywhere.
        // MT   Common Variable
        // CF   SO-1
        // BM   w,x,c,d                 # ⭐ Affects these fragment types
        // 

        var mod = Mods.AllRnaModsList.FirstOrDefault(m => m is BackboneModification);
        Assert.That(mod, Is.Not.Null, "Backbone modification should be loaded.");
        Assert.That(mod.ChemicalFormula, Is.EqualTo(ChemicalFormula.ParseFormula("SO-1")), "Backbone modification should have correct chemical formula.");
        Assert.That(mod.LocationRestriction, Is.EqualTo("Anywhere."), "Backbone modification should have correct location restriction.");
        Assert.That(mod.ModificationType, Is.EqualTo("Common Variable"), "Backbone modification should have correct modification type.");

        var backboneMod = mod as BackboneModification;
        Assert.That(backboneMod.ProductsContainingModMass, Is.Not.Null, "Backbone modification should have ProductsContainingModMass defined.");

        // Check that the correct product types are listed as containing the modification mass shift
        var expectedProductsContainingShift = new ProductType[] { ProductType.c, ProductType.d, ProductType.w, ProductType.x };

        foreach (var prod in expectedProductsContainingShift)
        {
            var types = prod.GetFragmentFamilyMembers();
            foreach (var type in types)
            {
                Assert.That(backboneMod.ProductsContainingModMass, Does.Contain(type), $"Backbone modification should contain mass shift for {type}.");
            }
        }
    }
}
