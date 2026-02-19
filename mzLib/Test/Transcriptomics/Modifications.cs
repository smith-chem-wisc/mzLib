using NUnit.Framework;
using Omics.Modifications;
using System.Linq;
using Chemistry;
using Omics.Fragmentation;
using Omics.Modifications.IO;
using System.Collections.Generic;
using System.IO;
using System;

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
        // BM   w,x,c,d  
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

    [Test]
    public static void BackboneModification_ToString_ContainsAllFields()
    {
        // Test that ToString properly includes the BM field
        var mod = Mods.AllRnaModsList.FirstOrDefault(m => m is BackboneModification) as BackboneModification;
        Assert.That(mod, Is.Not.Null);

        string modString = mod.ToString();
        
        // Check that all standard modification fields are present
        Assert.That(modString, Does.Contain("ID   "));
        Assert.That(modString, Does.Contain("MT   "));
        Assert.That(modString, Does.Contain("FT   "));
        Assert.That(modString, Does.Contain("TG   "));
        Assert.That(modString, Does.Contain("PP   "));
        Assert.That(modString, Does.Contain("CF   "));
        Assert.That(modString, Does.Contain("MM   "));
        
        // Check that backbone-specific BM field is present
        Assert.That(modString, Does.Contain("BM   "));
        
        // Check that the product types are listed
        Assert.That(modString, Does.Contain(ProductType.c.ToString()));
        Assert.That(modString, Does.Contain(ProductType.d.ToString()));
        Assert.That(modString, Does.Contain(ProductType.w.ToString()));
        Assert.That(modString, Does.Contain(ProductType.x.ToString()));
    }

    [Test]
    public static void BackboneModification_ToStringRoundTrip()
    {
        // Test that a BackboneModification can be written to string and read back successfully
        var originalMod = Mods.AllRnaModsList.FirstOrDefault(m => m is BackboneModification) as BackboneModification;
        Assert.That(originalMod, Is.Not.Null);

        // Convert to string
        string modString = originalMod.ToString();
        
        // Add the terminator that the reader expects
        modString += Environment.NewLine + "//" + Environment.NewLine;
        
        // Create a temporary file for testing
        string tempFile = Path.GetTempFileName();
        try
        {
            File.WriteAllText(tempFile, modString);
            
            // Read back using the modification loader
            var readMods = ModificationLoader.ReadModsFromFile(tempFile, new Dictionary<string, int>(), out var filteredMods);
            var readMod = readMods.FirstOrDefault() as BackboneModification;
            
            Assert.That(readMod, Is.Not.Null, "Should read back as BackboneModification");
            Assert.That(readMod.OriginalId, Is.EqualTo(originalMod.OriginalId));
            Assert.That(readMod.ModificationType, Is.EqualTo(originalMod.ModificationType));
            Assert.That(readMod.ChemicalFormula.Formula, Is.EqualTo(originalMod.ChemicalFormula.Formula));
            Assert.That(readMod.MonoisotopicMass, Is.EqualTo(originalMod.MonoisotopicMass).Within(0.0001));
            Assert.That(readMod.ProductsContainingModMass.Length, Is.EqualTo(originalMod.ProductsContainingModMass.Length));
            
            // Check that all product types match
            foreach (var prodType in originalMod.ProductsContainingModMass)
            {
                Assert.That(readMod.ProductsContainingModMass, Does.Contain(prodType));
            }
        }
        finally
        {
            // Clean up
            if (File.Exists(tempFile))
                File.Delete(tempFile);
        }
    }

    [Test]
    public static void BackboneModification_ValidModification_WithProductTypes()
    {
        // Test that a backbone modification with valid product types is considered valid
        ModificationMotif.TryGetMotif("X", out var motif);
        var validMod = new BackboneModification(
            _originalId: "TestMod",
            _modificationType: "TestType",
            _target: motif,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("H2O"),
            productsContainingShift: new[] { ProductType.c, ProductType.d, ProductType.w, ProductType.x });

        Assert.That(validMod.ValidModification, Is.True);
    }

    [Test]
    public static void BackboneModification_InvalidModification_WithEmptyProductTypes()
    {
        // Test that a backbone modification with empty product types array is invalid
        ModificationMotif.TryGetMotif("X", out var motif);
        var invalidMod = new BackboneModification(
            _originalId: "TestMod",
            _modificationType: "TestType",
            _target: motif,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("H2O"),
            productsContainingShift: Array.Empty<ProductType>());

        Assert.That(invalidMod.ValidModification, Is.False);
    }

    [Test]
    public static void BackboneModification_InvalidModification_WithForbiddenPairs()
    {
        // Test that backbone modification is invalid if it contains both sides of forbidden pairs
        // Forbidden pairs: (a,w), (b,x), (c,y), (d,z)
        ModificationMotif.TryGetMotif("X", out var motif);
        
        // Test a-w pair
        var invalidModAW = new BackboneModification(
            _originalId: "TestMod",
            _modificationType: "TestType",
            _target: motif,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("H2O"),
            productsContainingShift: new[] { ProductType.a, ProductType.w });
        Assert.That(invalidModAW.ValidModification, Is.False, "Mod with a and w should be invalid");

        // Test b-x pair
        var invalidModBX = new BackboneModification(
            _originalId: "TestMod",
            _modificationType: "TestType",
            _target: motif,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("H2O"),
            productsContainingShift: new[] { ProductType.b, ProductType.x });
        Assert.That(invalidModBX.ValidModification, Is.False, "Mod with b and x should be invalid");

        // Test c-y pair
        var invalidModCY = new BackboneModification(
            _originalId: "TestMod",
            _modificationType: "TestType",
            _target: motif,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("H2O"),
            productsContainingShift: new[] { ProductType.c, ProductType.y });
        Assert.That(invalidModCY.ValidModification, Is.False, "Mod with c and y should be invalid");

        // Test d-z pair
        var invalidModDZ = new BackboneModification(
            _originalId: "TestMod",
            _modificationType: "TestType",
            _target: motif,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("H2O"),
            productsContainingShift: new[] { ProductType.d, ProductType.z });
        Assert.That(invalidModDZ.ValidModification, Is.False, "Mod with d and z should be invalid");
    }

    [Test]
    public static void BackboneModification_ValidModification_WithoutForbiddenPairs()
    {
        // Test that backbone modification is valid when it doesn't contain forbidden pairs
        ModificationMotif.TryGetMotif("X", out var motif);
        
        var validMod = new BackboneModification(
            _originalId: "TestMod",
            _modificationType: "TestType",
            _target: motif,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("H2O"),
            productsContainingShift: new[] { ProductType.c, ProductType.d, ProductType.w, ProductType.x });
        
        Assert.That(validMod.ValidModification, Is.True);
    }

    [Test]
    public static void BackboneModification_InvalidModification_MissingBaseFields()
    {
        // Test that backbone modification respects base class validation
        ModificationMotif.TryGetMotif("X", out var motif);
        
        // Missing modification type
        var invalidMod1 = new BackboneModification(
            _originalId: "TestMod",
            _modificationType: null,
            _target: motif,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("H2O"),
            productsContainingShift: new[] { ProductType.c, ProductType.d });
        Assert.That(invalidMod1.ValidModification, Is.False);

        // Missing target
        var invalidMod2 = new BackboneModification(
            _originalId: "TestMod",
            _modificationType: "TestType",
            _target: null,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("H2O"),
            productsContainingShift: new[] { ProductType.c, ProductType.d });
        Assert.That(invalidMod2.ValidModification, Is.False);
    }

    [Test]
    public static void BackboneModification_ProductsContainingModMass_IsSorted()
    {
        // Test that the product types are sorted after construction
        ModificationMotif.TryGetMotif("X", out var motif);
        var mod = new BackboneModification(
            _originalId: "TestMod",
            _modificationType: "TestType",
            _target: motif,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("H2O"),
            productsContainingShift: new[] { ProductType.x, ProductType.c, ProductType.w, ProductType.d });

        // Check that array is sorted
        for (int i = 0; i < mod.ProductsContainingModMass.Length - 1; i++)
        {
            Assert.That(mod.ProductsContainingModMass[i], Is.LessThan(mod.ProductsContainingModMass[i + 1]));
        }
    }
}
