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
    #region BackboneModification Tests

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

    #endregion

    #region BaseModification Tests

    [Test]
    public static void BaseModification_LoadsCorrectly_Suppressed()
    {
        // Test 2'-O-Methyladenosine which should suppress base loss
        var mod = Mods.AllRnaModsList.FirstOrDefault(m => 
            m is BaseModification && 
            m.OriginalId != null && 
            m.OriginalId.Contains("2'-O-Methyl"));
        
        Assert.That(mod, Is.Not.Null, "2'-O-Methyl modification should be loaded as BaseModification.");
        
        var baseMod = mod as BaseModification;
        Assert.That(baseMod.BaseLossType, Is.EqualTo(BaseLossBehavior.Suppressed), 
            "2'-O-Methyl should have Suppressed base loss behavior.");
        Assert.That(baseMod.ChemicalFormula, Is.EqualTo(ChemicalFormula.ParseFormula("C1H2")), 
            "2'-O-Methyl should have correct chemical formula.");
    }

    [Test]
    [TestCase("2'-O-Methyl", "CH2", null, BaseLossBehavior.Suppressed)]
    [TestCase("N6-methyladenosine", "CH2", "CH2", BaseLossBehavior.Modified)]
    [TestCase("N6,2'-O-dimethyladenosine", "C2H4", "CH2", BaseLossBehavior.Suppressed)]

    public static void BaseModification_LoadsCorrectly(string modNameStart, string expectedCf, string expectedBaseLossCf, BaseLossBehavior expectedBehavior)
    {
        var mod = Mods.AllRnaModsList.FirstOrDefault(m => 
            m is BaseModification && 
            m.OriginalId != null && 
            m.OriginalId.StartsWith(modNameStart)) as BaseModification;
        
        Assert.That(mod, Is.Not.Null, $"{modNameStart} modification should be loaded as BaseModification.");
        
        Assert.That(mod.ChemicalFormula.Formula, Is.EqualTo(expectedCf), 
            $"{modNameStart} should have correct chemical formula.");
        Assert.That(mod.BaseLossType, Is.EqualTo(expectedBehavior), 
            $"{modNameStart} should have correct base loss behavior.");
        Assert.That(mod.BaseLossModification?.Formula, Is.EqualTo(expectedBaseLossCf), 
            $"{modNameStart} should have correct base loss modification formula.");
    }

    [Test]
    public static void BaseModification_LoadsCorrectly_Modified()
    {
        // Test m6A which should modify base loss mass
        var mod = Mods.AllRnaModsList.FirstOrDefault(m => 
            m is BaseModification && 
            m.OriginalId != null && 
            m.OriginalId.Contains("N6-methyl"));
        
        Assert.That(mod, Is.Not.Null, "N6-methyladenosine should be loaded as BaseModification.");
        
        var baseMod = mod as BaseModification;
        Assert.That(baseMod.BaseLossType, Is.EqualTo(BaseLossBehavior.Modified), 
            "N6-methyladenosine should have Modified base loss behavior.");
        Assert.That(baseMod.BaseLossModification, Is.Not.Null, 
            "N6-methyladenosine should have BaseLossModification formula defined.");
        Assert.That(baseMod.BaseLossModification.Formula, Is.EqualTo(ChemicalFormula.ParseFormula("C1H2").Formula),
            "N6-methyladenosine base loss modification should be CH2.");
    }

    [Test]
    public static void BaseModification_ToString_ContainsAllFields()
    {
        // Test that ToString properly includes the BL field
        var mod = Mods.AllRnaModsList.FirstOrDefault(m => 
            m is BaseModification && 
            m.OriginalId != null && 
            m.OriginalId.Contains("2'-O-Methyl")) as BaseModification;
        Assert.That(mod, Is.Not.Null);

        string modString = mod.ToString();
        
        // Check that all standard modification fields are present
        Assert.That(modString, Does.Contain("ID   "));
        Assert.That(modString, Does.Contain("MT   "));
        Assert.That(modString, Does.Contain("TG   "));
        Assert.That(modString, Does.Contain("PP   "));
        Assert.That(modString, Does.Contain("CF   "));
        
        // Check that base-specific BL field is present
        Assert.That(modString, Does.Contain("BL   "));
        Assert.That(modString, Does.Contain("Suppressed"));
    }

    [Test]
    public static void BaseModification_ToString_ModifiedWithFormula()
    {
        var mod = Mods.AllRnaModsList.FirstOrDefault(m => 
            m is BaseModification && 
            m.OriginalId != null && 
            m.OriginalId.Contains("N6-methyl")) as BaseModification;
        Assert.That(mod, Is.Not.Null);

        string modString = mod.ToString();
        
        // Check that BL field includes the formula for Modified type
        Assert.That(modString, Does.Contain("BL   "));
        Assert.That(modString, Does.Contain("Modified"));
        Assert.That(modString, Does.Contain(":"));
    }

    [Test]
    public static void BaseModification_ToStringRoundTrip_Suppressed()
    {
        // Test that a BaseModification can be written to string and read back successfully
        var originalMod = Mods.AllRnaModsList.FirstOrDefault(m => 
            m is BaseModification && 
            m.OriginalId != null && 
            m.OriginalId.Contains("2'-O-Methyl")) as BaseModification;
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
            var readMod = readMods.FirstOrDefault() as BaseModification;
            
            Assert.That(readMod, Is.Not.Null, "Should read back as BaseModification");
            Assert.That(readMod.OriginalId, Is.EqualTo(originalMod.OriginalId));
            Assert.That(readMod.ModificationType, Is.EqualTo(originalMod.ModificationType));
            Assert.That(readMod.ChemicalFormula.Formula, Is.EqualTo(originalMod.ChemicalFormula.Formula));
            Assert.That(readMod.MonoisotopicMass, Is.EqualTo(originalMod.MonoisotopicMass).Within(0.0001));
            Assert.That(readMod.BaseLossType, Is.EqualTo(originalMod.BaseLossType));
        }
        finally
        {
            // Clean up
            if (File.Exists(tempFile))
                File.Delete(tempFile);
        }
    }

    [Test]
    public static void BaseModification_ToStringRoundTrip_Modified()
    {
        var originalMod = Mods.AllRnaModsList.FirstOrDefault(m => 
            m is BaseModification && 
            m.OriginalId != null && 
            m.OriginalId.Contains("N6-methyl")) as BaseModification;
        Assert.That(originalMod, Is.Not.Null);

        string modString = originalMod.ToString();
        modString += Environment.NewLine + "//" + Environment.NewLine;
        
        string tempFile = Path.GetTempFileName();
        try
        {
            File.WriteAllText(tempFile, modString);
            
            var readMods = ModificationLoader.ReadModsFromFile(tempFile, new Dictionary<string, int>(), out var filteredMods);
            var readMod = readMods.FirstOrDefault() as BaseModification;
            
            Assert.That(readMod, Is.Not.Null, "Should read back as BaseModification");
            Assert.That(readMod.OriginalId, Is.EqualTo(originalMod.OriginalId));
            Assert.That(readMod.BaseLossType, Is.EqualTo(originalMod.BaseLossType));
            Assert.That(readMod.BaseLossModification, Is.Not.Null);
            Assert.That(readMod.BaseLossModification.Formula, Is.EqualTo(originalMod.BaseLossModification.Formula));
        }
        finally
        {
            if (File.Exists(tempFile))
                File.Delete(tempFile);
        }
    }

    [Test]
    public static void BaseModification_ValidModification_Default()
    {
        // Test that a base modification with Default behavior is valid
        ModificationMotif.TryGetMotif("A", out var motif);
        var validMod = new BaseModification(
            _originalId: "TestMod",
            _modificationType: "TestType",
            _target: motif,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("H2O"),
            baseLossType: BaseLossBehavior.Default);

        Assert.That(validMod.ValidModification, Is.True);
    }

    [Test]
    public static void BaseModification_ValidModification_Suppressed()
    {
        // Test that a base modification with Suppressed behavior is valid
        ModificationMotif.TryGetMotif("U", out var motif);
        var validMod = new BaseModification(
            _originalId: "TestMod",
            _modificationType: "TestType",
            _target: motif,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("CH2"),
            baseLossType: BaseLossBehavior.Suppressed);

        Assert.That(validMod.ValidModification, Is.True);
    }

    [Test]
    public static void BaseModification_ValidModification_ModifiedWithFormula()
    {
        // Test that a base modification with Modified behavior and formula is valid
        ModificationMotif.TryGetMotif("A", out var motif);
        var validMod = new BaseModification(
            _originalId: "TestMod",
            _modificationType: "TestType",
            _target: motif,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("CH2"),
            baseLossType: BaseLossBehavior.Modified,
            baseLossModification: ChemicalFormula.ParseFormula("CH2"));

        Assert.That(validMod.ValidModification, Is.True);
    }

    [Test]
    public static void BaseModification_InvalidModification_ModifiedWithoutFormula()
    {
        // Test that a base modification with Modified behavior but no formula is invalid
        ModificationMotif.TryGetMotif("A", out var motif);
        var invalidMod = new BaseModification(
            _originalId: "TestMod",
            _modificationType: "TestType",
            _target: motif,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("CH2"),
            baseLossType: BaseLossBehavior.Modified,
            baseLossModification: null);

        Assert.That(invalidMod.ValidModification, Is.False, 
            "Modified BaseLossBehavior without BaseLossModification formula should be invalid");
    }

    [Test]
    public static void BaseModification_InvalidModification_MissingBaseFields()
    {
        // Test that base modification respects base class validation
        ModificationMotif.TryGetMotif("A", out var motif);
        
        // Missing modification type
        var invalidMod1 = new BaseModification(
            _originalId: "TestMod",
            _modificationType: null,
            _target: motif,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("H2O"),
            baseLossType: BaseLossBehavior.Suppressed);
        Assert.That(invalidMod1.ValidModification, Is.False);

        // Missing target
        var invalidMod2 = new BaseModification(
            _originalId: "TestMod",
            _modificationType: "TestType",
            _target: null,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("H2O"),
            baseLossType: BaseLossBehavior.Suppressed);
        Assert.That(invalidMod2.ValidModification, Is.False);
    }

    #endregion
}
