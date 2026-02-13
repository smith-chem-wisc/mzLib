using System;
using Chemistry;
using NUnit.Framework;
using Omics.Modifications;
using System.Collections.Generic;
using System.Linq;
using static TorchSharp.torch.optim.lr_scheduler.impl.CyclicLR;

namespace Test.Omics;

[TestFixture]
public static class ModificationTest
{
    public class ModificationTestCase(Modification modification, ModificationNamingConvention convention, bool isProteinMod)
    {
        public Modification Modification { get; set; } = modification;
        public ModificationNamingConvention Convention { get; set; } = convention;
        public bool IsProteinMod { get; set; } = isProteinMod;
    }

    public static List<ModificationTestCase> GetAllModificationsTestCases()
    {
        List<ModificationTestCase> testCases =
        [
            new (Mods.AllKnownProteinModsDictionary["DVFQQQTGG (SUMO-2/3 Site human) on D"], ModificationNamingConvention.MetaMorpheus, true),
            new (Mods.AllKnownProteinModsDictionary["Phosphorylation on T"], ModificationNamingConvention.MetaMorpheus, true),
            new (Mods.AllKnownRnaModsDictionary["MethoxyEthoxylation on G"], ModificationNamingConvention.MetaMorpheus_Rna, false),
            new (Mods.AllKnownProteinModsDictionary["(3S)-3-hydroxyaspartate on D"], ModificationNamingConvention.UniProt, true),
            new (Mods.AllKnownProteinModsDictionary["ICAT-D:2H(8) on C"], ModificationNamingConvention.Unimod, true),
        ];

        return testCases;
    }   

    public static Modification CustomDummyMod { get; private set; }
    static ModificationTest()
    {
        ModificationMotif.TryGetMotif("M", out var motif);
        CustomDummyMod = new Modification("TestMod", "acc", "custom", "custom", motif, "Unassigned.", ChemicalFormula.ParseFormula("CF3"));
    }


    [Test]
    public static void UniprotModsAreAllUniprot()
    {
        var mods = Mods.UniprotModifications;
        foreach (var mod in mods)
        {
            Assert.That(mod.ModificationType, Is.EqualTo("UniProt"));
        }
    }

    [Test]
    public static void UnimodModsAreAllUnimod()
    {
        var mods = Mods.UnimodModifications;
        foreach (var mod in mods)
        {
            Assert.That(mod.ModificationType, Is.EqualTo("Unimod"));
            Assert.That(mod.DatabaseReference.Count, Is.EqualTo(1));
            Assert.That(mod.DatabaseReference.First().Key, Is.EqualTo("Unimod"));
        }
    }

    [Test]
    public static void MetaMorpheusModsAreNotUniprotOrUnimod()
    {
        var mods = Mods.MetaMorpheusModifications;
        foreach (var mod in mods)
        {
            Assert.That(mod.ModificationType, Is.Not.EqualTo("UniProt"));
            Assert.That(mod.ModificationType, Is.Not.EqualTo("Unimod"));
        }
    }

    [Test]
    [TestCaseSource(nameof(GetAllModificationsTestCases))]
    public static void GetModification_Global(ModificationTestCase testCase)
    {
        var mod = Mods.GetModification(testCase.Modification.IdWithMotif);
        Assert.That(mod, Is.EqualTo(testCase.Modification));

        var modFromProteinMods = Mods.GetModification(testCase.Modification.IdWithMotif, true, false);
        var modFromRnaMods = Mods.GetModification(testCase.Modification.IdWithMotif, false, true);

        if (testCase.IsProteinMod)
        {
            Assert.That(modFromProteinMods, Is.EqualTo(mod));
            Assert.That(modFromRnaMods, Is.Null);
        }
        else
        {
            Assert.That(modFromProteinMods, Is.Null);
            Assert.That(modFromRnaMods, Is.EqualTo(mod));
        }
    }

    [Test]
    [TestCaseSource(nameof(GetAllModificationsTestCases))]
    public static void GetModification_ByConvention(ModificationTestCase testCase)
    {
        var mod = Mods.GetModification(testCase.Modification.IdWithMotif, testCase.Convention);
        Assert.That(mod, Is.EqualTo(testCase.Modification));

        ModificationNamingConvention wrongConvention;
        if (testCase.Convention == ModificationNamingConvention.MetaMorpheus)
            wrongConvention = ModificationNamingConvention.UniProt;
        else
            wrongConvention = ModificationNamingConvention.MetaMorpheus;

        var modFromWrongConvention = Mods.GetModification(testCase.Modification.IdWithMotif, wrongConvention);
        Assert.That(modFromWrongConvention, Is.Null);
    }

    [Test]
    public static void GetMods_Throws()
    {
        try
        {
            var mod = Mods.GetModification("Oxidation on M", false, false);
            Assert.Fail("Expected an exception to be thrown for a non-existent modification.");
        }
        catch (ArgumentException)
        {
            Assert.Pass("ArgumentException was thrown as expected for searching neither proteins nor rna mods");
        }
    }

    [Test]
    public static void GetMods_NullOnNonExistingConvention()
    {
        var convention = (ModificationNamingConvention)999; // Invalid convention
        var mod = Mods.GetModification("Oxidation on M", convention);
        Assert.That(mod, Is.Null, "Expected null to be returned for a non-existent convention.");
    }

    [Test]
    public static void GetModifications_Protein()
    {
        var mods = Mods.GetModifications(p => p.IdWithMotif.EndsWith("on M"), true, false);

        foreach (var mod in mods)
        {
            Assert.That(mod.IdWithMotif.EndsWith("on M"), $"Expected modification ID to end with 'on M', but got '{mod.IdWithMotif}'");
            Assert.That(mod.Target.ToString(), Is.EqualTo("M"), $"Expected modification target to be 'M', but got '{mod.Target}'");
        }   
    }

    [Test]
    public static void GetModifications_Rna()
    {
        var mods = Mods.GetModifications(p => p.IdWithMotif.EndsWith("on G"), false, true);
        foreach (var mod in mods)
        {
            Assert.That(mod.IdWithMotif.EndsWith("on G"), $"Expected modification ID to end with 'on G', but got '{mod.IdWithMotif}'");
            Assert.That(mod.Target.ToString(), Is.EqualTo("G"), $"Expected modification target to be 'G', but got '{mod.Target}'");
        }
    }

    [Test]
    public static void GetModifications_Both()
    {
        var mods = Mods.GetModifications(p => p.IdWithMotif.EndsWith("on T")).ToList();
        foreach (var mod in mods)
        {
            Assert.That(mod.IdWithMotif.EndsWith("on T"), $"Expected modification ID to end with 'on T', but got '{mod.IdWithMotif}'");
            Assert.That(mod.Target.ToString(), Is.EqualTo("T"), $"Expected modification target to be 'T', but got '{mod.Target}'");
        }


        bool hasProteinMod = Mods.AllKnownProteinModsDictionary.Values.Any(m => mods.Contains(m));
        bool hasRnaMod = Mods.AllKnownRnaModsDictionary.Values.Any(m => mods.Contains(m));

        Assert.That(hasProteinMod, Is.True, "Expected at least one protein modification.");
        Assert.That(hasRnaMod, Is.True, "Expected at least one RNA modification.");
    }

    [Test]
    public static void AddOrUpdate_Protein()
    {
        Mods.AddOrUpdateModification(CustomDummyMod, false);
        var retrievedMod = Mods.GetModification("TestMod on M", true, false);
        Assert.That(retrievedMod, Is.EqualTo(CustomDummyMod));

        var newMod = new Modification(CustomDummyMod.OriginalId, CustomDummyMod.Accession, CustomDummyMod.ModificationType, CustomDummyMod.FeatureType, CustomDummyMod.Target, CustomDummyMod.LocationRestriction, ChemicalFormula.ParseFormula("C3H7P6"));

        Mods.AddOrUpdateModification(newMod, false);
        var retrievedUpdatedMod = Mods.GetModification("TestMod on M", true, false);
        Assert.That(retrievedUpdatedMod, Is.EqualTo(newMod));
    }

    [Test]
    public static void AddOrUpdate_RNA()
    {
        Mods.AddOrUpdateModification(CustomDummyMod, true);
        var retrievedMod2 = Mods.GetModification("TestMod on M", false, true);
        Assert.That(retrievedMod2, Is.EqualTo(CustomDummyMod));

        var newMod = new Modification(CustomDummyMod.OriginalId, CustomDummyMod.Accession, CustomDummyMod.ModificationType, CustomDummyMod.FeatureType, CustomDummyMod.Target, CustomDummyMod.LocationRestriction, ChemicalFormula.ParseFormula("C3H7P6"));

        Mods.AddOrUpdateModification(newMod, true);
        var retrievedUpdatedMod = Mods.GetModification("TestMod on M", false, true);
        Assert.That(retrievedUpdatedMod, Is.EqualTo(newMod));
    }
}
