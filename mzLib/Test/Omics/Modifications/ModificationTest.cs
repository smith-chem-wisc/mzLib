using System;
using Chemistry;
using NUnit.Framework;
using Omics.Modifications;
using Omics.Modifications.IO;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Linq;

namespace Test.Omics.Modifications;

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
            new (Mods.AllKnownProteinModsDictionary["Phosphorylation on T"], ModificationNamingConvention.MetaMorpheus_Protein, true),
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
            Assert.That(mod.DatabaseReference.Count, Is.EqualTo(1));
            Assert.That(mod.DatabaseReference.First().Key, Is.EqualTo("Unimod"));
        }
    }

    [Test]
    public static void MetaMorpheusModsAreNotUniprotOrUnimod()
    {
        var mods = Mods.MetaMorpheusProteinModifications;
        foreach (var mod in mods)
        {
            Assert.That(mod.ModificationType, Is.Not.EqualTo("UniProt"));
            Assert.That(mod.ModificationType, Is.Not.EqualTo("Unimod"));
        }
    }

    /// <summary>
    /// Lactylation adds lactic acid (C3H6O3) less one water, so C3H4O2. The Mods.txt entry once read
    /// C3H3O2, a hydrogen short of Unimod 2114 (72.021129), which shares its name.
    /// </summary>
    [Test]
    public static void LactylationAddsLacticAcidLessWater()
    {
        var lactylation = Mods.MetaMorpheusProteinModifications.Single(m => m.IdWithMotif == "Lactylation on K");

        Assert.That(lactylation.ChemicalFormula, Is.EqualTo(ChemicalFormula.ParseFormula("C3H4O2")));
        Assert.That(lactylation.MonoisotopicMass, Is.EqualTo(72.021129).Within(1e-5));
        Assert.That(lactylation.DatabaseReference["Unimod"], Does.Contain("2114"));
    }

    [Test]
    public static void GetModification_Nonsense()
    {
        var mod = Mods.GetModification("This modification does not exist");
        Assert.That(mod, Is.Null, "Expected null to be returned for a non-existent modification.");
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
        if (testCase.Convention.ToString().Contains("MetaMor", StringComparison.OrdinalIgnoreCase))
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

    /// <summary>
    /// The "DR   Unimod; N." lines in the shipped mod files are hand-written, not derived from
    /// unimod.xml. Where an entry declares its own chemical formula, that formula must equal the
    /// composition of the record it cites. Entries with no formula are skipped: the isobaric
    /// labels in TMT.txt give an MM averaged across label channels on purpose. A cited record that
    /// is absent from the shipped unimod.xml is also skipped, since that is a stale-XML problem
    /// rather than a wrong accession.
    /// </summary>
    [Test]
    public static void ShippedUnimodAccessionsMatchTheCitedRecordsComposition()
    {
        var unimodByAccession = Mods.UnimodModifications
            .Where(m => m.ModificationType == "Unimod"
                        && m.DatabaseReference != null
                        && m.DatabaseReference.ContainsKey("Unimod"))
            .GroupBy(m => m.DatabaseReference["Unimod"].First())
            .ToDictionary(g => g.Key, g => g.First());

        var shippedMods = Mods.MetaMorpheusProteinModifications
            .Concat(Mods.MetaMorpheusRnaModifications)
            .Concat(Mods.IsobaricLabelModifications)
            .Concat(ProteaseDictionary.LoadEmbeddedProteaseMods());

        var mismatches = new List<string>();

        foreach (var mod in shippedMods)
        {
            if (mod.ChemicalFormula == null
                || mod.DatabaseReference == null
                || !mod.DatabaseReference.TryGetValue("Unimod", out var accessions))
            {
                continue;
            }

            string accession = accessions.First();
            if (!unimodByAccession.TryGetValue(accession, out var cited))
                continue;

            if (!mod.ChemicalFormula.Equals(cited.ChemicalFormula))
            {
                mismatches.Add($"'{mod.IdWithMotif}' ({mod.ChemicalFormula.Formula}, {mod.MonoisotopicMass:F6}) cites " +
                               $"UNIMOD:{accession} '{cited.OriginalId}' ({cited.ChemicalFormula.Formula}, {cited.MonoisotopicMass:F6})");
            }
        }

        Assert.That(mismatches, Is.Empty, string.Join(Environment.NewLine, mismatches));
    }

    /// <summary>
    /// A loaded modification's monoisotopic mass must be the mass of its own formula. UniProt's
    /// ptmlist release 2026_01 had the MM and MA lines swapped on 29 complex N-glycans, so each
    /// carried its average mass (about 1 Da high), and the 2014 PSI-MOD charge list lacked
    /// N,N,N-trimethylglycine, so its formula kept one hydrogen too many. Both break this equality.
    /// </summary>
    [Test]
    public static void EveryLoadedModificationsMassIsItsFormulasMass()
    {
        var mismatches = Mods.UniprotModifications
            .Concat(Mods.MetaMorpheusProteinModifications)
            .Concat(Mods.IsobaricLabelModifications)
            .Concat(Mods.UnimodModifications)
            .Where(m => m.ChemicalFormula != null && m.MonoisotopicMass != null
                        && Math.Abs(m.MonoisotopicMass.Value - m.ChemicalFormula.MonoisotopicMass) > 0.001)
            .Select(m => $"'{m.IdWithMotif}' ({m.ModificationType}): mass {m.MonoisotopicMass:F5}, " +
                         $"{m.ChemicalFormula.Formula} {m.ChemicalFormula.MonoisotopicMass:F5}")
            .ToList();

        Assert.That(mismatches, Is.Empty, string.Join(Environment.NewLine, mismatches));
    }

    /// <summary>
    /// The embedded formal-charge table is generated from the current PSI-MOD.obo. It must keep every
    /// charge the 2014 PSI-MOD.obo.xml gave (the test fixture copy), with the same sign and size, and it
    /// adds N,N,N-trimethylglycine (MOD:01982), which UniProt cites and the 2014 file did not charge.
    /// </summary>
    [Test]
    public static void EmbeddedFormalChargesKeepEveryChargeOfThePsiModXml()
    {
        var assembly = typeof(Mods).Assembly;
        using var stream = assembly.GetManifestResourceStream($"{assembly.GetName().Name}.Resources.PsiModFormalCharges.tsv");
        using var reader = new System.IO.StreamReader(stream!);
        var embedded = ModificationLoader.ReadFormalChargesDictionary(reader);
        var fromXml = ModificationLoader.GetFormalChargesDictionary(ModificationLoader.LoadPsiMod(TestOntologies.PsiModXml));

        foreach (var (accession, charge) in fromXml)
        {
            Assert.That(embedded.TryGetValue(accession, out int embeddedCharge), $"{accession} is missing from the embedded table");
            Assert.That(embeddedCharge, Is.EqualTo(charge), accession);
        }
        Assert.That(embedded["PSI-MOD; MOD:01982"], Is.EqualTo(1));
        Assert.That(embedded.Values.Any(c => c < 0), "negative charges must keep their sign");
    }

    /// <summary>
    /// UniProt writes a charged modification's formula for the charged species. N,N,N-trimethylglycine
    /// must lose one hydrogen on load, exactly as N6,N6,N6-trimethyllysine always has, and get the same mass.
    /// </summary>
    [Test]
    public static void TrimethylglycineIsChargeCorrectedLikeTrimethyllysine()
    {
        var glycine = Mods.UniprotModifications.Single(m => m.IdWithMotif == "N,N,N-trimethylglycine on G");
        var lysine = Mods.UniprotModifications.Single(m => m.IdWithMotif == "N6,N6,N6-trimethyllysine on K");

        Assert.That(glycine.ChemicalFormula, Is.EqualTo(ChemicalFormula.ParseFormula("C3H6")));
        Assert.That(glycine.ChemicalFormula, Is.EqualTo(lysine.ChemicalFormula));
        Assert.That(glycine.MonoisotopicMass, Is.EqualTo(lysine.MonoisotopicMass).Within(1e-9));
        Assert.That(glycine.MonoisotopicMass, Is.EqualTo(glycine.ChemicalFormula.MonoisotopicMass).Within(1e-9),
            "UniProt writes this entry's MM as the neutral formula mass, not the cation's; the mass must come from the corrected formula");
    }

    [Test]
    public static void ReadFormalChargesDictionarySkipsCommentsKeepsSignsAndRefusesBadLines()
    {
        var charges = ModificationLoader.ReadFormalChargesDictionary(
            new System.IO.StringReader("# header\n\nMOD:00083\t1\nMOD:00147\t-3\n"));

        Assert.That(charges, Has.Count.EqualTo(2));
        Assert.That(charges["PSI-MOD; MOD:00083"], Is.EqualTo(1));
        Assert.That(charges["PSI-MOD; MOD:00147"], Is.EqualTo(-3));

        Assert.Throws<MzLibUtil.MzLibException>(() =>
            ModificationLoader.ReadFormalChargesDictionary(new System.IO.StringReader("MOD:00083 1\n")));
        Assert.Throws<MzLibUtil.MzLibException>(() =>
            ModificationLoader.ReadFormalChargesDictionary(new System.IO.StringReader("MOD:00083\t1+\n")));
    }

    /// <summary>
    /// Every curated UniProt-to-Unimod cross-reference must name an entry mzLib loads, and a Unimod record
    /// of the same composition (against the formula mzLib holds, after the formal-charge correction) that
    /// lists the entry's residue among its sites. Sameness of chemistry beyond that is a review decision,
    /// recorded in the table's header.
    /// </summary>
    [Test]
    public static void UniprotUnimodCrossReferencesCiteTheSameCompositionOnTheSameResidue()
    {
        var unimodRecords = Mods.UnimodModifications
            .Where(m => m.ModificationType == "Unimod" && m.DatabaseReference != null && m.DatabaseReference.ContainsKey("Unimod"))
            .ToLookup(m => m.DatabaseReference["Unimod"].First());

        var problems = new List<string>();
        foreach (var (accession, unimodId) in ModificationLoader.UniprotUnimodCrossReferences)
        {
            var entries = Mods.UniprotModifications.Where(m => m.Accession == accession).ToList();
            if (entries.Count == 0)
            {
                problems.Add($"{accession}: no loaded UniProt modification has this accession");
                continue;
            }
            var cited = unimodRecords[unimodId].ToList();
            if (cited.Count == 0)
            {
                problems.Add($"{accession}: UNIMOD:{unimodId} is not in the embedded unimod.xml");
                continue;
            }
            foreach (var entry in entries)
            {
                if (!cited.Any(u => u.ChemicalFormula != null && u.ChemicalFormula.Equals(entry.ChemicalFormula)))
                    problems.Add($"'{entry.IdWithMotif}' ({entry.ChemicalFormula?.Formula}) cites UNIMOD:{unimodId} " +
                                 $"({cited.First().ChemicalFormula?.Formula}), a different composition");
                if (!cited.Any(u => u.Target?.ToString() == entry.Target?.ToString()))
                    problems.Add($"'{entry.IdWithMotif}' cites UNIMOD:{unimodId}, which does not list {entry.Target} as a site");
                if (!entry.DatabaseReference["Unimod"].Contains(unimodId))
                    problems.Add($"'{entry.IdWithMotif}' did not receive UNIMOD:{unimodId}");
            }
        }

        Assert.That(ModificationLoader.UniprotUnimodCrossReferences, Is.Not.Empty);
        Assert.That(problems, Is.Empty, string.Join(Environment.NewLine, problems));
    }

    /// <summary>
    /// The curated table only fills a gap. An entry that carries its own Unimod reference keeps it, and an
    /// entry with no database references at all gets one.
    /// </summary>
    [Test]
    public static void UniprotUnimodCrossReferenceIsAddedOnlyWhenUniProtHasNone()
    {
        static Modification Read(string dr) => ModificationLoader.ReadModsFromFile(
            new System.IO.StreamReader(new System.IO.MemoryStream(System.Text.Encoding.UTF8.GetBytes(
                "ID   N6-succinyllysine\r\nAC   PTM-0438\r\nFT   MOD_RES\r\nTG   Lysine.\r\nPP   Anywhere.\r\n" +
                "CF   C4 H4 O3\r\nMM   100.016044\r\n" + dr + "//"))),
            new Dictionary<string, int>(), out _).Single();

        Assert.That(Read("").DatabaseReference["Unimod"], Is.EqualTo(new[] { "64" }));
        Assert.That(Read("DR   PSI-MOD; MOD:01819.\r\n").DatabaseReference["Unimod"], Is.EqualTo(new[] { "64" }));
        Assert.That(Read("DR   Unimod; 999.\r\n").DatabaseReference["Unimod"], Is.EqualTo(new[] { "999" }),
            "UniProt's own Unimod reference wins over the curated one");
    }

    /// <summary>
    /// The point of the table: a MetaMorpheus full-sequence name for a UniProt modification resolves to its
    /// Unimod accession.
    /// </summary>
    [Test]
    public static void SuccinyllysineResolvesToItsUnimodAccession()
    {
        var resolved = global::Omics.SequenceConversion.GlobalModificationLookup.Instance.TryResolve("UniProt:N6-succinyllysine on K");

        Assert.That(resolved, Is.Not.Null);
        Assert.That(resolved!.Value.UnimodId, Is.EqualTo(64));
    }

    [Test]
    public static void ReadUnimodCrossReferencesSkipsCommentsAndRefusesBadLines()
    {
        var crossReferences = ModificationLoader.ReadUnimodCrossReferences(
            new System.IO.StringReader("# header\n\nPTM-0438\t64\tN6-succinyllysine\tSuccinyl\n"));

        Assert.That(crossReferences, Has.Count.EqualTo(1));
        Assert.That(crossReferences["PTM-0438"], Is.EqualTo("64"));
        Assert.Throws<MzLibUtil.MzLibException>(() =>
            ModificationLoader.ReadUnimodCrossReferences(new System.IO.StringReader("PTM-0438 64\n")));
        Assert.Throws<MzLibUtil.MzLibException>(() =>
            ModificationLoader.ReadUnimodCrossReferences(new System.IO.StringReader("PTM-0438\tUNIMOD:64\n")));
    }

    [Test]
    public static void ReadFormalChargesDictionaryReadsAFile()
    {
        var path = System.IO.Path.Combine(TestContext.CurrentContext.WorkDirectory, $"{nameof(ReadFormalChargesDictionaryReadsAFile)}.tsv");
        System.IO.File.WriteAllText(path, "# accession\tcharge\nMOD:00083\t1\nMOD:00147\t-3\n");
        try
        {
            var charges = ModificationLoader.ReadFormalChargesDictionary(path);

            Assert.That(charges, Has.Count.EqualTo(2));
            Assert.That(charges["PSI-MOD; MOD:00083"], Is.EqualTo(1));
            Assert.That(charges["PSI-MOD; MOD:00147"], Is.EqualTo(-3));
        }
        finally
        {
            System.IO.File.Delete(path);
        }
    }
}
