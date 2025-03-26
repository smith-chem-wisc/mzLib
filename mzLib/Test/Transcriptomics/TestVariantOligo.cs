using NUnit.Framework;
using Omics.BioPolymer;
using Proteomics.ProteolyticDigestion;
using Proteomics;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using Omics.Modifications;
using Transcriptomics;
using Transcriptomics.Digestion;
using UsefulProteomicsDatabases;
using UsefulProteomicsDatabases.Transcriptomics;

namespace Test.Transcriptomics;

[TestFixture]
[ExcludeFromCodeCoverage]
public class TestVariantOligo
{
    static List<Modification> AllKnownMods;

    [OneTimeSetUp]
    public static void SetUpModifications()
    {
        var modPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "RnaMods.txt");
        AllKnownMods = PtmListLoader.ReadModsFromFile(modPath, out _).ToList();
    }

    [Test]
    public static void VariantRna()
    {
        RNA p = new RNA("CAAA","accession");
        RNA v = new RNA("CAUA", p, new[] { new SequenceVariation(3, "A", "U", "desc", null) }, null, null, null);
        Assert.That(v.NonVariant, Is.EqualTo(p));
    }

    [Test]
    public void VariantXml()
    {
        string file = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "SeqVar.xml");
        List<RNA> variantProteins = RnaDbLoader.LoadRnaXML(file, true, DecoyType.None, false, AllKnownMods, [], out _);

        Assert.That(variantProteins.First().NonVariant.SequenceVariations.Count(), Is.EqualTo(5));
        Assert.That(variantProteins.Count, Is.EqualTo(1)); // there is only one unique amino acid change
        Assert.That(variantProteins.First().NonVariant.BaseSequence, Is.Not.EqualTo(variantProteins.First().BaseSequence));
        Assert.That(variantProteins.First().NonVariant.BaseSequence[116], Is.EqualTo('C'));
        Assert.That(variantProteins.First().BaseSequence[116], Is.EqualTo('G'));
        Assert.That(variantProteins.First().NonVariant.Name, Is.Not.EqualTo(variantProteins.First().Name));
        Assert.That(variantProteins.First().NonVariant.FullName, Is.Not.EqualTo(variantProteins.First().FullName));
        Assert.That(variantProteins.First().NonVariant.Accession, Is.Not.EqualTo(variantProteins.First().Accession));

        List<OligoWithSetMods> oligos = variantProteins.SelectMany(vp => vp.Digest(new RnaDigestionParams(), null, null)).ToList();
    }

    [Test]
    [TestCase("oblm1.xml", 1, 1)] // mod on first residue
    [TestCase("oblm2.xml", 3, 3)] // mod on central residue
    [TestCase("oblm3.xml", 6, 6)] // mod on last residue
    public static void LoadSeqVarModifications(string databaseName, int modIdx, int reversedModIdx)
    {
        string dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", databaseName);
        var rna = RnaDbLoader.LoadRnaXML(dbPath, true, DecoyType.Reverse, false, AllKnownMods, [], out var unknownModifications);
        var target = rna[0];
        Assert.That(target.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(1));
        Assert.That(target.OneBasedPossibleLocalizedModifications.Single().Key, Is.EqualTo(modIdx));
        Assert.That(target.AppliedSequenceVariations.Count(), Is.EqualTo(1));
        Assert.That(target.AppliedSequenceVariations.Single().OneBasedBeginPosition, Is.EqualTo(modIdx));
        Assert.That(target.SequenceVariations.Count(), Is.EqualTo(1));
        Assert.That(target.SequenceVariations.Single().OneBasedBeginPosition, Is.EqualTo(modIdx));
        Assert.That(target.SequenceVariations.Single().OneBasedModifications.Count, Is.EqualTo(1));
        Assert.That(target.SequenceVariations.Single().OneBasedModifications.Single().Key, Is.EqualTo(modIdx)); //PEP[mod]TID, MEP[mod]TID
        var decoy = rna[1];
        Assert.That(decoy.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(1));
        Assert.That(decoy.OneBasedPossibleLocalizedModifications.Single().Key, Is.EqualTo(reversedModIdx)); //DITP[mod]EP, MDITP[mod]E
        Assert.That(decoy.AppliedSequenceVariations.Count(), Is.EqualTo(1));
        Assert.That(decoy.AppliedSequenceVariations.Single().OneBasedBeginPosition, Is.EqualTo(reversedModIdx));
        Assert.That(decoy.SequenceVariations.Count(), Is.EqualTo(1));
        Assert.That(decoy.SequenceVariations.Single().OneBasedBeginPosition, Is.EqualTo(reversedModIdx));
        Assert.That(decoy.SequenceVariations.Single().OneBasedModifications.Count, Is.EqualTo(1));
        Assert.That(decoy.SequenceVariations.Single().OneBasedModifications.Single().Key, Is.EqualTo(reversedModIdx));

        string rewriteDbName = $"{Path.GetFileNameWithoutExtension(databaseName)}rewrite.xml";
        ProteinDbWriter.WriteXmlDatabase(null, rna.Where(p => !p.IsDecoy).ToList(), Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", rewriteDbName));
        rna = RnaDbLoader.LoadRnaXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", rewriteDbName), true,
            DecoyType.Reverse, false, AllKnownMods, [], out unknownModifications);
        target = rna[0];
        Assert.That(target.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(1));
        Assert.That(target.OneBasedPossibleLocalizedModifications.Single().Key, Is.EqualTo(modIdx));
        Assert.That(target.AppliedSequenceVariations.Count(), Is.EqualTo(1));
        Assert.That(target.AppliedSequenceVariations.Single().OneBasedBeginPosition, Is.EqualTo(modIdx));
        Assert.That(target.SequenceVariations.Count(), Is.EqualTo(1));
        Assert.That(target.SequenceVariations.Single().OneBasedBeginPosition, Is.EqualTo(modIdx));
        Assert.That(target.SequenceVariations.Single().OneBasedModifications.Count, Is.EqualTo(1));
        Assert.That(target.SequenceVariations.Single().OneBasedModifications.Single().Key, Is.EqualTo(modIdx));
        decoy = rna[1];
        Assert.That(decoy.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(1));
        Assert.That(decoy.OneBasedPossibleLocalizedModifications.Single().Key, Is.EqualTo(reversedModIdx));
        Assert.That(decoy.AppliedSequenceVariations.Count(), Is.EqualTo(1));
        Assert.That(decoy.AppliedSequenceVariations.Single().OneBasedBeginPosition, Is.EqualTo(reversedModIdx));
        Assert.That(decoy.SequenceVariations.Count(), Is.EqualTo(1));
        Assert.That(decoy.SequenceVariations.Single().OneBasedBeginPosition, Is.EqualTo(reversedModIdx));
        Assert.That(decoy.SequenceVariations.Single().OneBasedModifications.Count, Is.EqualTo(1));
        Assert.That(decoy.SequenceVariations.Single().OneBasedModifications.Single().Key, Is.EqualTo(reversedModIdx));
    }
}