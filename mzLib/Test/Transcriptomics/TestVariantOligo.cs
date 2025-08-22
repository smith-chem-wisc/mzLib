using System;
using NUnit.Framework;
using Omics.BioPolymer;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using Omics.Modifications;
using Transcriptomics;
using Transcriptomics.Digestion;
using UsefulProteomicsDatabases;
using UsefulProteomicsDatabases.Transcriptomics;
using Omics;
using Proteomics;

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
        Assert.That(v.ConsensusVariant, Is.EqualTo(p));
    }

    [Test]
    public void VariantXml()
    {
        string file = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "SeqVar.xml");
        List<RNA> variantProteins = RnaDbLoader.LoadRnaXML(file, true, DecoyType.None, false, AllKnownMods, [], out _);

        Assert.That(variantProteins.First().ConsensusVariant.SequenceVariations.Count(), Is.EqualTo(5));
        Assert.That(variantProteins.Count, Is.EqualTo(1)); // there is only one unique amino acid change
        Assert.That(variantProteins.First().ConsensusVariant.BaseSequence, Is.Not.EqualTo(variantProteins.First().BaseSequence));
        Assert.That(variantProteins.First().ConsensusVariant.BaseSequence[116], Is.EqualTo('C'));
        Assert.That(variantProteins.First().BaseSequence[116], Is.EqualTo('G'));
        Assert.That(variantProteins.First().ConsensusVariant.Name, Is.Not.EqualTo(variantProteins.First().Name));
        Assert.That(variantProteins.First().ConsensusVariant.FullName, Is.Not.EqualTo(variantProteins.First().FullName));
        Assert.That(variantProteins.First().ConsensusVariant.Accession, Is.Not.EqualTo(variantProteins.First().Accession));

        List<OligoWithSetMods> oligos = variantProteins.SelectMany(vp => vp.Digest(new RnaDigestionParams(), null, null)).ToList();
    }

    [Test]
    [TestCase("oblm1.xml", 1, 6)] // mod on first residue
    [TestCase("oblm2.xml", 3, 4)] // mod on central residue
    [TestCase("oblm3.xml", 6, 1)] // mod on last residue
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
        ProteinDbWriter.WriteXmlDatabase(rna.Where(p => !p.IsDecoy).ToList(), Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", rewriteDbName));
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

    [TestCase("ranges1.xml", 1, 2, 5, 6)] // trunc excludes natural 3'
    [TestCase("ranges2.xml", 2, 1, 6, 5)] // trunc includes natural 3'
    public static void ReverseDecoyProteolysisProducts(string databaseName, int beginIdx, int reversedBeginIdx, int endIdx, int reversedEndIdx)
    {
        string dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", databaseName);
        var rna = RnaDbLoader.LoadRnaXML(dbPath, true, DecoyType.Reverse, false, AllKnownMods, [], out var unknownModifications);
        var target = rna[0];
        Assert.That(target.TruncationProducts.Count(), Is.EqualTo(1));
        Assert.That(target.TruncationProducts.Single().OneBasedBeginPosition, Is.EqualTo(beginIdx)); //[start]GUACU[end]G, G[start]UACUG[end]
        Assert.That(target.TruncationProducts.Single().OneBasedEndPosition, Is.EqualTo(endIdx));
        var decoy = rna[1];
        Assert.That(decoy.TruncationProducts.Count(), Is.EqualTo(1));
        Assert.That(decoy.TruncationProducts.Single().OneBasedBeginPosition, Is.EqualTo(reversedBeginIdx)); //DI[start]TPEP[end], M[start]DITP[end]E
        Assert.That(decoy.TruncationProducts.Single().OneBasedEndPosition, Is.EqualTo(reversedEndIdx));

        string rewriteDbName = $"{Path.GetFileNameWithoutExtension(databaseName)}rewrite.xml";
        ProteinDbWriter.WriteXmlDatabase(rna.Where(p => !p.IsDecoy).ToList(), Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", rewriteDbName));
        rna = RnaDbLoader.LoadRnaXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", rewriteDbName), true,
            DecoyType.Reverse, false, AllKnownMods, [], out unknownModifications);
        target = rna[0];
        Assert.That(target.TruncationProducts.Count(), Is.EqualTo(1));
        Assert.That(target.TruncationProducts.Single().OneBasedBeginPosition, Is.EqualTo(beginIdx));
        Assert.That(target.TruncationProducts.Single().OneBasedEndPosition, Is.EqualTo(endIdx));
        decoy = rna[1];
        Assert.That(decoy.TruncationProducts.Count(), Is.EqualTo(1));
        Assert.That(decoy.TruncationProducts.Single().OneBasedBeginPosition, Is.EqualTo(reversedBeginIdx));
        Assert.That(decoy.TruncationProducts.Single().OneBasedEndPosition, Is.EqualTo(reversedEndIdx));
    }

    [Test]
    [TestCase("HomozygousHLA.xml", 1, 18)]
    [TestCase("HomozygousHLA.xml", 10, 17)]
    public static void HomozygousVariantsAtVariedDepths(string filename, int minVariantDepth, int appliedCount)
    {
        string dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", filename);
        var rna = RnaDbLoader.LoadRnaXML(dbPath, true, DecoyType.None, false, AllKnownMods, [], out var unknownModifications, minAlleleDepth: minVariantDepth);
        Assert.That(rna.Count, Is.EqualTo(1));
        Assert.That(rna[0].SequenceVariations.Count(), Is.EqualTo(18)); // some redundant
        Assert.That(rna[0].SequenceVariations.Select(v => v.SimpleString()).Distinct().Count(), Is.EqualTo(18)); // unique changes
        Assert.That(rna[0].AppliedSequenceVariations.Count(), Is.EqualTo(appliedCount)); // some redundant
        Assert.That(rna[0].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count(), Is.EqualTo(appliedCount)); // unique changes
        Assert.That(rna[0].GetVariantBioPolymers().Count, Is.EqualTo(1));
        var variantProteins = rna[0].GetVariantBioPolymers();
        List<OligoWithSetMods> peptides = rna.SelectMany(vp => vp.Digest(new RnaDigestionParams(), null, null)).ToList();
    }

    [Test]
    public static void AppliedVariants()
    {
        ModificationMotif.TryGetMotif("C", out ModificationMotif motifP);
        Modification mp = new Modification("mod", null, "type", null, motifP, "Anywhere.", null, 42.01, new Dictionary<string, IList<string>>(), null, null, null, null, null);

        List<RNA> proteinsWithSeqVars = new List<RNA>
            {
                new RNA("GUACUGUA", "protein1", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "C", "U", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new RNA("GUACUGUA", "protein2", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 5, "CU", "AU", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new RNA("GUACUGUA", "protein3", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "C", "CCC", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new RNA("GUACCCUGUA", "protein4", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 6, "CCC", "C", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new RNA("GUACUGUA", "protein5", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "C", "CCC", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", new Dictionary<int, List<Modification>> {{ 5, new[] { mp }.ToList() } }) }),
             };
        var proteinsWithAppliedVariants = proteinsWithSeqVars.SelectMany(p => p.GetVariantBioPolymers()).ToList();
        var proteinsWithAppliedVariants2 = proteinsWithSeqVars.SelectMany(p => p.GetVariantBioPolymers()).ToList(); // should be stable
        string xml = Path.Combine(TestContext.CurrentContext.TestDirectory, "AppliedVariants.xml");
        ProteinDbWriter.WriteXmlDatabase(proteinsWithSeqVars, xml);
        var proteinsWithAppliedVariants3 = RnaDbLoader.LoadRnaXML(xml, true, DecoyType.None, false, AllKnownMods, null, out var un);

        var listArray = new List<RNA>[] { proteinsWithAppliedVariants, proteinsWithAppliedVariants2, proteinsWithAppliedVariants3 };
        for (int dbIdx = 0; dbIdx < listArray.Length; dbIdx++)
        {
            // sequences
            Assert.That(listArray[dbIdx][0].BaseSequence, Is.EqualTo("GUAUUGUA"));
            Assert.That(listArray[dbIdx][1].BaseSequence, Is.EqualTo("GUAAUGUA"));
            Assert.That(listArray[dbIdx][2].BaseSequence, Is.EqualTo("GUACCCUGUA"));
            Assert.That(listArray[dbIdx][3].BaseSequence, Is.EqualTo("GUACUGUA"));
            Assert.That(listArray[dbIdx][4].BaseSequence, Is.EqualTo("GUACCCUGUA"));
            Assert.That(listArray[dbIdx][4].OneBasedPossibleLocalizedModifications.Single().Key, Is.EqualTo(5));

            // SAV
            Assert.That(listArray[dbIdx][0].AppliedSequenceVariations.Single().OneBasedBeginPosition, Is.EqualTo(4));
            Assert.That(listArray[dbIdx][0].AppliedSequenceVariations.Single().OneBasedEndPosition, Is.EqualTo(4));

            // MNV
            Assert.That(listArray[dbIdx][1].AppliedSequenceVariations.Single().OneBasedBeginPosition, Is.EqualTo(4));
            Assert.That(listArray[dbIdx][1].AppliedSequenceVariations.Single().OneBasedEndPosition, Is.EqualTo(5));

            // insertion
            Assert.That(listArray[dbIdx][2].AppliedSequenceVariations.Single().OneBasedBeginPosition, Is.EqualTo(4));
            Assert.That(listArray[dbIdx][2].AppliedSequenceVariations.Single().OneBasedEndPosition, Is.EqualTo(6));

            // deletion
            Assert.That(listArray[dbIdx][3].AppliedSequenceVariations.Single().OneBasedBeginPosition, Is.EqualTo(4));
            Assert.That(listArray[dbIdx][3].AppliedSequenceVariations.Single().OneBasedEndPosition, Is.EqualTo(4));
        }
    }

    [Test]
    public static void AppliedVariants_AsBioPolymer()
    {
        ModificationMotif.TryGetMotif("C", out ModificationMotif motifP);
        Modification mp = new Modification("mod", null, "type", null, motifP, "Anywhere.", null, 42.01, new Dictionary<string, IList<string>>(), null, null, null, null, null);

        List<IBioPolymer> proteinsWithSeqVars = new List<IBioPolymer>
            {
                new RNA("GUACUGUA", "protein1", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "C", "U", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new RNA("GUACUGUA", "protein2", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 5, "CU", "AU", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new RNA("GUACUGUA", "protein3", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "C", "CCC", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new RNA("GUACCCUGUA", "protein4", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 6, "CCC", "C", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new RNA("GUACUGUA", "protein5", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "C", "CCC", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", new Dictionary<int, List<Modification>> {{ 5, new[] { mp }.ToList() } }) }),
             };
        var proteinsWithAppliedVariants = proteinsWithSeqVars.SelectMany(p => p.GetVariantBioPolymers()).ToList();
        var proteinsWithAppliedVariants2 = proteinsWithSeqVars.SelectMany(p => p.GetVariantBioPolymers()).ToList(); // should be stable
        string xml = Path.Combine(TestContext.CurrentContext.TestDirectory, "AppliedVariants.xml");
        ProteinDbWriter.WriteXmlDatabase(proteinsWithSeqVars, xml);
        var proteinsWithAppliedVariants3 = RnaDbLoader.LoadRnaXML(xml, true, DecoyType.None, false, AllKnownMods, null, out var un).Cast<IBioPolymer>().ToList();

        var listArray = new List<IBioPolymer>[] { proteinsWithAppliedVariants, proteinsWithAppliedVariants2, proteinsWithAppliedVariants3 };
        for (int dbIdx = 0; dbIdx < listArray.Length; dbIdx++)
        {
            // sequences
            Assert.That(listArray[dbIdx][0].BaseSequence, Is.EqualTo("GUAUUGUA"));
            Assert.That(listArray[dbIdx][1].BaseSequence, Is.EqualTo("GUAAUGUA"));
            Assert.That(listArray[dbIdx][2].BaseSequence, Is.EqualTo("GUACCCUGUA"));
            Assert.That(listArray[dbIdx][3].BaseSequence, Is.EqualTo("GUACUGUA"));
            Assert.That(listArray[dbIdx][4].BaseSequence, Is.EqualTo("GUACCCUGUA"));
            Assert.That(listArray[dbIdx][4].OneBasedPossibleLocalizedModifications.Single().Key, Is.EqualTo(5));

            // SAV
            Assert.That(listArray[dbIdx][0].AppliedSequenceVariations.Single().OneBasedBeginPosition, Is.EqualTo(4));
            Assert.That(listArray[dbIdx][0].AppliedSequenceVariations.Single().OneBasedEndPosition, Is.EqualTo(4));

            // MNV
            Assert.That(listArray[dbIdx][1].AppliedSequenceVariations.Single().OneBasedBeginPosition, Is.EqualTo(4));
            Assert.That(listArray[dbIdx][1].AppliedSequenceVariations.Single().OneBasedEndPosition, Is.EqualTo(5));

            // insertion
            Assert.That(listArray[dbIdx][2].AppliedSequenceVariations.Single().OneBasedBeginPosition, Is.EqualTo(4));
            Assert.That(listArray[dbIdx][2].AppliedSequenceVariations.Single().OneBasedEndPosition, Is.EqualTo(6));

            // deletion
            Assert.That(listArray[dbIdx][3].AppliedSequenceVariations.Single().OneBasedBeginPosition, Is.EqualTo(4));
            Assert.That(listArray[dbIdx][3].AppliedSequenceVariations.Single().OneBasedEndPosition, Is.EqualTo(4));
        }
    }

    [Test]
    public static void StopGained()
    {

        string dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "StopGained.xml");
        var rna = RnaDbLoader.LoadRnaXML(dbPath, true, DecoyType.None, false, AllKnownMods, [], out var unknownModifications);

        Assert.That(rna.Count, Is.EqualTo(2));
        Assert.That(rna[0].SequenceVariations.Count(), Is.EqualTo(1)); // some redundant
        Assert.That(rna[0].SequenceVariations.Select(v => v.SimpleString()).Distinct().Count(), Is.EqualTo(1)); // unique changes
        Assert.That(rna[0].AppliedSequenceVariations.Count(), Is.EqualTo(0)); // some redundant
        Assert.That(rna[0].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count(), Is.EqualTo(0)); // unique changes
        Assert.That(rna[1].AppliedSequenceVariations.Count(), Is.EqualTo(1)); // some redundant
        Assert.That(rna[1].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count(), Is.EqualTo(1)); // unique changes
        Assert.That(rna[0].Length, Is.EqualTo(191));
        Assert.That(rna[0][161 - 1], Is.EqualTo('G'));
        Assert.That(rna[1].Length, Is.EqualTo(161 - 1));
        Assert.That(rna[0].Length, Is.Not.EqualTo(rna[1].Length));

        rna = RnaDbLoader.LoadRnaXML(dbPath, true, DecoyType.None, false, AllKnownMods, [], out unknownModifications, minAlleleDepth: 400);

        Assert.That(rna.Count, Is.EqualTo(1));
        Assert.That(rna[0].AppliedSequenceVariations.Count(), Is.EqualTo(1)); // some redundant
        Assert.That(rna[0].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count(), Is.EqualTo(1)); // unique changes
        Assert.That(rna[0].Length, Is.EqualTo(161 - 1));
    }

    [Test]
    public static void MultipleAlternateAlleles()
    {
        string dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "MultipleAlternateAlleles.xml");
        var rna = RnaDbLoader.LoadRnaXML(dbPath, true, DecoyType.None, false, AllKnownMods, [], out var unknownModifications);
        Assert.That(rna.Count, Is.EqualTo(2));
        Assert.That(rna[0].SequenceVariations.Count(), Is.EqualTo(2)); // some redundant
        Assert.That(rna[0].SequenceVariations.Select(v => v.SimpleString()).Distinct().Count(), Is.EqualTo(2)); // unique changes

        Assert.That(rna[0].SequenceVariations.All(v => v.OneBasedBeginPosition == 63), Is.True); // there are two alternate alleles (1 and 2), but only 2 is in the genotype, so only that's applied
        Assert.That(rna[1].AppliedSequenceVariations.Count(), Is.EqualTo(1)); // some redundant
        Assert.That(rna[1].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count(), Is.EqualTo(1)); // unique changes
        Assert.That(rna[0].Length, Is.EqualTo(72));
        Assert.That(rna[1].Length, Is.EqualTo(72));
        Assert.That(rna[0][63 - 1], Is.EqualTo('G'));
        Assert.That(rna[1][63 - 1], Is.EqualTo('A'));

        rna = RnaDbLoader.LoadRnaXML(dbPath, true, DecoyType.None, false, AllKnownMods, [], out unknownModifications, minAlleleDepth: 10);

        Assert.That(rna.Count, Is.EqualTo(1));
        Assert.That(rna[0].AppliedSequenceVariations.Count(), Is.EqualTo(0)); // some redundant
        Assert.That(rna[0].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count(), Is.EqualTo(0)); // unique changes
        Assert.That(rna[0][63 - 1], Is.EqualTo('G')); // reference only
    }

    [Test]
    public static void CrashOnCreateVariantFromProtein()
    {
        var rnas = RnaDbLoader.LoadRnaXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "HomozygousHLA.xml"), true,
            DecoyType.None, false, null, null, out var unknownModifications);

        var protein = new Protein("PEPTIDE", "accession");
        NUnit.Framework.Assert.Throws<ArgumentException>(() =>
        {
            rnas[0].CreateVariant(rnas[0].BaseSequence, protein, [], [], new Dictionary<int, List<Modification>>(), "");
        });
    }

    [Test]
    public void IndelDecoyVariants()
    {
        string dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "DecoyVariants.xml");
        var variantRna = RnaDbLoader.LoadRnaXML(dbPath, true, DecoyType.Reverse, false, AllKnownMods, [], out var unknownModifications);

        Assert.That(variantRna.Count, Is.EqualTo(4));
        var homoTarget = variantRna[0];
        Assert.That(homoTarget.IsDecoy, Is.False);
        Assert.That(homoTarget.AppliedSequenceVariations.Count, Is.EqualTo(3));
        Assert.That(homoTarget.AppliedSequenceVariations[0].OriginalSequence, Is.EqualTo("C"));
        Assert.That(homoTarget.AppliedSequenceVariations[0].OneBasedBeginPosition, Is.EqualTo(1222));
        Assert.That(homoTarget.AppliedSequenceVariations[0].VariantSequence, Is.EqualTo("A"));
        Assert.That(homoTarget.AppliedSequenceVariations[1].OriginalSequence, Is.EqualTo("C"));
        Assert.That(homoTarget.AppliedSequenceVariations[1].OneBasedBeginPosition, Is.EqualTo(1488));
        Assert.That(homoTarget.AppliedSequenceVariations[1].VariantSequence, Is.EqualTo("G"));
        Assert.That(homoTarget.AppliedSequenceVariations[2].OriginalSequence, Is.EqualTo("C"));
        Assert.That(homoTarget.AppliedSequenceVariations[2].OneBasedBeginPosition, Is.EqualTo(1646));
        Assert.That(homoTarget.AppliedSequenceVariations[2].VariantSequence, Is.EqualTo("A"));

        var plusOneHeteroTarget = variantRna[1];
        Assert.That(plusOneHeteroTarget.IsDecoy, Is.False);
        Assert.That(plusOneHeteroTarget.AppliedSequenceVariations.Count, Is.EqualTo(4)); 
        Assert.That(plusOneHeteroTarget.AppliedSequenceVariations[0].OriginalSequence, Is.EqualTo("A"));
        Assert.That(plusOneHeteroTarget.AppliedSequenceVariations[0].OneBasedBeginPosition, Is.EqualTo(409));
        Assert.That(plusOneHeteroTarget.AppliedSequenceVariations[0].VariantSequence, Is.EqualTo("U"));
        Assert.That(plusOneHeteroTarget.AppliedSequenceVariations[1].OriginalSequence, Is.EqualTo("C"));
        Assert.That(plusOneHeteroTarget.AppliedSequenceVariations[1].OneBasedBeginPosition, Is.EqualTo(1222));
        Assert.That(plusOneHeteroTarget.AppliedSequenceVariations[1].VariantSequence, Is.EqualTo("A"));
        Assert.That(plusOneHeteroTarget.AppliedSequenceVariations[2].OriginalSequence, Is.EqualTo("C"));
        Assert.That(plusOneHeteroTarget.AppliedSequenceVariations[2].OneBasedBeginPosition, Is.EqualTo(1488));
        Assert.That(plusOneHeteroTarget.AppliedSequenceVariations[2].VariantSequence, Is.EqualTo("G"));
        Assert.That(plusOneHeteroTarget.AppliedSequenceVariations[3].OriginalSequence, Is.EqualTo("C"));
        Assert.That(plusOneHeteroTarget.AppliedSequenceVariations[3].OneBasedBeginPosition, Is.EqualTo(1646));
        Assert.That(plusOneHeteroTarget.AppliedSequenceVariations[3].VariantSequence, Is.EqualTo("A"));

        var homoDecoy = variantRna[2];
        Assert.That(homoDecoy.IsDecoy, Is.True);
        Assert.That(homoDecoy.AppliedSequenceVariations.Count, Is.EqualTo(3));
        Assert.That(homoDecoy.AppliedSequenceVariations[0].OriginalSequence, Is.EqualTo("C"));
        Assert.That(homoDecoy.AppliedSequenceVariations[0].OneBasedBeginPosition, Is.EqualTo(homoTarget.Length - 1646 + 1));
        Assert.That(homoDecoy.AppliedSequenceVariations[0].VariantSequence, Is.EqualTo("A"));
        Assert.That(homoDecoy.AppliedSequenceVariations[1].OriginalSequence, Is.EqualTo("C"));
        Assert.That(homoDecoy.AppliedSequenceVariations[1].OneBasedBeginPosition, Is.EqualTo(homoTarget.Length - 1488 + 1));
        Assert.That(homoDecoy.AppliedSequenceVariations[1].VariantSequence, Is.EqualTo("G"));
        Assert.That(homoDecoy.AppliedSequenceVariations[2].OriginalSequence, Is.EqualTo("C"));
        Assert.That(homoDecoy.AppliedSequenceVariations[2].OneBasedBeginPosition, Is.EqualTo(homoTarget.Length - 1222 + 1));
        Assert.That(homoDecoy.AppliedSequenceVariations[2].VariantSequence, Is.EqualTo("A"));

        var plusOneHeteroDecoy = variantRna[3];
        Assert.That(plusOneHeteroDecoy.IsDecoy, Is.True);
        Assert.That(plusOneHeteroDecoy.AppliedSequenceVariations.Count, Is.EqualTo(4));
        Assert.That(plusOneHeteroDecoy.AppliedSequenceVariations[0].OriginalSequence, Is.EqualTo("C"));
        Assert.That(plusOneHeteroDecoy.AppliedSequenceVariations[0].OneBasedBeginPosition, Is.EqualTo(plusOneHeteroTarget.Length - 1646 + 1));
        Assert.That(plusOneHeteroDecoy.AppliedSequenceVariations[0].VariantSequence, Is.EqualTo("A"));
        Assert.That(plusOneHeteroDecoy.AppliedSequenceVariations[1].OriginalSequence, Is.EqualTo("C"));
        Assert.That(plusOneHeteroDecoy.AppliedSequenceVariations[1].OneBasedBeginPosition, Is.EqualTo(plusOneHeteroTarget.Length - 1488 + 1));
        Assert.That(plusOneHeteroDecoy.AppliedSequenceVariations[1].VariantSequence, Is.EqualTo("G"));
        Assert.That(plusOneHeteroDecoy.AppliedSequenceVariations[2].OriginalSequence, Is.EqualTo("C"));
        Assert.That(plusOneHeteroDecoy.AppliedSequenceVariations[2].OneBasedBeginPosition, Is.EqualTo(plusOneHeteroTarget.Length - 1222 + 1));
        Assert.That(plusOneHeteroDecoy.AppliedSequenceVariations[2].VariantSequence, Is.EqualTo("A"));
        Assert.That(plusOneHeteroDecoy.AppliedSequenceVariations[3].OriginalSequence, Is.EqualTo("A"));
        Assert.That(plusOneHeteroDecoy.AppliedSequenceVariations[3].OneBasedBeginPosition, Is.EqualTo(plusOneHeteroTarget.Length - 409 + 1));
        Assert.That(plusOneHeteroDecoy.AppliedSequenceVariations[3].VariantSequence, Is.EqualTo("U"));
    }

    [Test]
    public void VariantModificationTest()
    {
        // This creates a heterozygous variant with 2 possible mods.
        // One of the mod residues is removed by the variant. 
        string dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "VariantModsGPTMD.xml");
        List<RNA> rna = RnaDbLoader.LoadRnaXML(dbPath, true, DecoyType.Reverse, false, AllKnownMods, [], out var unknownModifications);
        Assert.That(rna.All(p => p.SequenceVariations.Count == 1));

        List<RNA> targets = rna.Where(p => p.IsDecoy == false).ToList();
        RNA variantTarget = targets.First(p => p.AppliedSequenceVariations.Count >= 1);
        RNA nonVariantTarget = targets.First(p => p.AppliedSequenceVariations.Count == 0);

        Assert.That(variantTarget.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(1));
        Assert.That(nonVariantTarget.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(2));

        List<RNA> decoys = rna.Where(p => p.IsDecoy).ToList();
        RNA variantDecoy  = decoys.First(p => p.AppliedSequenceVariations.Count >= 1);
        RNA nonVariantDecoy = decoys.First(p => p.AppliedSequenceVariations.Count == 0);

        Assert.That(variantDecoy.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(1));
        Assert.That(nonVariantDecoy.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(2));

        var digestionParams = new RnaDigestionParams("top-down");
        List<OligoWithSetMods> oligos = rna.SelectMany(p => p.Digest(digestionParams, [], [])).ToList();

        string[] targetDigestedSequences = new[]
        {
            // Target Base Sequence and after application of 2 mods in database
            "GUACUGUAGCCUA", "GUA[Biological:Methylation on A]CUGUAGCCUA",
            "GUACUGUAGCCU[Biological:Methylation on U]A", "GUA[Biological:Methylation on A]CUGUAGCCU[Biological:Methylation on U]A",

            // Decoy Base Sequence and after application of 2 mods in database
            "AUCCGAUGUCAUG", "AUCCGAUGUCA[Biological:Methylation on A]UG",
            "AU[Biological:Methylation on U]CCGAUGUCAUG", "AU[Biological:Methylation on U]CCGAUGUCA[Biological:Methylation on A]UG", 

            // Target With Sequence Variant A3->U
            "GUUCUGUAGCCUA",
            "GUUCUGUAGCCU[Biological:Methylation on U]A",

            // Decoy With Sequence Variant A3->U
            "AUCCGAUGUCUUG",
            "AU[Biological:Methylation on U]CCGAUGUCUUG",
        };

        Assert.That(oligos.Count, Is.EqualTo(targetDigestedSequences.Length));
        for (int i = 0; i < oligos.Count; i++)
        {
            Assert.That(targetDigestedSequences.Contains(oligos[i].FullSequence));
        }
    }

    [Test]
    public void TwoTruncationsAndSequenceVariant_DbLoading()
    {
        string dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "TruncationAndVariantMods.xml");
        List<RNA> rna = RnaDbLoader.LoadRnaXML(dbPath, true, DecoyType.Reverse, false, AllKnownMods, [], out var unknownModifications);

        Assert.That(rna.All(p => p.SequenceVariations.Count == 1));
        Assert.That(rna.All(p => p.OriginalNonVariantModifications.Count == 2));

        List<RNA> targets = rna.Where(p => p.IsDecoy == false).ToList();
        RNA variantTarget = targets.First(p => p.AppliedSequenceVariations.Count >= 1);
        RNA nonVariantTarget = targets.First(p => p.AppliedSequenceVariations.Count == 0);

        Assert.That(variantTarget.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(1));
        Assert.That(nonVariantTarget.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(2));

        List<RNA> decoys = rna.Where(p => p.IsDecoy).ToList();
        RNA variantDecoy = decoys.First(p => p.AppliedSequenceVariations.Count >= 1);
        RNA nonVariantDecoy = decoys.First(p => p.AppliedSequenceVariations.Count == 0);

        Assert.That(variantDecoy.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(1));
        Assert.That(nonVariantDecoy.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(2));
    }

    [Test]
    [TestCase("NonVariantTarget", "GUACUGUAGCCUA", 0, new[] { "UACUG", "UAG", "CCUA", "UA[Biological:Methylation on A]CUG", "CCU[Biological:Methylation on U]A", "CUG" } )]
    [TestCase("VariantTarget", "GUUCUGUAGCCUA", 0, new[] { "UUCUG", "UAG", "CCUA", "CCU[Biological:Methylation on U]A", "CUG" } )]
    [TestCase("NonVariantDecoy", "AUCCGAUGUCAUG", 0, new[] { "AUCCG", "AUG", "UCAUG", "UCA[Biological:Methylation on A]UG", "AU[Biological:Methylation on U]CCG", "UG", "UC" } )]
    [TestCase("VariantDecoy", "AUCCGAUGUCUUG", 0, new[] { "AUCCG", "AUG", "UCUUG", "AU[Biological:Methylation on U]CCG", "UC", "UG" } )]
    [TestCase("NonVariantTarget", "GUACUGUAGCCUA", 1, new[] { "UACUG", "UAG", "CCUA", "UA[Biological:Methylation on A]CUG", "CCU[Biological:Methylation on U]A", "CUG", "GUACUG", "UACUGUAG", "GUA[Biological:Methylation on A]CUG", "UA[Biological:Methylation on A]CUGUAG", "UAGCCUA", "UAGCCU[Biological:Methylation on U]A", "UACUGU", "UA[Biological:Methylation on A]CUGU", "CUGUAG" } )]
    [TestCase("VariantTarget", "GUUCUGUAGCCUA", 1, new[] { "UUCUG", "UAG", "CCUA", "CCU[Biological:Methylation on U]A", "CUG", "GUUCUG", "UUCUGUAG", "UAGCCUA", "UAGCCU[Biological:Methylation on U]A", "CUGUAG", "UUCUGU" } )]
    [TestCase("NonVariantDecoy", "AUCCGAUGUCAUG", 1, new[] { "AUCCG", "AUG", "UCAUG", "UCA[Biological:Methylation on A]UG", "AU[Biological:Methylation on U]CCG", "UG", "UC", "AUCCGAUG", "AU[Biological:Methylation on U]CCGAUG", "AUGUCAUG", "AUGUCA[Biological:Methylation on A]UG", "AUGUC", "UGUCAUG", "UGUCA[Biological:Methylation on A]UG" } )]
    [TestCase("VariantDecoy", "AUCCGAUGUCUUG", 1, new[] { "AUCCG", "AUG", "UCUUG", "AU[Biological:Methylation on U]CCG", "UC", "UG", "AUCCGAUG", "AU[Biological:Methylation on U]CCGAUG", "AUGUCUUG", "AUGUC", "UGUCUUG" } )]
    public void TwoTruncationsAndSequenceVariant_Digestion(string testCase, string baseSequence, int missedCleavages, string[] expectedSequences)
    {
        string dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "TruncationAndVariantMods.xml");
        List<RNA> rna = RnaDbLoader.LoadRnaXML(dbPath, true, DecoyType.Reverse, false, AllKnownMods, [], out var unknownModifications);
        RnaDigestionParams digestionParams = new RnaDigestionParams("RNase T1", missedCleavages, 2);

        Assert.That(rna.All(p => p.SequenceVariations.Count == 1));
        Assert.That(rna.All(p => p.OriginalNonVariantModifications.Count == 2));
        Assert.That(rna.All(p => p.TruncationProducts.Count == 2));

        RNA toDigest = testCase switch
        {
            "NonVariantTarget" => rna[0],
            "VariantTarget" => rna[1],
            "NonVariantDecoy" => rna[2],
            "VariantDecoy" => rna[3],
            _ => throw new ArgumentException("Invalid test case")
        };

        var (truncation1, truncation2, expectedModCount) = testCase switch
        {
            "NonVariantTarget" => ((4, 13), (1, 7), 2),
            "VariantTarget" => ((4, 13), (1, 7), 1),
            "NonVariantDecoy" => ((1, 10), (7, 13), 2),
            "VariantDecoy" => ((1, 10), (7, 13), 1),
            _ => throw new ArgumentException("Invalid test case")
        };

        Assert.That(toDigest.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(expectedModCount));
        Assert.That(toDigest.TruncationProducts[0].OneBasedBeginPosition, Is.EqualTo(truncation1.Item1));
        Assert.That(toDigest.TruncationProducts[0].OneBasedEndPosition, Is.EqualTo(truncation1.Item2));
        Assert.That(toDigest.TruncationProducts[1].OneBasedBeginPosition, Is.EqualTo(truncation2.Item1));
        Assert.That(toDigest.TruncationProducts[1].OneBasedEndPosition, Is.EqualTo(truncation2.Item2));

        var oligos = toDigest.Digest(digestionParams, [], []).ToList();
        Assert.That(oligos.Count, Is.EqualTo(expectedSequences.Length));
        foreach (var oligo in oligos)
        {
            Assert.That(expectedSequences.Contains(oligo.FullSequence));
        }
    }
}
