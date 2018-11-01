using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public class TestVariantProtein
    {
        [Test]
        public static void VariantProtein()
        {
            Protein p = new Protein("MAAA", "accession");
            Protein v = new Protein("MAVA", p, new[] { new SequenceVariation(3, "A", "V", "desc", null) }, null, null, null);
            Assert.AreEqual(p, v.NonVariantProtein);
        }

        [Test]
        public void VariantXml()
        {
            string file = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "SeqVar.xml");
            List<Protein> variantProteins = ProteinDbLoader.LoadProteinXML(file, true, DecoyType.None, null, false, null, out var un);

            Assert.AreEqual(5, variantProteins.First().NonVariantProtein.SequenceVariations.Count());
            Assert.AreEqual(1, variantProteins.Count); // there is only one unique amino acid change
            Assert.AreNotEqual(variantProteins.First().NonVariantProtein.BaseSequence, variantProteins.First().BaseSequence);
            Assert.AreEqual('C', variantProteins.First().NonVariantProtein.BaseSequence[116]);
            Assert.AreEqual('Y', variantProteins.First().BaseSequence[116]);
            Assert.AreNotEqual(variantProteins.First().NonVariantProtein.Name, variantProteins.First().Name);
            Assert.AreNotEqual(variantProteins.First().NonVariantProtein.FullName, variantProteins.First().FullName);
            Assert.AreNotEqual(variantProteins.First().NonVariantProtein.Accession, variantProteins.First().Accession);

            List<PeptideWithSetModifications> peptides = variantProteins.SelectMany(vp => vp.Digest(new DigestionParams(), null, null)).ToList();
        }

        [Test]
        [TestCase("oblm2.xml", 3, 4)] // without starting methionine
        [TestCase("oblm3.xml", 3, 5)] // with starting methionine
        public static void LoadSeqVarModifications(string databaseName, int modIdx, int reversedModIdx)
        {
            var proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", databaseName), true,
                DecoyType.Reverse, null, false, null, out var unknownModifications);
            var target = proteins[0];
            Assert.AreEqual(1, target.OneBasedPossibleLocalizedModifications.Count);
            Assert.AreEqual(1, target.AppliedSequenceVariations.Count());
            Assert.AreEqual(1, target.SequenceVariations.Count());
            Assert.AreEqual(1, target.SequenceVariations.Single().OneBasedModifications.Count);
            Assert.AreEqual(modIdx, target.SequenceVariations.Single().OneBasedModifications.Single().Key); //PEP[mod]TID, MEP[mod]TID
            var decoy = proteins[1];
            Assert.AreEqual(1, decoy.OneBasedPossibleLocalizedModifications.Count);
            Assert.AreEqual(1, decoy.AppliedSequenceVariations.Count());
            Assert.AreEqual(1, decoy.SequenceVariations.Count());
            Assert.AreEqual(1, decoy.SequenceVariations.Single().OneBasedModifications.Count);
            Assert.AreEqual(reversedModIdx, decoy.SequenceVariations.Single().OneBasedModifications.Single().Key); //DITP[mod]EP, MDITP[mod]E

            string rewriteDbName = $"{Path.GetFileNameWithoutExtension(databaseName)}rewrite.xml";
            ProteinDbWriter.WriteXmlDatabase(null, proteins.Where(p => !p.IsDecoy).ToList(), Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", rewriteDbName));
            proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", rewriteDbName), true,
                DecoyType.Reverse, null, false, null, out unknownModifications);
            target = proteins[0];
            Assert.AreEqual(1, target.OneBasedPossibleLocalizedModifications.Count);
            Assert.AreEqual(1, target.AppliedSequenceVariations.Count());
            Assert.AreEqual(1, target.SequenceVariations.Count());
            Assert.AreEqual(1, target.SequenceVariations.Single().OneBasedModifications.Count);
            Assert.AreEqual(modIdx, target.SequenceVariations.Single().OneBasedModifications.Single().Key);
            decoy = proteins[1];
            Assert.AreEqual(1, decoy.OneBasedPossibleLocalizedModifications.Count);
            Assert.AreEqual(1, decoy.AppliedSequenceVariations.Count());
            Assert.AreEqual(1, decoy.SequenceVariations.Count());
            Assert.AreEqual(1, decoy.SequenceVariations.Single().OneBasedModifications.Count);
            Assert.AreEqual(reversedModIdx, decoy.SequenceVariations.Single().OneBasedModifications.Single().Key);
        }

        [TestCase("ranges1.xml", 1, 2, 5, 6)] // without starting methionine
        [TestCase("ranges2.xml", 1, 1, 5, 5)] // with starting methionine
        public static void ReverseDecoyProteolysisProducts(string databaseName, int beginIdx, int reversedBeginIdx, int endIdx, int reversedEndIdx)
        {
            var proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", databaseName), true,
                DecoyType.Reverse, null, false, null, out var unknownModifications);
            var target = proteins[0];
            Assert.AreEqual(1, target.ProteolysisProducts.Count());
            Assert.AreEqual(beginIdx, target.ProteolysisProducts.Single().OneBasedBeginPosition); //P[start]EPTI[end]D, M[start]EPTI[end]D
            Assert.AreEqual(endIdx, target.ProteolysisProducts.Single().OneBasedEndPosition);
            var decoy = proteins[1];
            Assert.AreEqual(1, decoy.ProteolysisProducts.Count());
            Assert.AreEqual(reversedBeginIdx, decoy.ProteolysisProducts.Single().OneBasedBeginPosition); //DI[start]TPEP[end], M[start]DITP[end]E
            Assert.AreEqual(reversedEndIdx, decoy.ProteolysisProducts.Single().OneBasedEndPosition);

            string rewriteDbName = $"{Path.GetFileNameWithoutExtension(databaseName)}rewrite.xml";
            ProteinDbWriter.WriteXmlDatabase(null, proteins.Where(p => !p.IsDecoy).ToList(), Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", rewriteDbName));
            proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", rewriteDbName), true,
                DecoyType.Reverse, null, false, null, out unknownModifications);
            target = proteins[0];
            Assert.AreEqual(1, target.ProteolysisProducts.Count());
            Assert.AreEqual(beginIdx, target.ProteolysisProducts.Single().OneBasedBeginPosition);
            Assert.AreEqual(endIdx, target.ProteolysisProducts.Single().OneBasedEndPosition);
            decoy = proteins[1];
            Assert.AreEqual(1, decoy.ProteolysisProducts.Count());
            Assert.AreEqual(reversedBeginIdx, decoy.ProteolysisProducts.Single().OneBasedBeginPosition);
            Assert.AreEqual(reversedEndIdx, decoy.ProteolysisProducts.Single().OneBasedEndPosition);
        }

        [TestCase("bonds1.xml", 2, 3, "DICPCP", 4, 5)] // without starting methionine
        [TestCase("bonds2.xml", 2, 4, "MDICPC", 4, 6)] // with starting methionine
        public static void ReverseDecoyDisulfideBonds(string databaseName, int beginIdx, int reversedBeginIdx, string reversedSequence, int endIdx, int reversedEndIdx)
        {
            var proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", databaseName), true,
                DecoyType.Reverse, null, false, null, out var unknownModifications);
            var target = proteins[0];
            Assert.AreEqual(1, target.DisulfideBonds.Count());
            Assert.AreEqual(beginIdx, target.DisulfideBonds.Single().OneBasedBeginPosition); //PC[start]PC[end]ID, MC[start]PC[end]ID
            Assert.AreEqual(endIdx, target.DisulfideBonds.Single().OneBasedEndPosition);
            var decoy = proteins[1];
            Assert.AreEqual(1, decoy.DisulfideBonds.Count());
            Assert.AreEqual(reversedSequence, decoy.BaseSequence);
            Assert.AreEqual(reversedBeginIdx, decoy.DisulfideBonds.Single().OneBasedBeginPosition); //DIC[start]PC[end]P, MDIC[start]PC[end]
            Assert.AreEqual(reversedEndIdx, decoy.DisulfideBonds.Single().OneBasedEndPosition);

            string rewriteDbName = $"{Path.GetFileNameWithoutExtension(databaseName)}rewrite.xml";
            ProteinDbWriter.WriteXmlDatabase(null, proteins.Where(p => !p.IsDecoy).ToList(), Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", rewriteDbName));
            proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", rewriteDbName), true,
                DecoyType.Reverse, null, false, null, out unknownModifications);
            target = proteins[0];
            Assert.AreEqual(1, target.DisulfideBonds.Count());
            Assert.AreEqual(beginIdx, target.DisulfideBonds.Single().OneBasedBeginPosition);
            Assert.AreEqual(endIdx, target.DisulfideBonds.Single().OneBasedEndPosition);
            decoy = proteins[1];
            Assert.AreEqual(1, decoy.DisulfideBonds.Count());
            Assert.AreEqual(reversedBeginIdx, decoy.DisulfideBonds.Single().OneBasedBeginPosition);
            Assert.AreEqual(reversedEndIdx, decoy.DisulfideBonds.Single().OneBasedEndPosition);
        }

        [Test]
        [TestCase("splices1.xml", 2, 4, 3, 5)] // range without starting methionine
        [TestCase("splices2.xml", 2, 5, 3, 6)] // range with starting methionine
        [TestCase("splices3.xml", 2, 5, 2, 5)] // site without starting methionine
        [TestCase("splices4.xml", 2, 6, 2, 6)] // site with starting methionine
        [TestCase("splices5.xml", 1, 6, 1, 6)] // start site without starting methionine
        [TestCase("splices6.xml", 1, 1, 1, 1)] // start site with starting methionine
        [TestCase("splices7.xml", 1, 5, 2, 6)] // range with start without starting methionine
        [TestCase("splices8.xml", 1, 5, 2, 6)] // range with start with starting methionine
        public static void ReverseDecoySpliceSites(string databaseName, int beginIdx, int reversedBeginIdx, int endIdx, int reversedEndIdx)
        {
            var proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", databaseName), true,
                DecoyType.Reverse, null, false, null, out var unknownModifications);
            var target = proteins[0];
            Assert.AreEqual(1, target.SpliceSites.Count());
            Assert.AreEqual(beginIdx, target.SpliceSites.Single().OneBasedBeginPosition); //PE[start]P[end]TID, ME[start]P[start]TID, PE[site]PTID, ME[site]PTID, P[site]EPTID, M[site]EPTID
            Assert.AreEqual(endIdx, target.SpliceSites.Single().OneBasedEndPosition);
            var decoy = proteins[1];
            Assert.AreEqual(1, decoy.SpliceSites.Count());
            Assert.AreEqual(reversedBeginIdx, decoy.SpliceSites.Single().OneBasedBeginPosition); //DITP[start]E[end]P, MDITP[start]E[end], DITPE[site]P, MDITPE[site], DITPEP[site], M[site]DITPE
            Assert.AreEqual(reversedEndIdx, decoy.SpliceSites.Single().OneBasedEndPosition);

            string rewriteDbName = $"{Path.GetFileNameWithoutExtension(databaseName)}rewrite.xml";
            ProteinDbWriter.WriteXmlDatabase(null, proteins.Where(p => !p.IsDecoy).ToList(), Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", rewriteDbName));
            proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", rewriteDbName), true,
                DecoyType.Reverse, null, false, null, out unknownModifications);
            target = proteins[0];
            Assert.AreEqual(1, target.SpliceSites.Count());
            Assert.AreEqual(beginIdx, target.SpliceSites.Single().OneBasedBeginPosition);
            Assert.AreEqual(endIdx, target.SpliceSites.Single().OneBasedEndPosition);
            decoy = proteins[1];
            Assert.AreEqual(1, decoy.SpliceSites.Count());
            Assert.AreEqual(reversedBeginIdx, decoy.SpliceSites.Single().OneBasedBeginPosition);
            Assert.AreEqual(reversedEndIdx, decoy.SpliceSites.Single().OneBasedEndPosition);
        }

        [Test]
        [TestCase("HomozygousHLA.xml", 1, 18)]
        [TestCase("HomozygousHLA.xml", 10, 17)]
        public static void HomozygousVariantsAtVariedDepths(string filename, int minVariantDepth, int appliedCount)
        {
            var proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", filename), true,
                DecoyType.None, null, false, null, out var unknownModifications, minAlleleDepth: minVariantDepth);
            Assert.AreEqual(1, proteins.Count);
            Assert.AreEqual(18, proteins[0].SequenceVariations.Count()); // some redundant
            Assert.AreEqual(18, proteins[0].SequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes
            Assert.AreEqual(appliedCount, proteins[0].AppliedSequenceVariations.Count()); // some redundant
            Assert.AreEqual(appliedCount, proteins[0].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes
            Assert.AreEqual(1, proteins[0].GetVariantProteins().Count);
            var variantProteins = proteins[0].GetVariantProteins();
        }

        [Test]
        public static void StopGained()
        {
            var proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "StopGained.xml"), true,
                DecoyType.None, null, false, null, out var unknownModifications);
            Assert.AreEqual(2, proteins.Count);
            Assert.AreEqual(1, proteins[0].SequenceVariations.Count()); // some redundant
            Assert.AreEqual(1, proteins[0].SequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes
            Assert.AreEqual(0, proteins[0].AppliedSequenceVariations.Count()); // some redundant
            Assert.AreEqual(0, proteins[0].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes
            Assert.AreEqual(1, proteins[1].AppliedSequenceVariations.Count()); // some redundant
            Assert.AreEqual(1, proteins[1].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes
            Assert.AreEqual(191, proteins[0].Length);
            Assert.AreEqual('Q', proteins[0][161 - 1]);
            Assert.AreEqual(161 - 1, proteins[1].Length);
            Assert.AreNotEqual(proteins[0].Length, proteins[1].Length);

            proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "StopGained.xml"), true,
                DecoyType.None, null, false, null, out unknownModifications, minAlleleDepth: 400);
            Assert.AreEqual(1, proteins.Count);
            Assert.AreEqual(1, proteins[0].AppliedSequenceVariations.Count()); // some redundant
            Assert.AreEqual(1, proteins[0].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes
            Assert.AreEqual(161 - 1, proteins[0].Length);
        }

        [Test]
        public static void MultipleAlternateAlleles()
        {
            var proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "MultipleAlternateAlleles.xml"), true,
                DecoyType.None, null, false, null, out var unknownModifications);
            Assert.AreEqual(2, proteins.Count);
            Assert.AreEqual(2, proteins[0].SequenceVariations.Count()); // some redundant
            Assert.AreEqual(2, proteins[0].SequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes

            Assert.IsTrue(proteins[0].SequenceVariations.All(v => v.OneBasedBeginPosition == 63)); // there are two alternate alleles (1 and 2), but only 2 is in the genotype, so only that's applied
            Assert.AreEqual(1, proteins[1].AppliedSequenceVariations.Count()); // some redundant
            Assert.AreEqual(1, proteins[1].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes
            Assert.AreEqual(72, proteins[0].Length);
            Assert.AreEqual(72, proteins[1].Length);
            Assert.AreEqual('K', proteins[0][63 - 1]);
            Assert.AreEqual('R', proteins[1][63 - 1]);

            proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "MultipleAlternateAlleles.xml"), true,
                DecoyType.None, null, false, null, out unknownModifications, minAlleleDepth: 10);
            Assert.AreEqual(1, proteins.Count);
            Assert.AreEqual(0, proteins[0].AppliedSequenceVariations.Count()); // some redundant
            Assert.AreEqual(0, proteins[0].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes
            Assert.AreEqual('K', proteins[0][63 - 1]); // reference only
        }

        [Test]
        public static void MultipleAlternateFrameshifts()
        {
            var proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "MultipleAlternateFrameshifts.xml"), true,
                DecoyType.None, null, false, null, out var unknownModifications);
            Assert.AreEqual(2, proteins.Count);
            Assert.AreEqual(3, proteins[0].SequenceVariations.Count()); // some redundant
            Assert.AreEqual(3, proteins[0].SequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes

            Assert.IsTrue(proteins[0].SequenceVariations.All(v => v.OneBasedBeginPosition == 471)); // there are two alternate alleles (1 and 2), but only 2 is in the genotype, so only that's applied
            Assert.AreEqual(1, proteins[1].AppliedSequenceVariations.Count()); // some redundant
            var applied = proteins[1].AppliedSequenceVariations.Single();
            Assert.AreEqual("KDKRATGRIKS", applied.VariantSequence);
            Assert.AreEqual(403 - 11, applied.OriginalSequence.Length - applied.VariantSequence.Length);
            Assert.AreEqual(1, proteins[1].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes
            Assert.AreEqual(873, proteins[0].Length);
            Assert.AreEqual(873 - 403 + 11, proteins[1].Length);
        }

        [Test]
        public void VariantSymbolWeirdnessXml()
        {
            string file = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "SeqVarSymbolWeirdness.xml");
            List<Protein> variantProteins = ProteinDbLoader.LoadProteinXML(file, true, DecoyType.None, null, false, null, out var un);
            Assert.AreEqual(12, variantProteins.First().NonVariantProtein.SequenceVariations.Count());
            Assert.AreEqual(2, variantProteins.First().NonVariantProtein.SequenceVariations.Count(v => v.Description.Heterozygous.Any(kv => kv.Value)));

            Assert.AreEqual(1, variantProteins.Count); // Should be 2^2 from combinitorics of heterozygous, but the giant indels overwrite them
            Assert.AreEqual(0, variantProteins.Where(v => v.BaseSequence == variantProteins.First().NonVariantProtein.BaseSequence).Count()); // Homozygous variations are included
            Assert.AreNotEqual(variantProteins.First().NonVariantProtein.Name, variantProteins.First().Name);
            Assert.AreNotEqual(variantProteins.First().NonVariantProtein.FullName, variantProteins.First().FullName);
            Assert.AreNotEqual(variantProteins.First().NonVariantProtein.Accession, variantProteins.First().Accession);

            List<PeptideWithSetModifications> peptides = variantProteins.SelectMany(vp => vp.Digest(new DigestionParams(), null, null)).ToList();
        }

        [Test]
        public void VariantSymbolWeirdness2Xml()
        {
            string file = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "SeqVarSymbolWeirdness2.xml");
            List<Protein> variantProteins = ProteinDbLoader.LoadProteinXML(file, true, DecoyType.None, null, false, null, out var un);

            Assert.AreEqual(1, variantProteins.First().NonVariantProtein.SequenceVariations.Count());
            Assert.AreEqual(2, variantProteins.Count); // there is only one unique amino acid change
            Assert.AreEqual(1, variantProteins.Where(v => v.BaseSequence == variantProteins.First().NonVariantProtein.BaseSequence).Count());
            var variantProteinRef = variantProteins.First();
            var variantProteinAlt = variantProteins.Last();
            Assert.AreEqual('R', variantProteins.First().NonVariantProtein.BaseSequence[2386]);
            Assert.AreEqual('R', variantProteinRef.BaseSequence[2386]);
            Assert.AreEqual('H', variantProteinAlt.BaseSequence[2386]);
            Assert.AreEqual(variantProteins.First().NonVariantProtein.Name, variantProteinRef.Name);
            Assert.AreNotEqual(variantProteins.First().NonVariantProtein.Name, variantProteinAlt.Name);
            Assert.AreEqual(variantProteins.First().NonVariantProtein.FullName, variantProteinRef.FullName);
            Assert.AreNotEqual(variantProteins.First().NonVariantProtein.FullName, variantProteinAlt.FullName);
            Assert.AreEqual(variantProteins.First().NonVariantProtein.Accession, variantProteinRef.Accession);
            Assert.AreNotEqual(variantProteins.First().NonVariantProtein.Accession, variantProteinAlt.Accession);
            List<PeptideWithSetModifications> peptides = variantProteins.SelectMany(vp => vp.Digest(new DigestionParams(), null, null)).ToList();
        }
    }
}