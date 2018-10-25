using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;
using System.Data;

namespace Test
{
    [TestFixture]
    public class VariantProteinTests
    {
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
        public static void LoadSeqVarModifications()
        {
            var proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "oblm2.xml"), true,
                DecoyType.None, null, false, null, out var unknownModifications);
            Assert.AreEqual(1, proteins[0].OneBasedPossibleLocalizedModifications.Count);
            Assert.AreEqual(1, proteins[0].AppliedSequenceVariations.Count());
            Assert.AreEqual(1, proteins[0].SequenceVariations.Count());
            Assert.AreEqual(1, proteins[0].SequenceVariations.First().OneBasedModifications.Count);

            ProteinDbWriter.WriteXmlDatabase(null, proteins, Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "oblm2rewrite.xml"));
            proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "oblm2rewrite.xml"), true,
                DecoyType.None, null, false, null, out unknownModifications);
            Assert.AreEqual(1, proteins[0].OneBasedPossibleLocalizedModifications.Count);
            Assert.AreEqual(1, proteins[0].SequenceVariations.Count());
            Assert.AreEqual(1, proteins[0].SequenceVariations.First().OneBasedModifications.Count);
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
            Assert.AreEqual(403-11, applied.OriginalSequence.Length - applied.VariantSequence.Length);
            Assert.AreEqual(1, proteins[1].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes
            Assert.AreEqual(873, proteins[0].Length);
            Assert.AreEqual(873-403+11, proteins[1].Length);
        }

        //[Test]
        //public static void HomoHetero()
        //{
        //    //var proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"F:\ProjectsActive\Spritz\jeko\SRX277277_1-trimmed-pair1Aligned.sortedByCoord.outProcessed.out.fixedQuals.split.NoIndels.snpEffAnnotated.protein.withmods.xml"), true,
        //    //    DecoyType.None, null, false, null, out var unknownModifications);
        //    var proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"E:\ProjectsActive\MCF7PacBio\IsoSeq_MCF7_2015edition_polished.unimapped.withcds.protein.withmods.xml"), true,
        //        DecoyType.None, null, false, null, out var unknownModifications);
        //    var variantProteins = proteins.SelectMany(p => p.GetVariantProteins()).ToList();
        //    proteins[0].GetHashCode();
        //    variantProteins[0].GetHashCode();
        //    //DataTable table = new DataTable();
        //    //table.Columns.Add("protein");
        //    //table.Columns.Add("heteroVariantCount");
        //    //table.Columns.Add("homoVariantCount");
        //    //table.Columns.Add("size_da");
        //    //foreach (var p in proteins)
        //    //{
        //    //    DataRow row = table.NewRow();
        //    //    row[0] = p.Accession;
        //    //    int homo = p.SequenceVariations.Count(x => x.Description.Genotypes.First().Value.Distinct().Count() == 1);
        //    //    int hetero = p.SequenceVariations.Count(x => x.Description.Genotypes.First().Value.Distinct().Count() > 1);
        //    //    row[1] = hetero.ToString();
        //    //    row[2] = homo.ToString();
        //    //    row[3] = new Proteomics.AminoAcidPolymer.Peptide(p.BaseSequence).MonoisotopicMass.ToString();
        //    //    table.Rows.Add(row);
        //    //}
        //    //var builder = new System.Text.StringBuilder();
        //    //foreach (DataRow row in table.Rows)
        //    //{
        //    //    builder.AppendLine(string.Join("\t", row.ItemArray));
        //    //}
        //    //File.WriteAllText(@"E:\ProjectsActive\JurkatProteogenomics\180831.1WithFixedSeqVarAndTranscriptIsoforms\homoheterotable.txt", builder.ToString());
        //}

        [Test]
        public void VariantLongDeletionXml()
        {
            string file = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "SeqVarLongDeletion.xml");
            List<Protein> variantProteins = ProteinDbLoader.LoadProteinXML(file, true, DecoyType.None, null, false, null, out var un);

            Assert.AreEqual(2, variantProteins.First().NonVariantProtein.SequenceVariations.Count());
            Assert.AreEqual(1, variantProteins.Count); // there is only one unique amino acid change
            Assert.AreNotEqual(variantProteins.First().NonVariantProtein.BaseSequence, variantProteins.First().BaseSequence);
            Assert.AreEqual('A', variantProteins.First().NonVariantProtein.BaseSequence[226]);
            Assert.AreNotEqual('A', variantProteins.First().BaseSequence[226]);
            Assert.AreNotEqual(variantProteins.First().NonVariantProtein.Name, variantProteins.First().Name);
            Assert.AreNotEqual(variantProteins.First().NonVariantProtein.FullName, variantProteins.First().FullName);
            Assert.AreNotEqual(variantProteins.First().NonVariantProtein.Accession, variantProteins.First().Accession);

            List<PeptideWithSetModifications> peptides = variantProteins.SelectMany(vp => vp.Digest(new DigestionParams(), null, null)).ToList();
        }

        //[Test]
        //public void test()
        //{
        //    string file = @"E:\ProjectsActive\JurkatProteogenomics\180413\combined_1-trimmed-pair1Aligned.sortedByCoord.outProcessed.out.fixedQuals.split.genotyped.NoIndels.snpEffAnnotated.protein.withmods.xml";
        //    List<Protein> variantProteins = ProteinDbLoader.LoadProteinXML(file, true, DecoyType.Reverse, null, false, null, out var un);

        //    List<PeptideWithSetModifications> peptides = variantProteins.SelectMany(vp => vp.Digest(new DigestionParams(), null, null)).ToList();
        //}

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