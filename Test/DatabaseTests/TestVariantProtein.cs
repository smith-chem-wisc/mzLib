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
    public class TestVariantProtein
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

        [Test]
        //[TestCase("trypsin")]
        //[TestCase("Lys-C (don't cleave before proline)")]
        //[TestCase("Glu-C")]
        [TestCase(@"E:\ProjectsActive\JurkatProteogenomics\180413\combined_1-trimmed-pair1Aligned.sortedByCoord.outProcessed.out.fixedQuals.split.genotyped.NoIndels.snpEffAnnotated.protein.withmods.xml", @"E:\ProjectsActive\JurkatProteogenomics\InSilicoVariantProteinDigestionTargetOnly.tsv")]
        public void test(string inputFile, string outputFile)//string protease)
        {
            List<Protein> variantProteins = ProteinDbLoader.LoadProteinXML(inputFile, true, DecoyType.None, null, false, null, out var un);
            var proteins = variantProteins.Select(vp => vp.NonVariantProtein).Distinct().ToList();
            var variantProteinsWithVariants = variantProteins.Where(vp => vp.AppliedSequenceVariations.Count > 0).ToList();
            var variantProteinsWithoutVariants = variantProteins.Except(variantProteinsWithVariants).ToList();
            double avgVariantsPerProtein = proteins.Average(p => p.SequenceVariations.Count());
            double avgVariantsAppliedPerVariantProtein = variantProteinsWithVariants.Count == 0 ? -1 : variantProteinsWithVariants.Average(p => p.AppliedSequenceVariations.Count);
            int proteinsUnder50kDa = proteins.Count(p => new Proteomics.AminoAcidPolymer.Peptide(p.BaseSequence).MonoisotopicMass < 50000);
            int variantProteinsUnder50kDa = variantProteins.Count(p => new Proteomics.AminoAcidPolymer.Peptide(p.BaseSequence).MonoisotopicMass < 50000);

            DataTable table = new DataTable();
            table.Columns.Add("protease");
            table.Columns.Add("proteins");
            table.Columns.Add("proteins_under_50kDa");
            table.Columns.Add("variant_proteins");
            table.Columns.Add("variant_proteins_under_50kDa");
            table.Columns.Add("avgVariantsPerProtein");
            table.Columns.Add("avgVariantsAppliedPerVariantProtein");
            table.Columns.Add("peptides");
            table.Columns.Add("peptides_crossing_variants");
            table.Columns.Add("unique_peptides");
            table.Columns.Add("proteins_with_unique_peptides");
            table.Columns.Add("unique_variant_peptides");
            table.Columns.Add("variant_proteins_with_unique_peptides");
            table.Columns.Add("unique_peptides_from_stop_gains");

            string[] proteases = new[] 
            {
                "Arg-C",
                "Asp-N",
                "trypsin",
                "chymotrypsin (don't cleave before proline)",
                "Glu-C",
                "Lys-C (don't cleave before proline)",
                "multiprotease"
            };

            Dictionary<string, (int, PeptideWithSetModifications)> totalBaseSequenceDict = new Dictionary<string, (int, PeptideWithSetModifications)>();
            foreach (string protease in proteases)
            {
                DataRow row = table.NewRow();
                row[0] = protease;
                row[1] = proteins.Count.ToString();
                row[2] = proteinsUnder50kDa.ToString();
                row[3] = variantProteins.Count.ToString();
                row[4] = variantProteinsUnder50kDa.ToString();
                row[5] = avgVariantsPerProtein.ToString();
                row[6] = avgVariantsAppliedPerVariantProtein.ToString();

                Dictionary<string, (int, PeptideWithSetModifications)> baseSequenceDict = new Dictionary<string, (int, PeptideWithSetModifications)>();
                if (protease != proteases.Last())
                {
                    List<PeptideWithSetModifications> peptides = variantProteins.SelectMany(vp => vp.Digest(new DigestionParams(protease), null, null)).ToList();
                    foreach (var pep in peptides)
                    {
                        if (baseSequenceDict.TryGetValue(pep.BaseSequence, out var peps)) peps.Item1 += 1;
                        else baseSequenceDict.Add(pep.BaseSequence, (1, pep));

                        if (totalBaseSequenceDict.TryGetValue(pep.BaseSequence, out peps)) peps.Item1 += 1;
                        else totalBaseSequenceDict.Add(pep.BaseSequence, (1, pep));
                    }
                    var peptidesCrossingVariants = peptides.Where(pep => pep.Protein.AppliedSequenceVariations.Any(v => pep.OneBasedStartResidueInProtein <= v.OneBasedBeginPosition && pep.OneBasedEndResidueInProtein >= v.OneBasedEndPosition)).ToList();
                    row[7] = peptides.Count.ToString();
                    row[8] = peptidesCrossingVariants.Count.ToString();
                }
                else
                {
                    row[7] = totalBaseSequenceDict.Keys.Count;
                    row[8] = "-1";
                }

                // TODO: ask about splice junctions 
                // 1) what percent of splice junction peptides are unique to a protein?
                // 2) what percent of peptides cross splice junction events?
                // 3) what percent of unique peptides cross splice junctions?
                // 4) in the pacbio and stringtie databases, are these numbers any different? what about subsetting out the ones unique to the new proteins
                // TODO: ask about variant proteins and peptides
                // 1) how about capturing all pwsms for each base sequence; how many variant peptides aren't unique because they're covered by multiple peptides
                // 2) considering those validating peptides of a variant, how many variant peptides can be validated? Higher than the 20% seen without that in mind?
                // 3) how many proteins have variants annotated?
                // 4) how many variant proteins does one with variants produce on average?
                var uniquePeptides = (protease != proteases.Last() ? baseSequenceDict : totalBaseSequenceDict).Where(kv => kv.Value.Item1 == 1).ToList();
                var uniquePepProteins = uniquePeptides.Select(kv => kv.Value.Item2.Protein).Distinct().ToList();
                var uniqueVariantPeptides = uniquePeptides.Where(kv => kv.Value.Item2.Protein.AppliedSequenceVariations.Count > 0 && kv.Value.Item2.Protein.AppliedSequenceVariations.Any(v => kv.Value.Item2.OneBasedEndResidueInProtein >= v.OneBasedBeginPosition && kv.Value.Item2.OneBasedStartResidueInProtein <= v.OneBasedEndPosition)).ToList();
                var uniqueVariantPepProteins = uniqueVariantPeptides.Select(kv => kv.Value.Item2.Protein).Distinct().ToList();
                double avgVarLength = uniqueVariantPeptides.Count == 0 ? -1 : uniqueVariantPeptides.Average(kv => kv.Value.Item2.Protein.AppliedSequenceVariations.Average(v => v.VariantSequence.Length));
                var stopGainUniqueVariantPeptides = uniqueVariantPeptides.Where(kv => kv.Value.Item2.OneBasedEndResidueInProtein == kv.Value.Item2.Protein.Length && kv.Value.Item2.Protein.AppliedSequenceVariations.Any(v => v.VariantSequence.EndsWith("*"))).ToList();

                row[9] = uniquePeptides.Count.ToString();
                row[10] = uniquePepProteins.Count.ToString();
                row[11] = uniqueVariantPeptides.Count.ToString();
                row[12] = uniqueVariantPepProteins.Count.ToString();
                row[13] = stopGainUniqueVariantPeptides.Count.ToString();

                table.Rows.Add(row);
            }

            var builder = new System.Text.StringBuilder();
            foreach (DataColumn column in table.Columns)
            {
                builder.Append($"{column.ColumnName}\t");
            }
            builder.Append(System.Environment.NewLine);
            foreach (DataRow row in table.Rows)
            {
                builder.AppendLine(string.Join("\t", row.ItemArray));
            }
            File.WriteAllText(outputFile, builder.ToString());
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