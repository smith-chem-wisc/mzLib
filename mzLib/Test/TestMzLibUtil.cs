using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using MzLibUtil;
using Readers;
using System.Collections.Generic;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestMzLibUtil
    {
        [Test]
        [TestCase(@"C:\Users\bubba\Documents\Projects\K562\K562_2\20100730_Velos1_TaGe_SA_K565_4.raw", "20100730_Velos1_TaGe_SA_K565_4")]
        [TestCase("C:\\Users\bubba\\Documents\\Projects\\K562\\K562_2\\20100730_Velos1_TaGe_SA_K565_4.raw", "20100730_Velos1_TaGe_SA_K565_4")]
        [TestCase("20100730_Velos1_TaGe_SA_K565_4.raw", "20100730_Velos1_TaGe_SA_K565_4")]
        [TestCase("20100730_Velos1_TaGe_SA_K565_4", "20100730_Velos1_TaGe_SA_K565_4")]
        [TestCase(@"C:\Users\bubba\Documents\Projects\K562\K562_2\20100730_Velos1_TaGe_SA_K565_4", "20100730_Velos1_TaGe_SA_K565_4")]
        //test extra period in folder name of path
        [TestCase(@"C:\Users\bubba\Documents.docs\Projects\K562\K562_2\20100730_Velos1_TaGe_SA_K565_4.raw", "20100730_Velos1_TaGe_SA_K565_4")]
        //test extra period in filename
        [TestCase(@"C:\Users\bubba\Documents\Projects\K562\K562_2\20100730_Velos1_.TaGe_SA_K565_4.raw", "20100730_Velos1_.TaGe_SA_K565_4")]
        [TestCase("/home/seth/Pictures/penguin.jpg","penguin")]
        [TestCase("/home/seth/Pictures/penguin", "penguin")]
        [TestCase("penguin.jpg", "penguin")]
        [TestCase("penguin", "penguin")]
        [TestCase("penguin.jpg.gz", "penguin")]
        [TestCase("penguin.jpg.zip", "penguin")]
        [TestCase("penguin.jpg.mzXML", "penguin.jpg")]
        public static void TestPeriodTolerantFilenameWithoutExtension(string filenameAndOrPath, string expectedResult)
        {
            string result = PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(filenameAndOrPath);
            string extensionResult = filenameAndOrPath.GetPeriodTolerantFilenameWithoutExtension();
            Assert.AreEqual(expectedResult, result);
            Assert.AreEqual(expectedResult, extensionResult);
        }
        [Test]
        public static void TestParseModificationsOneMod()
        {
            string fullSeq = "DM[Common Variable:Oxidation on M]MELVQPSISGVDLDK";
            var mods = fullSeq.ParseModifications();
            Assert.That(mods.Count == 1);
            Assert.That(mods.ContainsKey(2));
            Assert.That(mods[2].Count == 1);
            Assert.That(mods[2].Contains("Common Variable:Oxidation on M"));
        }
        public static void TestParseModificationsTwoModsZeroIndexing()
        {
            // sequence with two terminal mods with indexed termini (zero-based indexing)
            string fullSeq = "[UniProt:N-acetylglutamate on E]EEEIAALVID[Metal:Calcium on D]NGSGMC[Common Fixed:Carbamidomethyl on C]";
            var mods = fullSeq.ParseModifications(true, true);
            Assert.That(mods.Count == 3);
            Assert.That(mods.ContainsKey(0));
            Assert.That(mods.ContainsKey(10));
            Assert.That(mods.ContainsKey(17));
            Assert.That(mods[0].Count == 1);
            Assert.That(mods[10].Count == 1);
            Assert.That(mods[17].Count == 1);
            Assert.That(mods[0].Contains("UniProt:N-acetylglutamate on E"));
            Assert.That(mods[10].Contains("Metal:Calcium on D"));
            Assert.That(mods[17].Contains("Common Fixed:Carbamidomethyl on C"));
        }
        public static void TestParseModificationsTwoModsOneIndexing()
        {
            // sequence with two terminal mods with termini indexed at amino acid positions (one-based indexing)
            string fullSeq = "[UniProt:N-acetylglutamate on E]EEEIAALVID[Metal:Calcium on D]NGSGMC[Common Fixed:Carbamidomethyl on C]";
            var mods = fullSeq.ParseModifications();
            Assert.That(mods.Count == 3);
            Assert.That(mods.ContainsKey(1));
            Assert.That(mods.ContainsKey(10));
            Assert.That(mods.ContainsKey(16));
            Assert.That(mods[1].Count == 1);
            Assert.That(mods[10].Count == 1);
            Assert.That(mods[16].Count == 1);
            Assert.That(mods[1].Contains("UniProt:N-acetylglutamate on E"));
            Assert.That(mods[10].Contains("Metal:Calcium on D"));
            Assert.That(mods[16].Contains("Common Fixed:Carbamidomethyl on C"));
        }
        public static void TestParseModificationsTwoModsSameTerminusZeroIndexing()
        {
            // sequence with two mods on same terminus with with indexed termini  (zero-based indexing)
            string fullSeq = "[UniProt:N-acetylglutamate on E]|[Common Artifact:Water Loss on E]EEEIAALVID[Metal:Calcium on D]NGSGMC[Common Fixed:Carbamidomethyl on C]";
            var mods = fullSeq.ParseModifications(true, true);
            Assert.That(mods.Count == 3);
            Assert.That(mods.ContainsKey(0));
            Assert.That(mods.ContainsKey(10));
            Assert.That(mods.ContainsKey(17));
            Assert.That(mods[0].Count == 2);
            Assert.That(mods[10].Count == 1);
            Assert.That(mods[17].Count == 1);
            Assert.That(mods[0].Contains("UniProt:N-acetylglutamate on E"));
            Assert.That(mods[0].Contains("Common Artifact:Water Loss on E"));
            Assert.That(mods[10].Contains("Metal:Calcium on D"));
            Assert.That(mods[17].Contains("Common Fixed:Carbamidomethyl on C"));
        }
        public static void TestParseModificationsTwoModsSameTerminusOneIndexing()
        {
            // sequence with two mods on same terminus with termini indexed at amino acid positions (one-based indexing)
            string fullSeq = "[UniProt:N-acetylglutamate on E]|[Common Artifact:Water Loss on E]EEEIAALVID[Metal:Calcium on D]NGSGMC[Common Fixed:Carbamidomethyl on C]";
            var mods = fullSeq.ParseModifications();
            Assert.That(mods.Count == 3);
            Assert.That(mods.ContainsKey(1));
            Assert.That(mods.ContainsKey(10));
            Assert.That(mods.ContainsKey(16));
            Assert.That(mods[1].Count == 2);
            Assert.That(mods[10].Count == 1);
            Assert.That(mods[16].Count == 1);
            Assert.That(mods[1].Contains("UniProt:N-acetylglutamate on E"));
            Assert.That(mods[1].Contains("Common Artifact:Water Loss on E"));
            Assert.That(mods[10].Contains("Metal:Calcium on D"));
            Assert.That(mods[16].Contains("Common Fixed:Carbamidomethyl on C"));
        }
        public static void TestParseModificationsTwoModsTerminusAndSideChain()
        {
            // sequence with mod on N terminus and mod on first amino acid side chain
            string fullSeq = "[UniProt:N-acetylglutamate on E]E[Metal:Sodium on E]EEIAALVID[Metal:Calcium on D]NGSGMC[Common Fixed:Carbamidomethyl on C]";
            var mods = fullSeq.ParseModifications();
            Assert.That(mods.Count == 3);
            Assert.That(mods.ContainsKey(1));
            Assert.That(mods.ContainsKey(10));
            Assert.That(mods.ContainsKey(16));
            Assert.That(mods[1].Count == 2);
            Assert.That(mods[10].Count == 1);
            Assert.That(mods[16].Count == 1);
            Assert.That(mods[1].Contains("UniProt:N-acetylglutamate on E"));
            Assert.That(mods[1].Contains("Metal:Sodium on E"));
            Assert.That(mods[10].Contains("Metal:Calcium on D"));
            Assert.That(mods[16].Contains("Common Fixed:Carbamidomethyl on C"));
        }

        [Test]
        public void TestPeptidePTMOccupancy()
        {
            List<string> sequences = new List<string>();
            sequences.Add("[UniProt: N - acetylglutamate on E]EEEIAALVID[Metal: Calcium on D]NGSGMC[Common Fixed: Carbamidomethyl on C]K");
            sequences.Add("[UniProt: N - acetylglutamate on E]EEEIAALVID[Metal: Sodium on D]NGSGMC[Common Fixed: Carbamidomethyl on C]K");
            sequences.Add("[UniProt: N - acetylglutamate on E]EEEIAALVIDN[Common Artifact: Ammonia loss on N]GSGMC[Common Fixed: Carbamidomethyl on C]K");
            sequences.Add("[UniProt: N - acetylglutamate on E]EEEIAALVIDN[Common Biological: Hydroxylation on N]GSGMC[Common Fixed: Carbamidomethyl on C]K");
            sequences.Add("[UniProt: N - acetylglutamate on E]EEEIAALVIDNGSGM[Common Variable: Oxidation on M]C[Common Fixed: Carbamidomethyl on C]K");
            sequences.Add("[UniProt: N - acetylglutamate on E]EEEIAALVIDNGSGMC[Common Fixed: Carbamidomethyl on C]K");

            string baseSeq = "EEEIAALVIDNGSGMCK"; 
            
            List<string> pgs = new List<string>();
            pgs.Add("pg1");
            pgs.Add("pg2|pg3");

            var peptides = new List<(string, string, List<string>, double)>();
            foreach (var seq in sequences)
            {
                peptides.Add((seq, baseSeq, pgs, 1.0));
            }

            PositionFrequencyAnalysis pfa = new PositionFrequencyAnalysis();
            pfa.ProteinGroupsOccupancyByPeptide(peptides);
            var occupancy = pfa.Occupancy;  
              
            Assert.That(6.0 == occupancy["pg1"].Proteins["pg1"].Peptides[baseSeq].ModifiedAminoAcidPositions[0]["UniProt: N - acetylglutamate on E"].Intensity);
            Assert.That(1.0 == occupancy["pg1"].Proteins["pg1"].Peptides[baseSeq].ModifiedAminoAcidPositions[10]["Metal: Calcium on D"].Intensity);
            Assert.That(1.0 == occupancy["pg1"].Proteins["pg1"].Peptides[baseSeq].ModifiedAminoAcidPositions[10]["Metal: Sodium on D"].Intensity);
            Assert.That(1.0 == occupancy["pg1"].Proteins["pg1"].Peptides[baseSeq].ModifiedAminoAcidPositions[11]["Common Artifact: Ammonia loss on N"].Intensity);
            Assert.That(1.0 == occupancy["pg1"].Proteins["pg1"].Peptides[baseSeq].ModifiedAminoAcidPositions[11]["Common Biological: Hydroxylation on N"].Intensity);
            Assert.That(1.0 == occupancy["pg1"].Proteins["pg1"].Peptides[baseSeq].ModifiedAminoAcidPositions[15]["Common Variable: Oxidation on M"].Intensity);
            Assert.That(6.0 == occupancy["pg1"].Proteins["pg1"].Peptides[baseSeq].ModifiedAminoAcidPositions[16]["Common Fixed: Carbamidomethyl on C"].Intensity);
            Assert.That(6.0 == occupancy["pg1"].Proteins["pg1"].Peptides[baseSeq].Intensity);

            Assert.That(6.0 == occupancy["pg2|pg3"].Proteins["pg2"].Peptides[baseSeq].ModifiedAminoAcidPositions[0]["UniProt: N - acetylglutamate on E"].Intensity);
            Assert.That(1.0 == occupancy["pg2|pg3"].Proteins["pg2"].Peptides[baseSeq].ModifiedAminoAcidPositions[10]["Metal: Calcium on D"].Intensity);
            Assert.That(1.0 == occupancy["pg2|pg3"].Proteins["pg2"].Peptides[baseSeq].ModifiedAminoAcidPositions[10]["Metal: Sodium on D"].Intensity);
            Assert.That(1.0 == occupancy["pg2|pg3"].Proteins["pg2"].Peptides[baseSeq].ModifiedAminoAcidPositions[11]["Common Artifact: Ammonia loss on N"].Intensity);
            Assert.That(1.0 == occupancy["pg2|pg3"].Proteins["pg2"].Peptides[baseSeq].ModifiedAminoAcidPositions[11]["Common Biological: Hydroxylation on N"].Intensity);
            Assert.That(1.0 == occupancy["pg2|pg3"].Proteins["pg2"].Peptides[baseSeq].ModifiedAminoAcidPositions[15]["Common Variable: Oxidation on M"].Intensity);
            Assert.That(6.0 == occupancy["pg2|pg3"].Proteins["pg2"].Peptides[baseSeq].ModifiedAminoAcidPositions[16]["Common Fixed: Carbamidomethyl on C"].Intensity);
            Assert.That(6.0 == occupancy["pg2|pg3"].Proteins["pg2"].Peptides[baseSeq].Intensity);
        }

        [Test]
        public static void TestToEnum()
        {
            Assert.IsTrue(0.ToEnum<TimsTofMsMsType>(out var result));
            Assert.AreEqual(TimsTofMsMsType.MS, result);

            Assert.IsTrue(2.ToEnum<TimsTofMsMsType>(out result));
            Assert.AreEqual(TimsTofMsMsType.MSMSFragment, result);

            Assert.IsTrue(8.ToEnum<TimsTofMsMsType>(out result));
            Assert.AreEqual(TimsTofMsMsType.PASEF, result);

            Assert.IsTrue(9.ToEnum<TimsTofMsMsType>(out result));
            Assert.AreEqual(TimsTofMsMsType.DIA, result);

            Assert.IsTrue(10.ToEnum<TimsTofMsMsType>(out result));
            Assert.AreEqual(TimsTofMsMsType.PRM, result);

            Assert.IsTrue(0.ToEnum<TimsTofAcquisitionMode>(out var result2));
            Assert.AreEqual(TimsTofAcquisitionMode.MS, result2);

            Assert.IsFalse(1.ToEnum<TimsTofMsMsType>(out result));
            Assert.IsFalse(11.ToEnum<TimsTofMsMsType>(out result));
            Assert.IsFalse(7.ToEnum<TimsTofMsMsType>(out result));
            
        }
    }
}
