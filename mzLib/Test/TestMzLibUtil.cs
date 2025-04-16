using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using MzLibUtil;
using Readers;
using System.Linq;
using System.IO;
using System.Text.RegularExpressions;
using System.Collections.Generic;
using System;
using System.Runtime.InteropServices;
using System.Windows.Media;

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
        public static void TestParseModificationsSideChainModOnly()
        {
            string fullSeq = "DM[Common Variable:Oxidation on M]MELVQPSISGVDLDK";
            var mods = fullSeq.ParseModifications(ignoreTerminusMod: false);
            Assert.That(mods.Count == 1);
            Assert.That(mods.ContainsKey(2));
            Assert.That(mods[2] == ("Common Variable:Oxidation on M"));
        }

        [Test]
        public static void TestParseModificationsSideChainAndTerminusMods()
        {
            string fullSeq = "[UniProt:N-acetylglutamate on E]EDM[Common Variable:Oxidation on M]MELVQPSISGVDLDK[Test Mod2: ModName2 on K]-[Test Mod: ModName on K C-Terminus]";
            var mods = fullSeq.ParseModifications(ignoreTerminusMod: false);
            Assert.That(mods.Count == 4);
            Assert.That(mods.ContainsKey(0));
            Assert.That(mods.ContainsKey(3));
            Assert.That(mods.ContainsKey(18));
            Assert.That(mods.ContainsKey(19));
            Assert.That(mods[0] == "UniProt:N-acetylglutamate on E");
            Assert.That(mods[3] == "Common Variable:Oxidation on M");
            Assert.That(mods[18] == "Test Mod2: ModName2 on K");
            Assert.That(mods[19] == "Test Mod: ModName on K C-Terminus");
        }

        [Test]
        public static void TestParseModificationsIgnoreTerminusMod()
        {
            string fullSeq = "[UniProt:N-acetylglutamate on E]EDM[Common Variable:Oxidation on M]MELVQPSISGVDLDK[Test Mod2: ModName2 on K]-[Test Mod: ModName on K C-Terminus]";
            var mods = fullSeq.ParseModifications(ignoreTerminusMod: true);
            Assert.That(mods.Count == 2);
            Assert.That(mods.ContainsKey(3));
            Assert.That(mods.ContainsKey(18));
            Assert.That(mods[3] == "Common Variable:Oxidation on M");
            Assert.That(mods[18] == "Test Mod2: ModName2 on K");
        }

        [Test]
        public static void TestParseModificationsWithTsvExamples()
        {

            var path = @"ModificationTests\ModifiedFullSequencesAndModificationsExamples.txt";
            var lines = File.ReadAllLines(path);
            var header = lines.First().Split('\t');
            foreach (var line in lines.Skip(1))
            {
                if (!line.Contains('|')) // Skip any ambiguous sequences
                {
                    var parts = line.Split('\t');
                    var fullSeq = parts[1];
                    Regex expectedModsPattern = new(@"(?<=on [A-Z])\s(?=[A-Z])");
                    var expectedMods = string.Join(' ', expectedModsPattern.Split(parts[2]).ToList().Order()); // Sort the mods for consitency with foundMods
                    var mods = fullSeq.ParseModifications();
                    var foundMods = string.Join(' ', mods.Values.Select(x=> x.Split(':')[1]).ToList().Order());

                    Assert.AreEqual(expectedMods, foundMods);
                }
            }
        }

        [Test]
        public static void TestParseBaseSequence()
        {
            string fullSeq = "[UniProt:N-acetylglutamate on E]EDM[Common Variable:Oxidation on M]MELVQPSISGVDLDK[Test Mod2: ModName2 on K]-[Test Mod: ModName on K C-Terminus]";
            var baseSeq = fullSeq.ParseBaseSequence();
            Assert.That(baseSeq == "EDMMELVQPSISGVDLDK");
        }

        [Test]
        public static void TestQuantifiedClasses()
        {

            string fullSeq1 = "[UniProt: N-acetylglutamate on E]EDM[Common Variable1: Oxidation on M]AAAAAAK[Test Mod1: ModName1 on K]-[Test Mod: ModName on K C-Terminus]";
            string fullSeq2 = "[UniProt: N-acetylglutamate on E]EDM[Common Variable2: Oxidation on M]AAAAAAK[Test Mod1: ModName1 on K]-[Test Mod: ModName on K C-Terminus]";
            string fullSeq3 = "[UniProt: N-acetylglutamate on E]EDM[Common Variable3: Oxidation on M]IIIIIIK[Test Mod3: ModName3 on K]-[Test Mod: ModName on K C-Terminus]";
            string fullSeq4 = "[UniProt: N-acetylglutamate on E]EDM[Common Variable4: Oxidation on M]QQQQQQK[Test Mod4: ModName4 on K]-[Test Mod: ModName on K C-Terminus]";

            var fullSequences = new List<string> { fullSeq1, fullSeq2, fullSeq3, fullSeq4};
            var baseSequences = fullSequences.Select(x => x.ParseBaseSequence()).ToList();
            var oneBasedPeptideIndexInProtein = new List<int> { 1, 1, 11, 21 };
            var proteinSequence = string.Join("", baseSequences.ToHashSet());
            string proteinName = "TestProtein";

            var peptides = new Dictionary<string, QuantifiedPeptide> ();

            for (int i = 0; i < oneBasedPeptideIndexInProtein.Count; i++)
            {
                if (!peptides.ContainsKey(baseSequences[i]))
                {
                    peptides.Add(baseSequences[i], new QuantifiedPeptide(fullSequences[i], oneBasedPeptideIndexInProtein[i], intensity: 2 * i + 1));
                    peptides[baseSequences[i]].AddFullSequence(baseSequences[i], intensity: 1);
                }
                else
                {
                    peptides[baseSequences[i]].AddFullSequence(fullSequences[i], intensity: 1);
                }
            }

            var protein = new QuantifiedProtein(proteinName, proteinSequence, peptides);
            var stoich = protein.GetModStoichiometryFromProteinMods();

            Assert.AreEqual(stoich.Count, 8);
            Assert.AreEqual(stoich.Keys.ToList(), new List<int>{0, 3, 10, 13, 20, 23, 30, 31});
            //Assert.AreEqual(stoich.Values.SelectMany(x => x.Keys).ToHashSet().Order(), fullSequences.Select(x => x.ParseModifications()).SelectMany(x => x.Values).ToHashSet().Order());
            Assert.AreEqual(stoich[0]["UniProt: N-acetylglutamate on E"].Intensity, 2.0 / 3.0);
            Assert.AreEqual(stoich[3]["Common Variable1: Oxidation on M"].Intensity, 1.0 / 3.0);
            Assert.AreEqual(stoich[3]["Common Variable2: Oxidation on M"].Intensity, 1.0 / 3.0);
            Assert.AreEqual(stoich[10]["Test Mod1: ModName1 on K"].Intensity, 2.0 / 3.0);
            Assert.AreEqual(stoich[13]["Common Variable3: Oxidation on M"].Intensity, 5.0 / 6.0);
            Assert.AreEqual(stoich[20]["Test Mod3: ModName3 on K"].Intensity, 5.0 / 6.0);
            Assert.AreEqual(stoich[23]["C ommon Variable4: Oxidation on M"].Intensity, 7.0 / 8.0);
            Assert.AreEqual(stoich[30]["Test Mod4: ModName4 on K"].Intensity, 7.0 / 8.0);
            Assert.AreEqual(stoich[31]["Test Mod: ModName on K C-Terminus"].Intensity, 7.0 / 8.0);
        }

        [Test]
        public void TestProteinGroupsOccupancyByPeptide()
        {
            var fullSeq1 = "[UniProt: N-acetylglutamate on E]EDM[Common Variable1: Oxidation on M]AAAAAAK[Test Mod1: ModName1 on K]-[Test Mod: ModName on K C-Terminus]";
            var fullSeq2 = "[UniProt: N-acetylglutamate on E]EDMAAAAAAK[Test Mod1: ModName1 on K]-[Test Mod: ModName on K C-Terminus]";
            var fullSeq3 = "[UniProt: N-acetylglutamate on E]EDM[Common Variable1: Oxidation on M]QQQQQQK[Test Mod1: ModName1 on K]-[Test Mod: ModName on K C-Terminus]";


            var peptides = new List<(string fullSeq, List<string> proteinGroup, double intensity)>
            {
                (fullSeq1, new List<string>{"Protein1|Protein2", "Protein3"}, 2),
                (fullSeq1.ParseBaseSequence(), new List<string>{"Protein1|Protein2", "Protein3"}, 1),
                (fullSeq2, new List<string>{"Protein3"}, 2),
                (fullSeq2.ParseBaseSequence(), new List<string>{"Protein3"}, 1),
                (fullSeq3, new List<string>{"Protein4"}, 2),
                (fullSeq3.ParseBaseSequence(), new List<string>{"Protein4"}, 1)
            };

            var stoich = new PositionFrequencyAnalysis();
            stoich.ProteinGroupsOccupancyByPeptide(peptides);

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
