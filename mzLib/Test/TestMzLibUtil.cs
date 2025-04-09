using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using MzLibUtil;
using Readers;
using System.Linq;
using System.IO;
using System.Text.RegularExpressions;
using System.Collections.Generic;
using System;

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

            string fullSeq = "[UniProt NonProtein:N-acetylglutamate on E]EDM[Common Variable:Oxidation on M]MELVQPSISGVAADLDK[Test Mod2: ModName2 on K]-[Test Mod: ModName on K C-Terminus]";
            string baseSeq = fullSeq.ParseBaseSequence();
            Assert.That(baseSeq.Length == 20);

            string proteinSeq = string.Join("", Enumerable.Repeat(baseSeq, 5));
            string proteinName = "TestProtein";

            int[] oneBasedPeptideIndexInProtein = { 1, 3, 5 };

            var mods = fullSeq.ParseModifications()
                             .Select(x => new KeyValuePair<int, Dictionary<string, QuantifiedModification>>(
                                 x.Key, new Dictionary<string, QuantifiedModification>
                                 { { x.Value, new QuantifiedModification(x.Value, x.Key, intensity:1) } }))
                             .ToDictionary(pair => pair.Key, pair => pair.Value);

            var peptides = new Dictionary<string, QuantifiedPeptide>();
            for (int i = 0; i < oneBasedPeptideIndexInProtein.Length; i++)
            {
                peptides[baseSeq] = new QuantifiedPeptide(baseSeq, new List<string> { fullSeq, baseSeq}, mods, oneBasedPeptideIndexInProtein[i], intensity:2);
            }
            var protein = new QuantifiedProtein(proteinName, proteinSeq, peptides);
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
