using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using MzLibUtil;
using Readers;

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
            Assert.That(mods[2].Count == 1);
            Assert.That(mods[2].Contains("Common Variable:Oxidation on M"));
        }

        [Test]
        public static void TestParseModificationsSideChainAndCTerminusMods()
        {
            string fullSeq = "DM[Common Variable:Oxidation on M]MELVQPSISGVDLDK-[Test Mod: ModName on K C-Terminus]";
            var mods = fullSeq.ParseModifications(ignoreTerminusMod: false);
            Assert.That(mods.Count == 2);
            Assert.That(mods.ContainsKey(2));
            Assert.That(mods.ContainsKey(18));
            Assert.That(mods[2].Count == 1);
            Assert.That(mods[18].Count == 1);
            Assert.That(mods[2].Contains("Common Variable:Oxidation on M"));
            Assert.That(mods[18].Contains("Test Mod: ModName on K C-Terminus"));
        }

        [Test]
        public static void TestParseModificationsSideChainAndNTerminusMods()
        {
            // sequence with two terminal mods
            string fullSeq = "[UniProt:N-acetylglutamate on E]EEEIAALVIDNGSGMC[Common Fixed:Carbamidomethyl on C]";
            var mods = fullSeq.ParseModifications(ignoreTerminusMod: false);
            Assert.That(mods.Count == 2);
            Assert.That(mods.ContainsKey(0));
            Assert.That(mods.ContainsKey(16));
            Assert.That(mods[0].Count == 1);
            Assert.That(mods[16].Count == 1);
            Assert.That(mods[0].Contains("UniProt:N-acetylglutamate on E"));
            Assert.That(mods[16].Contains("Common Fixed:Carbamidomethyl on C"));
        }

        [Test]
        public static void TestParseModificationsTwoModsSamePosition()
        {
            // sequence with two mods on same terminus
            string fullSeq = "[UniProt:N-acetylglutamate on E]|[Common Artifact:Water Loss on E]EEEIAALVID[Metal:Calcium on D]NGSGMC";
            var mods = fullSeq.ParseModifications(ignoreTerminusMod: false);
            Assert.That(mods.Count == 2);
            Assert.That(mods.ContainsKey(0));
            Assert.That(mods.ContainsKey(10));
            Assert.That(mods[0].Count == 2);
            Assert.That(mods[10].Count == 1);
            Assert.That(mods[0].Contains("UniProt:N-acetylglutamate on E"));
            Assert.That(mods[0].Contains("Common Artifact:Water Loss on E"));
            Assert.That(mods[10].Contains("Metal:Calcium on D"));
        }

        [Test]
        public static void TestParseModificationsIgnoreTerminusMod()
        {
            // sequence with mod on both termini and mod on first amino acid side chain
            string fullSeq = "[UniProt:N-acetylglutamate on E]|[Common Artifact:Water Loss on E]E[Metal:Sodium[I] on E]EEIAALVID[Metal:Calcium[II] on D]NGSGMC[Common Fixed:Carbamidomethyl on C]";
            var mods = fullSeq.ParseModifications(ignoreTerminusMod: true);
            Assert.That(mods.Count == 3);
            Assert.That(mods.ContainsKey(1));
            Assert.That(mods.ContainsKey(10));
            Assert.That(mods.ContainsKey(16));
            Assert.That(mods[1].Count == 1);
            Assert.That(mods[10].Count == 1);
            Assert.That(mods[16].Count == 1);
            Assert.That(mods[1].Contains("Metal:Sodium[I] on E"));
            Assert.That(mods[10].Contains("Metal:Calcium[II] on D"));
            Assert.That(mods[16].Contains("Common Fixed:Carbamidomethyl on C"));
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
