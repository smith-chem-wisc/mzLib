using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using MzLibUtil;

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
        public static void TestPeriodTolerantFilenameWithoutExtension(string filenameAndOrPath, string expectedResult)
        {
            string result = PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(filenameAndOrPath);
            Assert.AreEqual(expectedResult, result);
        }
    }
}
