using IO.Mgf;
using NUnit.Framework;
using System.IO;


namespace Test
{
    [TestFixture]
    public sealed class TestMgf
    {

        [Test]
        public void TestLoadMgf()
        {
            try
            {
                Mgf crash = Mgf.LoadAllStaticData(Path.Combine(Directory.GetCurrentDirectory(), "ThereIsNothingHerePleaseDoNotGenerateThisFile.mgf"));
                Assert.IsTrue(false);
            }
            catch
            {
                
            }
            Mgf a = Mgf.LoadAllStaticData(Path.Combine(Directory.GetCurrentDirectory(),"Test","tester.mgf"));
            var ya = a.GetOneBasedScan(1).MassSpectrum;
            Assert.AreEqual(192, ya.Size);
            var ya2 = a.GetOneBasedScan(3).MassSpectrum;
            Assert.AreEqual(165, ya2.Size);
            var ya3 = a.GetOneBasedScan(5).MassSpectrum;
            Assert.AreEqual(551, ya3.Size);
            var ya4 = a.GetOneBasedScan(975).MassSpectrum;
            Assert.AreEqual(190, ya4.Size);
        }
    }
}
