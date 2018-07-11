using IO.Mgf;
using NUnit.Framework;
using System.IO;


namespace Test
{
    [TestFixture]
    public sealed class TestMgf
    {

        [Test]
        public static void TestLoadMgf()
        {
            try
            {
                Mgf.LoadAllStaticData(Path.Combine(TestContext.CurrentContext.TestDirectory, "ThereIsNothingHerePleaseDoNotGenerateThisFile.mgf"));
                Assert.IsTrue(false);
            }
            catch
            {
                //woohoo, there was an exception!
            }
            Mgf a = Mgf.LoadAllStaticData(Path.Combine(TestContext.CurrentContext.TestDirectory, "tester.mgf"));
            var ya = a.GetOneBasedScan(14);
            Assert.AreEqual(192, ya.MassSpectrum.Size);
            Assert.AreEqual(2, ya.MsnOrder);
            Assert.AreEqual(14, ya.OneBasedScanNumber);
            Assert.AreEqual(MassSpectrometry.Polarity.Positive, ya.Polarity);
            Assert.AreEqual(0.26666666666666666, ya.RetentionTime);
            Assert.AreEqual(571.806916, ya.IsolationMz);
            Assert.AreEqual(571.806916, ya.SelectedIonMZ);
            Assert.AreEqual(2, ya.SelectedIonChargeStateGuess);
            Assert.AreEqual(571.806916, ya.SelectedIonMonoisotopicGuessMz);
            Assert.AreEqual(1294963.5999999996, ya.TotalIonCurrent);
            Assert.AreEqual(110.0719, ya.ScanWindowRange.Minimum);
            Assert.AreEqual(1038.8018, ya.ScanWindowRange.Maximum);
            var ya2 = a.GetOneBasedScan(20).MassSpectrum;
            Assert.AreEqual(165, ya2.Size);
            var ya3 = a.GetOneBasedScan(2).MassSpectrum;
            Assert.AreEqual(551, ya3.Size);
        }
    }
}
