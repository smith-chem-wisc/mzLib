using Chromatography;
using NUnit.Framework;
using System;

namespace Test
{
    [TestFixture]
    public sealed class TestCZE
    {
        [Test]
        public static void TestElectrophoreticMobilityPredictions()
        {
            CZE testCZE = new CZE();

            double expMu = Math.Round(testCZE.ExperimentalElectrophoreticMobility(1, 1, 1), 0);
            Assert.AreEqual(expMu, 16666667);

            double predMu = Math.Round(testCZE.PredictedElectrophoreticMobility("ATPATEESTVPATQSSALPAAK", 2127.0708), 0);
            Assert.AreEqual(predMu, 12);
        }
    }
}