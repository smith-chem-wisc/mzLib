using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chromatography;

using NUnit.Framework;

namespace Test
{
    [TestFixture]
    public sealed class TestCZE
    {
        #region Public Methods


        [Test]
        public static void TestElectrophoreticMobilityPredictions()
        {
            CZE testCZE = new CZE();

            double expMu = Math.Round(testCZE.ExperimentalElectrophoreticMobility(1, 1, 1),0);

            Assert.AreEqual(expMu, 16666667);



            double predMu = Math.Round(testCZE.PredictedElectrophoreticMobility("ATPATEESTVPATQSSALPAAK", 2127.0708), 0);

            Assert.AreEqual(predMu, 12);

        }

        #endregion Public Methods
    }
}
