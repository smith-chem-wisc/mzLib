using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using MzLibUtil;
using NUnit.Framework;
using IO.Thermo;

namespace Test
{
    [TestFixture]
    public sealed class TestDecon
    {
        string filePath = @"C:\Data\Decon\deconTest.raw";

        [Test]
        public void TestDeconvolution()
        {
            var file = ThermoStaticData.LoadAllStaticData(filePath);

            var ms1Scan = file.Where(p => p.MsnOrder == 1).First();
            var ms2Scan = file.Where(p => p.MsnOrder == 2).First();

            //ms1Scan.MassSpectrum.Deconvolute(ms1Scan.ScanWindowRange, 30, , 5);
        }
    }
}
