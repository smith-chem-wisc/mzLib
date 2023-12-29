using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestMzSpectrum
    {

        [Test]
        public void Bubba()
        {
            double[] mzs = new double[3] { 1,1.0000001, 1.0000002 };
            double[] intensities = new double[3] { 1,2,3 };
            bool shouldCopy = false;
            MzSpectrum s = new MzSpectrum(mzs,intensities,shouldCopy);

            (mzs,intensities) = s.Yoyo(mzs.ToList(), intensities.ToList(), new PpmTolerance(100));
            Assert.IsTrue(false);
        }
    }
}
