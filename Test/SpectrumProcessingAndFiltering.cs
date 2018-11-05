using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MassSpectrometry;
using NUnit.Framework;

namespace Test
{
    [TestFixture]
    public sealed class SpectrumProcessingAndFiltering
    {
        [Test]
        [TestCase(10, null, 100)]
        [TestCase(10, 1, 100)]
        [TestCase(10, 5, 100)]
        public static void TestFilteringPeaksTopN_MultipleWindows(int peaksToKeep, int? windows, int peakCount)
        {
            double[] mzArray = new double[100];
            double[] intArray = new double[100];

            Random r = new Random();  

            for (int i = 0; i < 100; i++)
            {
                double randomMz = r.Next(100, 2000);
                double randomInst = r.Next(0, 100);

                mzArray[i] = randomMz;
                intArray[i] = randomInst;
            }

            FilteringParams f = new FilteringParams(peaksToKeep, null, windows, false, false);

            MsDataFile.WindowModeHelper(ref intArray, ref mzArray,f,100,2000);


            if (windows.HasValue)
            {
                Assert.LessOrEqual((decimal)mzArray.Count(), (decimal)peaksToKeep * (decimal)windows);
            }
            else
            {
                Assert.LessOrEqual((decimal)mzArray.Count(), (decimal)peaksToKeep * (decimal)1.0);
            }
        }


    }
}
