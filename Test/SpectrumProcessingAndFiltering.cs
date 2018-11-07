using MassSpectrometry;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;

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

            MsDataFile.WindowModeHelper(ref intArray, ref mzArray, f, 100, 2000);

            if (windows.HasValue)
            {
                Assert.LessOrEqual((decimal)mzArray.Count(), (decimal)peaksToKeep * (decimal)windows);
            }
            else
            {
                Assert.LessOrEqual((decimal)mzArray.Count(), (decimal)peaksToKeep * (decimal)1.0);
            }
        }

        [Test]
        public static void TestXcorrFilteringPeaksTopN_MultipleWindows()
        {
            List<double> masses = new List<double>();
            List<double> intensities = new List<double>();

            int startMass = 1;
            int incrementMass = 2;

            double startIntensity = 0.25;
            double incrementIntensity = 0.25;

            for (int i = 0; i < 100; i++)
            {
                masses.Add(startMass);
                intensities.Add(startIntensity);
                startMass = startMass + incrementMass;
                startIntensity = startIntensity + incrementIntensity * startMass;
            }

            double[] mzArray = masses.ToArray();
            double[] intArray = intensities.ToArray();

            MsDataFile.XCorrPrePreprocessing(ref intArray, ref mzArray, 1, 200, 197);

            //first mz rounded to nearest discrete mass bin 1.0005079
            Assert.AreEqual(1.0005079, mzArray[0]);

            //last mz rounded to nearest discrete mass bin 1.0005079
            Assert.AreEqual(Math.Round(200.10158,5), Math.Round(mzArray[199],5));
            
            //peaks within 1.5 thomson of precursor are absent
            Assert.AreEqual(0, intArray[195] + intArray[196] + intArray[197]);

            //not zero intensities
            Assert.AreEqual(94, intArray.Where(i => i > 0).ToList().Count);

            //zero intensities
            Assert.AreEqual(106, intArray.Where(i => i == 0).ToList().Count);

            //first peak with intensity
            Assert.AreEqual(Math.Round(11.60174419,5), Math.Round(intArray[10],5));

            //middle peak with intensity
            Assert.AreEqual(Math.Round(19.98981657,5), Math.Round(intArray[100],5));

            //middle peak with 0 intensity
            Assert.AreEqual(0, intArray[101]);

            //last peak with intensity
            Assert.AreEqual(Math.Round(27.20941558,5), Math.Round(intArray[198],5));
        }
    }
}