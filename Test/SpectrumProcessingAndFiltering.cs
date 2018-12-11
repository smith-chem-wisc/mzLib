﻿using IO.MzML;
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
        [TestCase(10, null, 100, null)]
        [TestCase(10, 1, 100, null)]
        [TestCase(10, 5, 100, null)]
        [TestCase(10, 5, 100, 500)]
        public static void TestFilteringPeaksTopN_MultipleWindows(int peaksToKeep, int? windows, int peakCount, double? normalizationMax)
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

            MsDataFile.WindowModeHelper(ref intArray, ref mzArray, f, 100, 2000, normalizationMax);

            if (windows.HasValue)
            {
                Assert.LessOrEqual((decimal)mzArray.Count(), (decimal)peaksToKeep * (decimal)windows);
            }
            else
            {
                Assert.LessOrEqual((decimal)mzArray.Count(), (decimal)peaksToKeep * (decimal)1.0);
            }
            if (normalizationMax.HasValue)
            {
                Assert.AreEqual(normalizationMax.Value, intArray.Max());
            }
        }

        [Test]
        public static void TestFilterLowIntensity()
        {
            double[] mzArray = new double[100];
            double[] intArray = new double[100];

            for (int i = 0; i < 100; i++)
            {
                mzArray[i] = (double)i;
                intArray[i] = (double)i;
            }

            FilteringParams f = new FilteringParams(null, 0.05, null, false, false);

            MsDataFile.WindowModeHelper(ref intArray, ref mzArray, f, mzArray.Min(), mzArray.Max());

            //The first five intensities are below 5% and therefore removed.
            Assert.AreEqual(95, intArray.Count());
            Assert.AreEqual(95, mzArray.Count());
        }

        [Test]
        public static void TestFilterNumberPeaksPerWindowOneWindow()
        {
            double[] mzArray = new double[200];
            double[] intArray = new double[200];

            for (int i = 0; i < 200; i++)
            {
                mzArray[i] = (double)(i + 1);
                intArray[i] = (double)(i + 1);
            }

            FilteringParams f = new FilteringParams(100, null, null, false, false);

            MsDataFile.WindowModeHelper(ref intArray, ref mzArray, f, mzArray.Min(), mzArray.Max());

            Assert.AreEqual(100, intArray.Count());
            Assert.AreEqual(100, mzArray.Count());
            Assert.AreEqual(101, intArray.Min());
            Assert.AreEqual(200, intArray.Max());
            Assert.AreEqual(101, mzArray.Min());
            Assert.AreEqual(200, mzArray.Max());
        }

        [Test]
        public static void TestFilterNumberPeaksPerWindowTenWindows()
        {
            double[] mzArray = new double[200];
            double[] intArray = new double[200];

            for (int i = 0; i < 200; i++)
            {
                mzArray[i] = (double)(i + 1);
                intArray[i] = (double)(i + 1);
            }

            FilteringParams f = new FilteringParams(10, null, 10, false, false);

            MsDataFile.WindowModeHelper(ref intArray, ref mzArray, f, mzArray.Min(), mzArray.Max());

            Assert.AreEqual(100, intArray.Count());
            Assert.AreEqual(100, mzArray.Count());
            Assert.AreEqual(11, intArray.Min());
            Assert.AreEqual(200, intArray.Max());
            Assert.AreEqual(11, mzArray.Min());
            Assert.AreEqual(200, mzArray.Max());
        }

        [Test]
        public static void TestXcorrFilteringPeaksTopN_MultipleWindows()
        {
            List<double> masses = new List<double>();
            List<double> intensities = new List<double>();

            int startMass = 1;
            int incrementMass = 5;

            double startIntensity = 0.5;
            double incrementIntensity = 5;

            Random r = new Random();

            while (startMass < 1969)
            {
                double smallDouble = (double)r.Next(0, 10000) / 100000d;
                masses.Add(startMass * 1.0005079 + smallDouble);
                intensities.Add(startIntensity);
                startMass = startMass + incrementMass;
                startIntensity = startIntensity + incrementIntensity * startMass;
            }

            double[] mzArray = masses.ToArray();
            double[] intArray = intensities.ToArray();
            Array.Sort(mzArray, intArray);

            var spectrum = new MzSpectrum(mzArray, intArray, false);
            spectrum.XCorrPrePreprocessing(mzArray.Min(), mzArray.Max(), 241.122);
  

            //first mz rounded to nearest discrete mass bin 1.0005079
            Assert.AreEqual(Math.Round(96.0487584, 5), Math.Round(spectrum.XArray.Min(), 5));

            //last mz rounded to nearest discrete mass bin 1.0005079
            Assert.AreEqual(Math.Round(1966.998531, 5), Math.Round(spectrum.XArray.Max(), 5));

            //peaks within 1.5 thomson of precursor 241.122 are absent
            double precursorIntensity = 0;
            for (int i = 0; i < spectrum.XArray.Length; i++)
            {
                if(spectrum.XArray[i] > (241.122 - 1.5) && spectrum.XArray[i] < (241.122 + 1.5))
                {
                    precursorIntensity += spectrum.YArray[i];
                }
            }
            Assert.AreEqual(0, precursorIntensity);

            //not zero intensities
            Assert.AreEqual(374, spectrum.YArray.Where(i => i > 0).ToList().Count);

            //zero intensities. there shouldn't be any.
            Assert.AreEqual(0, spectrum.YArray.Where(i => i == 0).ToList().Count);

            //first peak with intensity
            Assert.AreEqual(Math.Round(21.170981, 5), Math.Round(spectrum.YArray[0], 5));

            //last peak with intensity
            Assert.AreEqual(Math.Round(39.674212, 5), Math.Round(spectrum.YArray[373], 5));
        }

        [Test]
        public static void ProcessXcorrInMzSpectrum()
        {
            Dictionary<string, MsDataFile> MyMsDataFiles = new Dictionary<string, MsDataFile>();
            string origDataFile = System.IO.Path.Combine(TestContext.CurrentContext.TestDirectory, "BinGenerationTest.mzML");
            FilteringParams filter = new FilteringParams(200, 0.01, 1, false, true);

            MyMsDataFiles[origDataFile] = Mzml.LoadAllStaticData(origDataFile, filter, 1);

            var scans = MyMsDataFiles[origDataFile].GetAllScansList();

            foreach (MsDataScan scan in scans.Where(s => s.MsnOrder > 1))
            {
                scan.MassSpectrum.XCorrPrePreprocessing(0, 1968 * 1.0005079, scan.IsolationMz.Value);
            }

            Assert.AreEqual(6, scans[0].MassSpectrum.XArray.Count());
            Assert.AreEqual(20, scans[1].MassSpectrum.XArray.Count());
        }
    }
}