﻿// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work Copyright 2016 Stefan Solntsev
//
// This file (TestDataFile.cs) is part of MassSpectrometry.Tests.
//
// MassSpectrometry.Tests is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MassSpectrometry.Tests is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with MassSpectrometry.Tests. If not, see <http://www.gnu.org/licenses/>.

using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Generic;
using System.Linq;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test
{
    [TestFixture]
    public sealed class TestMsDataFile
    {
        private MzSpectrum _mzSpectrumA;
        private FakeMsDataFile myMsDataFile;
        private static Stopwatch Stopwatch { get; set; }

        [OneTimeSetUp]
        public void Setup()
        {
            Environment.CurrentDirectory = TestContext.CurrentContext.TestDirectory;

            UsefulProteomicsDatabases.Loaders.LoadElements();

            double[] mz = { 328.73795, 329.23935, 447.73849, 448.23987, 482.23792, 482.57089, 482.90393, 500.95358, 501.28732, 501.62131, 611.99377, 612.32806, 612.66187, 722.85217, 723.35345 };
            double[] intensities = { 81007096.0, 28604418.0, 78353512.0, 39291696.0, 122781408.0, 94147520.0, 44238040.0, 71198680.0, 54184096.0, 21975364.0, 44514172.0, 43061628.0, 23599424.0, 56022696.0, 41019144.0 };

            _mzSpectrumA = new MzSpectrum(mz, intensities, false);

            var peptide = new Peptide("KQEEQMETEQQNKDEGK");

            MzSpectrum MS1 = CreateSpectrum(peptide.GetChemicalFormula(), 300, 2000, 1);
            MzSpectrum MS2 = CreateMS2spectrum(peptide.Fragment(FragmentTypes.b | FragmentTypes.y, true), 100, 1500);

            MsDataScan[] Scans = new MsDataScan[2];
            Scans[0] = new MsDataScan(MS1, 1, 1, false, Polarity.Positive, 1.0, new MzRange(300, 2000), "first spectrum", MZAnalyzerType.Unknown, MS1.SumOfAllY, 1, null, "scan=1");

            Scans[1] = new MsDataScan(MS2, 2, 2, false, Polarity.Positive, 2.0, new MzRange(100, 1500), "second spectrum", MZAnalyzerType.Unknown, MS2.SumOfAllY, 1, null, "scan=2", 693.9892, 3, .3872, 693.99, 1, DissociationType.Unknown, 1, 693.6550);

            myMsDataFile = new FakeMsDataFile(Scans);
        }

        [SetUp]
        public static void Setuppp()
        {
            Stopwatch = new Stopwatch();
            Stopwatch.Start();
        }

        [TearDown]
        public static void TearDown()
        {
            Console.WriteLine($"Analysis time: {Stopwatch.Elapsed.Hours}h {Stopwatch.Elapsed.Minutes}m {Stopwatch.Elapsed.Seconds}s");
        }

        [Test]
        public void SpectrumCount()
        {
            Assert.AreEqual(15, _mzSpectrumA.Size);
        }

        [Test]
        public void SpectrumFirstMZ()
        {
            Assert.AreEqual(328.73795, _mzSpectrumA.FirstX);
        }

        [Test]
        public void SpectrumLastMZ()
        {
            Assert.AreEqual(723.35345, _mzSpectrumA.LastX);
        }

        [Test]
        public void DataFileTest()
        {
            MsDataScan theSpectrum = new MsDataScan(_mzSpectrumA, 1, 1, true, Polarity.Positive, 1, new MzRange(300, 1000), "fake scan filter", MZAnalyzerType.Unknown, _mzSpectrumA.SumOfAllY, 1, null, "scan=1");

            MsDataScan[] theList = new MsDataScan[1];

            theList[0] = theSpectrum;

            FakeMsDataFile thefile = new FakeMsDataFile(theList);

            var theOneBasedScan = thefile.GetOneBasedScan(1);

            Assert.AreEqual("Scan #1", theOneBasedScan.ToString());

            Assert.AreEqual(15, theOneBasedScan.MassSpectrum.Size);
            Assert.AreEqual(15, theOneBasedScan.MassSpectrum.Size);

            Assert.AreEqual(1, thefile.NumSpectra);
            Assert.AreEqual(1, thefile.NumSpectra);

            Assert.IsTrue(thefile.GetOneBasedScan(1).IsCentroid);

            foreach (var ok in thefile)
            {
                Assert.AreEqual(300, ok.ScanWindowRange.Minimum, 1e-9);
                Assert.AreEqual(1000, ok.ScanWindowRange.Maximum, 1e-9);
            }

            int ok1 = 0;
            foreach (var i in thefile.GetMsScansInTimeRange(0, 2))
                ok1 += 1;

            Assert.AreEqual(1, ok1);

            int ok2 = 0;
            foreach (var i in thefile.GetMsScansInTimeRange(2, 4))
                ok2 += 1;

            Assert.AreEqual(0, ok2);

            int ok3 = 0;
            foreach (var i in thefile.GetMsScansInTimeRange(-4, -2))
                ok3 += 1;

            Assert.AreEqual(0, ok3);
        }

        [Test]
        public void TestAMoreRealFile()
        {
            var theScan = myMsDataFile.GetOneBasedScan(2);
            Assert.AreEqual(1, theScan.IsolationRange.Width);
            Assert.AreEqual(DissociationType.Unknown, theScan.DissociationType);
            Assert.AreEqual(693.99, theScan.IsolationMz);
            Assert.AreEqual(1, theScan.IsolationRange.Maximum - theScan.IsolationRange.Minimum);
            Assert.AreEqual(1, theScan.OneBasedPrecursorScanNumber);
            Assert.AreEqual(3, theScan.SelectedIonChargeStateGuess.Value);

            var precursorScan = myMsDataFile.GetOneBasedScan(theScan.OneBasedPrecursorScanNumber.Value);
            theScan.RefineSelectedMzAndIntensity(precursorScan.MassSpectrum);
            Assert.AreEqual(.32872, theScan.SelectedIonIntensity, 0.01);
            Assert.AreEqual(693.9892, theScan.SelectedIonMZ, 0.01);
            Assert.AreEqual(693.655, theScan.SelectedIonMonoisotopicGuessMz, 0.001);

            Assert.AreNotEqual(0, myMsDataFile.GetOneBasedScan(2).MassSpectrum.FirstX);

            Assert.AreEqual(myMsDataFile.GetOneBasedScan(2).MassSpectrum.CopyTo2DArray()[0, 0], myMsDataFile.GetOneBasedScan(2).MassSpectrum.FirstX);

            Assert.AreNotEqual(0, myMsDataFile.GetOneBasedScan(2).MassSpectrum.LastX);

            theScan.ComputeMonoisotopicPeakIntensity(precursorScan.MassSpectrum);

            theScan.TransformMzs(b => 0, b => 0);

            Assert.AreEqual("Scan #2", myMsDataFile.GetOneBasedScan(2).ToString());

            Assert.AreEqual(0, myMsDataFile.GetOneBasedScan(2).MassSpectrum.FirstX);
            Assert.AreEqual(0, myMsDataFile.GetOneBasedScan(2).MassSpectrum.LastX);
            Assert.AreEqual(0, theScan.SelectedIonMZ);

            List<MsDataScan> a = myMsDataFile.GetAllScansList();
            foreach (var b in a)
                Assert.IsFalse(b.IsCentroid);
            foreach (var b in myMsDataFile)
                Assert.AreEqual(Polarity.Positive, b.Polarity);
        }

        private MzSpectrum CreateMS2spectrum(IEnumerable<Fragment> fragments, int v1, int v2)
        {
            List<double> allMasses = new List<double>();
            List<double> allIntensities = new List<double>();
            foreach (ChemicalFormulaFragment f in fragments)
            {
                var spec = CreateSpectrum(f.ThisChemicalFormula, v1, v2, 2);
                for (int i = 0; i < spec.Size; i++)
                {
                    allMasses.Add(spec.XArray[i]);
                    allIntensities.Add(spec.YArray[i]);
                }
            }
            var allMassesArray = allMasses.ToArray();
            var allIntensitiessArray = allIntensities.ToArray();

            Array.Sort(allMassesArray, allIntensitiessArray);
            return new MzSpectrum(allMassesArray, allIntensitiessArray, false);
        }

        private MzSpectrum CreateSpectrum(ChemicalFormula f, double lowerBound, double upperBound, int minCharge)
        {
            IsotopicDistribution isodist = IsotopicDistribution.GetDistribution(f, 0.1, 0.001);
            MzSpectrum notActuallyMzS = new MzSpectrum(isodist.Masses.ToArray(), isodist.Intensities.ToArray(), false);

            notActuallyMzS.ReplaceXbyApplyingFunction(b => b.Mz.ToMz(1));

            List<double> allMasses = new List<double>();
            List<double> allIntensitiess = new List<double>();

            while (notActuallyMzS.FirstX > lowerBound)
            {
                for (int i = 0; i < notActuallyMzS.Size; i++)
                {
                    if (notActuallyMzS.XArray[i] > lowerBound && notActuallyMzS.XArray[i] < upperBound)
                    {
                        allMasses.Add(notActuallyMzS.XArray[i]);
                        allIntensitiess.Add(notActuallyMzS.YArray[i]);
                    }
                }
                minCharge += 1;
                notActuallyMzS = new MzSpectrum(isodist.Masses.ToArray(), isodist.Intensities.ToArray(), false);
                notActuallyMzS.ReplaceXbyApplyingFunction(s => s.Mz.ToMz(minCharge));
            }

            var allMassesArray = allMasses.ToArray();
            var allIntensitiessArray = allIntensitiess.ToArray();

            Array.Sort(allMassesArray, allIntensitiessArray);

            return new MzSpectrum(allMassesArray, allIntensitiessArray, false);
        }
    }
}