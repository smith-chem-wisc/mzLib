// Copyright 2012, 2013, 2014 Derek J. Bailey
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
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public sealed class TestDataFile
    {
        #region Private Fields

        private MzmlMzSpectrum _mzSpectrumA;

        private FakeMsDataFile myMsDataFile;

        #endregion Private Fields

        #region Public Methods

        [OneTimeSetUp]
        public void Setup()
        {
            Environment.CurrentDirectory = TestContext.CurrentContext.TestDirectory;

            UsefulProteomicsDatabases.Loaders.LoadElements(@"elements.dat");

            double[] mz = { 328.73795, 329.23935, 447.73849, 448.23987, 482.23792, 482.57089, 482.90393, 500.95358, 501.28732, 501.62131, 611.99377, 612.32806, 612.66187, 722.85217, 723.35345 };
            double[] intensities = { 81007096.0, 28604418.0, 78353512.0, 39291696.0, 122781408.0, 94147520.0, 44238040.0, 71198680.0, 54184096.0, 21975364.0, 44514172.0, 43061628.0, 23599424.0, 56022696.0, 41019144.0 };

            _mzSpectrumA = new MzmlMzSpectrum(mz, intensities, false);

            var peptide = new Peptide("KQEEQMETEQQNKDEGK");

            MzmlMzSpectrum MS1 = createSpectrum(peptide.GetChemicalFormula(), 300, 2000, 1);
            MzmlMzSpectrum MS2 = createMS2spectrum(peptide.Fragment(FragmentTypes.b | FragmentTypes.y, true), 100, 1500);

            IMzmlScan[] Scans = new IMzmlScan[2];
            Scans[0] = new MzmlScan(1, MS1, "spectrum 1", 1, false, Polarity.Positive, 1.0, new MzRange(300, 2000), "first spectrum", MZAnalyzerType.Unknown, 1, MS1.SumOfAllY);

            Scans[1] = new MzmlScanWithPrecursor(2, MS2, "spectrum 2", 2, false, Polarity.Positive, 2.0, new MzRange(100, 1500), "second spectrum", MZAnalyzerType.Unknown, 1, MS2.SumOfAllY, "spectrum 1", 693.9892, 3, .3872, 693.99, 1, DissociationType.Unknown, 1, 0.32374, 693.6550);

            myMsDataFile = new FakeMsDataFile("myFakeFile", Scans);

            myMsDataFile.LoadAllScansInMemory();

            myMsDataFile.Open();
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
            MzmlScan theSpectrum = new MzmlScan(1, _mzSpectrumA, "first spectrum", 1, true, Polarity.Positive, 1, new MzRange(300, 1000), "fake scan filter", MZAnalyzerType.Unknown, 1, _mzSpectrumA.SumOfAllY);

            MzmlScan[] theList = new MzmlScan[1];

            theList[0] = theSpectrum;

            FakeMsDataFile thefile = new FakeMsDataFile("Somepath", theList);

            var theOneBasedScan = thefile.GetOneBasedScan(1);

            Assert.AreEqual("Scan #1", theOneBasedScan.ToString());

            Assert.AreEqual(15, theOneBasedScan.MassSpectrum.Size);
            Assert.AreEqual(15, theOneBasedScan.MassSpectrum.Size);

            Assert.AreEqual(1, thefile.NumSpectra);
            Assert.AreEqual(1, thefile.NumSpectra);

            Assert.IsTrue(thefile.GetOneBasedScan(1).IsCentroid);

            foreach (var ok in thefile.GetMsScans())
            {
                Assert.AreEqual(300, ok.ScanWindowRange.Minimum, 1e-9);
                Assert.AreEqual(1000, ok.ScanWindowRange.Maximum, 1e-9);
            }

            FakeMsDataFile okyee = thefile;

            Assert.AreEqual("Somepath (UnKnown)", okyee.ToString());

            Assert.AreEqual("Somepath", okyee.FilePath);

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

            thefile.Close();
        }

        [Test]
        public void TestAMoreRealFile()
        {
            var theScan = myMsDataFile.GetOneBasedScan(2) as IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>, IMzPeak>;
            Assert.AreEqual(1, theScan.IsolationRange.Width);
            Assert.AreEqual(DissociationType.Unknown, theScan.DissociationType);
            Assert.AreEqual(693.99, theScan.IsolationMz);
            Assert.AreEqual(1, theScan.IsolationWidth);
            Assert.AreEqual("spectrum 1", theScan.PrecursorID);
            Assert.AreEqual(1, theScan.OneBasedPrecursorScanNumber);
            Assert.AreEqual(3, theScan.SelectedIonGuessChargeStateGuess.Value);
            Assert.AreEqual(.3872, theScan.SelectedIonGuessIntensity);
            Assert.AreEqual(693.9892, theScan.SelectedIonGuessMZ);
            Assert.AreEqual(0.32374, theScan.SelectedIonGuessMonoisotopicIntensity);
            Assert.AreEqual(693.6550, theScan.SelectedIonGuessMonoisotopicMZ);

            Assert.AreNotEqual(0, myMsDataFile.GetOneBasedScan(2).MassSpectrum.FirstX);
            Assert.AreNotEqual(0, myMsDataFile.GetOneBasedScan(2).MassSpectrum.LastX);

            theScan.TranformByApplyingFunctionsToSpectraAndReplacingPrecursorMZs(b => 0, 0, 0);

            Assert.AreEqual("Scan #2", myMsDataFile.GetOneBasedScan(2).ToString());

            Assert.AreEqual(0, myMsDataFile.GetOneBasedScan(2).MassSpectrum.FirstX);
            Assert.AreEqual(0, myMsDataFile.GetOneBasedScan(2).MassSpectrum.LastX);
            Assert.AreEqual(0, theScan.SelectedIonGuessMZ);

            IEnumerable a = myMsDataFile;
            foreach (var b in a)
                Assert.IsFalse((b as IMsDataScan<IMzSpectrum<IMzPeak>, IMzPeak>).IsCentroid);
            foreach (var b in myMsDataFile)
                Assert.AreEqual(Polarity.Positive, b.Polarity);
        }

        #endregion Public Methods

        #region Private Methods

        private MzmlMzSpectrum createMS2spectrum(IEnumerable<Fragment> fragments, int v1, int v2)
        {
            List<double> allMasses = new List<double>();
            List<double> allIntensities = new List<double>();
            foreach (ChemicalFormulaFragment f in fragments)
            {
                foreach (var p in createSpectrum(f.ThisChemicalFormula, v1, v2, 2))
                {
                    allMasses.Add(p.Mz);
                    allIntensities.Add(p.Intensity);
                }
            }
            var allMassesArray = allMasses.ToArray();
            var allIntensitiessArray = allIntensities.ToArray();

            Array.Sort(allMassesArray, allIntensitiessArray);
            return new MzmlMzSpectrum(allMassesArray, allIntensitiessArray, false);
        }

        private MzmlMzSpectrum createSpectrum(ChemicalFormula f, double lowerBound, double upperBound, int minCharge)
        {
            IsotopicDistribution isodist = IsotopicDistribution.GetDistribution(f, 0.1, 0.001);
            MzmlMzSpectrum massSpectrum1 = new MzmlMzSpectrum(isodist.Masses.ToArray(), isodist.Intensities.ToArray(), false);

            return massSpectrum1;
            //var chargeToLookAt = minCharge;
            //var correctedSpectrum = massSpectrum1.NewSpectrumApplyFunctionToX(s => s.ToMz(chargeToLookAt));

            //List<double> allMasses = new List<double>();
            //List<double> allIntensitiess = new List<double>();

            //while (correctedSpectrum.FirstX > lowerBound)
            //{
            //    foreach (var thisPeak in correctedSpectrum)
            //    {
            //        if (thisPeak.Mz > lowerBound && thisPeak.Mz < upperBound)
            //        {
            //            allMasses.Add(thisPeak.Mz);
            //            allIntensitiess.Add(thisPeak.Intensity);
            //        }
            //    }
            //    chargeToLookAt += 1;
            //    correctedSpectrum = massSpectrum1.NewSpectrumApplyFunctionToX(s => s.ToMz(chargeToLookAt));
            //}

            //var allMassesArray = allMasses.ToArray();
            //var allIntensitiessArray = allIntensitiess.ToArray();

            //Array.Sort(allMassesArray, allIntensitiessArray);

            //return new MzmlMzSpectrum(allMassesArray, allIntensitiessArray, false);
        }

        #endregion Private Methods
    }
}