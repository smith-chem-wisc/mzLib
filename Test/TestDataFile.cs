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
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Spectra;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public sealed class TestDataFile
    {
        private DefaultMzSpectrum _mzSpectrumA;

        private FakeMsDataFile myMsDataFile;

        [OneTimeSetUp]
        public void Setup()
        {
            Environment.CurrentDirectory = TestContext.CurrentContext.TestDirectory;

            UsefulProteomicsDatabases.Loaders.LoadElements(@"elements.dat");

            double[] mz = { 328.73795, 329.23935, 447.73849, 448.23987, 482.23792, 482.57089, 482.90393, 500.95358, 501.28732, 501.62131, 611.99377, 612.32806, 612.66187, 722.85217, 723.35345 };
            double[] intensities = { 81007096.0, 28604418.0, 78353512.0, 39291696.0, 122781408.0, 94147520.0, 44238040.0, 71198680.0, 54184096.0, 21975364.0, 44514172.0, 43061628.0, 23599424.0, 56022696.0, 41019144.0 };

            _mzSpectrumA = new DefaultMzSpectrum(mz, intensities, false);

            var peptide = new Peptide("KQEEQMETEQQNKDEGK");

            DefaultMzSpectrum MS1 = createSpectrum(peptide.GetChemicalFormula(), 300, 2000, 1);
            DefaultMzSpectrum MS2 = createMS2spectrum(peptide.Fragment(FragmentTypes.b | FragmentTypes.y, true), 100, 1500);

            MsDataScan<IMzSpectrum<MzPeak>>[] Scans = new MsDataScan<IMzSpectrum<MzPeak>>[2];
            Scans[0] = new MsDataScan<IMzSpectrum<MzPeak>>(1, MS1.newSpectrumApplyFunctionToX(b => b + 0.00001 * b + 0.00001), "spectrum 1", 1, false, Polarity.Positive, 1.0, new MzRange(300, 2000), "first spectrum", MZAnalyzerType.Unknown, 1, MS1.SumOfAllY);

            Scans[1] = new MsDataScan<IMzSpectrum<MzPeak>>(2, MS2.newSpectrumApplyFunctionToX(b => b + 0.00001 * b + 0.00002), "spectrum 2", 2, false, Polarity.Positive, 2.0, new MzRange(100, 1500), "second spectrum", MZAnalyzerType.Unknown, 1, MS2.SumOfAllY, "spectrum 1", 693.9892, 3, .3872, 693.99, 1, DissociationType.Unknown, 1, 0.32374, 693.6550);

            myMsDataFile = new FakeMsDataFile("myFakeFile", Scans);

            myMsDataFile.LoadAllScansInMemory();

            myMsDataFile.Open();
        }



        private DefaultMzSpectrum createMS2spectrum(IEnumerable<Fragment> fragments, int v1, int v2)
        {
            List<double> allMasses = new List<double>();
            List<double> allIntensities = new List<double>();
            foreach (ChemicalFormulaFragment f in fragments)
            {
                foreach (var p in createSpectrum(f.ThisChemicalFormula, v1, v2, 2))
                {
                    allMasses.Add(p.MZ);
                    allIntensities.Add(p.Intensity);
                }
            }
            var allMassesArray = allMasses.ToArray();
            var allIntensitiessArray = allIntensities.ToArray();

            Array.Sort(allMassesArray, allIntensitiessArray);
            return new DefaultMzSpectrum(allMassesArray, allIntensitiessArray, false);
        }

        private DefaultMzSpectrum createSpectrum(ChemicalFormula f, double lowerBound, double upperBound, int minCharge)
        {

            IsotopicDistribution isodist = new IsotopicDistribution(f, 0.1, 0.001);
            DefaultMzSpectrum massSpectrum1 = new DefaultMzSpectrum(isodist.Masses.ToArray(), isodist.Intensities.ToArray(), false);

            var chargeToLookAt = minCharge;
            var correctedSpectrum = massSpectrum1.newSpectrumApplyFunctionToX(s => s.ToMassToChargeRatio(chargeToLookAt));

            List<double> allMasses = new List<double>();
            List<double> allIntensitiess = new List<double>();

            while (correctedSpectrum.FirstX > lowerBound)
            {
                foreach (var thisPeak in correctedSpectrum)
                {
                    if (thisPeak.MZ > lowerBound && thisPeak.MZ < upperBound)
                    {
                        allMasses.Add(thisPeak.MZ);
                        allIntensitiess.Add(thisPeak.Intensity);
                    }
                }
                chargeToLookAt += 1;
                correctedSpectrum = massSpectrum1.newSpectrumApplyFunctionToX(s => s.ToMassToChargeRatio(chargeToLookAt));
            }

            var allMassesArray = allMasses.ToArray();
            var allIntensitiessArray = allIntensitiess.ToArray();

            Array.Sort(allMassesArray, allIntensitiessArray);

            return new DefaultMzSpectrum(allMassesArray, allIntensitiessArray, false);
        }

        [Test]
        public void SpectrumCount()
        {
            Assert.AreEqual(15, _mzSpectrumA.Count);
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

            MsDataScan<IMzSpectrum<MzPeak>> theSpectrum = new MsDataScan<IMzSpectrum<MzPeak>>(1, _mzSpectrumA, "first spectrum", 1, true, Polarity.Positive, 1, new MzRange(300, 1000), "fake scan filter", MZAnalyzerType.Unknown, 1, _mzSpectrumA.SumOfAllY);

            MsDataScan<IMzSpectrum<MzPeak>>[] theList = new MsDataScan<IMzSpectrum<MzPeak>>[1];

            theList[0] = theSpectrum;

            FakeMsDataFile thefile = new FakeMsDataFile("Somepath", theList);

            Assert.AreEqual(15, thefile.GetScan(thefile.FirstSpectrumNumber).MassSpectrum.Count);
            Assert.AreEqual(15, thefile.GetScan(thefile.FirstSpectrumNumber).MassSpectrum.Count);

            Assert.AreEqual(1, thefile.LastSpectrumNumber);
            Assert.AreEqual(1, thefile.LastSpectrumNumber);


            Assert.IsTrue(thefile.GetScan(1).isCentroid);

            foreach (var ok in thefile.GetMsScans())
                Assert.AreEqual(new MzRange(300, 1000), ok.ScanWindowRange);


            IMsDataFile<IMzSpectrum<MzPeak>> okyee = thefile;

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

            MzRange yah;
            DissociationType d;
            double yahh;
            string s;
            int ja;
            Assert.IsFalse(thefile.GetScan(1).TryGetIsolationRange(out yah));
            Assert.IsFalse(thefile.GetScan(1).TryGetDissociationType(out d));
            Assert.IsFalse(thefile.GetScan(1).TryGetIsolationMZ(out yahh));
            Assert.IsFalse(thefile.GetScan(1).TryGetIsolationWidth(out yahh));
            Assert.IsFalse(thefile.GetScan(1).TryGetPrecursorID(out s));
            Assert.IsFalse(thefile.GetScan(1).TryGetPrecursorScanNumber(out ja));
            Assert.IsFalse(thefile.GetScan(1).TryGetSelectedIonGuessChargeStateGuess(out ja));
            Assert.IsFalse(thefile.GetScan(1).TryGetSelectedIonGuessIntensity(out yahh));
            Assert.IsFalse(thefile.GetScan(1).TryGetSelectedIonGuessMZ(out yahh));
            Assert.IsFalse(thefile.GetScan(1).TryGetSelectedIonGuessMonoisotopicIntensity(out yahh));
            Assert.IsFalse(thefile.GetScan(1).TryGetSelectedIonGuessMonoisotopicMZ(out yahh));

        }

        [Test]
        public void TestAMoreRealFile()
        {
            MzRange yah;
            DissociationType d;
            double yahh;
            string s;
            int ja;

            foreach (var aa in myMsDataFile.GetScan(1).MassSpectrum)
                Console.WriteLine(aa);

            Assert.IsTrue(myMsDataFile.GetScan(2).TryGetIsolationRange(out yah));
            Assert.AreEqual(1, yah.Width);
            Assert.IsTrue(myMsDataFile.GetScan(2).TryGetDissociationType(out d));
            Assert.AreEqual(DissociationType.Unknown, d);
            Assert.IsTrue(myMsDataFile.GetScan(2).TryGetIsolationMZ(out yahh));
            Assert.AreEqual(693.99, yahh);
            Assert.IsTrue(myMsDataFile.GetScan(2).TryGetIsolationWidth(out yahh));
            Assert.AreEqual(1, yahh);
            Assert.IsTrue(myMsDataFile.GetScan(2).TryGetPrecursorID(out s));
            Assert.AreEqual("spectrum 1", s);
            Assert.IsTrue(myMsDataFile.GetScan(2).TryGetPrecursorScanNumber(out ja));
            Assert.AreEqual(1, ja);
            Assert.IsTrue(myMsDataFile.GetScan(2).TryGetSelectedIonGuessChargeStateGuess(out ja));
            Assert.AreEqual(3, ja);
            Assert.IsTrue(myMsDataFile.GetScan(2).TryGetSelectedIonGuessIntensity(out yahh));
            Assert.AreEqual(.3872, yahh);
            Assert.IsTrue(myMsDataFile.GetScan(2).TryGetSelectedIonGuessMZ(out yahh));
            Assert.AreEqual(693.9892, yahh);
            Assert.IsTrue(myMsDataFile.GetScan(2).TryGetSelectedIonGuessMonoisotopicIntensity(out yahh));
            Assert.AreEqual(0.32374, yahh);
            Assert.IsTrue(myMsDataFile.GetScan(2).TryGetSelectedIonGuessMonoisotopicMZ(out yahh));
            Assert.AreEqual(693.6550, yahh);

            Assert.AreNotEqual(0, myMsDataFile.GetScan(2).MassSpectrum.FirstX);
            Assert.AreNotEqual(0, myMsDataFile.GetScan(2).MassSpectrum.LastX);
            double hehehe1;
            myMsDataFile.GetScan(2).TryGetSelectedIonGuessMZ(out hehehe1);
            Assert.AreNotEqual(0, hehehe1);

            myMsDataFile.GetScan(2).tranformByApplyingFunctionsToSpectraAndReplacingPrecursorMZs(b => 0, 0, 0);

            Assert.AreEqual("Scan #2", myMsDataFile.GetScan(2).ToString());

            Assert.AreEqual(0, myMsDataFile.GetScan(2).MassSpectrum.FirstX);
            Assert.AreEqual(0, myMsDataFile.GetScan(2).MassSpectrum.LastX);
            double hehehe;
            myMsDataFile.GetScan(2).TryGetSelectedIonGuessMZ(out hehehe);
            Assert.AreEqual(0, hehehe);

            IEnumerable a = myMsDataFile;
            foreach (var b in a)
                Assert.IsFalse((b as IMsDataScan<IMzSpectrum<MzPeak>>).isCentroid);
            foreach (var b in myMsDataFile)
                Assert.AreEqual(Polarity.Positive, b.Polarity);

        }
    }
}