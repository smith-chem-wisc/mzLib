// Copyright 2017 Stefan Solntsev
//
// This file (TestIsolation.cs) is part of Tests.
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
// License along with Tests. If not, see <http://www.gnu.org/licenses/>.

using Chemistry;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using System;
using System.Linq;

namespace Test
{
    [TestFixture]
    public sealed class TestIsolation
    {

        #region Public Methods

        [OneTimeSetUp]
        public void Setup()
        {
            Environment.CurrentDirectory = TestContext.CurrentContext.TestDirectory;

            UsefulProteomicsDatabases.Loaders.LoadElements(@"elements.dat");
        }

        [Test]
        public void TestCoIsolation()
        {
            Peptide pep1 = new Peptide("AAAAAA");
            Peptide pep2 = new Peptide("AAA[H]AAA");

            var dist1 = IsotopicDistribution.GetDistribution(pep1.GetChemicalFormula(), 0.1, 0.01);

            var dist2 = IsotopicDistribution.GetDistribution(pep2.GetChemicalFormula(), 0.1, 0.01);

            IMzmlScan[] Scans = new IMzmlScan[2];
            double[] ms1intensities = new double[] { 1, 1, 1, 1, 1, 1 };
            double[] ms1mzs = dist1.Masses.Concat(dist2.Masses).OrderBy(b => b).Select(b => b.ToMz(1)).ToArray();

            double selectedIonMz = ms1mzs[1];

            MzmlMzSpectrum MS1 = new MzmlMzSpectrum(ms1mzs, ms1intensities, false);

            Scans[0] = new MzmlScan(1, MS1, 1, false, Polarity.Positive, 1.0, new MzRange(300, 2000), "first spectrum", MZAnalyzerType.Unknown, MS1.SumOfAllY, null);

            // Horrible fragmentation, but we don't care about this!
            double[] ms2intensities = new double[] { 1000 };
            double[] ms2mzs = new double[] { 1000 };
            MzmlMzSpectrum MS2 = new MzmlMzSpectrum(ms2mzs, ms2intensities, false);
            double isolationMZ = selectedIonMz;
            Scans[1] = new MzmlScanWithPrecursor(2, MS2, 2, false, Polarity.Positive, 2.0, new MzRange(100, 1500), "second spectrum", MZAnalyzerType.Unknown, MS2.SumOfAllY, selectedIonMz, null, null, isolationMZ, 2.5, DissociationType.HCD, 1, null, null);

            var myMsDataFile = new FakeMsDataFile(Scans);

            var cool = myMsDataFile.Last() as IMsDataScanWithPrecursor<MzmlMzSpectrum>;

            int maxAssumedChargeState = 1;
            Tolerance massTolerance = new Tolerance("5 PPM");

            var isolatedMasses = cool.GetIsolatedMassesAndCharges(myMsDataFile.GetOneBasedScan(cool.OneBasedPrecursorScanNumber).MassSpectrum, maxAssumedChargeState, massTolerance);

            Assert.AreEqual(2, isolatedMasses.Count);
        }

        [Test]
        public void TestCoIsolationDifferentCharges()
        {
            Peptide pep1 = new Peptide("AAAAAA");
            Peptide pep2 = new Peptide("AAAAAAA[H21]AAAAA");

            var dist1 = IsotopicDistribution.GetDistribution(pep1.GetChemicalFormula(), 0.1, 0.01);

            var dist2 = IsotopicDistribution.GetDistribution(pep2.GetChemicalFormula(), 0.1, 0.01);

            IMzmlScan[] Scans = new IMzmlScan[2];
            double[] ms1intensities = new double[] { 1, 1, 1, 1, 1, 1 };
            double[] ms1mzs = dist1.Masses.Select(b => b.ToMz(1)).Concat(dist2.Masses.Select(b => b.ToMz(2))).OrderBy(b => b).ToArray();

            double selectedIonMz = ms1mzs[1];

            MzmlMzSpectrum MS1 = new MzmlMzSpectrum(ms1mzs, ms1intensities, false);

            Scans[0] = new MzmlScan(1, MS1, 1, false, Polarity.Positive, 1.0, new MzRange(300, 2000), "first spectrum", MZAnalyzerType.Unknown, MS1.SumOfAllY, null);

            // Horrible fragmentation, but we don't care about this!
            double[] ms2intensities = new double[] { 1000 };
            double[] ms2mzs = new double[] { 1000 };
            MzmlMzSpectrum MS2 = new MzmlMzSpectrum(ms2mzs, ms2intensities, false);
            double isolationMZ = selectedIonMz;
            Scans[1] = new MzmlScanWithPrecursor(2, MS2, 2, false, Polarity.Positive, 2.0, new MzRange(100, 1500), "second spectrum", MZAnalyzerType.Unknown, MS2.SumOfAllY, selectedIonMz, null, null, isolationMZ, 2.5, DissociationType.HCD, 1, null, null);

            var myMsDataFile = new FakeMsDataFile(Scans);

            var cool = myMsDataFile.Last() as IMsDataScanWithPrecursor<MzmlMzSpectrum>;

            int maxAssumedChargeState = 2;
            Tolerance massTolerance = new Tolerance("5 PPM");

            var isolatedMasses = cool.GetIsolatedMassesAndCharges(myMsDataFile.GetOneBasedScan(cool.OneBasedPrecursorScanNumber).MassSpectrum, maxAssumedChargeState, massTolerance);

            Assert.AreEqual(2, isolatedMasses.Count);
        }

        #endregion Public Methods

    }
}