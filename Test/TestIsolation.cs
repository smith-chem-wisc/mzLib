﻿// Copyright 2017 Stefan Solntsev
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
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics.AminoAcidPolymer;
using System;
using System.Linq;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test
{
    [TestFixture]
    public sealed class TestIsolation
    {
        private static Stopwatch Stopwatch { get; set; }

        [OneTimeSetUp]
        public void Setup()
        {
            Environment.CurrentDirectory = TestContext.CurrentContext.TestDirectory;

            UsefulProteomicsDatabases.Loaders.LoadElements();
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
        public void TestCoIsolation()
        {
            Peptide pep1 = new Peptide("AAAAAA");
            Peptide pep2 = new Peptide("AAA[H]AAA");

            var dist1 = IsotopicDistribution.GetDistribution(pep1.GetChemicalFormula(), 0.1, 0.01);

            var dist2 = IsotopicDistribution.GetDistribution(pep2.GetChemicalFormula(), 0.1, 0.01);

            MsDataScan[] Scans = new MsDataScan[2];
            double[] ms1intensities = new double[] { 0.8, 0.8, 0.2, 0.02, 0.2, 0.02 };
            double[] ms1mzs = dist1.Masses.Concat(dist2.Masses).OrderBy(b => b).Select(b => b.ToMz(1)).ToArray();

            double selectedIonMz = ms1mzs[1];

            MzSpectrum MS1 = new MzSpectrum(ms1mzs, ms1intensities, false);

            Scans[0] = new MsDataScan(MS1, 1, 1, false, Polarity.Positive, 1.0, new MzRange(300, 2000), "first spectrum", MZAnalyzerType.Unknown, MS1.SumOfAllY, null, null, null);

            // Horrible fragmentation, but we don't care about this!
            double[] ms2intensities = new double[] { 1000 };
            double[] ms2mzs = new double[] { 1000 };
            MzSpectrum MS2 = new MzSpectrum(ms2mzs, ms2intensities, false);
            double isolationMZ = selectedIonMz;
            Scans[1] = new MsDataScan(MS2, 2, 2, false, Polarity.Positive, 2.0, new MzRange(100, 1500), "second spectrum", MZAnalyzerType.Unknown, MS2.SumOfAllY, null, null, null, selectedIonMz, null, null, isolationMZ, 2.5, DissociationType.HCD, 1, null);

            var myMsDataFile = new FakeMsDataFile(Scans);

            var cool = myMsDataFile.GetAllScansList().Last();

            int maxAssumedChargeState = 1;
            Tolerance massTolerance = Tolerance.ParseToleranceString("10 PPM");

            var isolatedMasses = cool.GetIsolatedMassesAndCharges(myMsDataFile.GetOneBasedScan(cool.OneBasedPrecursorScanNumber.Value).MassSpectrum, 1, maxAssumedChargeState, 10, 5).ToList();

            Assert.AreEqual(2, isolatedMasses.Count);
            Assert.AreEqual(2, isolatedMasses.Count(b => b.charge == 1));
            Assert.AreEqual(pep1.MonoisotopicMass, isolatedMasses.Select(b => b.peaks.First().Item1.ToMass(b.charge)).Min(), 1e-9);
            Assert.AreEqual(pep2.MonoisotopicMass, isolatedMasses.Select(b => b.peaks.First().Item1.ToMass(b.charge)).Max(), 1e-9);
            Assert.AreEqual(pep1.MonoisotopicMass, isolatedMasses.Select(b => b.monoisotopicMass.ToMz(b.charge).ToMass(b.charge)).Min(), 1e-9);
            Assert.AreEqual(pep2.MonoisotopicMass, isolatedMasses.Select(b => b.monoisotopicMass.ToMz(b.charge).ToMass(b.charge)).Max(), 1e-9);
        }

        [Test]
        public void TestCoIsolationDifferentCharges()
        {
            Peptide pep1 = new Peptide("AAAAAA");
            Peptide pep2 = new Peptide("AAAAAAA[H21]AAAAA");

            var dist1 = IsotopicDistribution.GetDistribution(pep1.GetChemicalFormula(), 0.1, 0.01);

            var dist2 = IsotopicDistribution.GetDistribution(pep2.GetChemicalFormula(), 0.1, 0.01);

            MsDataScan[] Scans = new MsDataScan[2];
            double[] ms1intensities = new double[] { 0.8, 0.8, 0.2, 0.02, 0.2, 0.02 };
            double[] ms1mzs = dist1.Masses.Select(b => b.ToMz(1)).Concat(dist2.Masses.Select(b => b.ToMz(2))).OrderBy(b => b).ToArray();

            double selectedIonMz = ms1mzs[1];

            MzSpectrum MS1 = new MzSpectrum(ms1mzs, ms1intensities, false);

            Scans[0] = new MsDataScan(MS1, 1, 1, false, Polarity.Positive, 1.0, new MzRange(300, 2000), "first spectrum", MZAnalyzerType.Unknown, MS1.SumOfAllY, null, null, null);

            // Horrible fragmentation, but we don't care about this!
            double[] ms2intensities = new double[] { 1000 };
            double[] ms2mzs = new double[] { 1000 };
            MzSpectrum MS2 = new MzSpectrum(ms2mzs, ms2intensities, false);
            double isolationMZ = selectedIonMz;

            Scans[1] = new MsDataScan(MS2, 2, 2, false, Polarity.Positive, 2.0, new MzRange(100, 1500), "second spectrum", MZAnalyzerType.Unknown, MS2.SumOfAllY, null, null, null, selectedIonMz, null, null, isolationMZ, 2.5, DissociationType.HCD, 1, null);

            var myMsDataFile = new FakeMsDataFile(Scans);

            var cool = myMsDataFile.GetAllScansList().Last();

            int maxAssumedChargeState = 2;
            Tolerance massTolerance = Tolerance.ParseToleranceString("10 PPM");

            var isolatedMasses = cool.GetIsolatedMassesAndCharges(myMsDataFile.GetOneBasedScan(cool.OneBasedPrecursorScanNumber.Value).MassSpectrum, 1, maxAssumedChargeState, 10, 5).ToList();

            Assert.AreEqual(2, isolatedMasses.Count);
            Assert.AreEqual(1, isolatedMasses.Count(b => b.charge == 1));
            Assert.AreEqual(1, isolatedMasses.Count(b => b.charge == 2));
            Assert.AreEqual(pep1.MonoisotopicMass, isolatedMasses.Select(b => b.peaks.First().Item1.ToMass(b.charge)).Min(), 1e-9);
            Assert.AreEqual(pep2.MonoisotopicMass, isolatedMasses.Select(b => b.peaks.First().Item1.ToMass(b.charge)).Max(), 1e-9);
        }
    }
}