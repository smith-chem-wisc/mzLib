// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work Copyright 2016 Stefan Solntsev
//
// This file (TestSpectra.cs) is part of MassSpectrometry.Tests.
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

using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestSpectra
    {
        private MzSpectrum _mzSpectrumA;
        private static Stopwatch Stopwatch { get; set; }

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

        [SetUp]
        public void Setup()
        {
            double[] mz = { 328.73795, 329.23935, 447.73849, 448.23987, 482.23792, 482.57089, 482.90393, 500.95358, 501.28732, 501.62131, 611.99377, 612.32806, 612.66187, 722.85217, 723.35345 };
            double[] intensities = { 81007096.0, 28604418.0, 78353512.0, 39291696.0, 122781408.0, 94147520.0, 44238040.0, 71198680.0, 54184096.0, 21975364.0, 44514172.0, 43061628.0, 23599424.0, 56022696.0, 41019144.0 };

            _mzSpectrumA = new MzSpectrum(mz, intensities, false);
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
        public void SpectrumBasePeakIntensity()
        {
            double basePeakIntensity = _mzSpectrumA.YofPeakWithHighestY.Value;

            Assert.AreEqual(122781408.0, basePeakIntensity);
        }

        [Test]
        public void SpectrumTIC()
        {
            double tic = _mzSpectrumA.SumOfAllY;

            Assert.AreEqual(843998894.0, tic);
        }

        [Test]
        public void SpectrumGetIntensityFirst()
        {
            Assert.AreEqual(81007096.0, _mzSpectrumA.YArray[0]);
        }

        [Test]
        public void SpectrumGetIntensityRandom()
        {
            Assert.AreEqual(44238040.0, _mzSpectrumA.YArray[6]);
        }

        [Test]
        public void SpectrumGetMassFirst()
        {
            Assert.AreEqual(328.73795, _mzSpectrumA.FirstX);
        }

        [Test]
        public void SpectrumGetMassRandom()
        {
            Assert.AreEqual(482.90393, _mzSpectrumA.XArray[6]);
        }

        [Test]
        public void SpectrumContainsPeak()
        {
            Assert.IsTrue(_mzSpectrumA.Size > 0);
        }

        [Test]
        public void SpectrumContainsPeakInRange()
        {
            Assert.AreEqual(1, _mzSpectrumA.NumPeaksWithinRange(448.23987 - 0.001, 448.23987 + 0.001));
        }

        [Test]
        public void SpectrumContainsPeakInRangeEnd()
        {
            Assert.AreEqual(0, _mzSpectrumA.NumPeaksWithinRange(448.23987 - 0.001, 448.23987));
        }

        [Test]
        public void SpectrumContainsPeakInRangeStart()
        {
            Assert.AreEqual(1, _mzSpectrumA.NumPeaksWithinRange(448.23987, 448.23987 + 0.001));
        }

        [Test]
        public void SpectrumContainsPeakInRangeStartEnd()
        {
            Assert.AreEqual(0, _mzSpectrumA.NumPeaksWithinRange(448.23987, 448.23987));
        }

        [Test]
        public void SpectrumDoesntContainPeakInRange()
        {
            Assert.AreEqual(0, _mzSpectrumA.NumPeaksWithinRange(603.4243 - 0.001, 603.4243 + 0.001));
        }

        [Test]
        public void SpectrumMassRange()
        {
            MzRange range = new MzRange(328.73795, 723.35345);

            Assert.AreEqual(0, _mzSpectrumA.Range.Minimum - range.Minimum, 1e-9);
            Assert.AreEqual(0, _mzSpectrumA.Range.Maximum - range.Maximum, 1e-9);
        }

        [Test]
        public void SpectrumFilterCount()
        {
            var filteredMzSpectrum = _mzSpectrumA.FilterByY(28604417, 28604419);

            Assert.AreEqual(1, filteredMzSpectrum.Count());
        }

        [Test]
        public void FilterByNumberOfMostIntenseTest()
        {
            Assert.AreEqual(5, _mzSpectrumA.FilterByNumberOfMostIntense(5).Count());
        }

        [Test]
        public void FilterByNumberOfMostIntenseRobTest()
        {
            double[] x = new double[] { 50, 60, 70, 147.0764, 257.1244, 258.127, 275.135 };
            double[] y = new double[] { 1, 1, 1, 1, 1, 1, 1 };
            MzSpectrum spectrum = new MzSpectrum(x, y, false);
            Assert.AreEqual(7, spectrum.FilterByNumberOfMostIntense(200).Count());
        }

        [Test]
        public void GetBasePeak()
        {
            Assert.AreEqual(122781408.0, _mzSpectrumA.YofPeakWithHighestY);
        }

        [Test]
        public void GetClosestPeak()
        {
            Assert.AreEqual(448.23987, _mzSpectrumA.GetClosestPeakXvalue(448));
            Assert.AreEqual(447.73849, _mzSpectrumA.GetClosestPeakXvalue(447.9));
        }

        [Test]
        public void Extract()
        {
            Assert.AreEqual(3, _mzSpectrumA.Extract(500, 600).Count());
        }

        [Test]
        public void CorrectOrder()
        {
            _mzSpectrumA = new MzSpectrum(new double[] { 5, 6, 7 }, new double[] { 1, 2, 3 }, false);
            Assert.IsTrue(_mzSpectrumA.FilterByNumberOfMostIntense(2).First().Mz < _mzSpectrumA.FilterByNumberOfMostIntense(2).ToList()[1].Mz);
        }

        [Test]
        public void TestFunctionToX()
        {
            _mzSpectrumA.ReplaceXbyApplyingFunction(b => -1);
            Assert.AreEqual(-1, _mzSpectrumA.XArray[0]);
        }

        [Test]
        public void TestGetClosestPeakXValue()
        {
            Assert.AreEqual(447.73849, _mzSpectrumA.GetClosestPeakXvalue(447.73849));
            Assert.AreEqual(447.73849, _mzSpectrumA.GetClosestPeakXvalue(447));
            Assert.IsNull(new MzSpectrum(new double[0], new double[0], false).GetClosestPeakXvalue(1));
        }

        [Test]
        public void TestDotProduct()
        {
            double[] array1 = { 1 };
            double[] array2 = { 2 };
            double[] array3 = { 1, 2 };
            double[] array4 = { 1, 1 };

            MzSpectrum spec1 = new MzSpectrum(array1, array1, false);
            MzSpectrum spec2 = new MzSpectrum(array2, array1, false);
            MzSpectrum spec3 = new MzSpectrum(array3, array4, false);
            Tolerance tolerance = new PpmTolerance(10);

            Assert.AreEqual(spec1.CalculateDotProductSimilarity(spec3, tolerance), spec3.CalculateDotProductSimilarity(spec1, tolerance)); //comparison side shouldn't matter
            Assert.AreEqual(spec1.CalculateDotProductSimilarity(spec2, tolerance), 0); //orthogonal spectra give a score of zero
            Assert.AreEqual(spec2.CalculateDotProductSimilarity(spec2, tolerance), 1); //identical spectra give a score of 1
            Assert.IsTrue(tolerance.Within(spec3.CalculateDotProductSimilarity(spec2, tolerance), Math.Cos(Math.PI / 4)));
        }

        [Test]
        public void TestNumPeaksWithinRange()
        {
            double[] xArray = { 1, 2, 3, 4, 5, 6, 7 };
            double[] yArray = { 1, 2, 1, 5, 1, 2, 1 };

            var thisSpectrum = new MzSpectrum(xArray, yArray, false);

            Assert.AreEqual(7, thisSpectrum.NumPeaksWithinRange(double.MinValue, double.MaxValue));

            Assert.AreEqual(6, thisSpectrum.NumPeaksWithinRange(1, 7));

            Assert.AreEqual(0, thisSpectrum.NumPeaksWithinRange(1, 1));

            Assert.AreEqual(1, thisSpectrum.NumPeaksWithinRange(1, 2));

            Assert.AreEqual(2, thisSpectrum.NumPeaksWithinRange(0.001, 2.999));

            Assert.AreEqual(1, thisSpectrum.NumPeaksWithinRange(0, 1.5));

            Assert.AreEqual(1, thisSpectrum.NumPeaksWithinRange(6.5, 8));

            Assert.AreEqual(2, thisSpectrum.NumPeaksWithinRange(3, 5));

            Assert.AreEqual(2, thisSpectrum.NumPeaksWithinRange(3.5, 5.5));

            Assert.AreEqual(1, thisSpectrum.NumPeaksWithinRange(7, 8));

            Assert.AreEqual(0, thisSpectrum.NumPeaksWithinRange(8, 9));

            Assert.AreEqual(0, thisSpectrum.NumPeaksWithinRange(-2, -1));

            Assert.AreEqual("[1 to 7] m/z (Peaks 7)", thisSpectrum.ToString());

            //Assert.AreEqual(7, thisSpectrum.FilterByNumberOfMostIntense(7).Size);
            //Assert.AreEqual(1, thisSpectrum.FilterByNumberOfMostIntense(1).Size);
            //Assert.AreEqual(4, thisSpectrum.FilterByNumberOfMostIntense(1).FirstX);

            //Assert.AreEqual(2, thisSpectrum.FilterByNumberOfMostIntense(3).FirstX);

            //Assert.AreEqual(0, thisSpectrum.FilterByNumberOfMostIntense(0).Size);

            //Assert.AreEqual(2, thisSpectrum.WithRangeRemoved(2, 6).Size);
            //Assert.AreEqual(0, thisSpectrum.WithRangeRemoved(0, 100).Size);

            //Assert.AreEqual(6, thisSpectrum.WithRangeRemoved(7, 100).Size);

            //Assert.AreEqual(1, thisSpectrum.WithRangeRemoved(new DoubleRange(double.MinValue, 6)).Size);

            List<DoubleRange> xRanges = new List<DoubleRange>
            {
                new DoubleRange(2, 5),
                new DoubleRange(3, 6)
            };
            //Assert.AreEqual(2, thisSpectrum.WithRangesRemoved(xRanges).Size);

            //Assert.AreEqual(3, thisSpectrum.Extract(new DoubleRange(4.5, 10)).Size);

            //Assert.AreEqual(2, thisSpectrum.FilterByY(new DoubleRange(1.5, 2.5)).Size);

            //Assert.AreEqual(3, thisSpectrum.FilterByY(1.5, double.MaxValue).Size);

            //Assert.AreEqual(2, thisSpectrum.ApplyFunctionToX(b => b * 2).FirstX);

            Assert.AreEqual(1, thisSpectrum.GetClosestPeakXvalue(-100));

            Assert.AreEqual(7, thisSpectrum.GetClosestPeakXvalue(6.6));

            Assert.AreEqual(7, thisSpectrum.GetClosestPeakXvalue(7));

            Assert.AreEqual(7, thisSpectrum.GetClosestPeakXvalue(8));
        }

        [Test]
        public void TestEqualsAndHashCode()
        {
            // identical spectra, x and y arrays deep copied
            MzSpectrum identicalSpectrum = new(_mzSpectrumA.XArray, _mzSpectrumA.YArray, true);
            Assert.AreEqual(identicalSpectrum.GetHashCode(), _mzSpectrumA.GetHashCode());
            Assert.IsTrue(identicalSpectrum.Equals(_mzSpectrumA));
            Assert.IsTrue(identicalSpectrum.Equals((object)_mzSpectrumA));

            // changed x value
            identicalSpectrum.XArray[1] += 10;
            Assert.IsFalse(identicalSpectrum.Equals(_mzSpectrumA));
            Assert.IsFalse(identicalSpectrum.Equals((object)_mzSpectrumA));
            Assert.AreNotEqual(identicalSpectrum.GetHashCode(), _mzSpectrumA.GetHashCode());
            identicalSpectrum.XArray[1] -= 10;

            // changed y value
            identicalSpectrum.YArray[1] += 10;
            Assert.IsFalse(identicalSpectrum.Equals(_mzSpectrumA));
            Assert.IsFalse(identicalSpectrum.Equals((object)_mzSpectrumA));
            Assert.AreNotEqual(identicalSpectrum.GetHashCode(), _mzSpectrumA.GetHashCode());

            Assert.That(!_mzSpectrumA.Equals(null));
            Assert.That(!_mzSpectrumA.Equals((object)null));
            Assert.That(!_mzSpectrumA.Equals(2));
            Assert.That(!_mzSpectrumA.Equals((object)2));
        }
    }
}