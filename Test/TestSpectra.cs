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

using NUnit.Framework;
using Spectra;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public sealed class SpectrumTestFixture
    {

        #region Private Fields

        private DefaultMzSpectrum _mzSpectrumA;

        #endregion Private Fields

        #region Public Methods

        [SetUp]
        public void Setup()
        {
            double[] mz = { 328.73795, 329.23935, 447.73849, 448.23987, 482.23792, 482.57089, 482.90393, 500.95358, 501.28732, 501.62131, 611.99377, 612.32806, 612.66187, 722.85217, 723.35345 };
            double[] intensities = { 81007096.0, 28604418.0, 78353512.0, 39291696.0, 122781408.0, 94147520.0, 44238040.0, 71198680.0, 54184096.0, 21975364.0, 44514172.0, 43061628.0, 23599424.0, 56022696.0, 41019144.0 };

            _mzSpectrumA = new DefaultMzSpectrum(mz, intensities, false);
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
        public void SpectrumBasePeakIntensity()
        {
            double basePeakIntensity = _mzSpectrumA.YofPeakWithHighestY;

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
            double intensity = _mzSpectrumA.YArray[0];

            Assert.AreEqual(81007096.0, intensity);
        }

        [Test]
        public void SpectrumGetIntensityRandom()
        {
            double intensity = _mzSpectrumA.YArray[6];

            Assert.AreEqual(44238040.0, intensity);
        }

        [Test]
        public void SpectrumGetMassFirst()
        {
            double intensity = _mzSpectrumA.XArray[0];

            Assert.AreEqual(328.73795, intensity);
        }

        [Test]
        public void SpectrumGetMassRandom()
        {
            double intensity = _mzSpectrumA.XArray[6];

            Assert.AreEqual(482.90393, intensity);
        }

        [Test]
        public void SpectrumContainsPeak()
        {
            Assert.IsTrue(_mzSpectrumA.Count > 0);
        }

        [Test]
        public void SpectrumContainsPeakInRange()
        {
            Assert.AreEqual(1, _mzSpectrumA.NumPeaksWithinRange(448.23987 - 0.001, 448.23987 + 0.001));
        }

        [Test]
        public void SpectrumContainsPeakInRangeEnd()
        {
            Assert.AreEqual(1, _mzSpectrumA.NumPeaksWithinRange(448.23987 - 0.001, 448.23987));
        }

        [Test]
        public void SpectrumContainsPeakInRangeStart()
        {
            Assert.AreEqual(1, _mzSpectrumA.NumPeaksWithinRange(448.23987, 448.23987 + 0.001));
        }

        [Test]
        public void SpectrumContainsPeakInRangeStartEnd()
        {
            Assert.AreEqual(1, _mzSpectrumA.NumPeaksWithinRange(448.23987, 448.23987));
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
            var filteredMzSpectrum = _mzSpectrumA.NewSpectrumFilterByY(28604417, 28604419);

            Assert.AreEqual(1, filteredMzSpectrum.Count);
        }

        [Test]
        public void SpectrumSelect()
        {
            MzSpectrum<MzPeak> v2 = _mzSpectrumA;
            ISpectrum<Peak> v3 = v2;

            v3.Take(4);

            var v5 = v3.Select(b => b.X);
            Assert.AreEqual(328.73795, v5.First());

            var bn = v2[0];

            var bsrg = _mzSpectrumA[0];
        }

        [Test]
        public void FilterByNumberOfMostIntenseTest()
        {
            Assert.AreEqual(5, _mzSpectrumA.NewSpectrumFilterByNumberOfMostIntense(5).Count);
        }

        [Test]
        public void GetBasePeak()
        {
            Assert.AreEqual(122781408.0, _mzSpectrumA.PeakWithHighestY.Intensity);
        }

        [Test]
        public void GetClosestPeak()
        {
            Assert.AreEqual(448.23987, _mzSpectrumA.GetClosestPeak(448).Mz);
            Assert.AreEqual(447.73849, _mzSpectrumA.GetClosestPeak(447.9).Mz);
        }

        [Test]
        public void Extract()
        {
            Assert.AreEqual(3, _mzSpectrumA.NewSpectrumExtract(500, 600).Count);
        }

        [Test]
        public void CorrectOrder()
        {
            _mzSpectrumA = new DefaultMzSpectrum(new double[] { 5, 6, 7 }, new double[] { 1, 2, 3 }, false);
            Assert.IsTrue(_mzSpectrumA.NewSpectrumFilterByNumberOfMostIntense(2)[0].Mz < _mzSpectrumA.NewSpectrumFilterByNumberOfMostIntense(2)[1].Mz);
        }

        [Test]
        public void FilterByMZ()
        {
            List<DoubleRange> ok = new List<DoubleRange>();
            ok.Add(new DoubleRange(300, 400));
            ok.Add(new DoubleRange(700, 800));
            Assert.AreEqual(11, _mzSpectrumA.NewSpectrumWithRangesRemoved(ok).Count);

            Assert.AreEqual(10, _mzSpectrumA.NewSpectrumWithRangeRemoved(new DoubleRange(400, 500)).Count);
            Assert.AreEqual(10, _mzSpectrumA.NewSpectrumWithRangeRemoved(400, 500).Count);
        }

        [Test]
        public void Test2D()
        {
            double[,] array2D = { { 1, 2, 3, 4 }, { 5, 6, 7, 8 } };
            _mzSpectrumA = new DefaultMzSpectrum(array2D);
            Assert.AreEqual("1 - 4 m/z (Peaks 4)", _mzSpectrumA.ToString());
        }

        [Test]
        public void TestClone()
        {
            var ok = new DefaultSpectrum(_mzSpectrumA);
            Assert.AreNotEqual(ok[0], _mzSpectrumA[0]);
            Assert.AreEqual(ok[0].Y, _mzSpectrumA[0].Y);

            var ok2 = new DefaultSpectrum(_mzSpectrumA.XArray, _mzSpectrumA.YArray, true);

            ok2.XArray[0] = 0;
            Assert.AreNotEqual(ok2.XArray[0], _mzSpectrumA.XArray[0]);

            var ok3 = new DefaultSpectrum(_mzSpectrumA.XArray, _mzSpectrumA.YArray, false);

            ok3.XArray[0] = 0;
            Assert.AreEqual(ok3.XArray[0], _mzSpectrumA.XArray[0]);
        }

        [Test]
        public void TestFunctionToX()
        {
            var ok = _mzSpectrumA.NewSpectrumApplyFunctionToX(b => -1);
            Assert.AreEqual(-1, ok[0].X);
        }

        [Test]
        public void Test2dArray()
        {
            var ok = _mzSpectrumA.CopyTo2DArray();
            var ok2 = new DefaultMzSpectrum(ok);
            Assert.AreEqual(ok2[0].X, _mzSpectrumA[0].X);
        }

        [Test]
        public void TestFilterAndExtract()
        {
            Assert.AreEqual(447.73849, _mzSpectrumA.NewSpectrumFilterByY(new DoubleRange(78353510, 81007097)).NewSpectrumExtract(new DoubleRange(400, 500)).GetClosestPeak(327.5).X);
        }

        [Test]
        public void TestGetClosestPeakXValue()
        {
            Assert.AreEqual(447.73849, _mzSpectrumA.GetClosestPeakXvalue(447.73849));
            Assert.AreEqual(447.73849, _mzSpectrumA.GetClosestPeakXvalue(447));
            Assert.Throws<IndexOutOfRangeException>(() => { new DefaultMzSpectrum(new double[0], new double[0], false).GetClosestPeakXvalue(1); }, "No peaks in spectrum!");
        }

        [Test]
        public void IMzSpectrumTpeakTest()
        {
            //double[] mz = { 328.73795, 329.23935, 447.73849, 448.23987, 482.23792, 482.57089, 482.90393, 500.95358, 501.28732, 501.62131, 611.99377, 612.32806, 612.66187, 722.85217, 723.35345 };
            //double[] intensities = { 81007096.0, 28604418.0, 78353512.0, 39291696.0, 122781408.0, 94147520.0, 44238040.0, 71198680.0, 54184096.0, 21975364.0, 44514172.0, 43061628.0, 23599424.0, 56022696.0, 41019144.0 };

            IMzSpectrum<MzPeak> ok = _mzSpectrumA;
            IMzSpectrum<MzPeak> ok2 = new DefaultMzSpectrum(ok);

            Assert.Greater(100, ok2.NewSpectrumApplyFunctionToX(b => b / 10).LastX);

            Assert.AreEqual(1, ok2.NewSpectrumExtract(new DoubleRange(723, 1000)).Count);

            Assert.AreEqual(15, ok2.NewSpectrumExtract(double.MinValue, double.MaxValue).Count);

            Assert.AreEqual(0, ok2.NewSpectrumFilterByNumberOfMostIntense(0).Count);
            Assert.AreEqual(5, ok2.NewSpectrumFilterByNumberOfMostIntense(5).Count);
            Assert.AreEqual(15, ok2.NewSpectrumFilterByNumberOfMostIntense(15).Count);

            Assert.AreEqual(15, ok2.NewSpectrumFilterByY(new DoubleRange(double.MinValue, double.MaxValue)).Count);

            Assert.AreEqual(1, ok2.NewSpectrumFilterByY(39291695, 39291697).Count);

            Assert.AreEqual(2, ok2.NewSpectrumWithRangeRemoved(new DoubleRange(329, 723)).Count);

            Assert.AreEqual(15, ok2.NewSpectrumWithRangeRemoved(0, 1).Count);

            Assert.AreEqual(2, ok2.NewSpectrumWithRangesRemoved(new List<DoubleRange>{ new DoubleRange(329, 400), new DoubleRange(400, 723) }).Count);
        }

        [Test]
        public void TestDefaultSpectrum()
        {
            double[,] mzIntensities = new double[2, 3];
            mzIntensities[0, 0] = 1;
            mzIntensities[1, 0] = 2;
            mzIntensities[0, 1] = 3;
            mzIntensities[1, 1] = 4;
            mzIntensities[0, 2] = 5;
            mzIntensities[1, 2] = 6;

            DefaultMzSpectrum ok = new DefaultMzSpectrum(mzIntensities);

            Assert.AreEqual(6, ok.YofPeakWithHighestY);

            DefaultSpectrum ok2 = new DefaultSpectrum(mzIntensities);

            Assert.AreEqual(12, ok2.SumOfAllY);
            Assert.AreEqual(1, ok2.Range.Minimum);
            Assert.AreEqual(5, ok2.Range.Maximum);
        }

        [Test]
        public void TestNumPeaksWithinRange()
        {
            double[] xArray = { 1, 2, 3, 4, 5, 6, 7 };
            double[] yArray = { 1, 2, 1, 5, 1, 2, 1 };

            var thisSpectrum = new DefaultSpectrum(xArray, yArray, false);

            Assert.AreEqual(7, thisSpectrum.NumPeaksWithinRange(double.MinValue, double.MaxValue));

            Assert.AreEqual(7, thisSpectrum.NumPeaksWithinRange(1, 7));

            Assert.AreEqual(1, thisSpectrum.NumPeaksWithinRange(1, 1));

            Assert.AreEqual(2, thisSpectrum.NumPeaksWithinRange(1, 2));

            Assert.AreEqual(2, thisSpectrum.NumPeaksWithinRange(0.001, 2.999));

            Assert.AreEqual(1, thisSpectrum.NumPeaksWithinRange(0, 1.5));

            Assert.AreEqual(1, thisSpectrum.NumPeaksWithinRange(6.5, 8));

            Assert.AreEqual(3, thisSpectrum.NumPeaksWithinRange(3, 5));

            Assert.AreEqual(2, thisSpectrum.NumPeaksWithinRange(3.5, 5.5));

            Assert.AreEqual(1, thisSpectrum.NumPeaksWithinRange(7, 8));

            Assert.AreEqual(0, thisSpectrum.NumPeaksWithinRange(8, 9));

            Assert.AreEqual(0, thisSpectrum.NumPeaksWithinRange(-2, -1));

            Assert.AreEqual("[1 - 7] (Peaks 7)", thisSpectrum.ToString());

            Assert.AreEqual(7, thisSpectrum.NewSpectrumFilterByNumberOfMostIntense(7).Count);
            Assert.AreEqual(1, thisSpectrum.NewSpectrumFilterByNumberOfMostIntense(1).Count);
            Assert.AreEqual(4, thisSpectrum.NewSpectrumFilterByNumberOfMostIntense(1).FirstX);

            Assert.AreEqual(2, thisSpectrum.NewSpectrumFilterByNumberOfMostIntense(3).FirstX);

            Assert.AreEqual(0, thisSpectrum.NewSpectrumFilterByNumberOfMostIntense(0).Count);

            Assert.AreEqual(2, thisSpectrum.NewSpectrumWithRangeRemoved(2, 6).Count);
            Assert.AreEqual(0, thisSpectrum.NewSpectrumWithRangeRemoved(0, 100).Count);

            Assert.AreEqual(6, thisSpectrum.NewSpectrumWithRangeRemoved(7, 100).Count);

            Assert.AreEqual(1, thisSpectrum.NewSpectrumWithRangeRemoved(new DoubleRange(double.MinValue, 6)).Count);

            List<DoubleRange> xRanges = new List<DoubleRange>();
            xRanges.Add(new DoubleRange(2, 5));
            xRanges.Add(new DoubleRange(3, 6));
            Assert.AreEqual(2, thisSpectrum.NewSpectrumWithRangesRemoved(xRanges).Count);

            Assert.AreEqual(3, thisSpectrum.NewSpectrumExtract(new DoubleRange(4.5, 10)).Count);

            Assert.AreEqual(2, thisSpectrum.NewSpectrumFilterByY(new DoubleRange(1.5, 2.5)).Count);

            Assert.AreEqual(3, thisSpectrum.NewSpectrumFilterByY(1.5, double.MaxValue).Count);

            Assert.AreEqual(2, thisSpectrum.NewSpectrumApplyFunctionToX(b => b * 2).FirstX);

            Assert.AreEqual(7, thisSpectrum.GetClosestPeak(6.6).X);

            Assert.AreEqual(7, thisSpectrum.GetClosestPeak(7).X);

            Assert.AreEqual(7, thisSpectrum.GetClosestPeak(8).X);

            IEnumerable hnm = thisSpectrum;

            double dudu = 0;
            foreach (var ikik in hnm)
            {
                dudu += ((Peak)ikik).X;
            }
            Assert.AreEqual(1 + 2 + 3 + 4 + 5 + 6 + 7, dudu);
        }

        #endregion Public Methods

    }
}