﻿// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work Copyright 2016 Stefan Solntsev
//
// This file (TestRange.cs) is part of MassSpectrometry.Tests.
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

using MzLibUtil;
using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System;
using System.Linq;
using Stopwatch = System.Diagnostics.Stopwatch;
using Newtonsoft.Json.Linq;
using System.Collections.Generic;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestRangeAndTolerances
    {
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

        #region DoubleRangeTests

        [Test]
        public void RangeCompareToTest()
        {
            // In C#, the convention is that when writing CompareTo methods,
            // the instance is compared to the object (argument)
            // If the instance is less than the argument, a negative value is returned
            // If the instance is greater than the argument, a positive value is returned
            // An example of this is given below
            Double five = 5;
            double ten = 10;
            Assert.That(five.CompareTo(ten) < 0);

            // Range is greater than value
            var range = new DoubleRange(3, 10);
            int value = 1;

            int comparison = range.CompareTo(value);
            Assert.AreEqual(1, comparison);

            // Range contains value
            value = 5;

            comparison = range.CompareTo(value);
            Assert.AreEqual(0, comparison);
            
            //Range is less than value
            value = 12;

            comparison = range.CompareTo(value);
            Assert.AreEqual(-1, comparison);
        }

        [Test]
        public void RangeSubRange()
        {
            var a = new DoubleRange(5, 7);
            var b = new DoubleRange(0, 10);

            Assert.IsTrue(a.IsSubRange(b));
        }

        [Test]
        public void RangeSubRangeReverseIsFalse()
        {
            var a = new DoubleRange(5, 7);
            var b = new DoubleRange(0, 10);

            Assert.IsFalse(b.IsSubRange(a));
        }

        [Test]
        public void RangeOverlappingIsFalse()
        {
            var range1 = new DoubleRange(5, 10);
            var range2 = new DoubleRange(15, 20);

            Assert.IsFalse(range1.IsOverlapping(range2));
        }

        [Test]
        public void RangeOverlappingIsFalseReverse()
        {
            var range1 = new DoubleRange(5, 10);
            var range2 = new DoubleRange(15, 20);

            Assert.IsFalse(range2.IsOverlapping(range1));
        }

        [Test]
        public void RangeOverlappingIsTrue()
        {
            var range1 = new DoubleRange(5, 10);
            var range2 = new DoubleRange(7, 12);

            Assert.IsTrue(range1.IsOverlapping(range2));
        }

        [Test]
        public void RangeOverlappingIsTrueReverse()
        {
            var range1 = new DoubleRange(5, 10);
            var range2 = new DoubleRange(7, 12);

            Assert.IsTrue(range2.IsOverlapping(range1));
        }

        [Test]
        public void RangeOverlappingIsTrueLarger()
        {
            var range1 = new DoubleRange(0, 10);
            var range2 = new DoubleRange(3, 7);

            Assert.IsTrue(range1.IsOverlapping(range2));
        }

        [Test]
        public void RangeOverlappingIsTrueSmaller()
        {
            var range1 = new DoubleRange(0, 10);
            var range2 = new DoubleRange(3, 7);

            Assert.IsTrue(range2.IsOverlapping(range1));
        }

        [Test]
        public void RangeOverlappingIsTrueItSelf()
        {
            var range1 = new DoubleRange(0, 10);

            Assert.IsTrue(range1.IsOverlapping(range1));
        }

        [Test]
        public void RangeDoesContainItem()
        {
            var range1 = new DoubleRange(3, 10);

            Assert.IsTrue(range1.Contains(5));
        }

        [Test]
        public void RangeDoesnotContainItemHigher()
        {
            var range1 = new DoubleRange(3, 10);

            Assert.IsFalse(range1.Contains(12));
        }

        [Test]
        public void RangeDoesnotContainItemLower()
        {
            var range1 = new DoubleRange(3, 10);

            Assert.IsFalse(range1.Contains(1));
        }

        [Test]
        public void RangeDoesContainItemLowerBounds()
        {
            var range1 = new DoubleRange(3, 10);

            Assert.IsTrue(range1.Contains(3));
        }

        [Test]
        public void RangeDoesContainItemUpperBounds()
        {
            var range1 = new DoubleRange(3, 10);

            Assert.IsTrue(range1.Contains(10));
        }

        [Test]
        public void RangesAreEquivalentNotReference()
        {
            var range1 = new DoubleRange(3, 10);
            var range2 = new DoubleRange(3, 10);

            Assert.AreNotSame(range1, range2);
        }


        #endregion

        [Test]
        public void ToleranceParseAndWithin()
        {
            string tol = "4 Absolute";
            var range1 = new AbsoluteTolerance(4);
            var range2 = Tolerance.ParseToleranceString(tol);
            Assert.AreEqual(range1.Value, range2.Value);
            Assert.AreEqual(range1.Value, Tolerance.ParseToleranceString(range2.ToString()).Value);

            Assert.IsTrue(range1.Within(3.9, 0));
            Assert.IsTrue(range1.Within(4, 0));
            Assert.IsFalse(range1.Within(4.1, 0));

            Assert.IsTrue(range1.Within(0, 3.9));
            Assert.IsTrue(range1.Within(0, 4));
            Assert.IsFalse(range1.Within(0, 4.1));
        }

        [Test]
        public void MassRangeFromDAWidth()
        {
            var range1 = new AbsoluteTolerance(4).GetRange(10);

            Assert.AreEqual(8, range1.Width);
        }

        [Test]
        public void MassRangeFromDAMean()
        {
            var range1 = new AbsoluteTolerance(4).GetRange(10);

            Assert.AreEqual(10, range1.Mean);
        }

        [Test]
        public void MassRangeFromDAMin()
        {
            var range1 = new AbsoluteTolerance(4).GetRange(10);

            Assert.AreEqual(6, range1.Minimum);
        }

        [Test]
        public void MassRangeFromDAMax()
        {
            var range1 = new AbsoluteTolerance(4).GetRange(10);

            Assert.AreEqual(14, range1.Maximum);
        }

        [Test]
        public void MassRangeFromDANegative()
        {
            var range1 = new AbsoluteTolerance(4).GetRange(10);
            var range2 = new AbsoluteTolerance(-4).GetRange(10);

            Assert.AreEqual(0, range1.Minimum - range2.Minimum, 1e-9);
            Assert.AreEqual(0, range1.Maximum - range2.Maximum, 1e-9);
        }

        [Test]
        public void RangeFromRange()
        {
            var range1 = new AbsoluteTolerance(4).GetRange(10);
            var range2 = new DoubleRange(range1);
            Assert.AreEqual(0, range1.Minimum - range2.Minimum, 1e-9);
            Assert.AreEqual(0, range1.Maximum - range2.Maximum, 1e-9);
        }

        [Test]
        public void SuperRange()
        {
            var range1 = new AbsoluteTolerance(4).GetRange(10);
            var range2 = new AbsoluteTolerance(3).GetRange(10);
            Assert.IsTrue(range1.IsSuperRange(range2));
        }

        [Test]
        public static void TestMajorityWithin()
        {
            var testArr = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
            var testList = testArr.ToList();

            var goodRange = new DoubleRange(1, 6);
            var badRange = new DoubleRange(1, 4);
            Assert.IsTrue(goodRange.ContainsMajority(testList));
            Assert.IsFalse(badRange.ContainsMajority(testList));

            Assert.IsTrue(goodRange.ContainsMajority(testArr));
            Assert.IsFalse(badRange.ContainsMajority(testArr));
        }

        [Test]
        public void TestDoubleRangeStuff()
        {
            DoubleRange range1 = new DoubleRange(new DoubleRange(1000000 - 1, 1000000 + 1));
            DoubleRange range2 = new PpmTolerance(1).GetRange(1000000);

            Assert.AreEqual(0, range1.Minimum - range2.Minimum, 1e-9);
            Assert.AreEqual(0, range1.Maximum - range2.Maximum, 1e-9);
            Assert.AreEqual("[999999;1000001]", range1.ToString());
        }

        [Test]
        public void MassToleranceConstructorDaValue()
        {
            var tol = new AbsoluteTolerance(10);

            Assert.AreEqual(10, tol.Value);
        }

        [Test]
        public void MassToleranceImplicitValue()
        {
            var tol = Tolerance.ParseToleranceString("10 ppm");

            Assert.AreEqual(10, tol.Value);

            Assert.AreEqual(1 + 1e-5, tol.GetMaximumValue(1));
            Assert.AreEqual(1 - 1e-5, tol.GetMinimumValue(1));
        }

        [Test]
        public void ToleranceWithin1()
        {
            var tol = new PpmTolerance(10);

            Assert.IsTrue(tol.Within(500, 500.005));
        }

        [Test]
        public void ToleranceNewTest()
        {
            var tol = new AbsoluteTolerance(1);
            Assert.AreEqual(4, tol.GetRange(5).Minimum);
            Assert.AreEqual(6, tol.GetRange(5).Maximum);
        }

        [Test]
        public void ToleranceMinMaxTest()
        {
            var tol = new AbsoluteTolerance(1);
            Assert.AreEqual(2, tol.GetMaximumValue(1));
            Assert.AreEqual(0, tol.GetMinimumValue(1));
        }

        [Test]
        public void TolerancePPMGetRange()
        {
            var tol = new PpmTolerance(1);
            Assert.AreEqual(20, tol.GetRange(1e7).Width);

            Assert.AreEqual("±1.0000 PPM", tol.ToString());
        }

        [Test]
        public void PpmToleranceWithNotch_()
        {
            var tol = new PpmToleranceWithNotch(10, 2, 1);
            Assert.AreEqual(tol.GetMaximumValue(100), (100 + 2 * 1.00335483810) * (1 + (10 / 1e6)));
            Assert.AreEqual(tol.GetMinimumValue(100), (100 - 1 * 1.00335483810) * (1 - (10 / 1e6)));
            Assert.AreEqual(tol.GetRange(100).Maximum, (100 + 2 * 1.00335483810) * (1 + (10 / 1e6)));
            Assert.AreEqual(tol.GetRange(100).Minimum, (100 - 1 * 1.00335483810) * (1 - (10 / 1e6)));

            //notch +2
            Assert.IsTrue(tol.Within(102.0067, 100));
            //notch -1
            Assert.IsTrue(tol.Within(98.997, 100));
            //notch 0
            Assert.IsTrue(tol.Within(100.0009, 100));
            //between notches
            Assert.IsFalse(tol.Within(100.5, 100));
            //notch -2
            Assert.IsFalse(tol.Within(97.993, 100));
        }

        [Test]
        [TestCase(100, 10, 2, 1)]
        [TestCase(100, 10, 0, 0)]
        [TestCase(100, 10, 10, 10)]
        [TestCase(1e4, 10, 2, 1)]
        [TestCase(100, 1e-2, 2, 1)]
        [TestCase(100, 1e2, 2, 1)]
        [TestCase(100, 10, 10, 0)]
        [TestCase(100, 10, 0, 10)]
        [TestCase(100, 10, 1, 2)]
        [TestCase(100, 10, 0, 2)]
        [TestCase(100, 10, 2, 0)]
        [TestCase(100, 10, 50, 50)]
        [TestCase(1.10335483810, 10, 2, 1)]
        [TestCase(100, 10, 1, 1)]
        [TestCase(100, 10, 5, 5)]
        [TestCase(100, 10, 20, 20)]
        [TestCase(100, 10, 0, 0)]
        [TestCase(100, 10, 1, 0)]
        [TestCase(100, 10, 0, 1)]
        public void PpmToleranceWithNotchTest(double baseValue, double tolerance, int positiveNotches, int negativeNotches)
        {
            var tol = new PpmToleranceWithNotch(tolerance, positiveNotches, negativeNotches);
            double maxValue = (baseValue + positiveNotches * PpmToleranceWithNotch.NotchStep) * (1 + (tolerance / 1e6));
            double minValue = (baseValue - negativeNotches * PpmToleranceWithNotch.NotchStep) * (1 - (tolerance / 1e6));
            Assert.That(tol.GetMaximumValue(baseValue), Is.EqualTo(maxValue).Within(1e-6), "GetMaximumValue failed");
            Assert.That(tol.GetMinimumValue(baseValue), Is.EqualTo(minValue).Within(1e-6), "GetMinimumValue failed");
            Assert.That(tol.GetRange(baseValue).Maximum, Is.EqualTo(maxValue).Within(1e-6), "GetRange.Maximum failed");
            Assert.That(tol.GetRange(baseValue).Minimum, Is.EqualTo(minValue).Within(1e-6), "GetRange.Minimum failed");
            // Test values at notches and between
            double notchPlus = (baseValue + positiveNotches * PpmToleranceWithNotch.NotchStep) * (1 + (tolerance / 1e6)) - 1e-6;
            double notchMinus = (baseValue - negativeNotches * PpmToleranceWithNotch.NotchStep) * (1 - (tolerance / 1e6)) + 1e-6;
            double notchZero = baseValue * (1 + (tolerance / 1e6)) - 1e-6;
            double betweenNotches = baseValue + (positiveNotches + 1) * PpmToleranceWithNotch.NotchStep;
            double belowMin = minValue - 10;
            double aboveMax = maxValue + 10;
            Assert.That(tol.Within(notchPlus, baseValue), Is.True, "Within failed at notchPlus");
            Assert.That(tol.Within(notchMinus, baseValue), Is.True, "Within failed at notchMinus");
            Assert.That(tol.Within(notchZero, baseValue), Is.True, "Within failed at notchZero");
            if (tolerance <= 100) // Only test between notches if both notches exist
            {
                Assert.That(tol.Within(betweenNotches, baseValue), Is.False, "Within should not match at betweenNotches for high tolerance");
            }
            else
            {
                Assert.That(tol.Within(betweenNotches, baseValue), Is.True, "Within failed at betweenNotches");
            }
            Assert.That(tol.Within(belowMin, baseValue), Is.False, "Within failed at belowMin");
            Assert.That(tol.Within(aboveMax, baseValue), Is.False, "Within failed at aboveMax");
            // Edge: test all notches individually
            for (int i = -negativeNotches; i <= positiveNotches; i++)
            {
                double shifted = baseValue + i * PpmToleranceWithNotch.NotchStep;
                double testVal = shifted * (1 + (tolerance / 1e6)) - 1e-6;
                Assert.That(tol.Within(testVal, baseValue), Is.True, $"Within failed at notch {i} (testVal={testVal})");
            }
            // Edge: test just outside all notches
            for (int i = -negativeNotches; i <= positiveNotches; i++)
            {
                double shifted = baseValue + i * PpmToleranceWithNotch.NotchStep;
                double testVal = shifted * (1 + (tolerance / 1e6)) + 1e-3;
                Assert.That(tol.Within(testVal, baseValue), Is.False, $"Within failed at outside notch {i} (testVal={testVal})");
            }
            // Edge: test zero tolerance (should only match exact)
            if (tolerance == 0)
            {
                Assert.That(tol.Within(baseValue, baseValue), Is.True, "Zero tolerance should match exact");
                Assert.That(tol.Within(baseValue + 1e-6, baseValue), Is.False, "Zero tolerance should not match offset");
            }
        }
    }
}