// Copyright 2012, 2013, 2014 Derek J. Bailey
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
    }
}