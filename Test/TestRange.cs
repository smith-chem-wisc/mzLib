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

using NUnit.Framework;
using Spectra;
using System;
using System.Collections.Generic;

namespace Test
{
    [TestFixture]
    public sealed class RangeTest
    {

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
        public void RangeCompareToBelow()
        {
            var range = new DoubleRange(3, 10);
            int value = 1;

            int comp = range.CompareTo(value);

            Assert.AreEqual(-1, comp);
        }

        [Test]
        public void RangeCompareToWithin()
        {
            var range = new DoubleRange(3, 10);
            int value = 5;

            int comp = range.CompareTo(value);

            Assert.AreEqual(0, comp);
        }

        [Test]
        public void RangeCompareToAbove()
        {
            var range = new DoubleRange(3, 10);
            int value = 12;

            int comp = range.CompareTo(value);

            Assert.AreEqual(1, comp);
        }

        [Test]
        public void RangesAreEquivalent()
        {
            var range1 = new DoubleRange(3, 10);
            var range2 = new DoubleRange(3, 10);

            Assert.AreEqual(range1, range2);
        }

        [Test]
        public void RangesAreEquivalentNotReference()
        {
            var range1 = new DoubleRange(3, 10);
            var range2 = new DoubleRange(3, 10);

            Assert.AreNotSame(range1, range2);
        }

        [Test]
        public void RangeMinBiggerThanMax()
        {
            Assert.Throws<ArgumentException>(() => { new DoubleRange(10, 5); });
        }

        [Test]
        public void MassRangeFromDAWidth()
        {
            var range1 = new DoubleRange(10, new Tolerance(ToleranceUnit.Absolute, 4));

            Assert.AreEqual(8, range1.Width);
        }

        [Test]
        public void MassRangeFromDAMean()
        {
            var range1 = new DoubleRange(10, new Tolerance(ToleranceUnit.Absolute, 4));

            Assert.AreEqual(10, range1.Mean);
        }

        [Test]
        public void MassRangeFromDAMin()
        {
            var range1 = new DoubleRange(10, new Tolerance(ToleranceUnit.Absolute, 4));

            Assert.AreEqual(6, range1.Minimum);
        }

        [Test]
        public void MassRangeFromDAMax()
        {
            var range1 = new DoubleRange(10, new Tolerance(ToleranceUnit.Absolute, 4));

            Assert.AreEqual(14, range1.Maximum);
        }

        [Test]
        public void MassRangeFromDANullMean()
        {
            var range1 = new DoubleRange(10, null);

            Assert.AreEqual(10, range1.Mean);
        }

        [Test]
        public void MassRangeFromDANullWidth()
        {
            var range1 = new DoubleRange(10, null);

            Assert.AreEqual(0, range1.Width);
        }

        [Test]
        public void MassRangeFromDANullMin()
        {
            var range1 = new DoubleRange(10, null);

            Assert.AreEqual(10, range1.Minimum);
        }

        [Test]
        public void MassRangeFromDANullMax()
        {
            var range1 = new DoubleRange(10, null);

            Assert.AreEqual(10, range1.Maximum);
        }

        [Test]
        public void MassRangeFromDANegative()
        {
            var range1 = new DoubleRange(10, new Tolerance(ToleranceUnit.Absolute, 4));
            var range2 = new DoubleRange(10, new Tolerance(ToleranceUnit.Absolute, -4));

            Assert.AreEqual(range1, range2);
        }


        [Test]
        public void RangeFromRange()
        {
            var range1 = new DoubleRange(10, new Tolerance(ToleranceUnit.Absolute, 4));
            var range2 = new DoubleRange(range1);
            Assert.AreEqual(range1, range2);
        }

        [Test]
        public void SuperRange()
        {
            var range1 = new DoubleRange(10, new Tolerance(ToleranceUnit.Absolute, 4));
            var range2 = new DoubleRange(10, new Tolerance(ToleranceUnit.Absolute, 3));
            Assert.IsTrue(range1.IsSuperRange(range2));
        }

        [Test]
        public void TestDoubleRangeStuff()
        {
            DoubleRange range1 = new DoubleRange(new DoubleRange(1000000 - 1, 1000000 + 1));
            DoubleRange range2 = new DoubleRange(1000000, new Tolerance(ToleranceUnit.PPM, 1));

            Assert.IsTrue(range1.Equals(range2));
            Assert.AreEqual("[999999 - 1000001]", range1.ToString());
        }

        [Test]
        public void TestHashSet()
        {
            HashSet<DoubleRange> ok = new HashSet<DoubleRange>();
            ok.Add(new DoubleRange(1, 2));
            ok.Add(new DoubleRange(2, 3));
            ok.Add(new DoubleRange(1, 2));
            Assert.AreEqual(2, ok.Count);
        }

    }
}