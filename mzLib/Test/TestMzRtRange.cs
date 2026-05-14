using MzLibUtil;
using NUnit.Framework;
using System.Diagnostics.CodeAnalysis;
using Assert = NUnit.Framework.Legacy.ClassicAssert;

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestMzRtRange
    {
        [Test]
        public void Constructor_FourArgs_SetsAllFourBounds()
        {
            var range = new MzRtRange(minMZ: 100.0, maxMZ: 200.0, minRt: 5.0, maxRt: 10.0);
            Assert.AreEqual(100.0, range.Minimum);
            Assert.AreEqual(200.0, range.Maximum);
            Assert.AreEqual(5.0, range.MinimumRt);
            Assert.AreEqual(10.0, range.MaximumRt);
            Assert.AreEqual(100.0, range.MinimumMZ);
            Assert.AreEqual(200.0, range.MaximumMZ);
        }

        [Test]
        public void Constructor_FromMzRangeAndRetentionTime_AppliesTolerance()
        {
            var mzRange = new MzRange(100.0, 200.0);
            var range = new MzRtRange(mzRange, retentionTime: 7.5, rtTolerance: 0.25);
            Assert.AreEqual(100.0, range.Minimum);
            Assert.AreEqual(200.0, range.Maximum);
            Assert.AreEqual(7.25, range.MinimumRt, 1e-9);
            Assert.AreEqual(7.75, range.MaximumRt, 1e-9);
        }

        [Test]
        public void Constructor_FromMzRangeAndRetentionTime_DefaultToleranceIsTight()
        {
            var range = new MzRtRange(new MzRange(100.0, 200.0), retentionTime: 7.5);
            Assert.AreEqual(7.5 - 0.001, range.MinimumRt, 1e-9);
            Assert.AreEqual(7.5 + 0.001, range.MaximumRt, 1e-9);
        }

        [Test]
        public void DerivedProperties_Mean_And_Width_AreArithmetic()
        {
            var range = new MzRtRange(minMZ: 100.0, maxMZ: 200.0, minRt: 5.0, maxRt: 11.0);
            Assert.AreEqual(150.0, range.MeanMZ);
            Assert.AreEqual(8.0, range.MeanRt);
            Assert.AreEqual(100.0, range.MzRangeWidth);
            Assert.AreEqual(6.0, range.RtRangeWidth);
        }

        [Test]
        public void Contains_PointInside_True()
        {
            var range = new MzRtRange(100.0, 200.0, 5.0, 10.0);
            Assert.IsTrue(range.Contains(mz: 150.0, rt: 7.5));
        }

        [Test]
        public void Contains_MzBelowMinimum_False()
        {
            var range = new MzRtRange(100.0, 200.0, 5.0, 10.0);
            Assert.IsFalse(range.Contains(mz: 99.0, rt: 7.5));
        }

        [Test]
        public void Contains_MzAboveMaximum_False()
        {
            var range = new MzRtRange(100.0, 200.0, 5.0, 10.0);
            Assert.IsFalse(range.Contains(mz: 201.0, rt: 7.5));
        }

        [Test]
        public void Contains_RtBelowMinimum_False()
        {
            var range = new MzRtRange(100.0, 200.0, 5.0, 10.0);
            Assert.IsFalse(range.Contains(mz: 150.0, rt: 4.0));
        }

        [Test]
        public void Contains_RtAboveMaximum_False()
        {
            var range = new MzRtRange(100.0, 200.0, 5.0, 10.0);
            Assert.IsFalse(range.Contains(mz: 150.0, rt: 11.0));
        }

        [Test]
        public void Contains_BoundaryValues_Inclusive()
        {
            var range = new MzRtRange(100.0, 200.0, 5.0, 10.0);
            Assert.IsTrue(range.Contains(mz: 100.0, rt: 5.0));
            Assert.IsTrue(range.Contains(mz: 200.0, rt: 10.0));
        }

        [Test]
        public void CompareTo_PointInside_ReturnsZero()
        {
            var range = new MzRtRange(100.0, 200.0, 5.0, 10.0);
            Assert.AreEqual(0, range.CompareTo(mz: 150.0, rt: 7.5));
        }

        [Test]
        public void CompareTo_MzBelowMinimum_ReturnsPositive()
        {
            var range = new MzRtRange(100.0, 200.0, 5.0, 10.0);
            Assert.AreEqual(1, range.CompareTo(mz: 50.0, rt: 7.5));
        }

        [Test]
        public void CompareTo_MzAboveMaximum_ReturnsNegative()
        {
            var range = new MzRtRange(100.0, 200.0, 5.0, 10.0);
            Assert.AreEqual(-1, range.CompareTo(mz: 250.0, rt: 7.5));
        }

        [Test]
        public void CompareTo_RtBelowMinimum_ReturnsPositive()
        {
            var range = new MzRtRange(100.0, 200.0, 5.0, 10.0);
            Assert.AreEqual(1, range.CompareTo(mz: 150.0, rt: 1.0));
        }

        [Test]
        public void CompareTo_RtAboveMaximum_ReturnsNegative()
        {
            var range = new MzRtRange(100.0, 200.0, 5.0, 10.0);
            Assert.AreEqual(-1, range.CompareTo(mz: 150.0, rt: 99.0));
        }

        [Test]
        public void ToString_IncludesMzAndRtBoundsAndUnits()
        {
            var range = new MzRtRange(100.0, 200.0, 5.0, 10.0);
            string s = range.ToString("F1");
            // Doesn't pin exact format (culture-aware) but verifies both ranges appear.
            Assert.IsTrue(s.Contains("100.0") && s.Contains("200.0") && s.Contains("m/z"));
            Assert.IsTrue(s.Contains("5.0") && s.Contains("10.0") && s.Contains("RT"));
        }
    }
}
