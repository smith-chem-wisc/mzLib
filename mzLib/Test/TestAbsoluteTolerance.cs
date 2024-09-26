using NUnit.Framework;
using System.Diagnostics.CodeAnalysis;
using MzLibUtil;
using System.IO;
using System.Linq;
using System;
using System.Collections.Generic;

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestAbsoluteTolerance
    {
        private AbsoluteTolerance _tolerance = new AbsoluteTolerance(15); 
        [Test]
        public void TestAbsoluteToleranceInitialization()
        {
            Assert.AreEqual(_tolerance.Value, 15);
        }

        [Test]
        public void TestAbsoluteToleranceToString()
        {
            string output = "±15 Absolute"; 
            Assert.AreEqual(output, _tolerance.ToString());
        }

        [Test]
        public void TestDoubleRange()
        {
            double mean = 1000; 
            DoubleRange range = _tolerance.GetRange(mean);
            double expectedMax = mean + 15d;
            double expectedMin = mean - 15d;
            double width = expectedMax - expectedMin; 
            Assert.AreEqual(range.Maximum, expectedMax);
            Assert.AreEqual(range.Minimum, expectedMin);
            Assert.AreEqual(range.Width, Is.EqualTo(width).Within(0.001));
        }
    }
}