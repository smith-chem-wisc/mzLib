// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work Copyright 2016 Stefan Solntsev
//
// This file (TestTolerance.cs) is part of MassSpectrometry.Tests.
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

using ZMzLibUtil;
using NUnit.Framework;

namespace Test
{
    [TestFixture]
    public sealed class MassToleranceTestFixture
    {
        #region Public Methods

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

        #endregion Public Methods
    }
}