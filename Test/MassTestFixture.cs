// Copyright 2012, 2013, 2014 Derek J. Bailey
//
// This file (MassTestFixture.cs) is part of CSMSL.Tests.
//
// CSMSL.Tests is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// CSMSL.Tests is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with CSMSL.Tests. If not, see <http://www.gnu.org/licenses/>.

using Chemistry;
using NUnit.Framework;

namespace Test
{
    [TestFixture]
    public class MassTestFixture
    {
        [Test]
        public void MassToMzToMass()
        {
            ObjectWithMass1000 a = new ObjectWithMass1000();
            double mz = a.ToMz(2).ToMass(2);
            Assert.AreEqual(1000, mz);
        }

        [Test]
        public void MassToMzPositiveCharge()
        {
            ObjectWithMass1000 a = new ObjectWithMass1000();
            double mz = a.ToMz(2);
            Assert.AreEqual(501.00727646687898, mz);
        }

        [Test]
        public void MassToMzPositiveCharge_plus3()
        {
            ObjectWithMass1000 a = new ObjectWithMass1000();
            double mz = (double)ClassExtensions.RoundedDouble(a.ToMz(3),9);
            Assert.AreEqual(334.340609800, mz);
        }


        [Test]
        public void MassToMzNegativeCharge()
        {
            ObjectWithMass1000 a = new ObjectWithMass1000();
            double mz = a.ToMz(-2);
            Assert.AreEqual(498.99272353312102, mz);
        }

        [Test]
        public void MzToMassPostitiveCharge()
        {
            double a = 524.3;
            Assert.AreEqual(1046.5854470662418, a.ToMass(2));
        }


        [Test]
        public void MzToMassNegativeCharge()
        {
            double a = 524.3;
            Assert.AreEqual(1050.614552933758, a.ToMass(-2));
        }
    }

    internal class ObjectWithMass1000 : IHasMass
    {
        public double MonoisotopicMass
        {
            get
            {
                return 1000;
            }
        }
    }

    internal class ObjectWithMass100 : IHasMass
    {
        public double MonoisotopicMass
        {
            get
            {
                return 100;
            }
        }

        public override string ToString()
        {
            return "mass: 100";
        }
    }
}