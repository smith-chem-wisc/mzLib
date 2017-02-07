// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (ChemicalFormulaTestFixture.cs) is part of Chemistry Library.
//
// Chemistry Library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Chemistry Library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with Chemistry Library. If not, see <http://www.gnu.org/licenses/>

using Chemistry;
using NUnit.Framework;

namespace Test
{
    [TestFixture]
    public class ElementsAndIsotopesTest
    {
        #region Public Methods

        [Test]
        public void AddIsotopeWithExistingMassNumber()
        {
            var elementC = new Element("C", 6, 12.0106);
            elementC.AddIsotope(12, 12, 0.9893);
            elementC.AddIsotope(13, 13.00335483507, 0.0107);
            Isotope isotope = elementC[13];
            Assert.AreEqual("C{13}", isotope.ToString());
            Assert.AreEqual(6, isotope.Protons);
            Assert.AreEqual(7, isotope.Neutrons);
        }

        [Test]
        public void AddingExistingElementsTest()
        {
            var elementC = new Element("GGG", 127, 12.0106);
            PeriodicTable.Add(elementC);
        }

        #endregion Public Methods
    }
}