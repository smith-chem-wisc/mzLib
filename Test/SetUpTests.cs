// Copyright 2016 Stefan Solntsev
// 
// This file (SetUpTests.cs) is part of Proteomics.
// 
// Proteomics is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// Proteomics is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with Proteomics. If not, see <http://www.gnu.org/licenses/>.

using Chemistry;
using NUnit.Framework;

namespace Test
{

    [SetUpFixture]
    public class MySetUpClass
    {
        [OneTimeSetUp]
        public void RunBeforeAnyTests()
        {
            var elementH = new Element("H", 1, 1.007975);
            PeriodicTable.Add(elementH);
            elementH.AddIsotope(1, 1.00782503223, 0.999885);
            elementH.AddIsotope(2, 2.01410177812, 0.000115);

            var elementC = new Element("C", 6, 12.0106);
            PeriodicTable.Add(elementC);
            elementC.AddIsotope(12, 12, 0.9893);
            elementC.AddIsotope(13, 13.00335483507, 0.0107);


            var elementN = new Element("N", 7, 14.006855);
            PeriodicTable.Add(elementN);
            elementN.AddIsotope(14, 14.00307400443, 0.99636);
            elementN.AddIsotope(15, 15.00010889888, 0.00364);


            var elementO = new Element("O", 8, 15.9994);
            PeriodicTable.Add(elementO);
            elementO.AddIsotope(16, 15.99491461957, 0.99757);
            elementO.AddIsotope(17, 16.99913175650, 0.00038);
            elementO.AddIsotope(18, 17.99915961286, 0.00205);

            var elementFe = new Element("Fe", 26, 55.845);
            PeriodicTable.Add(elementFe);
            elementFe.AddIsotope(56, 55.93493633, 0.91754);

            var elementBr = new Element("Br", 35, 79.904);
            PeriodicTable.Add(elementBr);
            elementBr.AddIsotope(79, 78.9183376, 0.5069);

            var elementCa = new Element("Ca", 20, 40.078);
            PeriodicTable.Add(elementCa);
            elementCa.AddIsotope(40, 39.962590863, 0.96941);

            var elementS = new Element("S", 16, 32.0675);
            PeriodicTable.Add(elementS);
            elementS.AddIsotope(32, 31.9720711744, 0.9499);
            elementS.AddIsotope(33, 32.9714589098, 0.0075);
            elementS.AddIsotope(34, 33.967867004, 0.0425);
            elementS.AddIsotope(36, 35.96708071, 0.0001);


            var elementSe = new Element("Se", 34, 78.971);
            PeriodicTable.Add(elementSe);
            elementSe.AddIsotope(74, 73.922475934, 0.0089);


            var elementAl = new Element("Al", 13, 26.9815385);
            PeriodicTable.Add(elementAl);
            elementAl.AddIsotope(27, 26.98153853, 1);
        }
    }
}
