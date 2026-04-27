// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (FragmentInternalAndCopyTests.cs) is part of Proteomics.
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
using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Fragmentation.Peptide;
using Omics.Modifications;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class FragmentInternalAndCopyTests
    {
        private static Stopwatch Stopwatch { get; set; }

        [SetUp]
        public static void Setup()
        {
            Stopwatch = new Stopwatch();
            Stopwatch.Start();
        }

        [TearDown]
        public static void TearDown()
        {
            Console.WriteLine($"Analysis time: {Stopwatch.Elapsed.Hours}h {Stopwatch.Elapsed.Minutes}m {Stopwatch.Elapsed.Seconds}s");
        }

        [Test]
        public static void TestInternalFragments()
        {
            PeptideWithSetModifications pwsm = new PeptideWithSetModifications("PEPTIDE", null);
            List<Product> products = new List<Product>();

            //test with HCD
            pwsm.FragmentInternally(DissociationType.HCD, 3, products);


            List<Product> expectedProducts = new List<Product>
            {
                new Product(ProductType.y, FragmentationTerminus.None,327.14,2,3,0,ProductType.b,4), //EPT
                new Product(ProductType.y, FragmentationTerminus.None,440.23,2,4,0,ProductType.b,5), //EPTI
                new Product(ProductType.y, FragmentationTerminus.None,555.25,2,5,0,ProductType.b,6), //EPTID
                new Product(ProductType.y, FragmentationTerminus.None,311.18,3,3,0,ProductType.b,5), //PTI
                new Product(ProductType.y, FragmentationTerminus.None,426.21,3,4,0,ProductType.b,6), //PTID
                new Product(ProductType.y, FragmentationTerminus.None,329.16,4,3,0,ProductType.b,6), //TID
            };
            Assert.IsTrue(products.Count == expectedProducts.Count);
            for (int i = 0; i < products.Count; i++)
            {
                Assert.IsTrue(products[i].IsInternalFragment);
                Assert.IsTrue(products[i].Annotation.Equals(expectedProducts[i].Annotation));
                Assert.IsTrue(Math.Round(products[i].NeutralMass).Equals(Math.Round(expectedProducts[i].NeutralMass)));
            }

            //test with multiple different fragments (EThcD)
            pwsm.FragmentInternally(DissociationType.EThcD, 5, products);
            expectedProducts = new List<Product>
            {
                new Product(ProductType.y, FragmentationTerminus.None,555.25,2,5,0,ProductType.b,6), //EPTID, by
                new Product(ProductType.zDot, FragmentationTerminus.None,539.24,2,5,0,ProductType.b,6), //EPTID, bz
                new Product(ProductType.y, FragmentationTerminus.None,572.28,2,5,0,ProductType.c,6), //EPTID, cy
                new Product(ProductType.zDot, FragmentationTerminus.None,556.26,2,5,0,ProductType.c,6), //EPTID, cz
            };
            Assert.IsTrue(products.Count == expectedProducts.Count);
            for (int i = 0; i < products.Count; i++)
            {
                Assert.IsTrue(products[i].IsInternalFragment);
                Assert.IsTrue(products[i].Annotation.Equals(expectedProducts[i].Annotation));
                Assert.IsTrue(Math.Round(products[i].NeutralMass).Equals(Math.Round(expectedProducts[i].NeutralMass)));
            }

            //test string includes both termini
            Assert.IsTrue(products[0].ToString().Equals("yIb[2-6];555.25404-0"));

            //test that mods are incorporated
            ModificationMotif.TryGetMotif("T", out ModificationMotif target);
            Modification oxOnM = new Modification(_originalId: "Oxidation on M", _modificationType: "Common Variable", _target: target, _chemicalFormula: ChemicalFormula.ParseFormula("O"));
            pwsm = new PeptideWithSetModifications("PM[Common Variable:Oxidation on M]EPTIM[Common Variable:Oxidation on M]DE", new Dictionary<string, Modification> { { "Oxidation on M", oxOnM } });
            pwsm.FragmentInternally(DissociationType.EThcD, 5, products);
            Assert.IsTrue(Math.Round(products[0].NeutralMass).Equals(587)); //contains first ox M
            Assert.IsTrue(Math.Round(products[8].NeutralMass).Equals(849)); //contains both ox M

            //test that noncanonical amino acids are handled
            pwsm = new PeptideWithSetModifications("ZXCVZXCV", null);
            pwsm.FragmentInternally(DissociationType.ETD, 5, products); //shouldn't crash
            Assert.IsTrue(products[0].NeutralMass.Equals(double.NaN));
        }
    }
}
