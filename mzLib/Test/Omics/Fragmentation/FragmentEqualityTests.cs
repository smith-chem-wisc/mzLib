// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (FragmentEqualityTests.cs) is part of Proteomics.
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
using System;
using System.Diagnostics;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class FragmentEqualityTests
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
        public static void TestFragmentAnnotations()
        {
            Product p = new Product(ProductType.b, FragmentationTerminus.N, 505.505, 2, 3, 30.3);
            MatchedFragmentIon f = new MatchedFragmentIon(p, 400.0, 1000.0, 3);

            Assert.That(p.Annotation == "b2-30.30");
            Assert.That(f.Annotation == "(b2-30.30)+3");

            p = new Product(ProductType.b, FragmentationTerminus.N, 505.505, 2, 3, 0);
            f = new MatchedFragmentIon(p, 400.0, 1000.0, 3);

            Assert.That(p.Annotation == "b2");
            Assert.That(f.Annotation == "b2+3");
        }

        [Test]
        public static void TestFragmentErrors()
        {
            Product p = new Product(ProductType.b, FragmentationTerminus.N, 475.205, 2, 3, 30.3);
            MatchedFragmentIon f = new MatchedFragmentIon(p, 159.5, 1000.0, 3);

            double experMass = f.Mz.ToMass(f.Charge);
            double theorMass = p.NeutralMass;

            Assert.AreEqual(0.2732, Math.Round(experMass - theorMass, 4));
            Assert.AreEqual(574.85, Math.Round(f.MassErrorPpm, 2));
            Assert.AreEqual(0.2732, Math.Round(f.MassErrorDa, 4));
        }

        [Test]
        public static void TestFragmentEquality()
        {
            Product p1 = new Product(ProductType.b, FragmentationTerminus.N, 505.505, 2, 3, 30.3);
            MatchedFragmentIon f1 = new MatchedFragmentIon(p1, 400.0, 1000.0, 3);

            Product p2 = new Product(ProductType.b, FragmentationTerminus.N, 505.505, 2, 3, 30.3);
            MatchedFragmentIon f2 = new MatchedFragmentIon(p2, 400.0, 1000.0, 3);

            Assert.AreEqual(p1, p2);
            Assert.AreEqual(f1, f2);
        }
    }
}
