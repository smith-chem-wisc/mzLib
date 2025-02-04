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

using NUnit.Framework;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Threading;
using UsefulProteomicsDatabases;

namespace Test
{
    [SetUpFixture]
    [ExcludeFromCodeCoverage]
    public class MySetUpClass
    {
        [OneTimeSetUp]
        public void Setup()
        {
            Assert.That(Thread.CurrentThread.CurrentCulture == CultureInfo.InvariantCulture);
            Assert.That(Thread.CurrentThread.CurrentUICulture == CultureInfo.InvariantCulture);
        }
    }
}