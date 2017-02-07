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

using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public sealed class TestProteinReader
    {

        #region Public Methods

        [Test]
        public void XmlTest()
        {
            HashSet<Modification> nice = new HashSet<Modification>
            {
                new Modification("fayk",null, null,ModificationSites.A)
            };

            IDictionary<string, HashSet<Modification>> mods = new Dictionary<string, HashSet<Modification>>
            {
                {
                    "N6-acetyllysine",nice
                }
            };

            ProteinDbLoader.LoadProteinDb(@"xml.xml", true, mods, false).ToList();
        }

        [Test]
        public void FastaTest()
        {
            ProteinDbLoader.LoadProteinDb(@"fasta.fasta", true, null, false).ToList();
        }

        #endregion Public Methods

    }
}