// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (TestAminoAcids.cs) is part of Proteomics.
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
using Proteomics;
using System;

namespace Test
{
    [TestFixture]
    public sealed class TestAminoAcids
    {
        [Test]
        public void GetResidueNotInDictionary()
        {
            AminoAcid r;
            Assert.IsFalse(AminoAcid.TryGetResidue("?", out r));
            Assert.IsFalse(AminoAcid.TryGetResidue('?', out r));
        }

        [Test]
        public void GetResidueByCharacter()
        {
            AminoAcid aa = AminoAcid.GetResidue('A');

            Assert.AreEqual("Alanine", aa.Name);
        }

        [Test]
        public void GetNullResidue()
        {
            Assert.Throws<ArgumentNullException>(() => { AminoAcid.GetResidue(null); });
        }

        [Test]
        public void GetResidueByCharacterString()
        {
            AminoAcid aa = AminoAcid.GetResidue("A");

            Assert.AreEqual(aa.Name, "Alanine");
        }

        [Test]
        public void GetResidueByName()
        {
            AminoAcid aa = AminoAcid.GetResidue("Alanine");

            Assert.AreEqual("Alanine", aa.Name);
        }
    }
}