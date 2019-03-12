using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using NUnit.Framework;
using Proteomics.ProteolyticDigestion;

namespace Test
{
    [TestFixture]
    public static class TestMoreProtease
    {
        [Test]
        public static void SomeAdditionalProteaseTests()
        {
            Protease myProtease = new Protease("name", CleavageSpecificity.None, "psiAccessionNumber", "psiName", new List<DigestionMotif>());

            Assert.AreEqual(Proteomics.ProteolyticDigestion.CleavageSpecificity.None, myProtease.GetCleavageSpecificity("SEQUENCE",0,7));


        }
    }
}
