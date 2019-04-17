using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;

namespace Test
{
    [TestFixture]
    public static class TestMoreProtease
    {
        [Test]
        public static void SomeAdditionalProteaseTests()
        {
            Protease myProtease = new Protease("name", CleavageSpecificity.None, "psiAccessionNumber", "psiName", new List<DigestionMotif>());

            Assert.AreEqual(Proteomics.ProteolyticDigestion.CleavageSpecificity.None, myProtease.GetCleavageSpecificity("SEQUENCE", 0, 7));
        }
    }
}