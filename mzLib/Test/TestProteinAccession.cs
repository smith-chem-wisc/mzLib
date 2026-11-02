using System.Diagnostics.CodeAnalysis;
using MzLibUtil;
using NUnit.Framework;

namespace Test
{
    /// <summary>
    /// Parsing an accession into what can be said about it without reference data. Parse, never
    /// repair: a prefixed, lower-cased or padded id is a data problem to surface, not one to fix.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestProteinAccession
    {
        [Test]
        public void UniProtIsoformIsSplitFromTheEntry_NotEquatedWithTheCanonicalSequence()
        {
            // "-1" is an isoform number. UniProt's displayed (canonical) sequence can carry any number,
            // so EntryAccession means "same UniProt entry" and says nothing about sequence identity.
            var a = "P12345-1".ParseProteinAccession();

            Assert.That(a.Verbatim, Is.EqualTo("P12345-1"));
            Assert.That(a.EntryAccession, Is.EqualTo("P12345"));
            Assert.That(a.Isoform, Is.EqualTo(1));
            Assert.That(a.Namespace, Is.EqualTo(AccessionNamespace.UniProt));
            Assert.That("P12345".ParseProteinAccession().Isoform, Is.Null);
        }

        [Test]
        public void TenCharacterUniProtAccessionsParse()
        {
            var a = "A0A087X1C5".ParseProteinAccession();
            Assert.That((a.Namespace, a.EntryAccession), Is.EqualTo((AccessionNamespace.UniProt, "A0A087X1C5")));
        }

        [Test]
        public void RefSeqVersionIsSplitOff()
        {
            var a = "NP_000537.3".ParseProteinAccession();

            Assert.That(a.EntryAccession, Is.EqualTo("NP_000537"));
            Assert.That(a.Version, Is.EqualTo(3));
            Assert.That(a.Namespace, Is.EqualTo(AccessionNamespace.RefSeq));
        }

        [TestCase("CON__P02768")]
        [TestCase("p12345")]
        [TestCase(" P12345")]
        [TestCase("sp|P12345|X_HUMAN")]
        [TestCase("DECOY_P12345")]
        [TestCase("")]
        public void AnythingOutsideTheGrammarIsUnrecognized_AndKeptVerbatim(string accession)
        {
            var a = accession.ParseProteinAccession();

            Assert.That(a.Namespace, Is.EqualTo(AccessionNamespace.Unrecognized));
            Assert.That(a.EntryAccession, Is.EqualTo(accession));
            Assert.That(a.Isoform, Is.Null);
            Assert.That(a.Version, Is.Null);
        }

        [Test]
        public void NullParsesAsUnrecognizedEmpty()
        {
            var a = ((string)null).ParseProteinAccession();
            Assert.That((a.Namespace, a.Verbatim), Is.EqualTo((AccessionNamespace.Unrecognized, "")));
        }
    }
}
