using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Readers.ProForma;

namespace Test.FileReadingTests.ProForma
{
    /// <summary>
    /// Round-trip for ProForma constructs above a single proteoform term — chimeric (+), charge (/z),
    /// and multi-chain / branch (//) — handled by the mzLib multi-term facade
    /// (<see cref="ProFormaMultiTermReader"/> / <see cref="ProFormaMultiTermWriter"/>). These are the 13
    /// corpus forms previously pinned as "not yet supported" by the single-term parser.
    /// </summary>
    [TestFixture]
    internal class ProFormaMultiTermRoundTripTests
    {
        private static readonly HashSet<string> MultiTermIds = new()
        {
            "v2-4.2.3.2-01", "v2-4.2.3.2-02", "v2-4.2.3.3-03", // inter-chain crosslink //
            "v2-4.2.4-01", "v2-4.2.4-02",                       // branch //
            "v2-7.1-01", "v2-7.1-02", "v2-7.1-03", "v2-7.1-04", // charge /z [+adducts]
            "v2-7.1-05", "v2-7.1-06", "v2-7.1-07",
            "v2-7.2-01",                                         // chimeric +
        };

        private static IEnumerable<TestCaseData> MultiTermExamples()
        {
            foreach (var r in ProFormaTestCorpus.Load().Where(r => MultiTermIds.Contains(r.Id)))
                yield return new TestCaseData(r.ProformaString).SetName(r.Id);
        }

        [TestCaseSource(nameof(MultiTermExamples))]
        public void MultiTerm_RoundTrip_IsCanonicallyIdempotent(string proForma)
        {
            string first = ProFormaMultiTermWriter.Write(ProFormaMultiTermReader.Read(proForma));
            string second = ProFormaMultiTermWriter.Write(ProFormaMultiTermReader.Read(first));
            Assert.That(second, Is.EqualTo(first), $"non-idempotent multi-term form for '{proForma}' (first = '{first}')");
        }

        [Test]
        public void Charge_IsParsedAndWritten()
        {
            var mt = ProFormaMultiTermReader.Read("EMEVEESPEK/2");
            Assert.That(mt.Peptidoforms, Has.Count.EqualTo(1));
            Assert.That(mt.Peptidoforms[0].Chains, Has.Count.EqualTo(1));
            Assert.That(mt.Peptidoforms[0].Charge, Is.EqualTo(2));
            Assert.That(ProFormaMultiTermWriter.Write(mt), Is.EqualTo("EMEVEESPEK/2"));
        }

        [Test]
        public void Charge_WithAdducts_RoundTrips()
        {
            var mt = ProFormaMultiTermReader.Read("EMEVEESPEK/2[+2Na+,+H+]");
            Assert.That(mt.Peptidoforms[0].Charge, Is.EqualTo(2));
            Assert.That(mt.Peptidoforms[0].IonAdducts, Is.EqualTo("+2Na+,+H+"));
            Assert.That(ProFormaMultiTermWriter.Write(mt), Is.EqualTo("EMEVEESPEK/2[+2Na+,+H+]"));
        }

        [Test]
        public void Chimeric_SplitsPeptidoforms()
        {
            var mt = ProFormaMultiTermReader.Read("EMEVEESPEK/2+ELVISLIVER/3");
            Assert.That(mt.Peptidoforms, Has.Count.EqualTo(2));
            Assert.That(mt.Peptidoforms[0].Charge, Is.EqualTo(2));
            Assert.That(mt.Peptidoforms[1].Charge, Is.EqualTo(3));
            // Verify the split by content, not just cardinality: a mis-split would pass the counts above.
            Assert.That(mt.Peptidoforms[0].Chains[0].Sequence, Is.EqualTo("EMEVEESPEK"));
            Assert.That(mt.Peptidoforms[1].Chains[0].Sequence, Is.EqualTo("ELVISLIVER"));
        }

        [Test]
        public void MultiChain_SplitsOnDoubleSlash()
        {
            const string input = "SEK[XLMOD:02001#XL1]UENCE//EMEVTK[XLMOD:02001#XL1]SESPEK";
            var mt = ProFormaMultiTermReader.Read(input);
            Assert.That(mt.Peptidoforms, Has.Count.EqualTo(1));
            Assert.That(mt.Peptidoforms[0].Chains, Has.Count.EqualTo(2));
            Assert.That(mt.Peptidoforms[0].Chains[0].Sequence, Is.EqualTo("SEKUENCE"));
            Assert.That(mt.Peptidoforms[0].Chains[1].Sequence, Is.EqualTo("EMEVTKSESPEK"));
            // Exact input preservation: a writer that dropped a crosslink label or swapped the two
            // chains around // would still be idempotent, so assert the first write equals the input.
            Assert.That(ProFormaMultiTermWriter.Write(mt), Is.EqualTo(input));
        }

        [Test]
        public void Branch_RoundTrip_PreservesInput()
        {
            const string input = "ETFGD[MOD:00093#BRANCH]//R[#BRANCH]ATER";
            var mt = ProFormaMultiTermReader.Read(input);
            Assert.That(mt.Peptidoforms[0].Chains, Has.Count.EqualTo(2));
            Assert.That(ProFormaMultiTermWriter.Write(mt), Is.EqualTo(input));
        }
    }
}
