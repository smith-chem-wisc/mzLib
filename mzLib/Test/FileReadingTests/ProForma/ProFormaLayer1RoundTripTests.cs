using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Readers.ProForma;
using Tdp = TopDownProteomics.ProForma;

namespace Test.FileReadingTests.ProForma
{
    /// <summary>
    /// Phase 1-3: Layer-1 round-trip (string &lt;-&gt; ProFormaTerm), delegated to the wrapped
    /// TopDownProteomics parser/writer. Verifies parsing succeeds and that the writer is
    /// canonically idempotent — a practical proxy for AST equivalence:
    /// write(parse(s)) == write(parse(write(parse(s)))).
    /// Covers every valid corpus example except the documented exclusion sets below.
    /// The SDK fully covers a single proteoform term; the remaining gaps are documented in
    /// deliverables/known-limitations.md.
    /// </summary>
    [TestFixture]
    internal class ProFormaLayer1RoundTripTests
    {
        /// <summary>
        /// v1.0-only notations the v2-only SDK does not parse. Intentionally unsupported
        /// (v2.0 removed them; user decision 2026-05-24). Asserted by
        /// <see cref="V1_DefaultSourcePrefix_IsUnsupported"/>.
        /// </summary>
        private static readonly HashSet<string> KnownUnsupportedV1 = new()
        {
            // v1.0 Rule 6 / best-practice ii: "default-source" prefix [SOURCE]+SEQ, removed in v2.0
            "v1-rule6-01", "v1-rule6-02", "v1-rule6-03", "v1-bp-ii",
        };

        /// <summary>
        /// Valid v2.0 features that live above a single ProFormaTerm and so are not handled by the
        /// SDK's single-term ParseString: charge/adducts (/z), chimeric (+), and multi-chain (//).
        /// Supporting these requires facade-level splitting/charge handling (planned, README
        /// Phases 7-8). Pinned by <see cref="MultiTerm_ChargeChimericMultichain_NotYetSupported"/>.
        /// </summary>
        private static readonly HashSet<string> RequiresMultiTermFacade = new()
        {
            "v2-4.2.3.2-01", "v2-4.2.3.2-02", "v2-4.2.3.3-03", // inter-chain // (crosslink)
            "v2-4.2.4-01", "v2-4.2.4-02",                       // branch //
            "v2-7.1-01", "v2-7.1-02", "v2-7.1-03", "v2-7.1-04", // charge /z [+adducts]
            "v2-7.1-05", "v2-7.1-06", "v2-7.1-07",
            "v2-7.2-01",                                         // chimeric +
        };

        /// <summary>
        /// Examples the SDK parses but cannot write back. Each is pinned by a KnownGap_* test so a
        /// future TopDownProteomics release that fixes it trips an alert. See known-limitations.md.
        /// </summary>
        private static readonly HashSet<string> KnownSdkWriterGaps = new()
        {
            "v2-4.5-01", // spec 4.5: multiple modifications on one range — parses, writer throws
        };

        private static IEnumerable<TestCaseData> RoundTripExamples()
        {
            foreach (var r in ProFormaTestCorpus.Load()
                         .Where(r => r.Valid
                                     && !KnownUnsupportedV1.Contains(r.Id)
                                     && !RequiresMultiTermFacade.Contains(r.Id)
                                     && !KnownSdkWriterGaps.Contains(r.Id)))
                yield return new TestCaseData(r.ProformaString).SetName($"{r.ComplianceLevel}_{r.Id}");
        }

        [TestCaseSource(nameof(RoundTripExamples))]
        public void RoundTrip_IsCanonicallyIdempotent(string proForma)
        {
            string firstWrite = ProFormaWriter.Write(ProFormaReader.Read(proForma));
            string secondWrite = ProFormaWriter.Write(ProFormaReader.Read(firstWrite));
            Assert.That(secondWrite, Is.EqualTo(firstWrite),
                $"non-idempotent canonical form for '{proForma}' (first write = '{firstWrite}')");
        }

        /// <summary>
        /// Pins gap #2: the legacy v1.0 default-source prefix `[SOURCE]+SEQ` (dropped in v2.0) is
        /// not parsed by the v2-only SDK. Intentionally unsupported.
        /// </summary>
        [Test]
        public void V1_DefaultSourcePrefix_IsUnsupported()
        {
            var rows = ProFormaTestCorpus.Load().Where(r => KnownUnsupportedV1.Contains(r.Id)).ToList();
            Assert.That(rows, Has.Count.EqualTo(KnownUnsupportedV1.Count), "missing v1 default-source rows");
            foreach (var r in rows)
                Assert.Throws<Tdp.ProFormaParseException>(() => ProFormaReader.Read(r.ProformaString),
                    $"{r.Id} now parses — v1 default-source may be supported; revisit KnownUnsupportedV1.");
        }

        /// <summary>
        /// Pins gap #3: charge (/z), chimeric (+), and multi-chain (//) constructs are not handled
        /// by the SDK's single-term ParseString (it throws "/ is not an upper case letter"). When
        /// the multi-term facade is built (Phases 7-8), these move out of <see cref="RequiresMultiTermFacade"/>.
        /// </summary>
        [Test]
        public void MultiTerm_ChargeChimericMultichain_NotYetSupported()
        {
            var rows = ProFormaTestCorpus.Load().Where(r => RequiresMultiTermFacade.Contains(r.Id)).ToList();
            Assert.That(rows, Has.Count.EqualTo(RequiresMultiTermFacade.Count), "missing multi-term rows");
            foreach (var r in rows)
                Assert.Throws<Tdp.ProFormaParseException>(() => ProFormaReader.Read(r.ProformaString),
                    $"{r.Id} now parses via single-term Read — build the multi-term facade and move it out.");
        }

        /// <summary>
        /// Pins gap #1: a range bearing multiple modifications (spec 4.5) parses, but
        /// TopDownProteomics' ProFormaWriter throws when serializing it. If a future SDK release
        /// fixes this, this test fails — remove v2-4.5-01 from <see cref="KnownSdkWriterGaps"/>.
        /// </summary>
        [Test]
        public void KnownGap_MultiModRange_ParsesButSdkWriterThrows()
        {
            var row = ProFormaTestCorpus.Load().Single(r => r.Id == "v2-4.5-01");
            var term = ProFormaReader.Read(row.ProformaString);
            Assert.That(term, Is.Not.Null, "SDK should still parse a multi-mod range");
            Assert.Throws<Tdp.ProFormaParseException>(() => ProFormaWriter.Write(term),
                "SDK ProFormaWriter now serializes a multi-mod range — remove v2-4.5-01 from KnownSdkWriterGaps.");
        }
    }
}
