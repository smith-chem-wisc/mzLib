using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Readers.ProForma;

namespace Test.FileReadingTests.ProForma
{
    /// <summary>
    /// Phase 1: Layer-1 round-trip (string &lt;-&gt; ProFormaTerm), delegated to the wrapped
    /// TopDownProteomics parser/writer. Verifies parsing succeeds and that the writer is
    /// canonically idempotent — a practical proxy for AST equivalence:
    /// write(parse(s)) == write(parse(write(parse(s)))).
    /// Scoped to base-compliance v2 examples for Phase 1; later phases widen the filter.
    /// </summary>
    [TestFixture]
    internal class ProFormaLayer1RoundTripTests
    {
        /// <summary>
        /// Corpus ids the wrapped SDK cannot currently round-trip. Documented in
        /// deliverables/known-limitations.md; each has a dedicated KnownGap_* test below that
        /// pins the current SDK behavior so we are alerted if a future release fixes it.
        /// </summary>
        private static readonly HashSet<string> KnownSdkWriterGaps = new()
        {
            "v2-4.5-01", // spec 4.5: multiple modifications on one range — parses, writer throws
        };

        private static IEnumerable<TestCaseData> BaseExamples()
        {
            foreach (var r in ProFormaTestCorpus.Load()
                         .Where(r => r.Source == "v2-spec" && r.Valid && r.ComplianceLevel == "base"
                                     && !KnownSdkWriterGaps.Contains(r.Id)))
                yield return new TestCaseData(r.ProformaString).SetName($"Base_{r.Id}");
        }

        [TestCaseSource(nameof(BaseExamples))]
        public void RoundTrip_Base_IsCanonicallyIdempotent(string proForma)
        {
            string firstWrite = ProFormaWriter.Write(ProFormaReader.Read(proForma));
            string secondWrite = ProFormaWriter.Write(ProFormaReader.Read(firstWrite));
            Assert.That(secondWrite, Is.EqualTo(firstWrite),
                $"non-idempotent canonical form for '{proForma}' (first write = '{firstWrite}')");
        }

        /// <summary>
        /// Pins SDK gap: a range bearing multiple modifications (spec 4.5) parses, but
        /// TopDownProteomics' ProFormaWriter throws when serializing it. If a future SDK release
        /// fixes this, this test fails — remove v2-4.5-01 from <see cref="KnownSdkWriterGaps"/>.
        /// </summary>
        [Test]
        public void KnownGap_MultiModRange_ParsesButSdkWriterThrows()
        {
            var row = ProFormaTestCorpus.Load().Single(r => r.Id == "v2-4.5-01");
            var term = ProFormaReader.Read(row.ProformaString);
            Assert.That(term, Is.Not.Null, "SDK should still parse a multi-mod range");
            Assert.Throws<TopDownProteomics.ProForma.ProFormaParseException>(
                () => ProFormaWriter.Write(term),
                "SDK ProFormaWriter now serializes a multi-mod range — remove v2-4.5-01 from KnownSdkWriterGaps.");
        }
    }
}
