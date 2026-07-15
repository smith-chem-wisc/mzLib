using System.IO;
using System.Linq;
using NUnit.Framework;

namespace Test.FileReadingTests.ProForma
{
    /// <summary>
    /// Phase 0 scaffolding: verifies the manuscript-example corpus loads and is complete.
    /// The corpus is the authoritative test set for every later ProForma phase; these tests
    /// guard its integrity, not the parser. Source enumeration: corpus/EXAMPLE-INDEX.md.
    /// </summary>
    [TestFixture]
    internal class CorpusLoaderTests
    {
        [Test]
        public void Corpus_File_Exists()
        {
            Assert.That(File.Exists(ProFormaTestCorpus.TsvPath), Is.True,
                "Corpus file missing. Run Phase 0 to generate corpus/manuscript-examples.tsv.");
        }

        [Test]
        public void Corpus_Loads_With_Unique_Ids()
        {
            var rows = ProFormaTestCorpus.Load();
            Assert.That(rows, Is.Not.Empty);
            Assert.That(rows.Select(r => r.Id).Distinct().Count(), Is.EqualTo(rows.Count), "duplicate ids present");
        }

        [Test]
        public void Corpus_Has_All_V1_Paper_Examples()
        {
            var v1 = ProFormaTestCorpus.Load().Where(r => r.Source == "v1-paper").ToList();
            // EXAMPLE-INDEX.md enumerates 26 v1-paper examples (Rules 1-7 + best-practices i-iv).
            // NOTE: phase-00-setup.md's stub asserted 27 — that was an off-by-one vs the index.
            Assert.That(v1, Has.Count.EqualTo(26), "v1-paper example count mismatch vs EXAMPLE-INDEX.md");
        }

        [Test]
        public void Corpus_Has_All_V2_Spec_Examples()
        {
            var v2Good = ProFormaTestCorpus.Load().Where(r => r.Source == "v2-spec" && r.Valid).ToList();
            // Exact count (sections 4.2.1-7.2 of manuscript-examples.tsv), matching the precise-count
            // style of the sibling corpus tests so a deleted/mistyped valid row is caught.
            Assert.That(v2Good, Has.Count.EqualTo(89), "v2-spec valid example count mismatch vs manuscript-examples.tsv");
        }

        [Test]
        public void Corpus_Has_All_V2_Spec_CounterExamples()
        {
            var bad = ProFormaTestCorpus.Load().Where(r => r.Id.StartsWith("v2-bad-")).ToList();
            Assert.That(bad, Has.Count.EqualTo(8));
            Assert.That(bad.All(r => !r.Valid), Is.True, "counter-examples must be flagged valid=false");
        }

        [Test]
        public void Corpus_ComplianceLevels_AreKnown()
        {
            var known = new[] { "base", "level2", "top-down", "cross-linking", "glycans", "mass-spectrum" };
            Assert.That(ProFormaTestCorpus.Load().All(r => known.Contains(r.ComplianceLevel)), Is.True,
                "unexpected compliance_level value in corpus");
        }
    }
}
