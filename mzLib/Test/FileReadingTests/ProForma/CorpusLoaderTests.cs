using System.Collections.Generic;
using System.IO;
using System.Linq;
using NUnit.Framework;

namespace Test.FileReadingTests.ProForma
{
    /// <summary>
    /// Phase 0 scaffolding: verifies the manuscript-example corpus loads and is complete.
    /// The corpus is the authoritative test set for every later ProForma phase; these tests
    /// guard its integrity, not the (not-yet-written) parser.
    /// Source enumeration: corpus/EXAMPLE-INDEX.md.
    /// </summary>
    [TestFixture]
    internal class CorpusLoaderTests
    {
        private static string CorpusPath => Path.Combine(TestContext.CurrentContext.TestDirectory,
            @"FileReadingTests\ProForma\manuscript-examples.tsv");

        private record CorpusRow(string Id, string ProformaString, string Source, string Section,
            string ComplianceLevel, string Feature, string CanonicalForm, bool Valid, string Notes);

        private static List<CorpusRow> LoadCorpus()
        {
            var rows = new List<CorpusRow>();
            var lines = File.ReadAllLines(CorpusPath);
            for (int i = 1; i < lines.Length; i++) // skip header
            {
                if (lines[i].Length == 0) continue;
                var c = lines[i].Split('\t');
                Assert.That(c.Length, Is.EqualTo(9), $"Row {i + 1} does not have 9 tab-separated columns");
                rows.Add(new CorpusRow(c[0], c[1], c[2], c[3], c[4], c[5], c[6], bool.Parse(c[7]), c[8]));
            }
            return rows;
        }

        [Test]
        public void Corpus_File_Exists()
        {
            Assert.That(File.Exists(CorpusPath), Is.True,
                "Corpus file missing. Run Phase 0 to generate corpus/manuscript-examples.tsv.");
        }

        [Test]
        public void Corpus_Loads_With_Unique_Ids()
        {
            var rows = LoadCorpus();
            Assert.That(rows, Is.Not.Empty);
            Assert.That(rows.Select(r => r.Id).Distinct().Count(), Is.EqualTo(rows.Count), "duplicate ids present");
        }

        [Test]
        public void Corpus_Has_All_V1_Paper_Examples()
        {
            var rows = LoadCorpus();
            var v1 = rows.Where(r => r.Source == "v1-paper").ToList();
            // EXAMPLE-INDEX.md enumerates 26 v1-paper examples (Rules 1-7 + best-practices i-iv).
            // NOTE: phase-00-setup.md's stub asserted 27 — that was an off-by-one vs the index.
            Assert.That(v1, Has.Count.EqualTo(26), "v1-paper example count mismatch vs EXAMPLE-INDEX.md");
        }

        [Test]
        public void Corpus_Has_All_V2_Spec_Examples()
        {
            var rows = LoadCorpus();
            var v2Good = rows.Where(r => r.Source == "v2-spec" && r.Valid).ToList();
            Assert.That(v2Good.Count, Is.GreaterThan(70));
        }

        [Test]
        public void Corpus_Has_All_V2_Spec_CounterExamples()
        {
            var rows = LoadCorpus();
            var bad = rows.Where(r => r.Id.StartsWith("v2-bad-")).ToList();
            Assert.That(bad, Has.Count.EqualTo(8));
            Assert.That(bad.All(r => !r.Valid), Is.True, "counter-examples must be flagged valid=false");
        }

        [Test]
        public void Corpus_ComplianceLevels_AreKnown()
        {
            var known = new[] { "base", "level2", "top-down", "cross-linking", "glycans", "mass-spectrum" };
            var rows = LoadCorpus();
            Assert.That(rows.All(r => known.Contains(r.ComplianceLevel)), Is.True,
                "unexpected compliance_level value in corpus");
        }
    }
}
