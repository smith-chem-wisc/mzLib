using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Tests for the per-document sample verdict. Each fixture is a document every other SDRF
    /// instrument calls fine -- valid, no drift -- so what is under test is only whether the sample
    /// half says anything.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrfSampleInformativeness
    {
        private static readonly string[] Columns =
        {
            "source name", "characteristics[organism]", "characteristics[organism part]",
            "characteristics[biological replicate]", "assay name", "comment[data file]", "factor value[disease]"
        };

        private static SdrfDocument Doc(string[] columns, params string[][] rows)
        {
            var header = new SdrfHeader(columns);
            return new SdrfDocument(header, rows.Select(r => new SdrfRow(header, r)));
        }

        [Test]
        public void ADesignedStudyIsInformative()
        {
            var assessment = SdrfSampleInformativeness.Assess(Doc(Columns,
                new[] { "S1", "Homo sapiens", "liver", "1", "run 1", "a.raw", "normal" },
                new[] { "S2", "Homo sapiens", "liver", "2", "run 2", "b.raw", "hepatocellular carcinoma" }));

            Assert.That(assessment.Verdict, Is.EqualTo(SdrfSampleVerdict.Informative), assessment.ToString());
        }

        /// <summary>
        /// The shape a skeleton generator writes: one row per data file, every sample cell a reserved
        /// word, one replicate number throughout. Organism is filled -- a search can do that much
        /// unaided -- and must not rescue the verdict.
        /// </summary>
        [Test]
        public void ADataFileListWithReservedWordsIsASkeleton()
        {
            var assessment = SdrfSampleInformativeness.Assess(Doc(Columns,
                new[] { "GM1_a", "Homo sapiens", "not available", "1", "run 1", "GM1_a.raw", "not available" },
                new[] { "GM1_b", "Homo sapiens", "not available", "1", "run 2", "GM1_b.raw", "not available" },
                new[] { "GM8_a", "Homo sapiens", "Not Available", "1", "run 3", "GM8_a.raw", "not applicable" }));

            Assert.That(assessment.Verdict, Is.EqualTo(SdrfSampleVerdict.Skeleton), assessment.ToString());
            Assert.That(assessment.SampleCharacteristicColumns.Select(c => c.Column),
                Is.EqualTo(new[] { "characteristics[organism part]" }),
                "organism and biological replicate are not sample descriptions for this check");
        }

        /// <summary>
        /// A single-condition study legitimately has nothing to vary. It is Partial, not Skeleton:
        /// deciding whether that is good enough belongs to the caller.
        /// </summary>
        [Test]
        public void ASingleConditionStudyIsPartial()
        {
            var assessment = SdrfSampleInformativeness.Assess(Doc(Columns,
                new[] { "S1", "Homo sapiens", "liver", "1", "run 1", "a.raw", "normal" },
                new[] { "S1", "Homo sapiens", "liver", "1", "run 2", "b.raw", "normal" }));

            Assert.That(assessment.Verdict, Is.EqualTo(SdrfSampleVerdict.Partial));
            Assert.That(assessment.SampleIsDescribed, Is.True);
            Assert.That(assessment.FactorValueVaries, Is.False);
            Assert.That(assessment.BiologicalReplicateVaries, Is.False);
        }

        /// <summary>
        /// Every factor column counts, not only the first; and a factor that varies only by case or
        /// whitespace does not vary.
        /// </summary>
        [Test]
        public void EveryFactorColumnIsConsidered_AndCaseIsNotVariation()
        {
            string[] columns = { "source name", "characteristics[biological replicate]", "factor value[disease]", "Factor Value[age]" };

            var varies = SdrfSampleInformativeness.Assess(Doc(columns,
                new[] { "S1", "1", "normal", "58Y" },
                new[] { "S2", "1", "normal", "30Y" }));
            Assert.That(varies.FactorValueVaries, Is.True,
                "the second factor column varies, and its capitalised prefix is still a factor column");

            var caseOnly = SdrfSampleInformativeness.Assess(Doc(columns,
                new[] { "S1", "1", "Normal", "58Y" },
                new[] { "S2", "1", " normal ", "58Y" }));
            Assert.That(caseOnly.FactorValueVaries, Is.False);
        }

        [Test]
        public void ADocumentWithoutSampleColumnsIsASkeleton_AndHasNoReplicateColumn()
        {
            var assessment = SdrfSampleInformativeness.Assess(Doc(
                new[] { "source name", "assay name", "comment[data file]" },
                new[] { "S1", "run 1", "a.raw" }));

            Assert.That(assessment.Verdict, Is.EqualTo(SdrfSampleVerdict.Skeleton));
            Assert.That(assessment.BiologicalReplicate, Is.Null);
            Assert.That(assessment.FactorValueColumns, Is.Empty);
        }

        [Test]
        public void ANullDocumentIsRefused()
        {
            Assert.That(() => SdrfSampleInformativeness.Assess(null), Throws.ArgumentNullException);
        }

        /// <summary>
        /// Verdict counts over the curated corpus, for reference. [Explicit]; needs MZLIB_SDRF_CORPUS.
        /// </summary>
        [Test]
        [Explicit("Requires a local clone of bigbio/sdrf-annotated-datasets; set MZLIB_SDRF_CORPUS.")]
        public void CorpusVerdictReport()
        {
            string corpus = Environment.GetEnvironmentVariable("MZLIB_SDRF_CORPUS");
            if (string.IsNullOrWhiteSpace(corpus) || !Directory.Exists(corpus))
                Assert.Ignore($"MZLIB_SDRF_CORPUS not set or not found: '{corpus}'");

            var verdicts = Directory.GetFiles(corpus, "*.sdrf.tsv", SearchOption.AllDirectories)
                .Select(path =>
                {
                    var document = new SdrfDocument(path);
                    document.LoadResults();
                    return SdrfSampleInformativeness.Assess(document).Verdict;
                })
                .GroupBy(v => v)
                .OrderBy(g => g.Key)
                .ToList();

            foreach (var group in verdicts)
                TestContext.Progress.WriteLine($"{group.Key,-12} {group.Count()}");

            Assert.That(verdicts, Is.Not.Empty);
        }
    }
}
