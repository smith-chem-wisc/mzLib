using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Tests for <see cref="SdrfSample.Comments"/>: extension <c>comment[...]</c> columns written in the
    /// comment block AFTER the assay columns. The provenance columns (sdrf D29/D31) are the first user;
    /// smuggling them through the characteristics dictionaries would place them before <c>assay name</c>,
    /// which the reference validator rejects (mapping note N10).
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrfBuilderComments
    {
        private static SdrfRowInput Row(string source, string file, Dictionary<string, string> comments) => new(
            new SdrfSample { SourceName = source, Comments = comments },
            new SdrfAssay { DataFileName = file, AssayName = "run " + source });

        private static readonly SdrfBuilderOptions Lenient = new() { RequireSampleMetadata = false };

        [Test]
        public void CommentsAreWrittenAfterTheAssayColumnsAndBeforeTheFactors()
        {
            var rows = new[]
            {
                Row("s1", "a.raw", new() { ["comment[characteristics source]"] = "inferred" }) with
                    { Sample = new SdrfSample { SourceName = "s1", FactorValue = "WT", FactorValueColumn = "factor value[genotype]",
                      Comments = new Dictionary<string, string> { ["comment[characteristics source]"] = "inferred" } } },
            };

            var doc = SdrfBuilder.Build(rows, Lenient);

            var h = doc.Header.ToList();
            int comment = h.IndexOf("comment[characteristics source]");
            Assert.That(comment, Is.GreaterThan(h.IndexOf("assay name")));
            Assert.That(comment, Is.GreaterThan(h.IndexOf("comment[data file]")));
            Assert.That(comment, Is.LessThan(h.IndexOf("factor value[genotype]")));
            Assert.That(doc.Results[0]["comment[characteristics source]"], Is.EqualTo("inferred"));
        }

        [Test]
        public void EveryRowsKeysJoinOneUnionAndARowWithoutAKeyIsNotApplicable()
        {
            var rows = new[]
            {
                Row("s1", "a.raw", new() { ["comment[characteristics source]"] = "pride project record" }),
                Row("s2", "b.raw", new() { ["comment[characteristics source]"] = "pride project record", ["comment[disease source]"] = "inferred" }),
            };

            var doc = SdrfBuilder.Build(rows, Lenient);

            Assert.That(doc.Header.Count(h => h == "comment[disease source]"), Is.EqualTo(1));
            Assert.That(doc.Results.All(r => r.Header.Count == doc.Header.Count), "never ragged");
            Assert.That(doc.Results[0]["comment[disease source]"], Is.EqualTo("not applicable"), "no override on this row");
            Assert.That(doc.Results[1]["comment[disease source]"], Is.EqualTo("inferred"));
        }

        [Test]
        public void ACommentKeyMustBeACommentColumn()
        {
            var rows = new[] { Row("s1", "a.raw", new() { ["characteristics[age]"] = "40Y" }) };

            var e = Assert.Throws<ArgumentException>(() => SdrfBuilder.Build(rows, Lenient));
            Assert.That(e!.Message, Does.Contain("comment["));
        }

        [Test]
        public void ACommentThatCollidesWithABuiltInColumnThrows()
        {
            var rows = new[] { Row("s1", "a.raw", new() { ["comment[data file]"] = "x.raw" }) };

            Assert.Throws<ArgumentException>(() => SdrfBuilder.Build(rows, Lenient));
        }

        [Test]
        public void NoCommentsMeansTheSameDocumentAsBefore()
        {
            var plain = SdrfBuilder.Build(new[] { Row("s1", "a.raw", new()) }, Lenient);

            Assert.That(plain.Header.Any(h => h.StartsWith("comment[characteristics", StringComparison.Ordinal)), Is.False);
        }
    }
}
