using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Omics.BioPolymerGroup;

namespace Test.Omics.BioPolymerGroupTests
{
    /// <summary>
    /// Tests for the schema-driven writer itself, independent of any particular record type.
    ///
    /// The property the writer exists to guarantee is that a header and every row are the same width,
    /// because one schema produces both. These pin that on the writer's own terms -- a record that
    /// has no value for a column contributes an empty field rather than omitting it, which is what
    /// stops a ragged dataset shifting every later field on the affected rows.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TsvWriterTests
    {
        private sealed record Row(string Name, string? Value);

        private static IReadOnlyList<TsvColumn<Row>> Schema() =>
        [
            new TsvColumn<Row>("Name", r => r.Name),
            new TsvColumn<Row>("Value", r => r.Value!)
        ];

        /// <summary>
        /// Splits written output into lines, dropping only the break the final WriteLine leaves
        /// behind. Deliberately not TrimEnd: a row whose last column is empty ends in a tab, and
        /// trimming it away would delete the very field these tests exist to prove is present.
        /// </summary>
        private static string[] LinesOf(string written)
        {
            var lines = written.Split(Environment.NewLine, StringSplitOptions.None);
            return lines.Length > 0 && lines[^1].Length == 0 ? lines[..^1] : lines;
        }

        /// <summary>
        /// The whole-file entry point: one header, then one line per record, in the order given.
        /// </summary>
        [Test]
        public void WriteEmitsTheHeaderThenOneLinePerItem()
        {
            using var output = new StringWriter();

            TsvWriter.Write(output, Schema(), new[] { new Row("first", "1"), new Row("second", "2") });

            var lines = LinesOf(output.ToString());

            Assert.Multiple(() =>
            {
                Assert.That(lines, Has.Length.EqualTo(3), "header plus one line per record");
                Assert.That(lines[0], Is.EqualTo("Name\tValue"));
                Assert.That(lines[1], Is.EqualTo("first\t1"));
                Assert.That(lines[2], Is.EqualTo("second\t2"), "records keep the order they were given in");
            });
        }

        /// <summary>
        /// A file with no records is still a readable table: the header describes the columns that
        /// would have been there.
        /// </summary>
        [Test]
        public void WriteEmitsTheHeaderEvenWithNoItems()
        {
            using var output = new StringWriter();

            TsvWriter.Write(output, Schema(), Array.Empty<Row>());

            Assert.That(LinesOf(output.ToString()), Is.EqualTo(new[] { "Name\tValue" }));
        }

        /// <summary>
        /// Every row is exactly as wide as the header because both are rendered from the same schema.
        /// This is the ragged-row defect stated as the writer's own contract.
        /// </summary>
        [Test]
        public void EveryRowIsAsWideAsTheHeader()
        {
            var schema = Schema();
            using var output = new StringWriter();

            TsvWriter.Write(output, schema, new[] { new Row("a", "1"), new Row("b", null) });

            var lines = LinesOf(output.ToString());
            int headerWidth = lines[0].Split('\t').Length;

            Assert.That(lines.Skip(1).Select(l => l.Split('\t').Length),
                Has.All.EqualTo(headerWidth));
        }

        /// <summary>
        /// A column that cannot read a value contributes an empty field. Omitting it instead would
        /// shift every field after it on that row alone, which is the failure the schema removes.
        /// </summary>
        [Test]
        public void ANullValueBecomesAnEmptyFieldRatherThanAMissingOne()
        {
            var line = TsvWriter.RowLine(Schema(), new Row("a", null));

            Assert.Multiple(() =>
            {
                Assert.That(line, Is.EqualTo("a\t"));
                Assert.That(line.Split('\t'), Has.Length.EqualTo(2), "the field is present and empty");
            });
        }

        /// <summary>
        /// A column with no header or no accessor cannot render anything, and a schema that silently
        /// accepted one would fail later, per row, with nothing pointing at the column that caused it.
        /// </summary>
        [Test]
        public void AColumnRejectsAMissingHeaderOrAccessor()
        {
            Assert.Multiple(() =>
            {
                Assert.That(() => new TsvColumn<Row>(null!, r => r.Name),
                    Throws.ArgumentNullException.With.Property("ParamName").EqualTo("header"));
                Assert.That(() => new TsvColumn<Row>("Name", null!),
                    Throws.ArgumentNullException.With.Property("ParamName").EqualTo("getValue"));
            });
        }
    }
}
