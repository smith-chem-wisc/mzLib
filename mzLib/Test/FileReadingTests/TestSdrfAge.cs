using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Tests for reading characteristics[age] cells into years. The cases are the shapes the curated
    /// corpus actually contains; the refusals matter as much as the reads.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrfAge
    {
        private const double Tolerance = 1e-9;

        [TestCase("58Y", 58d)]
        [TestCase("30Y6M", 30.5d)]
        [TestCase("1Y6M15D", 1.5d + 15d / 365.25)]
        [TestCase("6M", 0.5d)]
        [TestCase("16W", 16d * 7d / 365.25)]
        [TestCase("3D", 3d / 365.25)]
        [TestCase("2.5Y", 2.5d)]
        [TestCase(" 54y ", 54d)]
        public void TheSpecificationsGrammarIsReadExactly(string cell, double years)
        {
            Assert.That(SdrfAge.TryParse(cell, out var age), Is.True);
            Assert.That(age!.Years, Is.EqualTo(years).Within(Tolerance));
            Assert.That(age.Precision, Is.EqualTo(SdrfAgePrecision.Exact));
            Assert.That(age.FollowsSpecification, Is.True);
            Assert.That(age.MinYears, Is.EqualTo(age.MaxYears));
        }

        /// <summary>
        /// Words the spec would not write but whose unit is not in doubt. Read, and marked as not
        /// following the specification, so a caller that trusts only the grammar can drop them.
        /// </summary>
        [TestCase("3 year", 3d)]
        [TestCase("51 years", 51d)]
        [TestCase("12 weeks old", 12d * 7d / 365.25)]
        [TestCase("1 week", 7d / 365.25)]
        [TestCase("4 hour", 4d / (365.25 * 24d))]
        [TestCase("2 Months", 2d / 12d)]
        public void UnambiguousWordsAreReadButNotCalledSpecification(string cell, double years)
        {
            Assert.That(SdrfAge.TryParse(cell, out var age), Is.True);
            Assert.That(age!.Years, Is.EqualTo(years).Within(Tolerance));
            Assert.That(age.Precision, Is.EqualTo(SdrfAgePrecision.Exact));
            Assert.That(age.FollowsSpecification, Is.False);
        }

        [TestCase("40Y-85Y", 40d, 85d, true)]
        [TestCase("8W-12W", 8d * 7d / 365.25, 12d * 7d / 365.25, true)]
        [TestCase("6-8 weeks", 6d * 7d / 365.25, 8d * 7d / 365.25, false)]
        [TestCase("0-2 hour", 0d, 2d / (365.25 * 24d), false)]
        [TestCase("40-85Y", 40d, 85d, false)]
        public void ARangeGivesBothEndsAndItsMidpoint(string cell, double min, double max, bool spec)
        {
            Assert.That(SdrfAge.TryParse(cell, out var age), Is.True);
            Assert.That(age!.Precision, Is.EqualTo(SdrfAgePrecision.Range));
            Assert.That(age.MinYears, Is.EqualTo(min).Within(Tolerance));
            Assert.That(age.MaxYears, Is.EqualTo(max).Within(Tolerance));
            Assert.That(age.Years, Is.EqualTo((min + max) / 2d).Within(Tolerance));
            Assert.That(age.FollowsSpecification, Is.EqualTo(spec));
        }

        /// <summary>
        /// A degenerate range pins one age, so it grades Exact rather than Range.
        ///
        /// This is not a cosmetic classification. The advice this library gives a caller who wants
        /// one age per sample is "filter on Precision == Exact"; grading 40Y-40Y as a Range makes
        /// that filter drop a cell that states an age exactly, and the drop looks identical to a
        /// refusal -- "this sample has no age". Raised by aging, thread 015.
        /// </summary>
        [TestCase("40Y-40Y", 40d, true)]
        [TestCase("8W-8W", 8d * 7d / 365.25, true)]
        [TestCase("6-6 weeks", 6d * 7d / 365.25, false)]
        [TestCase("40-40Y", 40d, false)]
        public void ADegenerateRangeIsExact(string cell, double years, bool spec)
        {
            Assert.That(SdrfAge.TryParse(cell, out var age), Is.True);
            Assert.That(age!.Precision, Is.EqualTo(SdrfAgePrecision.Exact),
                "A range whose ends are equal states one age, and Precision == Exact is the filter " +
                "every consumer is told to use for one.");
            Assert.That(age.Years, Is.EqualTo(years).Within(Tolerance));
            Assert.That(age.MinYears, Is.EqualTo(years).Within(Tolerance));
            Assert.That(age.MaxYears, Is.EqualTo(years).Within(Tolerance));
            Assert.That(age.FollowsSpecification, Is.EqualTo(spec),
                "Classifying it as Exact must not change what the cell's grammar was.");
        }

        /// <summary>
        /// A range with genuinely different ends is untouched by the degenerate case above.
        /// </summary>
        [Test]
        public void ARealRangeIsStillARange()
        {
            Assert.That(SdrfAge.TryParse("40Y-85Y", out var age), Is.True);
            Assert.That(age!.Precision, Is.EqualTo(SdrfAgePrecision.Range));
            Assert.That(age.Years, Is.EqualTo(62.5d).Within(Tolerance));
        }

        /// <summary>
        /// Every age carries the cell it was read from, so a normalised number can be audited back
        /// to the text a depositor wrote. Asked for by dataRepo, which stores the two side by side.
        /// </summary>
        [TestCase("58Y")]
        [TestCase(" 54y ")]
        [TestCase("4 hour")]
        [TestCase(">=90Y")]
        [TestCase("40Y-85Y")]
        [TestCase("40Y-40Y")]
        public void AnAgeKeepsTheCellItWasReadFrom(string cell)
        {
            Assert.That(SdrfAge.TryParse(cell, out var age), Is.True);
            Assert.That(age!.Cell, Is.EqualTo(cell.Trim()),
                "The parsed figure alone cannot be audited: 0.0109589 years does not say \"4 hour\".");
        }

                /// <summary>
        /// ">=90Y" is the de-identification cap: 750 cells in the corpus. The bound is the best single
        /// figure, and the open end is stated rather than invented.
        /// </summary>
        [Test]
        public void ABoundIsOpenAtOneEnd()
        {
            Assert.That(SdrfAge.TryParse(">=90Y", out var old), Is.True);
            Assert.That(old!.Precision, Is.EqualTo(SdrfAgePrecision.LowerBound));
            Assert.That(old.Years, Is.EqualTo(90d));
            Assert.That(old.MaxYears, Is.EqualTo(double.PositiveInfinity));

            Assert.That(SdrfAge.TryParse("<1Y", out var young), Is.True);
            Assert.That(young!.Precision, Is.EqualTo(SdrfAgePrecision.UpperBound));
            Assert.That(young.MinYears, Is.EqualTo(0d));
            Assert.That(young.MaxYears, Is.EqualTo(1d));
        }

        /// <summary>
        /// A bare number has no unit, and none can be recovered from the cell: 63 years and 63 days are
        /// both plausible in one study. Refusing is the whole point. So are reserved words, templates
        /// filled with one ("not availableY"), free text, and grammar out of order.
        /// </summary>
        [TestCase("63")]
        [TestCase("2.5")]
        [TestCase("not available")]
        [TestCase("Not Applicable")]
        [TestCase("not availableY")]
        [TestCase("not available years")]
        [TestCase("less than 1 day")]
        [TestCase("adult")]
        [TestCase("6M30Y")]
        [TestCase("30Y30Y")]
        [TestCase("85Y-40Y")]
        [TestCase("-5Y")]
        [TestCase("")]
        [TestCase("   ")]
        [TestCase(null)]
        public void WhatCannotBeReadWithoutGuessingIsRefused(string cell)
        {
            Assert.That(SdrfAge.TryParse(cell, out var age), Is.False);
            Assert.That(age, Is.Null);
        }

        /// <summary>
        /// How the curated corpus's age cells fall, for reference. [Explicit]; needs MZLIB_SDRF_CORPUS.
        /// </summary>
        [Test]
        [Explicit("Requires a local clone of bigbio/sdrf-annotated-datasets; set MZLIB_SDRF_CORPUS.")]
        public void CorpusAgeReport()
        {
            string corpus = Environment.GetEnvironmentVariable("MZLIB_SDRF_CORPUS");
            if (string.IsNullOrWhiteSpace(corpus) || !Directory.Exists(corpus))
                Assert.Ignore($"MZLIB_SDRF_CORPUS not set or not found: '{corpus}'");

            var cells = Directory.GetFiles(corpus, "*.sdrf.tsv", SearchOption.AllDirectories)
                .SelectMany(path =>
                {
                    var document = new SdrfDocument(path);
                    document.LoadResults();
                    return document.Results.SelectMany(r => r.All("characteristics[age]"));
                })
                .Where(c => !SdrfValidator.ReservedWords.Any(w => string.Equals(c?.Trim(), w, StringComparison.OrdinalIgnoreCase)))
                .ToList();

            var read = cells.Select(c => SdrfAge.TryParse(c, out var a) ? a : null).ToList();
            TestContext.Progress.WriteLine($"real age cells : {cells.Count}");
            TestContext.Progress.WriteLine($"refused        : {read.Count(a => a is null)}");
            foreach (var group in read.Where(a => a is not null)
                         .GroupBy(a => (a!.Precision, a.FollowsSpecification)).OrderBy(g => g.Key))
                TestContext.Progress.WriteLine($"{group.Key.Precision,-11} spec={group.Key.FollowsSpecification,-5} {group.Count()}");

            Assert.That(cells, Is.Not.Empty);
        }
    }
}
