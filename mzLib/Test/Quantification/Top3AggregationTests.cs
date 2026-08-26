using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Quantification;
using Quantification.Strategies;

namespace Test.Quantification
{
    /// <summary>
    /// Tests for Top3Aggregation and the Top3RollUp it drives -- the peptide-to-protein estimator
    /// that does not scale with how many peptides a search identified.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class Top3AggregationTests
    {
        private const string File = "file1.raw";

        private class TestExperimentalDesign : IExperimentalDesign
        {
            public Dictionary<string, ISampleInfo[]> FileNameSampleInfoDictionary { get; }

            public TestExperimentalDesign(Dictionary<string, ISampleInfo[]> dict)
            {
                FileNameSampleInfoDictionary = dict;
            }
        }

        private class Key : IEquatable<Key>
        {
            private readonly string _name;
            public Key(string name) { _name = name; }
            public bool Equals(Key other) => ReferenceEquals(this, other);
            public override bool Equals(object obj) => ReferenceEquals(this, obj);
            public override int GetHashCode() => System.Runtime.CompilerServices.RuntimeHelpers.GetHashCode(this);
            public override string ToString() => _name;
        }

        private static double Aggregate(params double[] values) =>
            new Top3Aggregation().Aggregate(new ReadOnlySpan<double>(values));

        [Test]
        public void AveragesTheThreeLargestValues()
        {
            // 100, 80 and 60 are the top three; the 5 and the 1 are ignored.
            Assert.That(Aggregate(5.0, 100.0, 60.0, 1.0, 80.0), Is.EqualTo(80.0));
        }

        [Test]
        public void OrderDoesNotMatter()
        {
            var ascending = Aggregate(1.0, 5.0, 60.0, 80.0, 100.0);
            var descending = Aggregate(100.0, 80.0, 60.0, 5.0, 1.0);

            Assert.That(ascending, Is.EqualTo(descending).Within(1e-12));
        }

        [Test]
        public void ExactlyThreeValuesAveragesAllOfThem()
        {
            Assert.That(Aggregate(10.0, 20.0, 30.0), Is.EqualTo(20.0));
        }

        [Test]
        public void FewerThanThreeAveragesWhatIsThere_NotDividedByThree()
        {
            Assert.Multiple(() =>
            {
                // Dividing by three here would report a protein seen twice as smaller than the same
                // protein seen three times -- the bias Top3 exists to avoid.
                Assert.That(Aggregate(10.0, 20.0), Is.EqualTo(15.0));
                Assert.That(Aggregate(10.0), Is.EqualTo(10.0));
            });
        }

        [Test]
        public void EmptyAggregatesToZero()
        {
            Assert.That(Aggregate(), Is.Zero);
        }

        [Test]
        public void RepeatedValuesAreCountedSeparately()
        {
            // Three peptides of equal intensity is a valid Top3, not one value seen three times.
            Assert.That(Aggregate(50.0, 50.0, 50.0, 10.0), Is.EqualTo(50.0));
        }

        [Test]
        public void IsNamedForParameterReadability()
        {
            Assert.Multiple(() =>
            {
                Assert.That(new Top3Aggregation().Name, Is.EqualTo("Top 3"));
                Assert.That(new Top3RollUp().Name, Is.EqualTo("Top 3 Roll-Up"));
            });
        }

        [Test]
        public void Top3RollUp_DoesNotScaleWithPeptideCount()
        {
            var columns = new List<ISampleInfo>
            {
                new IsobaricQuantSampleInfo(File, "Control", 0, 0, 0, 0, "126", 126.0, false)
            };
            var design = new TestExperimentalDesign(
                new Dictionary<string, ISampleInfo[]> { [File] = columns.ToArray() });

            // Six peptides, all at 100. Rows 0-2 are one protein, rows 0-5 the other.
            var rows = Enumerable.Range(0, 6).Select(i => new Key($"pep{i}")).ToList();
            var matrix = new QuantMatrix<Key>(rows, columns, design);
            foreach (var row in rows)
            {
                matrix.SetRow(row, new[] { 100.0 });
            }

            var threePeptides = new Key("threePeptideProtein");
            var sixPeptides = new Key("sixPeptideProtein");
            var map = new Dictionary<Key, List<int>>
            {
                [threePeptides] = new List<int> { 0, 1, 2 },
                [sixPeptides] = new List<int> { 0, 1, 2, 3, 4, 5 }
            };

            var top3 = new Top3RollUp().RollUp(matrix, map);
            var sum = new SumRollUp().RollUp(matrix, map);

            Assert.Multiple(() =>
            {
                // Same per-peptide abundance, so the same protein-level estimate.
                Assert.That(top3.GetRow(threePeptides)[0], Is.EqualTo(100.0));
                Assert.That(top3.GetRow(sixPeptides)[0], Is.EqualTo(100.0));

                // Summing, by contrast, reports the six-peptide protein as twice as abundant.
                Assert.That(sum.GetRow(threePeptides)[0], Is.EqualTo(300.0));
                Assert.That(sum.GetRow(sixPeptides)[0], Is.EqualTo(600.0));
            });
        }

        [Test]
        public void Top3RollUp_IgnoresUnobservedPeptides()
        {
            var columns = new List<ISampleInfo>
            {
                new IsobaricQuantSampleInfo(File, "Control", 0, 0, 0, 0, "126", 126.0, false)
            };
            var design = new TestExperimentalDesign(
                new Dictionary<string, ISampleInfo[]> { [File] = columns.ToArray() });

            var rows = new List<Key> { new Key("a"), new Key("b"), new Key("c"), new Key("d") };
            var matrix = new QuantMatrix<Key>(rows, columns, design);
            matrix.SetRow(rows[0], new[] { 30.0 });
            matrix.SetRow(rows[1], new[] { 0.0 });   // not observed
            matrix.SetRow(rows[2], new[] { 0.0 });   // not observed
            matrix.SetRow(rows[3], new[] { 50.0 });

            var protein = new Key("protein");
            var result = new Top3RollUp().RollUp(
                matrix, new Dictionary<Key, List<int>> { [protein] = new List<int> { 0, 1, 2, 3 } });

            // Two observations, so the mean of two -- the zeros are missing data, not low intensities.
            Assert.That(result.GetRow(protein)[0], Is.EqualTo(40.0));
        }
    }
}
