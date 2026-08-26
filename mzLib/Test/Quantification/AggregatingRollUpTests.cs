using System;
using System.Buffers;
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
    /// Tests for AggregatingRollUp, which SumRollUp and MedianRollUp now both derive from.
    /// The parity tests are the point: the two named strategies must keep producing exactly what
    /// they produced when each carried its own loop.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class AggregatingRollUpTests
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

        /// <summary>
        /// A roll-up that mzLib does not ship, standing in for one a caller would add. Two lines is the
        /// whole cost of a new estimator now, which is the point of the base class.
        /// </summary>
        private class MaxRollUp : AggregatingRollUp
        {
            public MaxRollUp() : base(new MaxAggregation()) { }
        }

        /// <summary>Exists only to reach the base constructor with null.</summary>
        private class NullAggregationRollUp : AggregatingRollUp
        {
            public NullAggregationRollUp() : base(null) { }
        }

        /// <summary>A row key that is just a name, so the tests are about arithmetic and nothing else.</summary>
        private class Key : IEquatable<Key>
        {
            private readonly string _name;
            public Key(string name) { _name = name; }
            public bool Equals(Key other) => ReferenceEquals(this, other);
            public override bool Equals(object obj) => ReferenceEquals(this, obj);
            public override int GetHashCode() => System.Runtime.CompilerServices.RuntimeHelpers.GetHashCode(this);
            public override string ToString() => _name;
        }

        /// <summary>Builds a matrix from rows of column values.</summary>
        private static (QuantMatrix<Key> matrix, List<Key> rows) Matrix(params double[][] rows)
        {
            int columnCount = rows[0].Length;
            var columns = Enumerable.Range(0, columnCount)
                .Select(i => (ISampleInfo)new IsobaricQuantSampleInfo(
                    File, "Control", 0, 0, 0, 0, $"12{i}", 126.0 + i, false))
                .ToList();

            var design = new TestExperimentalDesign(
                new Dictionary<string, ISampleInfo[]> { [File] = columns.ToArray() });

            var keys = rows.Select((_, i) => new Key($"row{i}")).ToList();
            var matrix = new QuantMatrix<Key>(keys, columns, design);

            for (int r = 0; r < rows.Length; r++)
            {
                matrix.SetRow(keys[r], rows[r]);
            }

            return (matrix, keys);
        }

        private static Dictionary<Key, List<int>> MapAllRowsTo(Key target, int rowCount) =>
            new Dictionary<Key, List<int>> { [target] = Enumerable.Range(0, rowCount).ToList() };

        [Test]
        public void SumRollUp_StillSums()
        {
            var (matrix, _) = Matrix(
                new[] { 10.0, 100.0 },
                new[] { 20.0, 200.0 },
                new[] { 30.0, 300.0 });

            var target = new Key("protein");
            var result = new SumRollUp().RollUp(matrix, MapAllRowsTo(target, 3));

            Assert.Multiple(() =>
            {
                Assert.That(result.GetRow(target)[0], Is.EqualTo(60.0));
                Assert.That(result.GetRow(target)[1], Is.EqualTo(600.0));
            });
        }

        [Test]
        public void MedianRollUp_StillTakesTheMedianOfObservedValues()
        {
            var (matrix, _) = Matrix(
                new[] { 10.0, 0.0 },
                new[] { 20.0, 50.0 },
                new[] { 90.0, 0.0 });

            var target = new Key("protein");
            var result = new MedianRollUp().RollUp(matrix, MapAllRowsTo(target, 3));

            Assert.Multiple(() =>
            {
                // Median of 10, 20, 90
                Assert.That(result.GetRow(target)[0], Is.EqualTo(20.0));
                // Zero means "not observed", so the median of the one observation is that observation
                // -- not the median of {0, 50, 0}.
                Assert.That(result.GetRow(target)[1], Is.EqualTo(50.0));
            });
        }

        [Test]
        public void ColumnWithNoObservations_RollsUpToZero()
        {
            var (matrix, _) = Matrix(
                new[] { 10.0, 0.0 },
                new[] { 20.0, 0.0 });

            var target = new Key("protein");

            Assert.Multiple(() =>
            {
                Assert.That(new MedianRollUp().RollUp(matrix, MapAllRowsTo(target, 2)).GetRow(target)[1], Is.Zero);
                Assert.That(new SumRollUp().RollUp(matrix, MapAllRowsTo(target, 2)).GetRow(target)[1], Is.Zero);
            });
        }

        [Test]
        public void ExcludingZerosDoesNotChangeASum()
        {
            var (matrix, _) = Matrix(
                new[] { 10.0 },
                new[] { 0.0 },
                new[] { 30.0 });

            var target = new Key("protein");
            var result = new SumRollUp().RollUp(matrix, MapAllRowsTo(target, 3));

            // The reason both named strategies keep their old results under one shared loop.
            Assert.That(result.GetRow(target)[0], Is.EqualTo(40.0));
        }

        [Test]
        public void NamesAreUnchanged()
        {
            Assert.Multiple(() =>
            {
                Assert.That(new SumRollUp().Name, Is.EqualTo("Sum Roll-Up"));
                Assert.That(new MedianRollUp().Name, Is.EqualTo("Median Roll-Up"));
                Assert.That(new MaxRollUp().Name, Is.EqualTo("Max Roll-Up"));
            });
        }

        [Test]
        public void ANewRollUpIsATwoLineSubclass()
        {
            var (matrix, _) = Matrix(
                new[] { 10.0 },
                new[] { 70.0 },
                new[] { 30.0 });

            var target = new Key("protein");
            var result = new MaxRollUp().RollUp(matrix, MapAllRowsTo(target, 3));

            // MaxRollUp is defined at the top of this file and is the whole implementation.
            // Previously this needed a full IRollUpStrategy with its own copy of the loop.
            Assert.That(result.GetRow(target)[0], Is.EqualTo(70.0));
        }

        [Test]
        public void TheBaseClassCannotBeUsedDirectly()
        {
            // A roll-up has to be named by a concrete class, so that QuantificationParameters records
            // what was actually done rather than "some aggregation".
            Assert.That(typeof(AggregatingRollUp).IsAbstract, Is.True);
        }

        [Test]
        public void SumIsCorrectEvenWhenTheSharedPoolHandsBackADirtyBuffer()
        {
            // SumRollUp used to accumulate with += into an array straight from ArrayPool.Shared,
            // which Rent does not zero. Aggregation strategies in the same pipeline return their
            // buffers without clearing, so a sum could pick up another strategy's leftovers.
            var pool = ArrayPool<double>.Shared;
            double[] poison = pool.Rent(2);
            for (int i = 0; i < poison.Length; i++)
            {
                poison[i] = 999_999.0;
            }
            pool.Return(poison); // no clearArray -- exactly what MedianAggregation does

            var (matrix, _) = Matrix(
                new[] { 10.0, 100.0 },
                new[] { 20.0, 200.0 });

            var target = new Key("protein");
            var result = new SumRollUp().RollUp(matrix, MapAllRowsTo(target, 2));

            Assert.Multiple(() =>
            {
                Assert.That(result.GetRow(target)[0], Is.EqualTo(30.0));
                Assert.That(result.GetRow(target)[1], Is.EqualTo(300.0));
            });
        }

        [Test]
        public void EmptyMapGivesAnEmptyMatrix()
        {
            var (matrix, _) = Matrix(new[] { 10.0 });

            var result = new SumRollUp().RollUp(matrix, new Dictionary<Key, List<int>>());

            Assert.That(result.RowKeys, Is.Empty);
        }

        [Test]
        public void GroupWithNoRowsRollsUpToZeros()
        {
            var (matrix, _) = Matrix(new[] { 10.0, 20.0 });

            var target = new Key("protein");
            var result = new SumRollUp().RollUp(
                matrix, new Dictionary<Key, List<int>> { [target] = new List<int>() });

            Assert.That(result.GetRow(target), Is.EqualTo(new[] { 0.0, 0.0 }));
        }

        [Test]
        public void NullAggregationIsRejected()
        {
            Assert.Throws<ArgumentNullException>(() => new NullAggregationRollUp());
        }
    }
}
