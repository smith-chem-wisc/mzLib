// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Collections.Generic;
using MassSpectrometry.Dia;
using NUnit.Framework;

namespace Test.Dia.Calibration
{
    [TestFixture]
    public class IrtLibraryIndexTests
    {
        /// <summary>
        /// Helper to create a list of LibraryPrecursorInput with specified iRT values.
        /// </summary>
        private static List<LibraryPrecursorInput> CreatePrecursors(params double?[] irtValues)
        {
            var precursors = new List<LibraryPrecursorInput>(irtValues.Length);
            for (int i = 0; i < irtValues.Length; i++)
            {
                precursors.Add(new LibraryPrecursorInput(
                    sequence: $"PEPTIDE{i}",
                    precursorMz: 500.0 + i,
                    chargeState: 2,
                    retentionTime: irtValues[i].HasValue ? irtValues[i].Value * 0.5 : null,
                    isDecoy: false,
                    fragmentMzs: new float[] { 200f, 300f },
                    fragmentIntensities: new float[] { 1000f, 500f },
                    irtValue: irtValues[i]
                ));
            }
            return precursors;
        }

        // ── Construction ─────────────────────────────────────────────────────

        [Test]
        public void Constructor_SortsbyIrt()
        {
            // Input iRTs in random order
            var precursors = CreatePrecursors(50.0, 10.0, 90.0, 30.0, 70.0);
            var index = new IrtLibraryIndex(precursors);

            Assert.That(index.Count, Is.EqualTo(5));
            Assert.That(index.MinIrt, Is.EqualTo(10.0));
            Assert.That(index.MaxIrt, Is.EqualTo(90.0));

            // Verify sorted order
            var sortedIrts = index.SortedIrts;
            for (int i = 1; i < sortedIrts.Length; i++)
            {
                Assert.That(sortedIrts[i], Is.GreaterThanOrEqualTo(sortedIrts[i - 1]));
            }
        }

        [Test]
        public void Constructor_SkipsNullIrt()
        {
            var precursors = CreatePrecursors(10.0, null, 30.0, null, 50.0);
            var index = new IrtLibraryIndex(precursors);

            Assert.That(index.Count, Is.EqualTo(3)); // Only 3 have iRT
        }

        [Test]
        public void Constructor_EmptyList_CreatesEmptyIndex()
        {
            var index = new IrtLibraryIndex(new List<LibraryPrecursorInput>());
            Assert.That(index.Count, Is.EqualTo(0));
            Assert.That(double.IsNaN(index.MinIrt));
            Assert.That(double.IsNaN(index.MaxIrt));
        }

        [Test]
        public void Constructor_AllNullIrt_CreatesEmptyIndex()
        {
            var precursors = CreatePrecursors(null, null, null);
            var index = new IrtLibraryIndex(precursors);
            Assert.That(index.Count, Is.EqualTo(0));
        }

        [Test]
        public void Constructor_NullList_Throws()
        {
            Assert.Throws<ArgumentNullException>(() => new IrtLibraryIndex(null));
        }

        // ── QueryRange ──────────────────────────────────────────────────────

        [Test]
        public void QueryRange_FullRange_ReturnsAll()
        {
            var precursors = CreatePrecursors(10.0, 20.0, 30.0, 40.0, 50.0);
            var index = new IrtLibraryIndex(precursors);

            index.QueryRange(0.0, 100.0, out int start, out int count);
            Assert.That(count, Is.EqualTo(5));
            Assert.That(start, Is.EqualTo(0));
        }

        [Test]
        public void QueryRange_SubRange_ReturnsCorrectSubset()
        {
            // iRTs: 10, 20, 30, 40, 50
            var precursors = CreatePrecursors(10.0, 20.0, 30.0, 40.0, 50.0);
            var index = new IrtLibraryIndex(precursors);

            // Query [15, 45] → should return 20, 30, 40
            index.QueryRange(15.0, 45.0, out int start, out int count);
            Assert.That(count, Is.EqualTo(3));

            // Verify the returned range has iRTs 20, 30, 40
            for (int i = 0; i < count; i++)
            {
                double irt = index.GetIrt(start + i);
                Assert.That(irt, Is.GreaterThanOrEqualTo(15.0));
                Assert.That(irt, Is.LessThanOrEqualTo(45.0));
            }
        }

        [Test]
        public void QueryRange_ExactBoundaryMatch_Inclusive()
        {
            var precursors = CreatePrecursors(10.0, 20.0, 30.0);
            var index = new IrtLibraryIndex(precursors);

            // Query [10, 30] — boundaries are inclusive
            index.QueryRange(10.0, 30.0, out _, out int count);
            Assert.That(count, Is.EqualTo(3));
        }

        [Test]
        public void QueryRange_NoMatch_ReturnsZero()
        {
            var precursors = CreatePrecursors(10.0, 20.0, 30.0);
            var index = new IrtLibraryIndex(precursors);

            index.QueryRange(50.0, 60.0, out _, out int count);
            Assert.That(count, Is.EqualTo(0));
        }

        [Test]
        public void QueryRange_BelowAll_ReturnsZero()
        {
            var precursors = CreatePrecursors(10.0, 20.0, 30.0);
            var index = new IrtLibraryIndex(precursors);

            index.QueryRange(-10.0, 5.0, out _, out int count);
            Assert.That(count, Is.EqualTo(0));
        }

        [Test]
        public void QueryRange_EmptyIndex_ReturnsZero()
        {
            var index = new IrtLibraryIndex(new List<LibraryPrecursorInput>());
            index.QueryRange(0.0, 100.0, out _, out int count);
            Assert.That(count, Is.EqualTo(0));
        }

        // ── GetOriginalIndices ──────────────────────────────────────────────

        [Test]
        public void GetOriginalIndices_MapsBackToInputList()
        {
            // Input order: iRT 50, 10, 90, 30, 70 → sorted: 10, 30, 50, 70, 90
            // Original indices:     1,  3,  0,  4,  2
            var precursors = CreatePrecursors(50.0, 10.0, 90.0, 30.0, 70.0);
            var index = new IrtLibraryIndex(precursors);

            // Query [25, 55] → iRTs 30, 50 → original indices 3, 0
            int[] originals = index.GetOriginalIndices(25.0, 55.0);

            Assert.That(originals.Length, Is.EqualTo(2));
            // Should contain original indices 3 (iRT=30) and 0 (iRT=50)
            Assert.That(originals, Does.Contain(3));
            Assert.That(originals, Does.Contain(0));

            // Verify the original precursors match
            foreach (int origIdx in originals)
            {
                Assert.That(precursors[origIdx].IrtValue.Value,
                    Is.GreaterThanOrEqualTo(25.0).And.LessThanOrEqualTo(55.0));
            }
        }

        [Test]
        public void GetOriginalIndices_NoMatch_ReturnsEmptyArray()
        {
            var precursors = CreatePrecursors(10.0, 20.0, 30.0);
            var index = new IrtLibraryIndex(precursors);

            int[] result = index.GetOriginalIndices(100.0, 200.0);
            Assert.That(result, Is.Empty);
        }

        [Test]
        public void GetOriginalIndex_ViaQueryRange_Consistent()
        {
            var precursors = CreatePrecursors(10.0, 20.0, 30.0, 40.0, 50.0);
            var index = new IrtLibraryIndex(precursors);

            index.QueryRange(15.0, 35.0, out int start, out int count);
            int[] fromMethod = index.GetOriginalIndices(15.0, 35.0);

            Assert.That(fromMethod.Length, Is.EqualTo(count));
            for (int i = 0; i < count; i++)
            {
                Assert.That(fromMethod[i], Is.EqualTo(index.GetOriginalIndex(start + i)));
            }
        }

        // ── Duplicate iRT values ────────────────────────────────────────────

        [Test]
        public void QueryRange_DuplicateIrts_ReturnsAll()
        {
            var precursors = CreatePrecursors(20.0, 20.0, 20.0, 40.0, 40.0);
            var index = new IrtLibraryIndex(precursors);

            index.QueryRange(19.0, 21.0, out _, out int count);
            Assert.That(count, Is.EqualTo(3)); // All three iRT=20 entries

            index.QueryRange(39.0, 41.0, out _, out int count2);
            Assert.That(count2, Is.EqualTo(2)); // Both iRT=40 entries
        }

        // ── Large-scale correctness ─────────────────────────────────────────

        [Test]
        public void QueryRange_LargeIndex_CorrectResults()
        {
            // Build 10K precursors with random iRTs in [0, 120]
            var rng = new Random(42);
            int n = 10_000;
            var precursors = new List<LibraryPrecursorInput>(n);

            for (int i = 0; i < n; i++)
            {
                double irt = rng.NextDouble() * 120.0;
                precursors.Add(new LibraryPrecursorInput(
                    $"PEP{i}", 500.0 + i * 0.1, 2, irt * 0.5, false,
                    new float[] { 300f }, new float[] { 1000f },
                    irtValue: irt));
            }

            var index = new IrtLibraryIndex(precursors);
            Assert.That(index.Count, Is.EqualTo(n));

            // Query a window and verify all results are in range
            double lower = 40.0, upper = 60.0;
            index.QueryRange(lower, upper, out int start, out int count);

            // Count expected by brute force
            int expectedCount = 0;
            for (int i = 0; i < n; i++)
            {
                if (precursors[i].IrtValue.Value >= lower && precursors[i].IrtValue.Value <= upper)
                    expectedCount++;
            }

            Assert.That(count, Is.EqualTo(expectedCount));

            // Verify each returned entry is in range
            for (int i = 0; i < count; i++)
            {
                double irt = index.GetIrt(start + i);
                Assert.That(irt, Is.GreaterThanOrEqualTo(lower));
                Assert.That(irt, Is.LessThanOrEqualTo(upper));
            }
        }

        // ── SortedIrts property ─────────────────────────────────────────────

        [Test]
        public void SortedIrts_IsStrictlyNonDecreasing()
        {
            var rng = new Random(99);
            var precursors = new List<LibraryPrecursorInput>(100);
            for (int i = 0; i < 100; i++)
            {
                precursors.Add(new LibraryPrecursorInput(
                    $"P{i}", 500, 2, null, false,
                    new float[] { 200f }, new float[] { 1000f },
                    irtValue: rng.NextDouble() * 100.0));
            }

            var index = new IrtLibraryIndex(precursors);
            var sorted = index.SortedIrts;

            for (int i = 1; i < sorted.Length; i++)
            {
                Assert.That(sorted[i], Is.GreaterThanOrEqualTo(sorted[i - 1]),
                    $"SortedIrts not non-decreasing at index {i}");
            }
        }

        // ── Single element ──────────────────────────────────────────────────

        [Test]
        public void QueryRange_SingleElement_FoundWhenInRange()
        {
            var precursors = CreatePrecursors(42.0);
            var index = new IrtLibraryIndex(precursors);

            index.QueryRange(40.0, 44.0, out _, out int count1);
            Assert.That(count1, Is.EqualTo(1));

            index.QueryRange(43.0, 50.0, out _, out int count2);
            Assert.That(count2, Is.EqualTo(0));
        }

        [Test]
        public void MinMaxIrt_SingleElement_Equal()
        {
            var precursors = CreatePrecursors(42.0);
            var index = new IrtLibraryIndex(precursors);

            Assert.That(index.MinIrt, Is.EqualTo(42.0));
            Assert.That(index.MaxIrt, Is.EqualTo(42.0));
        }
    }
}
