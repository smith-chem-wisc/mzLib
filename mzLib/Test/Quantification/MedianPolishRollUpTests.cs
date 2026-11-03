using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Quantification;
using Quantification.Strategies;
using FlashLfqResults = FlashLFQ.FlashLfqResults;
using FlashLfqDetectionType = FlashLFQ.DetectionType;
using FlashLfqIdentification = FlashLFQ.Identification;
using FlashLfqPeptide = FlashLFQ.Peptide;
using FlashLfqProteinGroup = FlashLFQ.ProteinGroup;

namespace Test.Quantification
{
    /// <summary>
    /// Tests for MedianPolishRollUp. The parity tests are the point: the roll-up wraps FlashLFQ's median
    /// polish, so for one sample per column it must give what FlashLFQ itself gives, except where FlashLFQ
    /// reports NaN.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class MedianPolishRollUpTests
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
        private static QuantMatrix<Key> Matrix(params double[][] rows)
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

            return matrix;
        }

        private static double[] RollUpAllRows(double[][] rows)
        {
            var matrix = Matrix(rows);
            var protein = new Key("protein");
            var map = new Dictionary<Key, List<int>> { [protein] = Enumerable.Range(0, rows.Length).ToList() };

            return new MedianPolishRollUp().RollUp(matrix, map).GetRow(protein);
        }

        /// <summary>
        /// Runs the same table through FlashLFQ's own protein quantification, one biological replicate
        /// per column, the way TestFlashLFQ does.
        /// </summary>
        private static double[] FlashLfqProteinIntensities(double[][] peptideIntensities)
        {
            int columnCount = peptideIntensities[0].Length;
            var files = Enumerable.Range(0, columnCount)
                .Select(col => new SpectraFileInfo("", "cond", col, 0, 0))
                .ToList();

            var res = new FlashLfqResults(files, new List<FlashLfqIdentification>());
            var protein = new FlashLfqProteinGroup("accession1", "gene1", "organism1");
            res.ProteinGroups.Add(protein.ProteinGroupName, protein);

            for (int row = 0; row < peptideIntensities.Length; row++)
            {
                var pep = new FlashLfqPeptide("PEPTIDE" + row, "PEPTIDE" + row, true,
                    new HashSet<FlashLfqProteinGroup> { protein });
                res.PeptideModifiedSequences.Add(pep.Sequence, pep);

                for (int col = 0; col < columnCount; col++)
                {
                    pep.SetIntensity(files[col], peptideIntensities[row][col]);
                    pep.SetDetectionType(files[col], FlashLfqDetectionType.MSMS);
                }
            }

            res.CalculateProteinResultsMedianPolish(useSharedPeptides: false);

            return files.Select(f => protein.GetIntensity(f)).ToArray();
        }

        public static IEnumerable<TestCaseData> FlashLfqTables()
        {
            yield return new TestCaseData((object)new[]
            {
                new[] { 1000.0, 1010, 2050, 2010 },
                new[] { 2000.0, 1900, 3900, 4100 },
            }).SetName("TrueChanger");

            yield return new TestCaseData((object)new[]
            {
                new[] { 19007964.63, 18208648.31, 0 },
                new[] { 14890408.31, 14411359.03, 14910408.31 },
                new[] { 27894671.25, 27384454.5, 27914671.25 },
                new[] { 20857567.75, 21912047.75, 20877567.75 },
                new[] { 9142974.708, 17019194.54, 9162974.708 },
                new[] { 29630634, 29207244.17, 0 },
                new[] { 4463091.969, 3933777.311, 4483091.969 },
                new[] { 18686402.29, 19290511.46, 18706402.29 },
                new[] { 4073149.652, 4451360.921, 4093149.652 },
            }).SetName("SimilarIntensityRanks");

            yield return new TestCaseData((object)new[]
            {
                new[] { 9796546, 9852023.625, 0 },
                new[] { 2193142.286, 2132802.28, 0 },
                new[] { 13807677.53, 12251660.38, 0 },
            }).SetName("MissingProteinValue");

            yield return new TestCaseData((object)new[]
            {
                new[] { 9796546.0, 9852023.625 },
                new[] { 22670965.15, 0 },
                new[] { 0, 0.0 },
                new[] { 0, 2121691.667 },
                new[] { 13807677.53, 12251660.38 },
            }).SetName("SinglyMeasuredPeptides");

            yield return new TestCaseData((object)new[]
            {
                new[] { 9796546.0, 0 },
            }).SetName("OneValidValue");
        }

        [Test, TestCaseSource(nameof(FlashLfqTables))]
        public void MatchesFlashLfq(double[][] peptideIntensities)
        {
            double[] expected = FlashLfqProteinIntensities(peptideIntensities);
            double[] actual = RollUpAllRows(peptideIntensities);

            Assert.That(expected.Any(double.IsNaN), Is.False, "Parity tables should not contain unquantifiable samples");
            Assert.That(actual, Is.EqualTo(expected));
        }

        /// <summary>
        /// No peptide is seen in more than one sample, so there is nothing to compare across samples. Each sample is
        /// quantified from its own peptides, as a protein seen only in that sample is, scaled to the group's number of
        /// peptides. Before, the fit set every column effect to 0 and the unquantifiable check then reported NaN in
        /// every sample, so adding a second peptide in a second sample took away the value a one-peptide group had.
        /// </summary>
        [Test]
        public void PeptidesEachSeenInOneDifferentSample_QuantifyEachSampleFromItsOwnPeptides()
        {
            var table = new[]
            {
                new[] { 0, 1000.0 },
                new[] { 1000.0, 0 },
            };

            Assert.That(FlashLfqProteinIntensities(table), Is.EqualTo(new[] { 2000.0, 2000.0 }).Within(1e-6));
            Assert.That(RollUpAllRows(table), Is.EqualTo(new[] { 2000.0, 2000.0 }).Within(1e-6));
        }

        /// <summary>
        /// Real groups from MetaMorpheus's islet search (master, 4 label-free files), each with every peptide seen in one
        /// file only. A2BDX3: two peptides in files 1 and 3; it had no value in any file. Q91V09: one peptide in file 1
        /// and three in file 3; it had the same number in both files, 5,883,728, which is neither file's peptides.
        /// </summary>
        [Test]
        public void IsletGroupsWithEveryPeptideInOneFile_GetEachFilesOwnValue()
        {
            var a2bdx3 = new[]
            {
                new[] { 2163521.0, 0, 0, 0 },
                new[] { 0, 0, 109879500.0, 0 },
            };
            Assert.That(FlashLfqProteinIntensities(a2bdx3), Is.EqualTo(new[] { 2 * 2163521.0, 0, 2 * 109879500.0, 0 }).Within(1e-3));

            var q91v09 = new[]
            {
                new[] { 0, 0, 193941.0, 0 },
                new[] { 779644.0, 0, 0, 0 },
                new[] { 0, 0, 2266565.0, 0 },
                new[] { 0, 0, 2613322.0, 0 },
            };
            // file 3 gets what its three peptides give as a protein seen only there, scaled from 3 peptides to the group's 4
            double file3Alone = MedianPolishRollUp.QuantifyGroup(new[] { new[] { 193941.0 }, new[] { 2266565.0 }, new[] { 2613322.0 } })[0];
            Assert.That(FlashLfqProteinIntensities(q91v09), Is.EqualTo(new[] { 4 * 779644.0, 0, file3Alone * 4 / 3, 0 }).Within(1e-3));
        }

        /// <summary>A protein seen in one sample only is not this case: its value comes from the usual fit, unchanged.</summary>
        [Test]
        public void ProteinSeenInOneSampleOnly_KeepsItsFittedValue()
        {
            var oneSample = new[] { new[] { 39282.17 }, new[] { 573731.04 }, new[] { 21413355.83 }, new[] { 137757.15 }, new[] { 588566.70 } };
            var oneOfFour = oneSample.Select(row => new[] { 0, 0, 0, row[0] }).ToArray();

            Assert.That(MedianPolishRollUp.QuantifyGroup(oneOfFour)[3], Is.EqualTo(MedianPolishRollUp.QuantifyGroup(oneSample)[0]).Within(1e-6));
        }

        /// <summary>
        /// The fix is only for a group with no peptide in two samples. One peptide seen in two samples links them, and
        /// the fit is the median polish it always was; a peptide seen in one sample is then set aside, as before.
        /// </summary>
        [Test]
        public void APeptideSeenInTwoSamples_KeepsTheMedianPolishFit()
        {
            var table = new[]
            {
                new[] { 1000.0, 2000.0, 0 },
                new[] { 0, 0, 5000.0 },
            };

            double[] protein = FlashLfqProteinIntensities(table);
            Assert.That(protein[1] / protein[0], Is.EqualTo(2.0).Within(1e-9), "the linked samples keep their fitted ratio");
        }

        [Test]
        public void ConstantPeptideRatios_GiveThatRatio()
        {
            // Two peptides that ionize 10x differently, both doubling between samples.
            double[] protein = RollUpAllRows(new[]
            {
                new[] { 100.0, 200.0 },
                new[] { 1000.0, 2000.0 },
            });

            Assert.That(protein[0], Is.GreaterThan(0));
            Assert.That(protein[1] / protein[0], Is.EqualTo(2.0).Within(1e-9));
        }

        [Test]
        public void MissingOrNegativeValues_AreNotObserved()
        {
            double[] protein = RollUpAllRows(new[]
            {
                new[] { 100.0, 200.0, 0 },
                new[] { 1000.0, 2000.0, -5 },
            });

            Assert.That(protein[2], Is.EqualTo(0));
        }

        [Test]
        public void EmptyOrNullGroup_RollsUpToZeros()
        {
            var matrix = Matrix(new[] { 100.0, 200.0 });
            var empty = new Key("empty");
            var nullGroup = new Key("null");
            var map = new Dictionary<Key, List<int>> { [empty] = new List<int>(), [nullGroup] = null };

            var result = new MedianPolishRollUp().RollUp(matrix, map);

            Assert.That(result.GetRow(empty), Is.EqualTo(new[] { 0.0, 0.0 }));
            Assert.That(result.GetRow(nullGroup), Is.EqualTo(new[] { 0.0, 0.0 }));
        }

        [Test]
        public void DuplicateIndices_CountOnce()
        {
            var matrix = Matrix(
                new[] { 100.0, 150.0, 90.0 },
                new[] { 1000.0, 2100.0, 800.0 },
                new[] { 50.0, 60.0, 70.0 });
            var once = new Key("once");
            var twice = new Key("twice");
            var map = new Dictionary<Key, List<int>>
            {
                [once] = new List<int> { 0, 1, 2 },
                [twice] = new List<int> { 0, 1, 1, 2, 2 },
            };

            var result = new MedianPolishRollUp().RollUp(matrix, map);

            Assert.That(result.GetRow(twice), Is.EqualTo(result.GetRow(once)));
        }

        [Test]
        public void EachGroup_UsesOnlyItsOwnRows()
        {
            var rowsA = new[] { new[] { 100.0, 200.0, 300.0 }, new[] { 1000.0, 1900.0, 3100.0 } };
            var rowsB = new[] { new[] { 5.0, 50.0, 7.0 }, new[] { 8.0, 90.0, 6.0 } };
            var matrix = Matrix(rowsA.Concat(rowsB).ToArray());
            var a = new Key("a");
            var b = new Key("b");
            var map = new Dictionary<Key, List<int>> { [a] = new List<int> { 0, 1 }, [b] = new List<int> { 2, 3 } };

            var result = new MedianPolishRollUp().RollUp(matrix, map);

            Assert.That(result.GetRow(a), Is.EqualTo(RollUpAllRows(rowsA)));
            Assert.That(result.GetRow(b), Is.EqualTo(RollUpAllRows(rowsB)));
        }

        [Test]
        public void Name_IsMedianPolish()
        {
            Assert.That(new MedianPolishRollUp().Name, Is.EqualTo("Median Polish Roll-Up"));
        }
    }
}
