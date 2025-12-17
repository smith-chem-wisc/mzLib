using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry.ExperimentalDesign;
using NUnit.Framework;

namespace Test.Omics
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class ISampleInfoTests
    {
        #region Test Implementation

        /// <summary>
        /// Test implementation of <see cref="ISampleInfo"/> for unit testing purposes.
        /// Implements equality based on all identity components: FullFilePathWithExtension,
        /// Condition, BiologicalReplicate, TechnicalReplicate, and Fraction.
        /// </summary>
        private class TestSampleInfo : ISampleInfo, IEquatable<TestSampleInfo>
        {
            public string FullFilePathWithExtension { get; }
            public string Condition { get; }
            public int BiologicalReplicate { get; }
            public int TechnicalReplicate { get; }
            public int Fraction { get; }

            public TestSampleInfo(
                string fullFilePathWithExtension,
                string condition,
                int biologicalReplicate,
                int technicalReplicate,
                int fraction)
            {
                FullFilePathWithExtension = fullFilePathWithExtension ?? string.Empty;
                Condition = condition ?? string.Empty;
                BiologicalReplicate = biologicalReplicate;
                TechnicalReplicate = technicalReplicate;
                Fraction = fraction;
            }

            public bool Equals(TestSampleInfo? other)
            {
                if (other is null) return false;
                if (ReferenceEquals(this, other)) return true;
                return string.Equals(FullFilePathWithExtension, other.FullFilePathWithExtension, StringComparison.Ordinal)
                    && string.Equals(Condition, other.Condition, StringComparison.Ordinal)
                    && BiologicalReplicate == other.BiologicalReplicate
                    && TechnicalReplicate == other.TechnicalReplicate
                    && Fraction == other.Fraction;
            }

            public bool Equals(ISampleInfo? other)
            {
                if (other is TestSampleInfo testSampleInfo)
                {
                    return Equals(testSampleInfo);
                }
                return false;
            }

            public override bool Equals(object? obj)
            {
                if (obj is null) return false;
                if (obj is TestSampleInfo testSampleInfo)
                {
                    return Equals(testSampleInfo);
                }
                if (obj is ISampleInfo sampleInfo)
                {
                    return Equals(sampleInfo);
                }
                return false;
            }

            public override int GetHashCode()
            {
                return HashCode.Combine(
                    StringComparer.Ordinal.GetHashCode(FullFilePathWithExtension),
                    StringComparer.Ordinal.GetHashCode(Condition),
                    BiologicalReplicate,
                    TechnicalReplicate,
                    Fraction);
            }

            public override string ToString()
            {
                return $"{Condition}_B{BiologicalReplicate}_T{TechnicalReplicate}_F{Fraction}";
            }

            public int CompareTo(ISampleInfo? other)
            {
                // Non-null comes before null
                if (other is null) return -1;

                // Compare by FullFilePathWithExtension (A before B)
                int comparison = string.Compare(FullFilePathWithExtension, other.FullFilePathWithExtension, StringComparison.Ordinal);
                if (comparison != 0) return comparison;

                // Compare by Condition (A before B)
                comparison = string.Compare(Condition, other.Condition, StringComparison.Ordinal);
                if (comparison != 0) return comparison;

                // Compare by BiologicalReplicate (1 before 2)
                comparison = BiologicalReplicate.CompareTo(other.BiologicalReplicate);
                if (comparison != 0) return comparison;

                // Compare by TechnicalReplicate (1 before 2)
                comparison = TechnicalReplicate.CompareTo(other.TechnicalReplicate);
                if (comparison != 0) return comparison;

                // Compare by Fraction (1 before 2)
                return Fraction.CompareTo(other.Fraction);
            }

            public static bool operator ==(TestSampleInfo? left, TestSampleInfo? right)
            {
                if (left is null) return right is null;
                return left.Equals(right);
            }

            public static bool operator !=(TestSampleInfo? left, TestSampleInfo? right)
            {
                return !(left == right);
            }
        }

        #endregion

        #region Test Data

        private TestSampleInfo _sample1 = null!;
        private TestSampleInfo _sample2 = null!;
        private TestSampleInfo _sampleIdentical = null!;
        private TestSampleInfo _sampleDifferentCondition = null!;
        private TestSampleInfo _sampleDifferentBioRep = null!;
        private TestSampleInfo _sampleDifferentTechRep = null!;
        private TestSampleInfo _sampleDifferentFraction = null!;
        private TestSampleInfo _sampleDifferentFilePath = null!;
        private TestSampleInfo _sampleEmptyStrings = null!;

        [SetUp]
        public void Setup()
        {
            _sample1 = new TestSampleInfo(
                @"C:\Data\sample1.raw",
                "Control",
                1,
                1,
                0);

            _sample2 = new TestSampleInfo(
                @"C:\Data\sample2.raw",
                "Treatment",
                2,
                2,
                1);

            _sampleIdentical = new TestSampleInfo(
                @"C:\Data\sample1.raw",
                "Control",
                1,
                1,
                0);

            _sampleDifferentCondition = new TestSampleInfo(
                @"C:\Data\sample1.raw",
                "Treatment",
                1,
                1,
                0);

            _sampleDifferentBioRep = new TestSampleInfo(
                @"C:\Data\sample1.raw",
                "Control",
                2,
                1,
                0);

            _sampleDifferentTechRep = new TestSampleInfo(
                @"C:\Data\sample1.raw",
                "Control",
                1,
                2,
                0);

            _sampleDifferentFraction = new TestSampleInfo(
                @"C:\Data\sample1.raw",
                "Control",
                1,
                1,
                1);

            _sampleDifferentFilePath = new TestSampleInfo(
                @"C:\Data\different.raw",
                "Control",
                1,
                1,
                0);

            _sampleEmptyStrings = new TestSampleInfo(
                string.Empty,
                string.Empty,
                0,
                0,
                0);
        }

        #endregion

        #region Constructor Tests

        [Test]
        public void Constructor_WithValidParameters_SetsAllProperties()
        {
            var sample = new TestSampleInfo(
                @"C:\Data\test.raw",
                "Control",
                1,
                2,
                3);

            Assert.That(sample.FullFilePathWithExtension, Is.EqualTo(@"C:\Data\test.raw"));
            Assert.That(sample.Condition, Is.EqualTo("Control"));
            Assert.That(sample.BiologicalReplicate, Is.EqualTo(1));
            Assert.That(sample.TechnicalReplicate, Is.EqualTo(2));
            Assert.That(sample.Fraction, Is.EqualTo(3));
        }

        [Test]
        public void Constructor_WithNullFilePath_SetsEmptyString()
        {
            var sample = new TestSampleInfo(null!, "Control", 1, 1, 0);
            Assert.That(sample.FullFilePathWithExtension, Is.EqualTo(string.Empty));
        }

        [Test]
        public void Constructor_WithNullCondition_SetsEmptyString()
        {
            var sample = new TestSampleInfo(@"C:\test.raw", null!, 1, 1, 0);
            Assert.That(sample.Condition, Is.EqualTo(string.Empty));
        }

        [Test]
        public void Constructor_WithNegativeValues_AcceptsValues()
        {
            var sample = new TestSampleInfo(@"C:\test.raw", "Control", -1, -2, -3);
            Assert.That(sample.BiologicalReplicate, Is.EqualTo(-1));
            Assert.That(sample.TechnicalReplicate, Is.EqualTo(-2));
            Assert.That(sample.Fraction, Is.EqualTo(-3));
        }

        [Test]
        public void Constructor_WithZeroValues_AcceptsValues()
        {
            var sample = new TestSampleInfo(@"C:\test.raw", "Control", 0, 0, 0);
            Assert.That(sample.BiologicalReplicate, Is.EqualTo(0));
            Assert.That(sample.TechnicalReplicate, Is.EqualTo(0));
            Assert.That(sample.Fraction, Is.EqualTo(0));
        }

        [Test]
        public void Constructor_WithMaxIntValues_AcceptsValues()
        {
            var sample = new TestSampleInfo(@"C:\test.raw", "Control", int.MaxValue, int.MaxValue, int.MaxValue);
            Assert.That(sample.BiologicalReplicate, Is.EqualTo(int.MaxValue));
            Assert.That(sample.TechnicalReplicate, Is.EqualTo(int.MaxValue));
            Assert.That(sample.Fraction, Is.EqualTo(int.MaxValue));
        }

        #endregion

        #region Equals Tests

        [Test]
        public void Equals_SameReference_ReturnsTrue()
        {
            Assert.That(_sample1.Equals(_sample1), Is.True);
        }

        [Test]
        public void Equals_IdenticalValues_ReturnsTrue()
        {
            Assert.That(_sample1.Equals(_sampleIdentical), Is.True);
        }

        [Test]
        public void Equals_Null_ReturnsFalse()
        {
            Assert.That(_sample1.Equals((TestSampleInfo?)null), Is.False);
            Assert.That(_sample1.Equals((ISampleInfo?)null), Is.False);
            Assert.That(_sample1.Equals((object?)null), Is.False);
        }

        [Test]
        public void Equals_DifferentFilePath_ReturnsFalse()
        {
            Assert.That(_sample1.Equals(_sampleDifferentFilePath), Is.False);
        }

        [Test]
        public void Equals_DifferentCondition_ReturnsFalse()
        {
            Assert.That(_sample1.Equals(_sampleDifferentCondition), Is.False);
        }

        [Test]
        public void Equals_DifferentBiologicalReplicate_ReturnsFalse()
        {
            Assert.That(_sample1.Equals(_sampleDifferentBioRep), Is.False);
        }

        [Test]
        public void Equals_DifferentTechnicalReplicate_ReturnsFalse()
        {
            Assert.That(_sample1.Equals(_sampleDifferentTechRep), Is.False);
        }

        [Test]
        public void Equals_DifferentFraction_ReturnsFalse()
        {
            Assert.That(_sample1.Equals(_sampleDifferentFraction), Is.False);
        }

        [Test]
        public void Equals_ObjectEquals_IdenticalValues_ReturnsTrue()
        {
            Assert.That(_sample1.Equals((object)_sampleIdentical), Is.True);
        }

        [Test]
        public void Equals_ObjectEquals_DifferentType_ReturnsFalse()
        {
            Assert.That(_sample1.Equals("not a sample"), Is.False);
        }

        [Test]
        public void Equals_CaseSensitive_Condition_ReturnsFalse()
        {
            var sample = new TestSampleInfo(@"C:\Data\sample1.raw", "control", 1, 1, 0);
            Assert.That(_sample1.Equals(sample), Is.False);
        }

        [Test]
        public void Equals_CaseSensitive_FilePath_ReturnsFalse()
        {
            var sample = new TestSampleInfo(@"C:\DATA\SAMPLE1.RAW", "Control", 1, 1, 0);
            Assert.That(_sample1.Equals(sample), Is.False);
        }

        [Test]
        public void Equals_EmptyStrings_BothEmpty_ReturnsTrue()
        {
            var sample1 = new TestSampleInfo(string.Empty, string.Empty, 0, 0, 0);
            var sample2 = new TestSampleInfo(string.Empty, string.Empty, 0, 0, 0);
            Assert.That(sample1.Equals(sample2), Is.True);
        }

        #endregion

        #region Equality Operator Tests

        [Test]
        public void EqualityOperator_IdenticalValues_ReturnsTrue()
        {
            Assert.That(_sample1 == _sampleIdentical, Is.True);
        }

        [Test]
        public void EqualityOperator_DifferentValues_ReturnsFalse()
        {
            Assert.That(_sample1 == _sample2, Is.False);
        }

        [Test]
        public void EqualityOperator_BothNull_ReturnsTrue()
        {
            TestSampleInfo? left = null;
            TestSampleInfo? right = null;
            Assert.That(left == right, Is.True);
        }

        [Test]
        public void EqualityOperator_LeftNull_ReturnsFalse()
        {
            TestSampleInfo? left = null;
            Assert.That(left == _sample1, Is.False);
        }

        [Test]
        public void EqualityOperator_RightNull_ReturnsFalse()
        {
            TestSampleInfo? right = null;
            Assert.That(_sample1 == right, Is.False);
        }

        [Test]
        public void InequalityOperator_DifferentValues_ReturnsTrue()
        {
            Assert.That(_sample1 != _sample2, Is.True);
        }

        [Test]
        public void InequalityOperator_IdenticalValues_ReturnsFalse()
        {
            Assert.That(_sample1 != _sampleIdentical, Is.False);
        }

        [Test]
        public void InequalityOperator_OneNull_ReturnsTrue()
        {
            TestSampleInfo? nullSample = null;
            Assert.That(_sample1 != nullSample, Is.True);
            Assert.That(nullSample != _sample1, Is.True);
        }

        #endregion

        #region GetHashCode Tests

        [Test]
        public void GetHashCode_IdenticalValues_ReturnsSameHashCode()
        {
            Assert.That(_sample1.GetHashCode(), Is.EqualTo(_sampleIdentical.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentFilePath_ReturnsDifferentHashCode()
        {
            Assert.That(_sample1.GetHashCode(), Is.Not.EqualTo(_sampleDifferentFilePath.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentCondition_ReturnsDifferentHashCode()
        {
            Assert.That(_sample1.GetHashCode(), Is.Not.EqualTo(_sampleDifferentCondition.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentBiologicalReplicate_ReturnsDifferentHashCode()
        {
            Assert.That(_sample1.GetHashCode(), Is.Not.EqualTo(_sampleDifferentBioRep.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentTechnicalReplicate_ReturnsDifferentHashCode()
        {
            Assert.That(_sample1.GetHashCode(), Is.Not.EqualTo(_sampleDifferentTechRep.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentFraction_ReturnsDifferentHashCode()
        {
            Assert.That(_sample1.GetHashCode(), Is.Not.EqualTo(_sampleDifferentFraction.GetHashCode()));
        }

        [Test]
        public void GetHashCode_EmptyStrings_ProducesValidHashCode()
        {
            var hash = _sampleEmptyStrings.GetHashCode();
            Assert.That(hash, Is.TypeOf<int>());
        }

        [Test]
        public void GetHashCode_Consistency_MultipleCalls_ReturnsSameValue()
        {
            var hash1 = _sample1.GetHashCode();
            var hash2 = _sample1.GetHashCode();
            var hash3 = _sample1.GetHashCode();
            Assert.That(hash1, Is.EqualTo(hash2));
            Assert.That(hash2, Is.EqualTo(hash3));
        }

        [Test]
        public void GetHashCode_WorksInHashSet()
        {
            var set = new HashSet<TestSampleInfo>
            {
                _sample1,
                _sampleIdentical, // Should not be added (duplicate)
                _sample2
            };

            Assert.That(set.Count, Is.EqualTo(2));
            Assert.That(set.Contains(_sample1), Is.True);
            Assert.That(set.Contains(_sample2), Is.True);
        }

        [Test]
        public void GetHashCode_WorksAsDictionaryKey()
        {
            var dict = new Dictionary<TestSampleInfo, string>
            {
                { _sample1, "First" },
                { _sample2, "Second" }
            };

            Assert.That(dict[_sample1], Is.EqualTo("First"));
            Assert.That(dict[_sampleIdentical], Is.EqualTo("First")); // Same key as _sample1
            Assert.That(dict[_sample2], Is.EqualTo("Second"));
        }

        [Test]
        public void GetHashCode_NegativeValues_ProducesValidHashCode()
        {
            var sample = new TestSampleInfo(@"C:\test.raw", "Control", -1, -2, -3);
            var hash = sample.GetHashCode();
            Assert.That(hash, Is.TypeOf<int>());
        }

        [Test]
        public void GetHashCode_MaxIntValues_ProducesValidHashCode()
        {
            var sample = new TestSampleInfo(@"C:\test.raw", "Control", int.MaxValue, int.MaxValue, int.MaxValue);
            var hash = sample.GetHashCode();
            Assert.That(hash, Is.TypeOf<int>());
        }

        [Test]
        public void GetHashCode_MinIntValues_ProducesValidHashCode()
        {
            var sample = new TestSampleInfo(@"C:\test.raw", "Control", int.MinValue, int.MinValue, int.MinValue);
            var hash = sample.GetHashCode();
            Assert.That(hash, Is.TypeOf<int>());
        }

        #endregion

        #region ToString Tests

        [Test]
        public void ToString_ReturnsExpectedFormat()
        {
            var result = _sample1.ToString();
            Assert.That(result, Is.EqualTo("Control_B1_T1_F0"));
        }

        [Test]
        public void ToString_WithDifferentValues_ReturnsCorrectFormat()
        {
            var result = _sample2.ToString();
            Assert.That(result, Is.EqualTo("Treatment_B2_T2_F1"));
        }

        [Test]
        public void ToString_EmptyCondition_ReturnsFormatWithEmptyCondition()
        {
            var result = _sampleEmptyStrings.ToString();
            Assert.That(result, Is.EqualTo("_B0_T0_F0"));
        }

        [Test]
        public void ToString_NegativeValues_ReturnsCorrectFormat()
        {
            var sample = new TestSampleInfo(@"C:\test.raw", "Control", -1, -2, -3);
            var result = sample.ToString();
            Assert.That(result, Is.EqualTo("Control_B-1_T-2_F-3"));
        }

        [Test]
        public void ToString_LargeValues_ReturnsCorrectFormat()
        {
            var sample = new TestSampleInfo(@"C:\test.raw", "Control", 999, 888, 777);
            var result = sample.ToString();
            Assert.That(result, Is.EqualTo("Control_B999_T888_F777"));
        }

        [Test]
        public void ToString_SpecialCharactersInCondition_ReturnsCorrectFormat()
        {
            var sample = new TestSampleInfo(@"C:\test.raw", "Control-Group_1", 1, 1, 0);
            var result = sample.ToString();
            Assert.That(result, Is.EqualTo("Control-Group_1_B1_T1_F0"));
        }

        [Test]
        public void ToString_WhitespaceCondition_ReturnsCorrectFormat()
        {
            var sample = new TestSampleInfo(@"C:\test.raw", "   ", 1, 1, 0);
            var result = sample.ToString();
            Assert.That(result, Is.EqualTo("   _B1_T1_F0"));
        }

        [Test]
        public void ToString_IsNotNull()
        {
            Assert.That(_sample1.ToString(), Is.Not.Null);
            Assert.That(_sampleEmptyStrings.ToString(), Is.Not.Null);
        }

        [Test]
        public void ToString_IsNotEmpty()
        {
            Assert.That(_sample1.ToString(), Is.Not.Empty);
            Assert.That(_sampleEmptyStrings.ToString(), Is.Not.Empty);
        }

        #endregion

        #region CompareTo Tests

        [Test]
        public void CompareTo_Null_ReturnsNegative()
        {
            Assert.That(_sample1.CompareTo(null), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_SameInstance_ReturnsZero()
        {
            Assert.That(_sample1.CompareTo(_sample1), Is.EqualTo(0));
        }

        [Test]
        public void CompareTo_IdenticalValues_ReturnsZero()
        {
            Assert.That(_sample1.CompareTo(_sampleIdentical), Is.EqualTo(0));
        }

        [Test]
        public void CompareTo_DifferentFilePath_AscendingOrder()
        {
            var sampleA = new TestSampleInfo(@"A:\file.raw", "Control", 1, 1, 0);
            var sampleB = new TestSampleInfo(@"B:\file.raw", "Control", 1, 1, 0);
            var sampleZ = new TestSampleInfo(@"Z:\file.raw", "Control", 1, 1, 0);

            Assert.That(sampleA.CompareTo(sampleB), Is.LessThan(0));
            Assert.That(sampleB.CompareTo(sampleA), Is.GreaterThan(0));
            Assert.That(sampleA.CompareTo(sampleZ), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_DifferentCondition_AscendingOrder()
        {
            var sampleA = new TestSampleInfo(@"C:\file.raw", "Alpha", 1, 1, 0);
            var sampleB = new TestSampleInfo(@"C:\file.raw", "Beta", 1, 1, 0);
            var sampleZ = new TestSampleInfo(@"C:\file.raw", "Zeta", 1, 1, 0);

            Assert.That(sampleA.CompareTo(sampleB), Is.LessThan(0));
            Assert.That(sampleB.CompareTo(sampleA), Is.GreaterThan(0));
            Assert.That(sampleA.CompareTo(sampleZ), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_DifferentBiologicalReplicate_AscendingOrder()
        {
            var sample1 = new TestSampleInfo(@"C:\file.raw", "Control", 1, 1, 0);
            var sample2 = new TestSampleInfo(@"C:\file.raw", "Control", 2, 1, 0);
            var sample10 = new TestSampleInfo(@"C:\file.raw", "Control", 10, 1, 0);

            Assert.That(sample1.CompareTo(sample2), Is.LessThan(0));
            Assert.That(sample2.CompareTo(sample1), Is.GreaterThan(0));
            Assert.That(sample1.CompareTo(sample10), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_DifferentTechnicalReplicate_AscendingOrder()
        {
            var sample1 = new TestSampleInfo(@"C:\file.raw", "Control", 1, 1, 0);
            var sample2 = new TestSampleInfo(@"C:\file.raw", "Control", 1, 2, 0);
            var sample10 = new TestSampleInfo(@"C:\file.raw", "Control", 1, 10, 0);

            Assert.That(sample1.CompareTo(sample2), Is.LessThan(0));
            Assert.That(sample2.CompareTo(sample1), Is.GreaterThan(0));
            Assert.That(sample1.CompareTo(sample10), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_DifferentFraction_AscendingOrder()
        {
            var sample0 = new TestSampleInfo(@"C:\file.raw", "Control", 1, 1, 0);
            var sample1 = new TestSampleInfo(@"C:\file.raw", "Control", 1, 1, 1);
            var sample10 = new TestSampleInfo(@"C:\file.raw", "Control", 1, 1, 10);

            Assert.That(sample0.CompareTo(sample1), Is.LessThan(0));
            Assert.That(sample1.CompareTo(sample0), Is.GreaterThan(0));
            Assert.That(sample0.CompareTo(sample10), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_ComparisonPriority_FilePathFirst()
        {
            var sampleAZZZ = new TestSampleInfo(@"A:\file.raw", "Zeta", 99, 99, 99);
            var sampleBAAA = new TestSampleInfo(@"B:\file.raw", "Alpha", 1, 1, 0);

            // A:\file.raw should come before B:\file.raw regardless of other values
            Assert.That(sampleAZZZ.CompareTo(sampleBAAA), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_ComparisonPriority_ConditionSecond()
        {
            var sampleAlpha = new TestSampleInfo(@"C:\file.raw", "Alpha", 99, 99, 99);
            var sampleBeta = new TestSampleInfo(@"C:\file.raw", "Beta", 1, 1, 0);

            // Same path, but Alpha should come before Beta
            Assert.That(sampleAlpha.CompareTo(sampleBeta), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_ComparisonPriority_BioRepThird()
        {
            var sample1 = new TestSampleInfo(@"C:\file.raw", "Control", 1, 99, 99);
            var sample2 = new TestSampleInfo(@"C:\file.raw", "Control", 2, 1, 0);

            // Same path and condition, but BioRep 1 should come before 2
            Assert.That(sample1.CompareTo(sample2), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_ComparisonPriority_TechRepFourth()
        {
            var sample1 = new TestSampleInfo(@"C:\file.raw", "Control", 1, 1, 99);
            var sample2 = new TestSampleInfo(@"C:\file.raw", "Control", 1, 2, 0);

            // Same path, condition, and BioRep, but TechRep 1 should come before 2
            Assert.That(sample1.CompareTo(sample2), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_ComparisonPriority_FractionLast()
        {
            var sample0 = new TestSampleInfo(@"C:\file.raw", "Control", 1, 1, 0);
            var sample1 = new TestSampleInfo(@"C:\file.raw", "Control", 1, 1, 1);

            // Only Fraction differs, 0 should come before 1
            Assert.That(sample0.CompareTo(sample1), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_NegativeValues_OrderedCorrectly()
        {
            var sampleNeg = new TestSampleInfo(@"C:\file.raw", "Control", -1, -1, -1);
            var sampleZero = new TestSampleInfo(@"C:\file.raw", "Control", 0, 0, 0);
            var samplePos = new TestSampleInfo(@"C:\file.raw", "Control", 1, 1, 1);

            Assert.That(sampleNeg.CompareTo(sampleZero), Is.LessThan(0));
            Assert.That(sampleZero.CompareTo(samplePos), Is.LessThan(0));
            Assert.That(sampleNeg.CompareTo(samplePos), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_EmptyStringCondition_ComesBeforeNonEmpty()
        {
            var sampleEmpty = new TestSampleInfo(@"C:\file.raw", string.Empty, 1, 1, 0);
            var sampleAlpha = new TestSampleInfo(@"C:\file.raw", "Alpha", 1, 1, 0);

            Assert.That(sampleEmpty.CompareTo(sampleAlpha), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_Sorting_ProducesCorrectOrder()
        {
            var samples = new List<TestSampleInfo>
            {
                new(@"C:\file.raw", "Control", 2, 1, 0),
                new(@"C:\file.raw", "Control", 1, 2, 0),
                new(@"C:\file.raw", "Control", 1, 1, 1),
                new(@"C:\file.raw", "Control", 1, 1, 0),
                new(@"B:\file.raw", "Control", 1, 1, 0),
                new(@"A:\file.raw", "Control", 1, 1, 0),
            };

            samples.Sort((a, b) => a.CompareTo(b));

            Assert.That(samples[0].FullFilePathWithExtension, Is.EqualTo(@"A:\file.raw"));
            Assert.That(samples[1].FullFilePathWithExtension, Is.EqualTo(@"B:\file.raw"));
            Assert.That(samples[2].BiologicalReplicate, Is.EqualTo(1));
            Assert.That(samples[2].TechnicalReplicate, Is.EqualTo(1));
            Assert.That(samples[2].Fraction, Is.EqualTo(0));
            Assert.That(samples[5].BiologicalReplicate, Is.EqualTo(2));
        }

        [Test]
        public void CompareTo_IsTransitive()
        {
            var a = new TestSampleInfo(@"A:\file.raw", "Alpha", 1, 1, 0);
            var b = new TestSampleInfo(@"B:\file.raw", "Beta", 2, 2, 1);
            var c = new TestSampleInfo(@"C:\file.raw", "Gamma", 3, 3, 2);

            // If a < b and b < c, then a < c
            Assert.That(a.CompareTo(b), Is.LessThan(0));
            Assert.That(b.CompareTo(c), Is.LessThan(0));
            Assert.That(a.CompareTo(c), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_IsAntiSymmetric()
        {
            var a = new TestSampleInfo(@"A:\file.raw", "Control", 1, 1, 0);
            var b = new TestSampleInfo(@"B:\file.raw", "Control", 1, 1, 0);

            // If a < b, then b > a
            var aToB = a.CompareTo(b);
            var bToA = b.CompareTo(a);
            Assert.That(aToB, Is.LessThan(0));
            Assert.That(bToA, Is.GreaterThan(0));
            Assert.That(aToB, Is.EqualTo(-bToA).Or.EqualTo(-1 * Math.Sign(bToA)));
        }

        [Test]
        public void CompareTo_IsReflexive()
        {
            Assert.That(_sample1.CompareTo(_sample1), Is.EqualTo(0));
            Assert.That(_sample2.CompareTo(_sample2), Is.EqualTo(0));
            Assert.That(_sampleEmptyStrings.CompareTo(_sampleEmptyStrings), Is.EqualTo(0));
        }

        [Test]
        public void CompareTo_CaseSensitiveCondition()
        {
            var sampleLower = new TestSampleInfo(@"C:\file.raw", "alpha", 1, 1, 0);
            var sampleUpper = new TestSampleInfo(@"C:\file.raw", "Alpha", 1, 1, 0);

            // In ordinal comparison, uppercase letters come before lowercase
            Assert.That(sampleUpper.CompareTo(sampleLower), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_CaseSensitiveFilePath()
        {
            var sampleLower = new TestSampleInfo(@"c:\file.raw", "Control", 1, 1, 0);
            var sampleUpper = new TestSampleInfo(@"C:\file.raw", "Control", 1, 1, 0);

            // In ordinal comparison, uppercase letters come before lowercase
            Assert.That(sampleUpper.CompareTo(sampleLower), Is.LessThan(0));
        }

        #endregion

        #region Interface Implementation Tests

        [Test]
        public void ImplementsISampleInfo()
        {
            Assert.That(_sample1, Is.InstanceOf<ISampleInfo>());
        }

        [Test]
        public void ImplementsIEquatableOfTestSampleInfo()
        {
            Assert.That(_sample1, Is.InstanceOf<IEquatable<TestSampleInfo>>());
        }

        [Test]
        public void ImplementsIEquatableOfISampleInfo()
        {
            Assert.That(_sample1, Is.InstanceOf<IEquatable<ISampleInfo>>());
        }

        [Test]
        public void ImplementsIComparableOfISampleInfo()
        {
            Assert.That(_sample1, Is.InstanceOf<IComparable<ISampleInfo>>());
        }

        [Test]
        public void ISampleInfo_Equals_WithTypedMethod()
        {
            ISampleInfo sample1 = _sample1;
            ISampleInfo sample2 = _sampleIdentical;
            Assert.That(sample1.Equals(sample2), Is.True);
        }

        [Test]
        public void ISampleInfo_CompareTo_WorksViaInterface()
        {
            ISampleInfo sample1 = _sample1;
            ISampleInfo sample2 = _sample2;
            Assert.That(sample1.CompareTo(sample2), Is.Not.EqualTo(0));
        }

        #endregion
    }
}