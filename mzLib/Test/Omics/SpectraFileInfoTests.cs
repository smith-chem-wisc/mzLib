using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;

namespace Test.Omics
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class SpectraFileInfoTests
    {
        #region Test Data

        private SpectraFileInfo _sample1 = null!;
        private SpectraFileInfo _sample2 = null!;
        private SpectraFileInfo _sampleIdenticalToSample1 = null!;
        private SpectraFileInfo _sampleDifferentFilePath = null!;
        private SpectraFileInfo _sampleDifferentCondition = null!;
        private SpectraFileInfo _sampleDifferentBioRep = null!;
        private SpectraFileInfo _sampleDifferentTechRep = null!;
        private SpectraFileInfo _sampleDifferentFraction = null!;
        private SpectraFileInfo _sampleEmptyStrings = null!;

        [SetUp]
        public void Setup()
        {
            _sample1 = new SpectraFileInfo(
                fullFilePathWithExtension: @"C:\Data\sample1.raw",
                condition: "Control",
                biorep: 1,
                techrep: 1,
                fraction: 0);

            _sample2 = new SpectraFileInfo(
                fullFilePathWithExtension: @"C:\Data\sample2.mzML",
                condition: "Treatment",
                biorep: 2,
                techrep: 2,
                fraction: 1);

            _sampleIdenticalToSample1 = new SpectraFileInfo(
                fullFilePathWithExtension: @"C:\Data\sample1.raw",
                condition: "Control",
                biorep: 1,
                techrep: 1,
                fraction: 0);

            _sampleDifferentFilePath = new SpectraFileInfo(
                fullFilePathWithExtension: @"C:\Data\different.raw",
                condition: "Control",
                biorep: 1,
                techrep: 1,
                fraction: 0);

            _sampleDifferentCondition = new SpectraFileInfo(
                fullFilePathWithExtension: @"C:\Data\sample1.raw",
                condition: "Treatment",
                biorep: 1,
                techrep: 1,
                fraction: 0);

            _sampleDifferentBioRep = new SpectraFileInfo(
                fullFilePathWithExtension: @"C:\Data\sample1.raw",
                condition: "Control",
                biorep: 2,
                techrep: 1,
                fraction: 0);

            _sampleDifferentTechRep = new SpectraFileInfo(
                fullFilePathWithExtension: @"C:\Data\sample1.raw",
                condition: "Control",
                biorep: 1,
                techrep: 2,
                fraction: 0);

            _sampleDifferentFraction = new SpectraFileInfo(
                fullFilePathWithExtension: @"C:\Data\sample1.raw",
                condition: "Control",
                biorep: 1,
                techrep: 1,
                fraction: 1);

            _sampleEmptyStrings = new SpectraFileInfo(
                fullFilePathWithExtension: string.Empty,
                condition: string.Empty,
                biorep: 0,
                techrep: 0,
                fraction: 0);
        }

        #endregion

        #region Constructor Tests

        [Test]
        public void Constructor_WithValidParameters_SetsAllProperties()
        {
            var sample = new SpectraFileInfo(
                fullFilePathWithExtension: @"C:\Data\Experiment\test_sample.raw",
                condition: "Control",
                biorep: 1,
                techrep: 2,
                fraction: 3);

            Assert.That(sample.FullFilePathWithExtension, Is.EqualTo(@"C:\Data\Experiment\test_sample.raw"));
            Assert.That(sample.FilenameWithoutExtension, Is.EqualTo("test_sample"));
            Assert.That(sample.Condition, Is.EqualTo("Control"));
            Assert.That(sample.BiologicalReplicate, Is.EqualTo(1));
            Assert.That(sample.TechnicalReplicate, Is.EqualTo(2));
            Assert.That(sample.Fraction, Is.EqualTo(3));
        }

        [Test]
        public void Constructor_FilenameWithoutExtension_ExtractsCorrectly()
        {
            var sample = new SpectraFileInfo(@"C:\Folder\SubFolder\MyFile.raw", "Control", 1, 1, 0);
            Assert.That(sample.FilenameWithoutExtension, Is.EqualTo("MyFile"));
        }

        [Test]
        public void Constructor_FilenameWithoutExtension_HandlesMultipleDots()
        {
            var sample = new SpectraFileInfo(@"C:\Data\sample.test.raw", "Control", 1, 1, 0);
            Assert.That(sample.FilenameWithoutExtension, Is.EqualTo("sample.test"));
        }

        [Test]
        public void Constructor_FilenameWithoutExtension_HandlesMzML()
        {
            var sample = new SpectraFileInfo(@"C:\Data\sample.mzML", "Control", 1, 1, 0);
            Assert.That(sample.FilenameWithoutExtension, Is.EqualTo("sample"));
        }

        [Test]
        public void Constructor_FilenameWithoutExtension_HandlesNoExtension()
        {
            var sample = new SpectraFileInfo(@"C:\Data\samplefile", "Control", 1, 1, 0);
            Assert.That(sample.FilenameWithoutExtension, Is.EqualTo("samplefile"));
        }

        [Test]
        public void Constructor_WithEmptyFilePath_SetsEmptyFilename()
        {
            var sample = new SpectraFileInfo(string.Empty, "Control", 1, 1, 0);
            Assert.That(sample.FullFilePathWithExtension, Is.EqualTo(string.Empty));
            Assert.That(sample.FilenameWithoutExtension, Is.EqualTo(string.Empty));
        }

        [Test]
        public void Constructor_WithNullFilePath_SetsNullValues()
        {
            // In .NET 8, Path.GetFileNameWithoutExtension(null) returns null, not an exception
            var sample = new SpectraFileInfo(null!, "Control", 1, 1, 0);
            Assert.That(sample.FullFilePathWithExtension, Is.Null);
            Assert.That(sample.FilenameWithoutExtension, Is.Null);
        }

        [Test]
        public void Constructor_WithNullCondition_SetsNullCondition()
        {
            var sample = new SpectraFileInfo(@"C:\Data\sample.raw", null!, 1, 1, 0);
            Assert.That(sample.Condition, Is.Null);
        }

        [Test]
        public void Constructor_WithNegativeValues_AcceptsValues()
        {
            var sample = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", -1, -2, -3);
            Assert.That(sample.BiologicalReplicate, Is.EqualTo(-1));
            Assert.That(sample.TechnicalReplicate, Is.EqualTo(-2));
            Assert.That(sample.Fraction, Is.EqualTo(-3));
        }

        [Test]
        public void Constructor_WithZeroValues_AcceptsValues()
        {
            var sample = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 0, 0, 0);
            Assert.That(sample.BiologicalReplicate, Is.EqualTo(0));
            Assert.That(sample.TechnicalReplicate, Is.EqualTo(0));
            Assert.That(sample.Fraction, Is.EqualTo(0));
        }

        [Test]
        public void Constructor_WithMaxIntValues_AcceptsValues()
        {
            var sample = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", int.MaxValue, int.MaxValue, int.MaxValue);
            Assert.That(sample.BiologicalReplicate, Is.EqualTo(int.MaxValue));
            Assert.That(sample.TechnicalReplicate, Is.EqualTo(int.MaxValue));
            Assert.That(sample.Fraction, Is.EqualTo(int.MaxValue));
        }

        [Test]
        public void Constructor_WithMinIntValues_AcceptsValues()
        {
            var sample = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", int.MinValue, int.MinValue, int.MinValue);
            Assert.That(sample.BiologicalReplicate, Is.EqualTo(int.MinValue));
            Assert.That(sample.TechnicalReplicate, Is.EqualTo(int.MinValue));
            Assert.That(sample.Fraction, Is.EqualTo(int.MinValue));
        }

        [Test]
        public void Constructor_WithUnixStylePath_HandlesCorrectly()
        {
            var sample = new SpectraFileInfo("/home/user/data/sample.raw", "Control", 1, 1, 0);
            Assert.That(sample.FilenameWithoutExtension, Is.EqualTo("sample"));
        }

        [Test]
        public void Constructor_WithRelativePath_HandlesCorrectly()
        {
            var sample = new SpectraFileInfo(@".\data\sample.raw", "Control", 1, 1, 0);
            Assert.That(sample.FilenameWithoutExtension, Is.EqualTo("sample"));
        }

        [Test]
        public void Constructor_WithNetworkPath_HandlesCorrectly()
        {
            var sample = new SpectraFileInfo(@"\\server\share\data\sample.raw", "Control", 1, 1, 0);
            Assert.That(sample.FilenameWithoutExtension, Is.EqualTo("sample"));
        }

        #endregion

        #region Equals(SpectraFileInfo) Tests

        [Test]
        public void Equals_Typed_SameReference_ReturnsTrue()
        {
            Assert.That(_sample1.Equals(_sample1), Is.True);
        }

        [Test]
        public void Equals_Typed_IdenticalValues_ReturnsTrue()
        {
            Assert.That(_sample1.Equals(_sampleIdenticalToSample1), Is.True);
        }

        [Test]
        public void Equals_Typed_Null_ReturnsFalse()
        {
            Assert.That(_sample1.Equals((SpectraFileInfo?)null), Is.False);
        }

        [Test]
        public void Equals_Typed_DifferentFilePath_ReturnsFalse()
        {
            Assert.That(_sample1.Equals(_sampleDifferentFilePath), Is.False);
        }

        [Test]
        public void Equals_Typed_DifferentCondition_ReturnsFalse()
        {
            Assert.That(_sample1.Equals(_sampleDifferentCondition), Is.False);
        }

        [Test]
        public void Equals_Typed_DifferentBiologicalReplicate_ReturnsFalse()
        {
            Assert.That(_sample1.Equals(_sampleDifferentBioRep), Is.False);
        }

        [Test]
        public void Equals_Typed_DifferentTechnicalReplicate_ReturnsFalse()
        {
            Assert.That(_sample1.Equals(_sampleDifferentTechRep), Is.False);
        }

        [Test]
        public void Equals_Typed_DifferentFraction_ReturnsFalse()
        {
            Assert.That(_sample1.Equals(_sampleDifferentFraction), Is.False);
        }

        [Test]
        public void Equals_Typed_CaseSensitiveFilePath()
        {
            var sampleLower = new SpectraFileInfo(@"c:\data\sample.raw", "Control", 1, 1, 0);
            var sampleUpper = new SpectraFileInfo(@"C:\Data\SAMPLE.raw", "Control", 1, 1, 0);
            Assert.That(sampleLower.Equals(sampleUpper), Is.False);
        }

        [Test]
        public void Equals_Typed_CaseSensitiveCondition()
        {
            var sampleLower = new SpectraFileInfo(@"C:\Data\sample.raw", "control", 1, 1, 0);
            var sampleUpper = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            Assert.That(sampleLower.Equals(sampleUpper), Is.False);
        }

        [Test]
        public void Equals_Typed_EmptyStrings_BothEmpty_ReturnsTrue()
        {
            var sample1 = new SpectraFileInfo(string.Empty, string.Empty, 0, 0, 0);
            var sample2 = new SpectraFileInfo(string.Empty, string.Empty, 0, 0, 0);
            Assert.That(sample1.Equals(sample2), Is.True);
        }

        #endregion

        #region Equals(ISampleInfo) Tests

        [Test]
        public void Equals_ISampleInfo_SameSpectraFileInfo_ReturnsTrue()
        {
            ISampleInfo other = _sampleIdenticalToSample1;
            Assert.That(_sample1.Equals(other), Is.True);
        }

        [Test]
        public void Equals_ISampleInfo_Null_ReturnsFalse()
        {
            Assert.That(_sample1.Equals((ISampleInfo?)null), Is.False);
        }

        [Test]
        public void Equals_ISampleInfo_DifferentSpectraFileInfo_ReturnsFalse()
        {
            ISampleInfo other = _sample2;
            Assert.That(_sample1.Equals(other), Is.False);
        }

        [Test]
        public void Equals_ISampleInfo_IsobaricQuantSampleInfo_ReturnsFalse()
        {
            var isobaric = new IsobaricQuantSampleInfo(
                @"C:\Data\sample1.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            Assert.That(_sample1.Equals(isobaric), Is.False);
        }

        #endregion

        #region Equals(object) Tests

        [Test]
        public void Equals_Object_SameReference_ReturnsTrue()
        {
            Assert.That(_sample1.Equals((object)_sample1), Is.True);
        }

        [Test]
        public void Equals_Object_IdenticalValues_ReturnsTrue()
        {
            Assert.That(_sample1.Equals((object)_sampleIdenticalToSample1), Is.True);
        }

        [Test]
        public void Equals_Object_Null_ReturnsFalse()
        {
            Assert.That(_sample1.Equals((object?)null), Is.False);
        }

        [Test]
        public void Equals_Object_DifferentType_ReturnsFalse()
        {
            Assert.That(_sample1.Equals("not a SpectraFileInfo"), Is.False);
        }

        [Test]
        public void Equals_Object_IntegerType_ReturnsFalse()
        {
            Assert.That(_sample1.Equals(123), Is.False);
        }

        [Test]
        public void Equals_Object_IsobaricQuantSampleInfo_ReturnsFalse()
        {
            var isobaric = new IsobaricQuantSampleInfo(
                @"C:\Data\sample1.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            Assert.That(_sample1.Equals((object)isobaric), Is.False);
        }

        #endregion

        #region GetHashCode Tests

        [Test]
        public void GetHashCode_IdenticalValues_ReturnsSameHashCode()
        {
            Assert.That(_sample1.GetHashCode(), Is.EqualTo(_sampleIdenticalToSample1.GetHashCode()));
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
            var set = new HashSet<SpectraFileInfo>
            {
                _sample1,
                _sampleIdenticalToSample1, // Should not be added (duplicate)
                _sample2
            };

            Assert.That(set.Count, Is.EqualTo(2));
            Assert.That(set.Contains(_sample1), Is.True);
            Assert.That(set.Contains(_sample2), Is.True);
            Assert.That(set.Contains(_sampleIdenticalToSample1), Is.True);
        }

        [Test]
        public void GetHashCode_WorksAsDictionaryKey()
        {
            var dict = new Dictionary<SpectraFileInfo, string>
            {
                { _sample1, "First" },
                { _sample2, "Second" }
            };

            Assert.That(dict[_sample1], Is.EqualTo("First"));
            Assert.That(dict[_sampleIdenticalToSample1], Is.EqualTo("First"));
            Assert.That(dict[_sample2], Is.EqualTo("Second"));
        }

        [Test]
        public void GetHashCode_NegativeValues_ProducesValidHashCode()
        {
            var sample = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", -1, -2, -3);
            var hash = sample.GetHashCode();
            Assert.That(hash, Is.TypeOf<int>());
        }

        [Test]
        public void GetHashCode_MaxIntValues_ProducesValidHashCode()
        {
            var sample = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", int.MaxValue, int.MaxValue, int.MaxValue);
            var hash = sample.GetHashCode();
            Assert.That(hash, Is.TypeOf<int>());
        }

        [Test]
        public void GetHashCode_MinIntValues_ProducesValidHashCode()
        {
            var sample = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", int.MinValue, int.MinValue, int.MinValue);
            var hash = sample.GetHashCode();
            Assert.That(hash, Is.TypeOf<int>());
        }

        [Test]
        public void GetHashCode_NullCondition_ProducesValidHashCode()
        {
            var sample = new SpectraFileInfo(@"C:\Data\sample.raw", null!, 1, 1, 0);
            var hash = sample.GetHashCode();
            Assert.That(hash, Is.TypeOf<int>());
        }

        #endregion

        #region ToString Tests

        [Test]
        public void ToString_ReturnsFileNameWithExtension()
        {
            Assert.That(_sample1.ToString(), Is.EqualTo("sample1.raw"));
        }

        [Test]
        public void ToString_WithMzMLExtension_ReturnsCorrectFileName()
        {
            Assert.That(_sample2.ToString(), Is.EqualTo("sample2.mzML"));
        }

        [Test]
        public void ToString_WithEmptyPath_ReturnsEmptyString()
        {
            Assert.That(_sampleEmptyStrings.ToString(), Is.EqualTo(string.Empty));
        }

        [Test]
        public void ToString_WithDeepPath_ReturnsOnlyFileName()
        {
            var sample = new SpectraFileInfo(@"C:\Very\Deep\Folder\Structure\file.raw", "Control", 1, 1, 0);
            Assert.That(sample.ToString(), Is.EqualTo("file.raw"));
        }

        [Test]
        public void ToString_WithNoExtension_ReturnsFileName()
        {
            var sample = new SpectraFileInfo(@"C:\Data\filename", "Control", 1, 1, 0);
            Assert.That(sample.ToString(), Is.EqualTo("filename"));
        }

        [Test]
        public void ToString_WithMultipleDots_ReturnsFullFileName()
        {
            var sample = new SpectraFileInfo(@"C:\Data\sample.test.data.raw", "Control", 1, 1, 0);
            Assert.That(sample.ToString(), Is.EqualTo("sample.test.data.raw"));
        }

        [Test]
        public void ToString_WithUnixPath_ReturnsFileName()
        {
            var sample = new SpectraFileInfo("/home/user/data/sample.raw", "Control", 1, 1, 0);
            Assert.That(sample.ToString(), Is.EqualTo("sample.raw"));
        }

        [Test]
        public void ToString_WithNetworkPath_ReturnsFileName()
        {
            var sample = new SpectraFileInfo(@"\\server\share\data\sample.raw", "Control", 1, 1, 0);
            Assert.That(sample.ToString(), Is.EqualTo("sample.raw"));
        }

        [Test]
        public void ToString_IsNotNull()
        {
            Assert.That(_sample1.ToString(), Is.Not.Null);
        }

        [Test]
        public void ToString_WithSpecialCharactersInFileName_ReturnsCorrectly()
        {
            var sample = new SpectraFileInfo(@"C:\Data\sample-test_01 (copy).raw", "Control", 1, 1, 0);
            Assert.That(sample.ToString(), Is.EqualTo("sample-test_01 (copy).raw"));
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
            Assert.That(_sample1.CompareTo(_sampleIdenticalToSample1), Is.EqualTo(0));
        }

        [Test]
        public void CompareTo_DifferentFilePath_AscendingOrder()
        {
            var sampleA = new SpectraFileInfo(@"A:\file.raw", "Control", 1, 1, 0);
            var sampleB = new SpectraFileInfo(@"B:\file.raw", "Control", 1, 1, 0);
            var sampleZ = new SpectraFileInfo(@"Z:\file.raw", "Control", 1, 1, 0);

            Assert.That(sampleA.CompareTo(sampleB), Is.LessThan(0));
            Assert.That(sampleB.CompareTo(sampleA), Is.GreaterThan(0));
            Assert.That(sampleA.CompareTo(sampleZ), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_DifferentCondition_AscendingOrder()
        {
            var sampleAlpha = new SpectraFileInfo(@"C:\file.raw", "Alpha", 1, 1, 0);
            var sampleBeta = new SpectraFileInfo(@"C:\file.raw", "Beta", 1, 1, 0);
            var sampleZeta = new SpectraFileInfo(@"C:\file.raw", "Zeta", 1, 1, 0);

            Assert.That(sampleAlpha.CompareTo(sampleBeta), Is.LessThan(0));
            Assert.That(sampleBeta.CompareTo(sampleAlpha), Is.GreaterThan(0));
            Assert.That(sampleAlpha.CompareTo(sampleZeta), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_DifferentBiologicalReplicate_AscendingOrder()
        {
            var sample1 = new SpectraFileInfo(@"C:\file.raw", "Control", 1, 1, 0);
            var sample2 = new SpectraFileInfo(@"C:\file.raw", "Control", 2, 1, 0);
            var sample10 = new SpectraFileInfo(@"C:\file.raw", "Control", 10, 1, 0);

            Assert.That(sample1.CompareTo(sample2), Is.LessThan(0));
            Assert.That(sample2.CompareTo(sample1), Is.GreaterThan(0));
            Assert.That(sample1.CompareTo(sample10), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_DifferentFraction_AscendingOrder()
        {
            var sample0 = new SpectraFileInfo(@"C:\file.raw", "Control", 1, 1, 0);
            var sample1 = new SpectraFileInfo(@"C:\file.raw", "Control", 1, 1, 1);
            var sample10 = new SpectraFileInfo(@"C:\file.raw", "Control", 1, 1, 10);

            Assert.That(sample0.CompareTo(sample1), Is.LessThan(0));
            Assert.That(sample1.CompareTo(sample0), Is.GreaterThan(0));
            Assert.That(sample0.CompareTo(sample10), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_DifferentTechnicalReplicate_AscendingOrder()
        {
            var sample1 = new SpectraFileInfo(@"C:\file.raw", "Control", 1, 1, 0);
            var sample2 = new SpectraFileInfo(@"C:\file.raw", "Control", 1, 2, 0);
            var sample10 = new SpectraFileInfo(@"C:\file.raw", "Control", 1, 10, 0);

            Assert.That(sample1.CompareTo(sample2), Is.LessThan(0));
            Assert.That(sample2.CompareTo(sample1), Is.GreaterThan(0));
            Assert.That(sample1.CompareTo(sample10), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_ComparisonPriority_FilePathFirst()
        {
            var sampleAZZZ = new SpectraFileInfo(@"A:\file.raw", "Zeta", 99, 99, 99);
            var sampleBAAA = new SpectraFileInfo(@"B:\file.raw", "Alpha", 1, 1, 0);

            Assert.That(sampleAZZZ.CompareTo(sampleBAAA), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_ComparisonPriority_ConditionSecond()
        {
            var sampleAlpha = new SpectraFileInfo(@"C:\file.raw", "Alpha", 99, 99, 99);
            var sampleBeta = new SpectraFileInfo(@"C:\file.raw", "Beta", 1, 1, 0);

            Assert.That(sampleAlpha.CompareTo(sampleBeta), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_ComparisonPriority_BioRepThird()
        {
            var sample1 = new SpectraFileInfo(@"C:\file.raw", "Control", 1, 99, 99);
            var sample2 = new SpectraFileInfo(@"C:\file.raw", "Control", 2, 1, 0);

            Assert.That(sample1.CompareTo(sample2), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_ComparisonPriority_FractionFourth()
        {
            var sample0 = new SpectraFileInfo(@"C:\file.raw", "Control", 1, 99, 0);
            var sample1 = new SpectraFileInfo(@"C:\file.raw", "Control", 1, 1, 1);

            Assert.That(sample0.CompareTo(sample1), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_ComparisonPriority_TechRepLast()
        {
            var sample1 = new SpectraFileInfo(@"C:\file.raw", "Control", 1, 1, 0);
            var sample2 = new SpectraFileInfo(@"C:\file.raw", "Control", 1, 2, 0);

            Assert.That(sample1.CompareTo(sample2), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_NegativeValues_OrderedCorrectly()
        {
            var sampleNeg = new SpectraFileInfo(@"C:\file.raw", "Control", -1, -1, -1);
            var sampleZero = new SpectraFileInfo(@"C:\file.raw", "Control", 0, 0, 0);
            var samplePos = new SpectraFileInfo(@"C:\file.raw", "Control", 1, 1, 1);

            Assert.That(sampleNeg.CompareTo(sampleZero), Is.LessThan(0));
            Assert.That(sampleZero.CompareTo(samplePos), Is.LessThan(0));
            Assert.That(sampleNeg.CompareTo(samplePos), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_EmptyStringCondition_ComesBeforeNonEmpty()
        {
            var sampleEmpty = new SpectraFileInfo(@"C:\file.raw", string.Empty, 1, 1, 0);
            var sampleAlpha = new SpectraFileInfo(@"C:\file.raw", "Alpha", 1, 1, 0);

            Assert.That(sampleEmpty.CompareTo(sampleAlpha), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_Sorting_ProducesCorrectOrder()
        {
            var samples = new List<SpectraFileInfo>
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
            Assert.That(samples[2].Fraction, Is.EqualTo(0));
            Assert.That(samples[2].TechnicalReplicate, Is.EqualTo(1));
            Assert.That(samples[5].BiologicalReplicate, Is.EqualTo(2));
        }

        [Test]
        public void CompareTo_IsTransitive()
        {
            var a = new SpectraFileInfo(@"A:\file.raw", "Alpha", 1, 1, 0);
            var b = new SpectraFileInfo(@"B:\file.raw", "Beta", 2, 2, 1);
            var c = new SpectraFileInfo(@"C:\file.raw", "Gamma", 3, 3, 2);

            Assert.That(a.CompareTo(b), Is.LessThan(0));
            Assert.That(b.CompareTo(c), Is.LessThan(0));
            Assert.That(a.CompareTo(c), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_IsAntiSymmetric()
        {
            var a = new SpectraFileInfo(@"A:\file.raw", "Control", 1, 1, 0);
            var b = new SpectraFileInfo(@"B:\file.raw", "Control", 1, 1, 0);

            var aToB = a.CompareTo(b);
            var bToA = b.CompareTo(a);
            Assert.That(aToB, Is.LessThan(0));
            Assert.That(bToA, Is.GreaterThan(0));
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
            var sampleLower = new SpectraFileInfo(@"C:\file.raw", "alpha", 1, 1, 0);
            var sampleUpper = new SpectraFileInfo(@"C:\file.raw", "Alpha", 1, 1, 0);

            // Uppercase comes before lowercase in ordinal comparison
            Assert.That(sampleUpper.CompareTo(sampleLower), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_CaseSensitiveFilePath()
        {
            var sampleLower = new SpectraFileInfo(@"c:\file.raw", "Control", 1, 1, 0);
            var sampleUpper = new SpectraFileInfo(@"C:\file.raw", "Control", 1, 1, 0);

            // Uppercase comes before lowercase in ordinal comparison
            Assert.That(sampleUpper.CompareTo(sampleLower), Is.LessThan(0));
        }
        [Test]
        public void CompareTo_WithIsobaricQuantSampleInfo_CrossTypeOrdering()
        {
            // When comparing SpectraFileInfo to a non-SpectraFileInfo ISampleInfo,
            // SpectraFileInfo.CompareTo returns positive (sorts after other types).
            // This is by design to provide consistent ordering within same-type collections.
            var isobaric = new IsobaricQuantSampleInfo(
                @"C:\Data\sample1.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);

            Assert.That(_sample1.CompareTo(isobaric), Is.GreaterThan(0));
        }

        [Test]
        public void CompareTo_WithIsobaricQuantSampleInfo_DifferentProperties_ReturnsPositive()
        {
            // Cross-type comparison always returns positive regardless of property values
            var isobaric = new IsobaricQuantSampleInfo(
                @"D:\Data\different.raw", "Treatment", 2, 2, 1, 1, "126", 126.0, false);

            Assert.That(_sample1.CompareTo(isobaric), Is.GreaterThan(0));
        }

        #endregion

        #region Interface Implementation Tests

        [Test]
        public void ImplementsISampleInfo()
        {
            Assert.That(_sample1, Is.InstanceOf<ISampleInfo>());
        }

        [Test]
        public void ImplementsIEquatableOfSpectraFileInfo()
        {
            Assert.That(_sample1, Is.InstanceOf<IEquatable<SpectraFileInfo>>());
        }

        [Test]
        public void ImplementsIComparableOfISampleInfo()
        {
            Assert.That(_sample1, Is.InstanceOf<IComparable<ISampleInfo>>());
        }

        [Test]
        public void ISampleInfo_Properties_AreAccessibleViaInterface()
        {
            ISampleInfo sample = _sample1;

            Assert.That(sample.FullFilePathWithExtension, Is.EqualTo(@"C:\Data\sample1.raw"));
            Assert.That(sample.Condition, Is.EqualTo("Control"));
            Assert.That(sample.BiologicalReplicate, Is.EqualTo(1));
            Assert.That(sample.TechnicalReplicate, Is.EqualTo(1));
            Assert.That(sample.Fraction, Is.EqualTo(0));
        }

        [Test]
        public void ISampleInfo_CompareTo_WorksViaInterface()
        {
            ISampleInfo sample1 = _sample1;
            ISampleInfo sample2 = _sample2;
            Assert.That(sample1.CompareTo(sample2), Is.Not.EqualTo(0));
        }

        [Test]
        public void ISampleInfo_Equals_WorksViaInterface()
        {
            ISampleInfo sample1 = _sample1;
            ISampleInfo identical = _sampleIdenticalToSample1;
            Assert.That(sample1.Equals(identical), Is.True);
        }

        #endregion

        #region Edge Case Tests

        [Test]
        public void EdgeCase_WhitespaceCondition()
        {
            var sample = new SpectraFileInfo(@"C:\Data\sample.raw", "   ", 1, 1, 0);
            Assert.That(sample.Condition, Is.EqualTo("   "));
        }

        [Test]
        public void EdgeCase_VeryLongFilePath()
        {
            var longPath = @"C:\" + new string('a', 200) + @"\file.raw";
            var sample = new SpectraFileInfo(longPath, "Control", 1, 1, 0);
            Assert.That(sample.FullFilePathWithExtension, Is.EqualTo(longPath));
            Assert.That(sample.ToString(), Is.EqualTo("file.raw"));
        }

        [Test]
        public void EdgeCase_SpecialCharactersInCondition()
        {
            var sample = new SpectraFileInfo(@"C:\Data\sample.raw", "Control-Group_1/A", 1, 1, 0);
            Assert.That(sample.Condition, Is.EqualTo("Control-Group_1/A"));
        }

        [Test]
        public void EdgeCase_UnicodeInFilePath()
        {
            var sample = new SpectraFileInfo(@"C:\Data\αβγ_sample.raw", "Control", 1, 1, 0);
            Assert.That(sample.FilenameWithoutExtension, Is.EqualTo("αβγ_sample"));
        }

        [Test]
        public void EdgeCase_UnicodeInCondition()
        {
            var sample = new SpectraFileInfo(@"C:\Data\sample.raw", "Контроль", 1, 1, 0);
            Assert.That(sample.Condition, Is.EqualTo("Контроль"));
        }

        [Test]
        public void EdgeCase_MultipleIdenticalSamplesInList()
        {
            var list = new List<SpectraFileInfo>
            {
                _sample1,
                _sampleIdenticalToSample1,
                _sample1
            };

            var distinctCount = list.Distinct().Count();
            Assert.That(distinctCount, Is.EqualTo(1));
        }

        [Test]
        public void EdgeCase_AllFieldsDifferent_NotEqual()
        {
            Assert.That(_sample1.Equals(_sample2), Is.False);
        }

        [Test]
        public void EdgeCase_OnlyOneFieldDifferent_NotEqual()
        {
            // Each of these differs from _sample1 in only one field
            Assert.That(_sample1.Equals(_sampleDifferentFilePath), Is.False);
            Assert.That(_sample1.Equals(_sampleDifferentCondition), Is.False);
            Assert.That(_sample1.Equals(_sampleDifferentBioRep), Is.False);
            Assert.That(_sample1.Equals(_sampleDifferentTechRep), Is.False);
            Assert.That(_sample1.Equals(_sampleDifferentFraction), Is.False);
        }

        [Test]
        public void EdgeCase_FilePathWithSpaces()
        {
            var sample = new SpectraFileInfo(@"C:\My Data\Sample Files\test sample.raw", "Control", 1, 1, 0);
            Assert.That(sample.FilenameWithoutExtension, Is.EqualTo("test sample"));
            Assert.That(sample.ToString(), Is.EqualTo("test sample.raw"));
        }

        [Test]
        public void EdgeCase_DifferentExtensions_SameBaseName()
        {
            var sampleRaw = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            var sampleMzML = new SpectraFileInfo(@"C:\Data\sample.mzML", "Control", 1, 1, 0);

            // Same filename without extension, but different full paths
            Assert.That(sampleRaw.FilenameWithoutExtension, Is.EqualTo(sampleMzML.FilenameWithoutExtension));
            Assert.That(sampleRaw.Equals(sampleMzML), Is.False); // Different FullFilePathWithExtension
        }

        #endregion

        #region Cross-Type Comparison Tests

        [Test]
        public void CrossType_SpectraFileInfo_NotEqualTo_IsobaricQuantSampleInfo()
        {
            var spectra = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            var isobaric = new IsobaricQuantSampleInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);

            Assert.That(spectra.Equals(isobaric), Is.False);
            Assert.That(spectra.Equals((ISampleInfo)isobaric), Is.False);
            Assert.That(spectra.Equals((object)isobaric), Is.False);
        }

        [Test]
        public void CrossType_InHashSet_TreatedAsDifferent()
        {
            var spectra = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            var isobaric = new IsobaricQuantSampleInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);

            var set = new HashSet<ISampleInfo> { spectra, isobaric };
            Assert.That(set.Count, Is.EqualTo(2));
        }

        [Test]
        public void CrossType_InDictionary_TreatedAsDifferent()
        {
            var spectra = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            var isobaric = new IsobaricQuantSampleInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);

            // This should work because they're different types with different hash codes
            var dict = new Dictionary<ISampleInfo, string>
            {
                { spectra, "Spectra" },
                { isobaric, "Isobaric" }
            };

            Assert.That(dict.Count, Is.EqualTo(2));
        }

        #endregion
    }
}