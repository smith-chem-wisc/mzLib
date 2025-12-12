using MassSpectrometry;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Test.Omics
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class IsobaricQuantSampleInfoTests
    {
        #region Test Data

        private IsobaricQuantSampleInfo _sample1 = null!;
        private IsobaricQuantSampleInfo _sample2 = null!;
        private IsobaricQuantSampleInfo _sampleEqualToSample1 = null!;

        [SetUp]
        public void Setup()
        {
            _sample1 = new IsobaricQuantSampleInfo(
                fullFilePathWithExtension: @"C:\Data\sample1.raw",
                condition: "Control",
                biologicalReplicate: 1,
                technicalReplicate: 1,
                fraction: 0,
                plexId: 1,
                channelLabel: "126",
                reporterIonMz: 126.127726,
                isReferenceChannel: false);

            _sample2 = new IsobaricQuantSampleInfo(
                fullFilePathWithExtension: @"C:\Data\sample2.raw",
                condition: "Treatment",
                biologicalReplicate: 2,
                technicalReplicate: 2,
                fraction: 1,
                plexId: 2,
                channelLabel: "127N",
                reporterIonMz: 127.124761,
                isReferenceChannel: true);

            // Same FullFilePathWithExtension and ChannelLabel as _sample1 (should be equal per Equals)
            _sampleEqualToSample1 = new IsobaricQuantSampleInfo(
                fullFilePathWithExtension: @"C:\Data\sample1.raw",
                condition: "Different",
                biologicalReplicate: 99,
                technicalReplicate: 99,
                fraction: 99,
                plexId: 99,
                channelLabel: "126",
                reporterIonMz: 999.999,
                isReferenceChannel: true);
        }

        #endregion

        #region Constructor Tests

        [Test]
        public void Constructor_WithValidParameters_SetsAllProperties()
        {
            var sample = new IsobaricQuantSampleInfo(
                fullFilePathWithExtension: @"C:\Data\Experiment\test.raw",
                condition: "Control",
                biologicalReplicate: 1,
                technicalReplicate: 2,
                fraction: 3,
                plexId: 4,
                channelLabel: "127N",
                reporterIonMz: 127.124761,
                isReferenceChannel: true);

            Assert.That(sample.FullFilePathWithExtension, Is.EqualTo(@"C:\Data\Experiment\test.raw"));
            Assert.That(sample.Condition, Is.EqualTo("Control"));
            Assert.That(sample.BiologicalReplicate, Is.EqualTo(1));
            Assert.That(sample.TechnicalReplicate, Is.EqualTo(2));
            Assert.That(sample.Fraction, Is.EqualTo(3));
            Assert.That(sample.PlexId, Is.EqualTo(4));
            Assert.That(sample.ChannelLabel, Is.EqualTo("127N"));
            Assert.That(sample.ReporterIonMz, Is.EqualTo(127.124761));
            Assert.That(sample.IsReferenceChannel, Is.True);
        }

        [Test]
        public void Constructor_UniqueIdentifier_ComputedFromFilePathAndChannelLabel()
        {
            var sample1 = new IsobaricQuantSampleInfo(@"C:\a.raw", "A", 1, 1, 0, 1, "126", 126.0, false);
            var sample2 = new IsobaricQuantSampleInfo(@"C:\a.raw", "B", 99, 99, 99, 99, "126", 999.0, true);

            // Same FilePath and ChannelLabel should produce same UniqueIdentifier
            Assert.That(sample1.UniqueIdentifier, Is.EqualTo(sample2.UniqueIdentifier));
        }

        [Test]
        public void Constructor_UniqueIdentifier_DifferentFilePath_DifferentValue()
        {
            var sample1 = new IsobaricQuantSampleInfo(@"C:\a.raw", "A", 1, 1, 0, 1, "126", 126.0, false);
            var sample2 = new IsobaricQuantSampleInfo(@"C:\b.raw", "A", 1, 1, 0, 1, "126", 126.0, false);

            Assert.That(sample1.UniqueIdentifier, Is.Not.EqualTo(sample2.UniqueIdentifier));
        }

        [Test]
        public void Constructor_UniqueIdentifier_DifferentChannelLabel_DifferentValue()
        {
            var sample1 = new IsobaricQuantSampleInfo(@"C:\a.raw", "A", 1, 1, 0, 1, "126", 126.0, false);
            var sample2 = new IsobaricQuantSampleInfo(@"C:\a.raw", "A", 1, 1, 0, 1, "127N", 126.0, false);

            Assert.That(sample1.UniqueIdentifier, Is.Not.EqualTo(sample2.UniqueIdentifier));
        }

        [Test]
        public void Constructor_WithEmptyStrings_AcceptsValues()
        {
            var sample = new IsobaricQuantSampleInfo(
                string.Empty, string.Empty, 0, 0, 0, 0, string.Empty, 0.0, false);

            Assert.That(sample.FullFilePathWithExtension, Is.EqualTo(string.Empty));
            Assert.That(sample.Condition, Is.EqualTo(string.Empty));
            Assert.That(sample.ChannelLabel, Is.EqualTo(string.Empty));
        }

        [Test]
        public void Constructor_WithNegativeValues_AcceptsValues()
        {
            var sample = new IsobaricQuantSampleInfo(
                @"C:\test.raw", "Control", -1, -2, -3, -4, "126", -126.0, false);

            Assert.That(sample.BiologicalReplicate, Is.EqualTo(-1));
            Assert.That(sample.TechnicalReplicate, Is.EqualTo(-2));
            Assert.That(sample.Fraction, Is.EqualTo(-3));
            Assert.That(sample.PlexId, Is.EqualTo(-4));
            Assert.That(sample.ReporterIonMz, Is.EqualTo(-126.0));
        }

        [Test]
        public void Constructor_WithMaxIntValues_AcceptsValues()
        {
            var sample = new IsobaricQuantSampleInfo(
                @"C:\test.raw", "Control", int.MaxValue, int.MaxValue, int.MaxValue, int.MaxValue,
                "126", double.MaxValue, false);

            Assert.That(sample.BiologicalReplicate, Is.EqualTo(int.MaxValue));
            Assert.That(sample.TechnicalReplicate, Is.EqualTo(int.MaxValue));
            Assert.That(sample.Fraction, Is.EqualTo(int.MaxValue));
            Assert.That(sample.PlexId, Is.EqualTo(int.MaxValue));
            Assert.That(sample.ReporterIonMz, Is.EqualTo(double.MaxValue));
        }

        #endregion

        #region Equals(IsobaricQuantSampleInfo) Tests

        [Test]
        public void Equals_Typed_SameReference_ReturnsTrue()
        {
            Assert.That(_sample1.Equals(_sample1), Is.True);
        }

        [Test]
        public void Equals_Typed_SameFilePathAndChannelLabel_ReturnsTrue()
        {
            Assert.That(_sample1.Equals(_sampleEqualToSample1), Is.True);
        }

        [Test]
        public void Equals_Typed_Null_ReturnsFalse()
        {
            Assert.That(_sample1.Equals((IsobaricQuantSampleInfo?)null), Is.False);
        }

        [Test]
        public void Equals_Typed_DifferentFilePath_ReturnsFalse()
        {
            var other = new IsobaricQuantSampleInfo(@"C:\different.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            Assert.That(_sample1.Equals(other), Is.False);
        }

        [Test]
        public void Equals_Typed_DifferentChannelLabel_ReturnsFalse()
        {
            var other = new IsobaricQuantSampleInfo(@"C:\Data\sample1.raw", "Control", 1, 1, 0, 1, "127N", 126.0, false);
            Assert.That(_sample1.Equals(other), Is.False);
        }

        [Test]
        public void Equals_Typed_CaseSensitiveChannelLabel()
        {
            var lower = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "127n", 127.0, false);
            var upper = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "127N", 127.0, false);

            Assert.That(lower.Equals(upper), Is.False);
        }

        [Test]
        public void Equals_Typed_CaseSensitiveFilePath()
        {
            var lower = new IsobaricQuantSampleInfo(@"c:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            var upper = new IsobaricQuantSampleInfo(@"C:\Test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);

            Assert.That(lower.Equals(upper), Is.False);
        }

        [Test]
        public void Equals_Typed_DifferentCondition_SameIdentity_ReturnsTrue()
        {
            var other = new IsobaricQuantSampleInfo(@"C:\Data\sample1.raw", "Treatment", 1, 1, 0, 1, "126", 126.0, false);
            Assert.That(_sample1.Equals(other), Is.True);
        }

        [Test]
        public void Equals_Typed_DifferentPlexId_SameIdentity_ReturnsTrue()
        {
            var other = new IsobaricQuantSampleInfo(@"C:\Data\sample1.raw", "Control", 1, 1, 0, 99, "126", 126.0, false);
            Assert.That(_sample1.Equals(other), Is.True);
        }

        #endregion

        #region Equals(ISampleInfo) Tests

        [Test]
        public void Equals_ISampleInfo_SameIsobaricInstance_ReturnsTrue()
        {
            ISampleInfo other = _sampleEqualToSample1;
            Assert.That(_sample1.Equals(other), Is.True);
        }

        [Test]
        public void Equals_ISampleInfo_Null_ReturnsFalse()
        {
            Assert.That(_sample1.Equals((ISampleInfo?)null), Is.False);
        }

        [Test]
        public void Equals_ISampleInfo_DifferentIsobaric_ReturnsFalse()
        {
            ISampleInfo other = _sample2;
            Assert.That(_sample1.Equals(other), Is.False);
        }

        [Test]
        public void Equals_ISampleInfo_SpectraFileInfo_ReturnsFalse()
        {
            var spectraFile = new SpectraFileInfo(@"C:\Data\sample1.raw", "Control", 1, 1, 0);
            Assert.That(_sample1.Equals(spectraFile), Is.False);
        }

        #endregion

        #region Equals(object) Tests

        [Test]
        public void Equals_Object_SameReference_ReturnsTrue()
        {
            Assert.That(_sample1.Equals((object)_sample1), Is.True);
        }

        [Test]
        public void Equals_Object_EqualInstance_ReturnsTrue()
        {
            Assert.That(_sample1.Equals((object)_sampleEqualToSample1), Is.True);
        }

        [Test]
        public void Equals_Object_Null_ReturnsFalse()
        {
            Assert.That(_sample1.Equals((object?)null), Is.False);
        }

        [Test]
        public void Equals_Object_DifferentType_ReturnsFalse()
        {
            Assert.That(_sample1.Equals("not an IsobaricQuantSampleInfo"), Is.False);
            Assert.That(_sample1.Equals(123), Is.False);
        }

        [Test]
        public void Equals_Object_SpectraFileInfo_ReturnsFalse()
        {
            var spectraFile = new SpectraFileInfo(@"C:\Data\sample1.raw", "Control", 1, 1, 0);
            Assert.That(_sample1.Equals((object)spectraFile), Is.False);
        }

        #endregion

        #region Equality Operator Tests

        [Test]
        public void EqualityOperator_EqualInstances_ReturnsTrue()
        {
            Assert.That(_sample1 == _sampleEqualToSample1, Is.True);
        }

        [Test]
        public void EqualityOperator_DifferentInstances_ReturnsFalse()
        {
            Assert.That(_sample1 == _sample2, Is.False);
        }

        [Test]
        public void EqualityOperator_BothNull_ReturnsTrue()
        {
            IsobaricQuantSampleInfo? left = null;
            IsobaricQuantSampleInfo? right = null;
            Assert.That(left == right, Is.True);
        }

        [Test]
        public void EqualityOperator_LeftNull_ReturnsFalse()
        {
            IsobaricQuantSampleInfo? left = null;
            Assert.That(left == _sample1, Is.False);
        }

        [Test]
        public void EqualityOperator_RightNull_ReturnsFalse()
        {
            IsobaricQuantSampleInfo? right = null;
            Assert.That(_sample1 == right, Is.False);
        }

        #endregion

        #region Inequality Operator Tests

        [Test]
        public void InequalityOperator_DifferentInstances_ReturnsTrue()
        {
            Assert.That(_sample1 != _sample2, Is.True);
        }

        [Test]
        public void InequalityOperator_EqualInstances_ReturnsFalse()
        {
            Assert.That(_sample1 != _sampleEqualToSample1, Is.False);
        }

        [Test]
        public void InequalityOperator_BothNull_ReturnsFalse()
        {
            IsobaricQuantSampleInfo? left = null;
            IsobaricQuantSampleInfo? right = null;
            Assert.That(left != right, Is.False);
        }

        [Test]
        public void InequalityOperator_OneNull_ReturnsTrue()
        {
            IsobaricQuantSampleInfo? nullSample = null;
            Assert.That(nullSample != _sample1, Is.True);
            Assert.That(_sample1 != nullSample, Is.True);
        }

        #endregion

        #region GetHashCode Tests

        [Test]
        public void GetHashCode_EqualInstances_SameHashCode()
        {
            Assert.That(_sample1.GetHashCode(), Is.EqualTo(_sampleEqualToSample1.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentFilePath_DifferentHashCode()
        {
            var other = new IsobaricQuantSampleInfo(@"C:\different.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            Assert.That(_sample1.GetHashCode(), Is.Not.EqualTo(other.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentChannelLabel_DifferentHashCode()
        {
            var other = new IsobaricQuantSampleInfo(@"C:\Data\sample1.raw", "Control", 1, 1, 0, 1, "127N", 126.0, false);
            Assert.That(_sample1.GetHashCode(), Is.Not.EqualTo(other.GetHashCode()));
        }

        [Test]
        public void GetHashCode_MatchesUniqueIdentifier()
        {
            Assert.That(_sample1.GetHashCode(), Is.EqualTo(_sample1.UniqueIdentifier));
        }

        [Test]
        public void GetHashCode_Consistency_MultipleCalls()
        {
            var hash1 = _sample1.GetHashCode();
            var hash2 = _sample1.GetHashCode();
            Assert.That(hash1, Is.EqualTo(hash2));
        }

        [Test]
        public void GetHashCode_WorksInHashSet()
        {
            var set = new HashSet<IsobaricQuantSampleInfo> { _sample1, _sampleEqualToSample1, _sample2 };

            Assert.That(set.Count, Is.EqualTo(2)); // _sample1 and _sampleEqualToSample1 are equal
            Assert.That(set.Contains(_sample1), Is.True);
            Assert.That(set.Contains(_sample2), Is.True);
        }

        [Test]
        public void GetHashCode_WorksAsDictionaryKey()
        {
            var dict = new Dictionary<IsobaricQuantSampleInfo, string>
            {
                { _sample1, "First" },
                { _sample2, "Second" }
            };

            Assert.That(dict[_sample1], Is.EqualTo("First"));
            Assert.That(dict[_sampleEqualToSample1], Is.EqualTo("First")); // Same key as _sample1
            Assert.That(dict[_sample2], Is.EqualTo("Second"));
        }

        #endregion

        #region ToString Tests

        [Test]
        public void ToString_ReturnsFileNameAndChannelLabel()
        {
            Assert.That(_sample1.ToString(), Is.EqualTo("sample1_126"));
        }

        [Test]
        public void ToString_WithDifferentFile_CorrectFormat()
        {
            Assert.That(_sample2.ToString(), Is.EqualTo("sample2_127N"));
        }

        [Test]
        public void ToString_WithEmptyFilePath_ReturnsUnderscoreAndLabel()
        {
            var sample = new IsobaricQuantSampleInfo(string.Empty, "Control", 1, 1, 0, 1, "126", 126.0, false);
            Assert.That(sample.ToString(), Is.EqualTo("_126"));
        }

        [Test]
        public void ToString_WithEmptyChannelLabel_ReturnsFileNameAndUnderscore()
        {
            var sample = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, string.Empty, 126.0, false);
            Assert.That(sample.ToString(), Is.EqualTo("test_"));
        }

        [Test]
        public void ToString_UsesFileNameWithoutExtension()
        {
            var sample = new IsobaricQuantSampleInfo(@"C:\Path\To\MyExperiment.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            Assert.That(sample.ToString(), Is.EqualTo("MyExperiment_126"));
        }

        #endregion

        #region CompareTo Tests - Basic

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
        public void CompareTo_EqualInstance_ReturnsZero()
        {
            // Same FilePath and ChannelLabel triggers Equals check, returns 0
            Assert.That(_sample1.CompareTo(_sampleEqualToSample1), Is.EqualTo(0));
        }

        [Test]
        public void CompareTo_NonIsobaricSampleInfo_ReturnsPositive()
        {
            var spectraFile = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            Assert.That(_sample1.CompareTo(spectraFile), Is.GreaterThan(0));
        }

        #endregion

        #region CompareTo Tests - Comparison Order (Condition → BioRep → Fraction → TechRep → FilePath)

        [Test]
        public void CompareTo_DifferentCondition_OrdersByCondition()
        {
            // Different FilePath to avoid Equals short-circuit
            var alpha = new IsobaricQuantSampleInfo(@"C:\a.raw", "Alpha", 1, 1, 0, 1, "126", 126.0, false);
            var beta = new IsobaricQuantSampleInfo(@"C:\b.raw", "Beta", 1, 1, 0, 1, "127N", 126.0, false);

            Assert.That(alpha.CompareTo(beta), Is.LessThan(0));
            Assert.That(beta.CompareTo(alpha), Is.GreaterThan(0));
        }

        [Test]
        public void CompareTo_SameCondition_DifferentBioRep_OrdersByBioRep()
        {
            // Different FilePath/ChannelLabel to avoid Equals short-circuit
            var bioRep1 = new IsobaricQuantSampleInfo(@"C:\a.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            var bioRep2 = new IsobaricQuantSampleInfo(@"C:\b.raw", "Control", 2, 1, 0, 1, "127N", 126.0, false);

            Assert.That(bioRep1.CompareTo(bioRep2), Is.LessThan(0));
            Assert.That(bioRep2.CompareTo(bioRep1), Is.GreaterThan(0));
        }

        [Test]
        public void CompareTo_SameConditionAndBioRep_DifferentFraction_OrdersByFraction()
        {
            var fraction0 = new IsobaricQuantSampleInfo(@"C:\a.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            var fraction1 = new IsobaricQuantSampleInfo(@"C:\b.raw", "Control", 1, 1, 1, 1, "127N", 126.0, false);

            Assert.That(fraction0.CompareTo(fraction1), Is.LessThan(0));
            Assert.That(fraction1.CompareTo(fraction0), Is.GreaterThan(0));
        }

        [Test]
        public void CompareTo_SameConditionBioRepFraction_DifferentTechRep_OrdersByTechRep()
        {
            var techRep1 = new IsobaricQuantSampleInfo(@"C:\a.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            var techRep2 = new IsobaricQuantSampleInfo(@"C:\b.raw", "Control", 1, 2, 0, 1, "127N", 126.0, false);

            Assert.That(techRep1.CompareTo(techRep2), Is.LessThan(0));
            Assert.That(techRep2.CompareTo(techRep1), Is.GreaterThan(0));
        }

        [Test]
        public void CompareTo_AllSameExceptFilePath_OrdersByFilePath()
        {
            // Different ChannelLabel to avoid Equals short-circuit
            var fileA = new IsobaricQuantSampleInfo(@"C:\a.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            var fileZ = new IsobaricQuantSampleInfo(@"C:\z.raw", "Control", 1, 1, 0, 1, "127N", 126.0, false);

            Assert.That(fileA.CompareTo(fileZ), Is.LessThan(0));
            Assert.That(fileZ.CompareTo(fileA), Is.GreaterThan(0));
        }

        [Test]
        public void CompareTo_ConditionTakesPriority_OverOtherProperties()
        {
            // Alpha condition with "worse" other values should still come first
            var alpha = new IsobaricQuantSampleInfo(@"C:\z.raw", "Alpha", 99, 99, 99, 1, "999", 126.0, false);
            var beta = new IsobaricQuantSampleInfo(@"C:\a.raw", "Beta", 1, 1, 0, 1, "126", 126.0, false);

            Assert.That(alpha.CompareTo(beta), Is.LessThan(0));
        }

        #endregion

        #region CompareTo Tests - Sorting

        [Test]
        public void CompareTo_Sorting_ProducesCorrectOrder()
        {
            var samples = new List<IsobaricQuantSampleInfo>
            {
                new(@"C:\b.raw", "Control", 1, 1, 0, 1, "126", 126.0, false),
                new(@"C:\a.raw", "Control", 1, 1, 0, 1, "131C", 131.0, false),
                new(@"C:\a.raw", "Control", 1, 1, 1, 1, "127N", 127.0, false),
                new(@"C:\a.raw", "Control", 2, 1, 0, 1, "128C", 128.0, false),
                new(@"C:\a.raw", "Alpha", 1, 1, 0, 1, "129N", 129.0, false),
            };

            samples.Sort((a, b) => a.CompareTo(b));

            // Expected order: Alpha first, then Control sorted by BioRep/Fraction/TechRep/FilePath
            Assert.That(samples[0].Condition, Is.EqualTo("Alpha"));
            Assert.That(samples[1].Condition, Is.EqualTo("Control"));
            Assert.That(samples[1].FullFilePathWithExtension, Is.EqualTo(@"C:\a.raw"));
            Assert.That(samples[1].ChannelLabel, Is.EqualTo("131C"));
            Assert.That(samples[2].FullFilePathWithExtension, Is.EqualTo(@"C:\b.raw"));
            Assert.That(samples[3].Fraction, Is.EqualTo(1));
            Assert.That(samples[4].BiologicalReplicate, Is.EqualTo(2));
        }

        [Test]
        public void CompareTo_IsTransitive()
        {
            var a = new IsobaricQuantSampleInfo(@"C:\a.raw", "Alpha", 1, 1, 0, 1, "126", 126.0, false);
            var b = new IsobaricQuantSampleInfo(@"C:\b.raw", "Beta", 1, 1, 0, 1, "127N", 126.0, false);
            var c = new IsobaricQuantSampleInfo(@"C:\c.raw", "Gamma", 1, 1, 0, 1, "128C", 126.0, false);

            Assert.That(a.CompareTo(b), Is.LessThan(0));
            Assert.That(b.CompareTo(c), Is.LessThan(0));
            Assert.That(a.CompareTo(c), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_IsAntiSymmetric()
        {
            var a = new IsobaricQuantSampleInfo(@"C:\a.raw", "Alpha", 1, 1, 0, 1, "126", 126.0, false);
            var b = new IsobaricQuantSampleInfo(@"C:\b.raw", "Beta", 1, 1, 0, 1, "127N", 126.0, false);

            Assert.That(a.CompareTo(b), Is.LessThan(0));
            Assert.That(b.CompareTo(a), Is.GreaterThan(0));
        }

        [Test]
        public void CompareTo_IsReflexive()
        {
            Assert.That(_sample1.CompareTo(_sample1), Is.EqualTo(0));
            Assert.That(_sample2.CompareTo(_sample2), Is.EqualTo(0));
        }

        #endregion

        #region CompareTo Tests - Case Sensitivity

        [Test]
        public void CompareTo_CaseSensitiveCondition()
        {
            // Different identity to avoid Equals short-circuit
            var upper = new IsobaricQuantSampleInfo(@"C:\a.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            var lower = new IsobaricQuantSampleInfo(@"C:\b.raw", "control", 1, 1, 0, 1, "127N", 126.0, false);

            // Uppercase 'C' comes before lowercase 'c' in ordinal comparison
            Assert.That(upper.CompareTo(lower), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_CaseSensitiveFilePath()
        {
            // Different ChannelLabel to avoid Equals short-circuit
            var upper = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            var lower = new IsobaricQuantSampleInfo(@"c:\test.raw", "Control", 1, 1, 0, 1, "127N", 126.0, false);

            // Uppercase 'C' comes before lowercase 'c'
            Assert.That(upper.CompareTo(lower), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_EmptyCondition_ComesBeforeNonEmpty()
        {
            var empty = new IsobaricQuantSampleInfo(@"C:\a.raw", string.Empty, 1, 1, 0, 1, "126", 126.0, false);
            var alpha = new IsobaricQuantSampleInfo(@"C:\b.raw", "Alpha", 1, 1, 0, 1, "127N", 126.0, false);

            Assert.That(empty.CompareTo(alpha), Is.LessThan(0));
        }

        #endregion

        #region Interface Implementation Tests

        [Test]
        public void ImplementsISampleInfo()
        {
            Assert.That(_sample1, Is.InstanceOf<ISampleInfo>());
        }

        [Test]
        public void ImplementsIEquatableOfIsobaricQuantSampleInfo()
        {
            Assert.That(_sample1, Is.InstanceOf<IEquatable<IsobaricQuantSampleInfo>>());
        }

        [Test]
        public void ImplementsIComparableOfISampleInfo()
        {
            Assert.That(_sample1, Is.InstanceOf<IComparable<ISampleInfo>>());
        }

        [Test]
        public void ISampleInfo_PropertiesAccessibleViaInterface()
        {
            ISampleInfo sample = _sample1;

            Assert.That(sample.FullFilePathWithExtension, Is.EqualTo(@"C:\Data\sample1.raw"));
            Assert.That(sample.Condition, Is.EqualTo("Control"));
            Assert.That(sample.BiologicalReplicate, Is.EqualTo(1));
            Assert.That(sample.TechnicalReplicate, Is.EqualTo(1));
            Assert.That(sample.Fraction, Is.EqualTo(0));
        }

        #endregion

        #region Cross-Type Tests

        [Test]
        public void CrossType_NotEqualTo_SpectraFileInfo()
        {
            var isobaric = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            var spectra = new SpectraFileInfo(@"C:\test.raw", "Control", 1, 1, 0);

            Assert.That(isobaric.Equals(spectra), Is.False);
            Assert.That(isobaric.Equals((ISampleInfo)spectra), Is.False);
            Assert.That(isobaric.Equals((object)spectra), Is.False);
        }

        [Test]
        public void CrossType_InHashSet_TreatedAsDifferent()
        {
            var isobaric = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            var spectra = new SpectraFileInfo(@"C:\test.raw", "Control", 1, 1, 0);

            var set = new HashSet<ISampleInfo> { isobaric, spectra };
            Assert.That(set.Count, Is.EqualTo(2));
        }

        #endregion

        #region Edge Case Tests

        [Test]
        public void EdgeCase_WhitespaceChannelLabel()
        {
            var sample = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "   ", 126.0, false);
            Assert.That(sample.ChannelLabel, Is.EqualTo("   "));
            Assert.That(sample.ToString(), Is.EqualTo("test_   "));
        }

        [Test]
        public void EdgeCase_SpecialCharactersInCondition()
        {
            var sample = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control-Group_1/A", 1, 1, 0, 1, "126", 126.0, false);
            Assert.That(sample.Condition, Is.EqualTo("Control-Group_1/A"));
        }

        [Test]
        public void EdgeCase_UnicodeInChannelLabel()
        {
            var sample = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "αβγ", 126.0, false);
            Assert.That(sample.ChannelLabel, Is.EqualTo("αβγ"));
            Assert.That(sample.ToString(), Is.EqualTo("test_αβγ"));
        }

        [Test]
        public void EdgeCase_DoubleNaN_ReporterMz()
        {
            var sample = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "126", double.NaN, false);
            Assert.That(double.IsNaN(sample.ReporterIonMz), Is.True);
        }

        [Test]
        public void EdgeCase_DoubleInfinity_ReporterMz()
        {
            var posInf = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "126", double.PositiveInfinity, false);
            var negInf = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "127N", double.NegativeInfinity, false);

            Assert.That(double.IsPositiveInfinity(posInf.ReporterIonMz), Is.True);
            Assert.That(double.IsNegativeInfinity(negInf.ReporterIonMz), Is.True);
        }

        [Test]
        public void EdgeCase_UnixStylePath()
        {
            var sample = new IsobaricQuantSampleInfo("/home/user/data/sample.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            Assert.That(sample.ToString(), Is.EqualTo("sample_126"));
        }

        [Test]
        public void EdgeCase_NetworkPath()
        {
            var sample = new IsobaricQuantSampleInfo(@"\\server\share\data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            Assert.That(sample.ToString(), Is.EqualTo("sample_126"));
        }

        [Test]
        public void EdgeCase_MultipleIdenticalInList_DistinctRemovesDuplicates()
        {
            var list = new List<IsobaricQuantSampleInfo> { _sample1, _sampleEqualToSample1, _sample1 };
            var distinctCount = list.Distinct().Count();
            Assert.That(distinctCount, Is.EqualTo(1));
        }

        [Test]
        public void EdgeCase_TMTChannelLabels()
        {
            var channels = new[] { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C" };
            var samples = channels.Select((ch, i) =>
                new IsobaricQuantSampleInfo($@"C:\test{i}.raw", "Control", 1, 1, 0, 1, ch, 126.0 + i, false)).ToList();

            Assert.That(samples.Count, Is.EqualTo(11));
            Assert.That(samples.Distinct().Count(), Is.EqualTo(11)); // All unique
        }

        #endregion
    }
}