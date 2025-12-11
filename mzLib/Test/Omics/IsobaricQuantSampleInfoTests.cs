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
        private IsobaricQuantSampleInfo _sampleIdenticalToSample1 = null!;
        private IsobaricQuantSampleInfo _sampleDifferentPlexId = null!;
        private IsobaricQuantSampleInfo _sampleDifferentChannelLabel = null!;
        private IsobaricQuantSampleInfo _sampleDifferentCondition = null!;
        private IsobaricQuantSampleInfo _sampleDifferentBioRep = null!;
        private IsobaricQuantSampleInfo _sampleDifferentTechRep = null!;
        private IsobaricQuantSampleInfo _sampleDifferentFraction = null!;
        private IsobaricQuantSampleInfo _sampleDifferentFilePath = null!;
        private IsobaricQuantSampleInfo _sampleDifferentReporterMz = null!;
        private IsobaricQuantSampleInfo _sampleReferenceChannel = null!;

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
                isReferenceChannel: false);

            // Same FullFilePathWithExtension and ChannelLabel as _sample1 (should be equal)
            _sampleIdenticalToSample1 = new IsobaricQuantSampleInfo(
                fullFilePathWithExtension: @"C:\Data\sample1.raw",
                condition: "Different",
                biologicalReplicate: 99,
                technicalReplicate: 99,
                fraction: 99,
                plexId: 99,
                channelLabel: "126",
                reporterIonMz: 999.999,
                isReferenceChannel: true);

            // Same FilePath and ChannelLabel as _sample1, different PlexId (should be equal)
            _sampleDifferentPlexId = new IsobaricQuantSampleInfo(
                fullFilePathWithExtension: @"C:\Data\sample1.raw",
                condition: "Control",
                biologicalReplicate: 1,
                technicalReplicate: 1,
                fraction: 0,
                plexId: 2,
                channelLabel: "126",
                reporterIonMz: 126.127726,
                isReferenceChannel: false);

            _sampleDifferentChannelLabel = new IsobaricQuantSampleInfo(
                fullFilePathWithExtension: @"C:\Data\sample1.raw",
                condition: "Control",
                biologicalReplicate: 1,
                technicalReplicate: 1,
                fraction: 0,
                plexId: 1,
                channelLabel: "127N",
                reporterIonMz: 126.127726,
                isReferenceChannel: false);

            _sampleDifferentCondition = new IsobaricQuantSampleInfo(
                fullFilePathWithExtension: @"C:\Data\sample1.raw",
                condition: "Treatment",
                biologicalReplicate: 1,
                technicalReplicate: 1,
                fraction: 0,
                plexId: 1,
                channelLabel: "126",
                reporterIonMz: 126.127726,
                isReferenceChannel: false);

            _sampleDifferentBioRep = new IsobaricQuantSampleInfo(
                fullFilePathWithExtension: @"C:\Data\sample1.raw",
                condition: "Control",
                biologicalReplicate: 2,
                technicalReplicate: 1,
                fraction: 0,
                plexId: 1,
                channelLabel: "126",
                reporterIonMz: 126.127726,
                isReferenceChannel: false);

            _sampleDifferentTechRep = new IsobaricQuantSampleInfo(
                fullFilePathWithExtension: @"C:\Data\sample1.raw",
                condition: "Control",
                biologicalReplicate: 1,
                technicalReplicate: 2,
                fraction: 0,
                plexId: 1,
                channelLabel: "126",
                reporterIonMz: 126.127726,
                isReferenceChannel: false);

            _sampleDifferentFraction = new IsobaricQuantSampleInfo(
                fullFilePathWithExtension: @"C:\Data\sample1.raw",
                condition: "Control",
                biologicalReplicate: 1,
                technicalReplicate: 1,
                fraction: 1,
                plexId: 1,
                channelLabel: "126",
                reporterIonMz: 126.127726,
                isReferenceChannel: false);

            _sampleDifferentFilePath = new IsobaricQuantSampleInfo(
                fullFilePathWithExtension: @"C:\Data\other.raw",
                condition: "Control",
                biologicalReplicate: 1,
                technicalReplicate: 1,
                fraction: 0,
                plexId: 1,
                channelLabel: "126",
                reporterIonMz: 126.127726,
                isReferenceChannel: false);

            _sampleDifferentReporterMz = new IsobaricQuantSampleInfo(
                fullFilePathWithExtension: @"C:\Data\sample1.raw",
                condition: "Control",
                biologicalReplicate: 1,
                technicalReplicate: 1,
                fraction: 0,
                plexId: 1,
                channelLabel: "126",
                reporterIonMz: 999.999,
                isReferenceChannel: false);

            _sampleReferenceChannel = new IsobaricQuantSampleInfo(
                fullFilePathWithExtension: @"C:\Data\sample1.raw",
                condition: "Control",
                biologicalReplicate: 1,
                technicalReplicate: 1,
                fraction: 0,
                plexId: 1,
                channelLabel: "126",
                reporterIonMz: 126.127726,
                isReferenceChannel: true);
        }

        #endregion

        #region Constructor Tests

        [Test]
        public void Constructor_SampleIdentifier_IsComputedFromFilePathAndChannelLabel()
        {
            var sample1 = new IsobaricQuantSampleInfo(
                @"C:\a.raw", "A", 1, 1, 0, 1, "126", 126.0, false);
            var sample2 = new IsobaricQuantSampleInfo(
                @"C:\a.raw", "B", 2, 2, 1, 99, "126", 127.0, true);

            // Same FullFilePathWithExtension and ChannelLabel should produce same UniqueIdentifier
            Assert.That(sample1.UniqueIdentifier, Is.EqualTo(sample2.UniqueIdentifier));
        }

        [Test]
        public void Constructor_SampleIdentifier_DifferentFilePath_ProducesDifferentValue()
        {
            var sample1 = new IsobaricQuantSampleInfo(
                @"C:\a.raw", "A", 1, 1, 0, 1, "126", 126.0, false);
            var sample2 = new IsobaricQuantSampleInfo(
                @"C:\b.raw", "A", 1, 1, 0, 1, "126", 126.0, false);

            Assert.That(sample1.UniqueIdentifier, Is.Not.EqualTo(sample2.UniqueIdentifier));
        }

        [Test]
        public void Constructor_SampleIdentifier_DifferentChannelLabel_ProducesDifferentValue()
        {
            var sample1 = new IsobaricQuantSampleInfo(
                @"C:\a.raw", "A", 1, 1, 0, 1, "126", 126.0, false);
            var sample2 = new IsobaricQuantSampleInfo(
                @"C:\a.raw", "A", 1, 1, 0, 1, "127N", 126.0, false);

            Assert.That(sample1.UniqueIdentifier, Is.Not.EqualTo(sample2.UniqueIdentifier));
        }

        [Test]
        public void Constructor_SampleIdentifier_DifferentPlexId_SameFilePathAndChannel_ProducesSameValue()
        {
            // PlexId is no longer part of identity - FilePath is
            var sample1 = new IsobaricQuantSampleInfo(
                @"C:\a.raw", "A", 1, 1, 0, 1, "126", 126.0, false);
            var sample2 = new IsobaricQuantSampleInfo(
                @"C:\a.raw", "A", 1, 1, 0, 99, "126", 126.0, false);

            Assert.That(sample1.UniqueIdentifier, Is.EqualTo(sample2.UniqueIdentifier));
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
            // _sampleIdenticalToSample1 has same FullFilePathWithExtension and ChannelLabel but different everything else
            Assert.That(_sample1.Equals(_sampleIdenticalToSample1), Is.True);
        }

        [Test]
        public void Equals_Typed_Null_ReturnsFalse()
        {
            Assert.That(_sample1.Equals((IsobaricQuantSampleInfo?)null), Is.False);
        }

        [Test]
        public void Equals_Typed_DifferentFilePath_ReturnsFalse()
        {
            Assert.That(_sample1.Equals(_sampleDifferentFilePath), Is.False);
        }

        [Test]
        public void Equals_Typed_DifferentChannelLabel_ReturnsFalse()
        {
            Assert.That(_sample1.Equals(_sampleDifferentChannelLabel), Is.False);
        }

        [Test]
        public void Equals_Typed_DifferentPlexId_SameFilePathAndChannel_ReturnsTrue()
        {
            // PlexId is NOT part of equality - only FullFilePathWithExtension and ChannelLabel
            Assert.That(_sample1.Equals(_sampleDifferentPlexId), Is.True);
        }

        [Test]
        public void Equals_Typed_DifferentCondition_SameFilePathAndChannel_ReturnsTrue()
        {
            Assert.That(_sample1.Equals(_sampleDifferentCondition), Is.True);
        }

        [Test]
        public void Equals_Typed_DifferentBioRep_SameFilePathAndChannel_ReturnsTrue()
        {
            Assert.That(_sample1.Equals(_sampleDifferentBioRep), Is.True);
        }

        [Test]
        public void Equals_Typed_DifferentTechRep_SameFilePathAndChannel_ReturnsTrue()
        {
            Assert.That(_sample1.Equals(_sampleDifferentTechRep), Is.True);
        }

        [Test]
        public void Equals_Typed_DifferentFraction_SameFilePathAndChannel_ReturnsTrue()
        {
            Assert.That(_sample1.Equals(_sampleDifferentFraction), Is.True);
        }

        [Test]
        public void Equals_Typed_DifferentReporterMz_SameFilePathAndChannel_ReturnsTrue()
        {
            Assert.That(_sample1.Equals(_sampleDifferentReporterMz), Is.True);
        }

        [Test]
        public void Equals_Typed_DifferentIsReferenceChannel_SameFilePathAndChannel_ReturnsTrue()
        {
            Assert.That(_sample1.Equals(_sampleReferenceChannel), Is.True);
        }

        [Test]
        public void Equals_Typed_CaseSensitiveChannelLabel()
        {
            var sampleLower = new IsobaricQuantSampleInfo(
                @"C:\test.raw", "Control", 1, 1, 0, 1, "127n", 127.0, false);
            var sampleUpper = new IsobaricQuantSampleInfo(
                @"C:\test.raw", "Control", 1, 1, 0, 1, "127N", 127.0, false);

            Assert.That(sampleLower.Equals(sampleUpper), Is.False);
        }

        [Test]
        public void Equals_Typed_CaseSensitiveFilePath()
        {
            var sampleLower = new IsobaricQuantSampleInfo(
                @"c:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            var sampleUpper = new IsobaricQuantSampleInfo(
                @"C:\Test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);

            Assert.That(sampleLower.Equals(sampleUpper), Is.False);
        }

        #endregion

        #region Equals(ISampleInfo) Tests

        [Test]
        public void Equals_ISampleInfo_SameIsobaricInstance_ReturnsTrue()
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
        public void Equals_ISampleInfo_DifferentIsobaricSample_ReturnsFalse()
        {
            ISampleInfo other = _sample2;
            Assert.That(_sample1.Equals(other), Is.False);
        }

        [Test]
        public void Equals_ISampleInfo_SpectraFileInfo_ReturnsFalse()
        {
            var spectraFileInfo = new SpectraFileInfo(@"C:\Data\sample1.raw", "Control", 1, 1, 0);
            Assert.That(_sample1.Equals(spectraFileInfo), Is.False);
        }

        #endregion

        #region Equals(object) Tests

        [Test]
        public void Equals_Object_SameReference_ReturnsTrue()
        {
            Assert.That(_sample1.Equals((object)_sample1), Is.True);
        }

        [Test]
        public void Equals_Object_IdenticalFilePathAndChannel_ReturnsTrue()
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
            Assert.That(_sample1.Equals("not an IsobaricQuantSampleInfo"), Is.False);
        }

        [Test]
        public void Equals_Object_IntegerType_ReturnsFalse()
        {
            Assert.That(_sample1.Equals(123), Is.False);
        }

        [Test]
        public void Equals_Object_SpectraFileInfo_ReturnsFalse()
        {
            var spectraFileInfo = new SpectraFileInfo(@"C:\Data\sample1.raw", "Control", 1, 1, 0);
            Assert.That(_sample1.Equals((object)spectraFileInfo), Is.False);
        }

        #endregion

        #region Equality Operator (==) Tests

        [Test]
        public void EqualityOperator_SameFilePathAndChannelLabel_ReturnsTrue()
        {
            Assert.That(_sample1 == _sampleIdenticalToSample1, Is.True);
        }

        [Test]
        public void EqualityOperator_DifferentValues_ReturnsFalse()
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
        public void EqualityOperator_LeftNull_RightNotNull_ReturnsFalse()
        {
            IsobaricQuantSampleInfo? left = null;
            Assert.That(left == _sample1, Is.False);
        }

        [Test]
        public void EqualityOperator_LeftNotNull_RightNull_ReturnsFalse()
        {
            IsobaricQuantSampleInfo? right = null;
            Assert.That(_sample1 == right, Is.False);
        }

        [Test]
        public void EqualityOperator_SameReference_ReturnsTrue()
        {
            var sample = _sample1;
            Assert.That(_sample1 == sample, Is.True);
        }

        [Test]
        public void EqualityOperator_DifferentFilePath_ReturnsFalse()
        {
            Assert.That(_sample1 == _sampleDifferentFilePath, Is.False);
        }

        [Test]
        public void EqualityOperator_DifferentChannelLabel_ReturnsFalse()
        {
            Assert.That(_sample1 == _sampleDifferentChannelLabel, Is.False);
        }

        [Test]
        public void EqualityOperator_DifferentPlexId_SameFilePathAndChannel_ReturnsTrue()
        {
            // PlexId is NOT part of equality
            Assert.That(_sample1 == _sampleDifferentPlexId, Is.True);
        }

        #endregion

        #region Inequality Operator (!=) Tests

        [Test]
        public void InequalityOperator_DifferentValues_ReturnsTrue()
        {
            Assert.That(_sample1 != _sample2, Is.True);
        }

        [Test]
        public void InequalityOperator_SameFilePathAndChannelLabel_ReturnsFalse()
        {
            Assert.That(_sample1 != _sampleIdenticalToSample1, Is.False);
        }

        [Test]
        public void InequalityOperator_BothNull_ReturnsFalse()
        {
            IsobaricQuantSampleInfo? left = null;
            IsobaricQuantSampleInfo? right = null;
            Assert.That(left != right, Is.False);
        }

        [Test]
        public void InequalityOperator_LeftNull_ReturnsTrue()
        {
            IsobaricQuantSampleInfo? left = null;
            Assert.That(left != _sample1, Is.True);
        }

        [Test]
        public void InequalityOperator_RightNull_ReturnsTrue()
        {
            IsobaricQuantSampleInfo? right = null;
            Assert.That(_sample1 != right, Is.True);
        }

        [Test]
        public void InequalityOperator_DifferentFilePath_ReturnsTrue()
        {
            Assert.That(_sample1 != _sampleDifferentFilePath, Is.True);
        }

        [Test]
        public void InequalityOperator_DifferentChannelLabel_ReturnsTrue()
        {
            Assert.That(_sample1 != _sampleDifferentChannelLabel, Is.True);
        }

        [Test]
        public void InequalityOperator_DifferentPlexId_SameFilePathAndChannel_ReturnsFalse()
        {
            // PlexId is NOT part of equality
            Assert.That(_sample1 != _sampleDifferentPlexId, Is.False);
        }

        #endregion

        #region GetHashCode Tests

        [Test]
        public void GetHashCode_SameFilePathAndChannelLabel_ReturnsSameHashCode()
        {
            Assert.That(_sample1.GetHashCode(), Is.EqualTo(_sampleIdenticalToSample1.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentFilePath_ReturnsDifferentHashCode()
        {
            Assert.That(_sample1.GetHashCode(), Is.Not.EqualTo(_sampleDifferentFilePath.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentChannelLabel_ReturnsDifferentHashCode()
        {
            Assert.That(_sample1.GetHashCode(), Is.Not.EqualTo(_sampleDifferentChannelLabel.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentPlexId_SameFilePathAndChannel_ReturnsSameHashCode()
        {
            // PlexId is NOT part of hash code
            Assert.That(_sample1.GetHashCode(), Is.EqualTo(_sampleDifferentPlexId.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentCondition_SameFilePathAndChannel_ReturnsSameHashCode()
        {
            Assert.That(_sample1.GetHashCode(), Is.EqualTo(_sampleDifferentCondition.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentBioRep_SameFilePathAndChannel_ReturnsSameHashCode()
        {
            Assert.That(_sample1.GetHashCode(), Is.EqualTo(_sampleDifferentBioRep.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentTechRep_SameFilePathAndChannel_ReturnsSameHashCode()
        {
            Assert.That(_sample1.GetHashCode(), Is.EqualTo(_sampleDifferentTechRep.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentFraction_SameFilePathAndChannel_ReturnsSameHashCode()
        {
            Assert.That(_sample1.GetHashCode(), Is.EqualTo(_sampleDifferentFraction.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentReporterMz_SameFilePathAndChannel_ReturnsSameHashCode()
        {
            Assert.That(_sample1.GetHashCode(), Is.EqualTo(_sampleDifferentReporterMz.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentIsReferenceChannel_SameFilePathAndChannel_ReturnsSameHashCode()
        {
            Assert.That(_sample1.GetHashCode(), Is.EqualTo(_sampleReferenceChannel.GetHashCode()));
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
        public void GetHashCode_MatchesSampleIdentifier()
        {
            // UniqueIdentifier is computed the same way as GetHashCode
            Assert.That(_sample1.GetHashCode(), Is.EqualTo(_sample1.UniqueIdentifier));
        }

        [Test]
        public void GetHashCode_WorksInHashSet()
        {
            var set = new HashSet<IsobaricQuantSampleInfo>
            {
                _sample1,
                _sampleIdenticalToSample1, // Should not be added (duplicate)
                _sample2
            };

            Assert.That(set.Count, Is.EqualTo(2));
            Assert.That(set.Contains(_sample1), Is.True);
            Assert.That(set.Contains(_sample2), Is.True);
            Assert.That(set.Contains(_sampleIdenticalToSample1), Is.True); // Same as _sample1
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
            Assert.That(dict[_sampleIdenticalToSample1], Is.EqualTo("First")); // Same key as _sample1
            Assert.That(dict[_sample2], Is.EqualTo("Second"));
        }

        [Test]
        public void GetHashCode_EmptyChannelLabel_ProducesValidHashCode()
        {
            var sample = new IsobaricQuantSampleInfo(
                @"C:\test.raw", "Control", 1, 1, 0, 1, string.Empty, 126.0, false);
            var hash = sample.GetHashCode();
            Assert.That(hash, Is.TypeOf<int>());
        }

        [Test]
        public void GetHashCode_NegativePlexId_ProducesValidHashCode()
        {
            var sample = new IsobaricQuantSampleInfo(
                @"C:\test.raw", "Control", 1, 1, 0, -1, "126", 126.0, false);
            var hash = sample.GetHashCode();
            Assert.That(hash, Is.TypeOf<int>());
        }

        [Test]
        public void GetHashCode_MaxPlexId_ProducesValidHashCode()
        {
            var sample = new IsobaricQuantSampleInfo(
                @"C:\test.raw", "Control", 1, 1, 0, int.MaxValue, "126", 126.0, false);
            var hash = sample.GetHashCode();
            Assert.That(hash, Is.TypeOf<int>());
        }

        #endregion

        #region ToString Tests

        [Test]
        public void ToString_ReturnsExpectedFormat()
        {
            Assert.That(_sample1.ToString(), Is.EqualTo("sample1_126"));
        }

        [Test]
        public void ToString_WithDifferentFile_ReturnsCorrectFormat()
        {
            Assert.That(_sample2.ToString(), Is.EqualTo("sample2_127N"));
        }

        [Test]
        public void ToString_WithEmptyFilePath_ReturnsCorrectFormat()
        {
            var sample = new IsobaricQuantSampleInfo(
                string.Empty, "Control", 1, 1, 0, 1, "126", 126.0, false);
            Assert.That(sample.ToString(), Is.EqualTo("_126"));
        }

        [Test]
        public void ToString_WithEmptyChannelLabel_ReturnsCorrectFormat()
        {
            var sample = new IsobaricQuantSampleInfo(
                @"C:\test.raw", "Control", 1, 1, 0, 1, string.Empty, 126.0, false);
            Assert.That(sample.ToString(), Is.EqualTo("test_"));
        }

        [Test]
        public void ToString_WithLongChannelLabel_ReturnsCorrectFormat()
        {
            var sample = new IsobaricQuantSampleInfo(
                @"C:\test.raw", "Control", 1, 1, 0, 1, "VeryLongChannelLabel", 126.0, false);
            Assert.That(sample.ToString(), Is.EqualTo("test_VeryLongChannelLabel"));
        }

        [Test]
        public void ToString_IsNotNull()
        {
            Assert.That(_sample1.ToString(), Is.Not.Null);
        }

        [Test]
        public void ToString_IsNotEmpty()
        {
            Assert.That(_sample1.ToString(), Is.Not.Empty);
        }

        [Test]
        public void ToString_UsesFileNameWithoutExtension()
        {
            var sample = new IsobaricQuantSampleInfo(
                @"C:\Data\MyExperiment.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            Assert.That(sample.ToString(), Is.EqualTo("MyExperiment_126"));
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
        public void CompareTo_IdenticalFilePathAndChannelLabel_ReturnsZero()
        {
            Assert.That(_sample1.CompareTo(_sampleIdenticalToSample1), Is.EqualTo(0));
        }

        [Test]
        public void CompareTo_DifferentFilePath_AscendingOrder()
        {
            var fileA = new IsobaricQuantSampleInfo(@"C:\a.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            var fileB = new IsobaricQuantSampleInfo(@"C:\b.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            var fileZ = new IsobaricQuantSampleInfo(@"C:\z.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);

            Assert.That(fileA.CompareTo(fileB), Is.LessThan(0));
            Assert.That(fileB.CompareTo(fileA), Is.GreaterThan(0));
            Assert.That(fileA.CompareTo(fileZ), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_SameFilePath_DifferentChannelLabel_AscendingOrder()
        {
            var channel126 = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            var channel127N = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "127N", 127.0, false);
            var channel131C = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "131C", 131.0, false);

            Assert.That(channel126.CompareTo(channel127N), Is.LessThan(0));
            Assert.That(channel127N.CompareTo(channel126), Is.GreaterThan(0));
            Assert.That(channel126.CompareTo(channel131C), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_FilePathTakesPriorityOverChannelLabel()
        {
            var fileAChannel131 = new IsobaricQuantSampleInfo(@"C:\a.raw", "Control", 1, 1, 0, 1, "131C", 131.0, false);
            var fileBChannel126 = new IsobaricQuantSampleInfo(@"C:\b.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);

            // File A should come before File B, even though 131C > 126 alphabetically
            Assert.That(fileAChannel131.CompareTo(fileBChannel126), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_DifferentPlexId_SameFilePathAndChannel_ReturnsZero()
        {
            // PlexId is NOT part of comparison - only FilePath and ChannelLabel
            Assert.That(_sample1.CompareTo(_sampleDifferentPlexId), Is.EqualTo(0));
        }

        [Test]
        public void CompareTo_NonIsobaricSampleInfo_ReturnsPositive()
        {
            var spectraFileInfo = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            Assert.That(_sample1.CompareTo(spectraFileInfo), Is.GreaterThan(0));
        }

        [Test]
        public void CompareTo_EmptyChannelLabel_ComesBeforeNonEmpty()
        {
            var channelEmpty = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, string.Empty, 126.0, false);
            var channel126 = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);

            Assert.That(channelEmpty.CompareTo(channel126), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_CaseSensitiveChannelLabel()
        {
            var channelLower = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "127n", 127.0, false);
            var channelUpper = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "127N", 127.0, false);

            // Uppercase comes before lowercase in ordinal comparison
            Assert.That(channelUpper.CompareTo(channelLower), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_CaseSensitiveFilePath()
        {
            var filePathLower = new IsobaricQuantSampleInfo(@"c:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            var filePathUpper = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);

            // Uppercase comes before lowercase in ordinal comparison
            Assert.That(filePathUpper.CompareTo(filePathLower), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_Sorting_ProducesCorrectOrder()
        {
            var samples = new List<IsobaricQuantSampleInfo>
            {
                new(@"C:\b.raw", "Control", 1, 1, 0, 1, "126", 126.0, false),
                new(@"C:\a.raw", "Control", 1, 1, 0, 1, "131C", 131.0, false),
                new(@"C:\a.raw", "Control", 1, 1, 0, 1, "127N", 127.0, false),
                new(@"C:\a.raw", "Control", 1, 1, 0, 1, "126", 126.0, false),
                new(@"C:\b.raw", "Control", 1, 1, 0, 1, "127N", 127.0, false),
            };

            samples.Sort((a, b) => a.CompareTo(b));

            Assert.That(samples[0].FullFilePathWithExtension, Is.EqualTo(@"C:\a.raw"));
            Assert.That(samples[0].ChannelLabel, Is.EqualTo("126"));
            Assert.That(samples[1].FullFilePathWithExtension, Is.EqualTo(@"C:\a.raw"));
            Assert.That(samples[1].ChannelLabel, Is.EqualTo("127N"));
            Assert.That(samples[2].FullFilePathWithExtension, Is.EqualTo(@"C:\a.raw"));
            Assert.That(samples[2].ChannelLabel, Is.EqualTo("131C"));
            Assert.That(samples[3].FullFilePathWithExtension, Is.EqualTo(@"C:\b.raw"));
            Assert.That(samples[3].ChannelLabel, Is.EqualTo("126"));
            Assert.That(samples[4].FullFilePathWithExtension, Is.EqualTo(@"C:\b.raw"));
            Assert.That(samples[4].ChannelLabel, Is.EqualTo("127N"));
        }

        [Test]
        public void CompareTo_IsTransitive()
        {
            var a = new IsobaricQuantSampleInfo(@"C:\a.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            var b = new IsobaricQuantSampleInfo(@"C:\a.raw", "Control", 1, 1, 0, 1, "127N", 127.0, false);
            var c = new IsobaricQuantSampleInfo(@"C:\b.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);

            // a < b (same file, 126 < 127N)
            // b < c (file a < file b)
            // Therefore a < c
            Assert.That(a.CompareTo(b), Is.LessThan(0));
            Assert.That(b.CompareTo(c), Is.LessThan(0));
            Assert.That(a.CompareTo(c), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_IsAntiSymmetric()
        {
            var a = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            var b = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "127N", 127.0, false);

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
        }

        #endregion

        #region Edge Case Tests

        [Test]
        public void EdgeCase_WhitespaceChannelLabel()
        {
            var sample = new IsobaricQuantSampleInfo(
                @"C:\test.raw", "Control", 1, 1, 0, 1, "   ", 126.0, false);
            Assert.That(sample.ChannelLabel, Is.EqualTo("   "));
            Assert.That(sample.ToString(), Is.EqualTo("test_   "));
        }

        [Test]
        public void EdgeCase_VeryLargePlexId()
        {
            var sample = new IsobaricQuantSampleInfo(
                @"C:\test.raw", "Control", 1, 1, 0, int.MaxValue, "126", 126.0, false);
            Assert.That(sample.PlexId, Is.EqualTo(int.MaxValue));
            // PlexId is no longer part of ToString - it uses filename
            Assert.That(sample.ToString(), Is.EqualTo("test_126"));
        }

        [Test]
        public void EdgeCase_SpecialCharactersInCondition()
        {
            var sample = new IsobaricQuantSampleInfo(
                @"C:\test.raw", "Control-Group_1/A", 1, 1, 0, 1, "126", 126.0, false);
            Assert.That(sample.Condition, Is.EqualTo("Control-Group_1/A"));
        }

        [Test]
        public void EdgeCase_UnicodeInChannelLabel()
        {
            var sample = new IsobaricQuantSampleInfo(
                @"C:\test.raw", "Control", 1, 1, 0, 1, "αβγ", 126.0, false);
            Assert.That(sample.ChannelLabel, Is.EqualTo("αβγ"));
            Assert.That(sample.ToString(), Is.EqualTo("test_αβγ"));
        }

        [Test]
        public void EdgeCase_DoubleNaN_ReporterMz()
        {
            var sample = new IsobaricQuantSampleInfo(
                @"C:\test.raw", "Control", 1, 1, 0, 1, "126", double.NaN, false);
            Assert.That(double.IsNaN(sample.ReporterIonMz), Is.True);
        }

        [Test]
        public void EdgeCase_DoubleInfinity_ReporterMz()
        {
            var samplePos = new IsobaricQuantSampleInfo(
                @"C:\test.raw", "Control", 1, 1, 0, 1, "126", double.PositiveInfinity, false);
            var sampleNeg = new IsobaricQuantSampleInfo(
                @"C:\test.raw", "Control", 1, 1, 0, 1, "127N", double.NegativeInfinity, false);

            Assert.That(double.IsPositiveInfinity(samplePos.ReporterIonMz), Is.True);
            Assert.That(double.IsNegativeInfinity(sampleNeg.ReporterIonMz), Is.True);
        }

        [Test]
        public void EdgeCase_MultipleIdenticalSamplesInList()
        {
            var list = new List<IsobaricQuantSampleInfo>
            {
                _sample1,
                _sampleIdenticalToSample1,
                _sample1
            };

            var distinctCount = list.Distinct().Count();
            Assert.That(distinctCount, Is.EqualTo(1));
        }

        #endregion
    }
}