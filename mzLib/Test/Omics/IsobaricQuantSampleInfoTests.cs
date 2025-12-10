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

            // Same PlexId and ChannelLabel as _sample1 (should be equal)
            _sampleIdenticalToSample1 = new IsobaricQuantSampleInfo(
                fullFilePathWithExtension: @"C:\Data\different.raw",
                condition: "Different",
                biologicalReplicate: 99,
                technicalReplicate: 99,
                fraction: 99,
                plexId: 1,
                channelLabel: "126",
                reporterIonMz: 999.999,
                isReferenceChannel: true);

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
        public void Constructor_WithValidParameters_SetsAllProperties()
        {
            var sample = new IsobaricQuantSampleInfo(
                fullFilePathWithExtension: @"C:\Data\test.raw",
                condition: "Control",
                biologicalReplicate: 1,
                technicalReplicate: 2,
                fraction: 3,
                plexId: 4,
                channelLabel: "127N",
                reporterIonMz: 127.124761,
                isReferenceChannel: true);

            Assert.That(sample.FullFilePathWithExtension, Is.EqualTo(@"C:\Data\test.raw"));
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
        public void Constructor_WithNullChannelLabel_ThrowsArgumentNullException()
        {
            Assert.Throws<ArgumentNullException>(() => new IsobaricQuantSampleInfo(
                fullFilePathWithExtension: @"C:\Data\test.raw",
                condition: "Control",
                biologicalReplicate: 1,
                technicalReplicate: 1,
                fraction: 0,
                plexId: 1,
                channelLabel: null!,
                reporterIonMz: 126.127726,
                isReferenceChannel: false));
        }

        [Test]
        public void Constructor_WithEmptyChannelLabel_AcceptsValue()
        {
            var sample = new IsobaricQuantSampleInfo(
                fullFilePathWithExtension: @"C:\Data\test.raw",
                condition: "Control",
                biologicalReplicate: 1,
                technicalReplicate: 1,
                fraction: 0,
                plexId: 1,
                channelLabel: string.Empty,
                reporterIonMz: 126.127726,
                isReferenceChannel: false);

            Assert.That(sample.ChannelLabel, Is.EqualTo(string.Empty));
        }

        [Test]
        public void Constructor_WithNegativeValues_AcceptsValues()
        {
            var sample = new IsobaricQuantSampleInfo(
                fullFilePathWithExtension: @"C:\Data\test.raw",
                condition: "Control",
                biologicalReplicate: -1,
                technicalReplicate: -2,
                fraction: -3,
                plexId: -4,
                channelLabel: "126",
                reporterIonMz: -126.0,
                isReferenceChannel: false);

            Assert.That(sample.BiologicalReplicate, Is.EqualTo(-1));
            Assert.That(sample.TechnicalReplicate, Is.EqualTo(-2));
            Assert.That(sample.Fraction, Is.EqualTo(-3));
            Assert.That(sample.PlexId, Is.EqualTo(-4));
            Assert.That(sample.ReporterIonMz, Is.EqualTo(-126.0));
        }

        [Test]
        public void Constructor_WithZeroValues_AcceptsValues()
        {
            var sample = new IsobaricQuantSampleInfo(
                fullFilePathWithExtension: string.Empty,
                condition: string.Empty,
                biologicalReplicate: 0,
                technicalReplicate: 0,
                fraction: 0,
                plexId: 0,
                channelLabel: "126",
                reporterIonMz: 0.0,
                isReferenceChannel: false);

            Assert.That(sample.BiologicalReplicate, Is.EqualTo(0));
            Assert.That(sample.TechnicalReplicate, Is.EqualTo(0));
            Assert.That(sample.Fraction, Is.EqualTo(0));
            Assert.That(sample.PlexId, Is.EqualTo(0));
            Assert.That(sample.ReporterIonMz, Is.EqualTo(0.0));
        }

        [Test]
        public void Constructor_WithMaxIntValues_AcceptsValues()
        {
            var sample = new IsobaricQuantSampleInfo(
                fullFilePathWithExtension: @"C:\Data\test.raw",
                condition: "Control",
                biologicalReplicate: int.MaxValue,
                technicalReplicate: int.MaxValue,
                fraction: int.MaxValue,
                plexId: int.MaxValue,
                channelLabel: "126",
                reporterIonMz: double.MaxValue,
                isReferenceChannel: false);

            Assert.That(sample.BiologicalReplicate, Is.EqualTo(int.MaxValue));
            Assert.That(sample.TechnicalReplicate, Is.EqualTo(int.MaxValue));
            Assert.That(sample.Fraction, Is.EqualTo(int.MaxValue));
            Assert.That(sample.PlexId, Is.EqualTo(int.MaxValue));
            Assert.That(sample.ReporterIonMz, Is.EqualTo(double.MaxValue));
        }

        [Test]
        public void Constructor_SampleIdentifier_IsComputedFromPlexIdAndChannelLabel()
        {
            var sample1 = new IsobaricQuantSampleInfo(
                @"C:\a.raw", "A", 1, 1, 0, 1, "126", 126.0, false);
            var sample2 = new IsobaricQuantSampleInfo(
                @"C:\b.raw", "B", 2, 2, 1, 1, "126", 127.0, true);

            // Same PlexId and ChannelLabel should produce same SampleIdentifier
            Assert.That(sample1.SampleIdentifier, Is.EqualTo(sample2.SampleIdentifier));
        }

        [Test]
        public void Constructor_SampleIdentifier_DifferentPlexId_ProducesDifferentValue()
        {
            var sample1 = new IsobaricQuantSampleInfo(
                @"C:\a.raw", "A", 1, 1, 0, 1, "126", 126.0, false);
            var sample2 = new IsobaricQuantSampleInfo(
                @"C:\a.raw", "A", 1, 1, 0, 2, "126", 126.0, false);

            Assert.That(sample1.SampleIdentifier, Is.Not.EqualTo(sample2.SampleIdentifier));
        }

        [Test]
        public void Constructor_SampleIdentifier_DifferentChannelLabel_ProducesDifferentValue()
        {
            var sample1 = new IsobaricQuantSampleInfo(
                @"C:\a.raw", "A", 1, 1, 0, 1, "126", 126.0, false);
            var sample2 = new IsobaricQuantSampleInfo(
                @"C:\a.raw", "A", 1, 1, 0, 1, "127N", 126.0, false);

            Assert.That(sample1.SampleIdentifier, Is.Not.EqualTo(sample2.SampleIdentifier));
        }

        [Test]
        public void Constructor_WithTmtChannelLabels_AcceptsAllFormats()
        {
            var channels = new[] { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C" };
            foreach (var channel in channels)
            {
                var sample = new IsobaricQuantSampleInfo(
                    @"C:\test.raw", "Control", 1, 1, 0, 1, channel, 126.0, false);
                Assert.That(sample.ChannelLabel, Is.EqualTo(channel));
            }
        }

        [Test]
        public void Constructor_WithItraqChannelLabels_AcceptsAllFormats()
        {
            var channels = new[] { "113", "114", "115", "116", "117", "118", "119", "121" };
            foreach (var channel in channels)
            {
                var sample = new IsobaricQuantSampleInfo(
                    @"C:\test.raw", "Control", 1, 1, 0, 1, channel, 113.0, false);
                Assert.That(sample.ChannelLabel, Is.EqualTo(channel));
            }
        }

        #endregion

        #region Equals(IsobaricQuantSampleInfo) Tests

        [Test]
        public void Equals_Typed_SameReference_ReturnsTrue()
        {
            Assert.That(_sample1.Equals(_sample1), Is.True);
        }

        [Test]
        public void Equals_Typed_SamePlexIdAndChannelLabel_ReturnsTrue()
        {
            // _sampleIdenticalToSample1 has same PlexId and ChannelLabel but different everything else
            Assert.That(_sample1.Equals(_sampleIdenticalToSample1), Is.True);
        }

        [Test]
        public void Equals_Typed_Null_ReturnsFalse()
        {
            Assert.That(_sample1.Equals((IsobaricQuantSampleInfo?)null), Is.False);
        }

        [Test]
        public void Equals_Typed_DifferentPlexId_ReturnsFalse()
        {
            Assert.That(_sample1.Equals(_sampleDifferentPlexId), Is.False);
        }

        [Test]
        public void Equals_Typed_DifferentChannelLabel_ReturnsFalse()
        {
            Assert.That(_sample1.Equals(_sampleDifferentChannelLabel), Is.False);
        }

        [Test]
        public void Equals_Typed_DifferentCondition_SamePlexAndChannel_ReturnsTrue()
        {
            // Condition is NOT part of equality - only PlexId and ChannelLabel
            Assert.That(_sample1.Equals(_sampleDifferentCondition), Is.True);
        }

        [Test]
        public void Equals_Typed_DifferentBioRep_SamePlexAndChannel_ReturnsTrue()
        {
            Assert.That(_sample1.Equals(_sampleDifferentBioRep), Is.True);
        }

        [Test]
        public void Equals_Typed_DifferentTechRep_SamePlexAndChannel_ReturnsTrue()
        {
            Assert.That(_sample1.Equals(_sampleDifferentTechRep), Is.True);
        }

        [Test]
        public void Equals_Typed_DifferentFraction_SamePlexAndChannel_ReturnsTrue()
        {
            Assert.That(_sample1.Equals(_sampleDifferentFraction), Is.True);
        }

        [Test]
        public void Equals_Typed_DifferentFilePath_SamePlexAndChannel_ReturnsTrue()
        {
            Assert.That(_sample1.Equals(_sampleDifferentFilePath), Is.True);
        }

        [Test]
        public void Equals_Typed_DifferentReporterMz_SamePlexAndChannel_ReturnsTrue()
        {
            Assert.That(_sample1.Equals(_sampleDifferentReporterMz), Is.True);
        }

        [Test]
        public void Equals_Typed_DifferentIsReferenceChannel_SamePlexAndChannel_ReturnsTrue()
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
        public void Equals_Object_IdenticalPlexAndChannel_ReturnsTrue()
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
        public void EqualityOperator_SamePlexIdAndChannelLabel_ReturnsTrue()
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
        public void EqualityOperator_DifferentPlexId_ReturnsFalse()
        {
            Assert.That(_sample1 == _sampleDifferentPlexId, Is.False);
        }

        [Test]
        public void EqualityOperator_DifferentChannelLabel_ReturnsFalse()
        {
            Assert.That(_sample1 == _sampleDifferentChannelLabel, Is.False);
        }

        #endregion

        #region Inequality Operator (!=) Tests

        [Test]
        public void InequalityOperator_DifferentValues_ReturnsTrue()
        {
            Assert.That(_sample1 != _sample2, Is.True);
        }

        [Test]
        public void InequalityOperator_SamePlexIdAndChannelLabel_ReturnsFalse()
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
        public void InequalityOperator_DifferentPlexId_ReturnsTrue()
        {
            Assert.That(_sample1 != _sampleDifferentPlexId, Is.True);
        }

        [Test]
        public void InequalityOperator_DifferentChannelLabel_ReturnsTrue()
        {
            Assert.That(_sample1 != _sampleDifferentChannelLabel, Is.True);
        }

        #endregion

        #region GetHashCode Tests

        [Test]
        public void GetHashCode_SamePlexIdAndChannelLabel_ReturnsSameHashCode()
        {
            Assert.That(_sample1.GetHashCode(), Is.EqualTo(_sampleIdenticalToSample1.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentPlexId_ReturnsDifferentHashCode()
        {
            Assert.That(_sample1.GetHashCode(), Is.Not.EqualTo(_sampleDifferentPlexId.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentChannelLabel_ReturnsDifferentHashCode()
        {
            Assert.That(_sample1.GetHashCode(), Is.Not.EqualTo(_sampleDifferentChannelLabel.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentCondition_SamePlexAndChannel_ReturnsSameHashCode()
        {
            Assert.That(_sample1.GetHashCode(), Is.EqualTo(_sampleDifferentCondition.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentBioRep_SamePlexAndChannel_ReturnsSameHashCode()
        {
            Assert.That(_sample1.GetHashCode(), Is.EqualTo(_sampleDifferentBioRep.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentTechRep_SamePlexAndChannel_ReturnsSameHashCode()
        {
            Assert.That(_sample1.GetHashCode(), Is.EqualTo(_sampleDifferentTechRep.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentFraction_SamePlexAndChannel_ReturnsSameHashCode()
        {
            Assert.That(_sample1.GetHashCode(), Is.EqualTo(_sampleDifferentFraction.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentFilePath_SamePlexAndChannel_ReturnsSameHashCode()
        {
            Assert.That(_sample1.GetHashCode(), Is.EqualTo(_sampleDifferentFilePath.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentReporterMz_SamePlexAndChannel_ReturnsSameHashCode()
        {
            Assert.That(_sample1.GetHashCode(), Is.EqualTo(_sampleDifferentReporterMz.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentIsReferenceChannel_SamePlexAndChannel_ReturnsSameHashCode()
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
            // SampleIdentifier is computed the same way as GetHashCode
            Assert.That(_sample1.GetHashCode(), Is.EqualTo(_sample1.SampleIdentifier));
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
            Assert.That(_sample1.ToString(), Is.EqualTo("Plex1_126"));
        }

        [Test]
        public void ToString_WithDifferentPlexId_ReturnsCorrectFormat()
        {
            Assert.That(_sample2.ToString(), Is.EqualTo("Plex2_127N"));
        }

        [Test]
        public void ToString_WithZeroPlexId_ReturnsCorrectFormat()
        {
            var sample = new IsobaricQuantSampleInfo(
                @"C:\test.raw", "Control", 1, 1, 0, 0, "126", 126.0, false);
            Assert.That(sample.ToString(), Is.EqualTo("Plex0_126"));
        }

        [Test]
        public void ToString_WithNegativePlexId_ReturnsCorrectFormat()
        {
            var sample = new IsobaricQuantSampleInfo(
                @"C:\test.raw", "Control", 1, 1, 0, -1, "126", 126.0, false);
            Assert.That(sample.ToString(), Is.EqualTo("Plex-1_126"));
        }

        [Test]
        public void ToString_WithEmptyChannelLabel_ReturnsCorrectFormat()
        {
            var sample = new IsobaricQuantSampleInfo(
                @"C:\test.raw", "Control", 1, 1, 0, 1, string.Empty, 126.0, false);
            Assert.That(sample.ToString(), Is.EqualTo("Plex1_"));
        }

        [Test]
        public void ToString_WithLongChannelLabel_ReturnsCorrectFormat()
        {
            var sample = new IsobaricQuantSampleInfo(
                @"C:\test.raw", "Control", 1, 1, 0, 1, "VeryLongChannelLabel", 126.0, false);
            Assert.That(sample.ToString(), Is.EqualTo("Plex1_VeryLongChannelLabel"));
        }

        [Test]
        public void ToString_WithSpecialCharactersInChannelLabel_ReturnsCorrectFormat()
        {
            var sample = new IsobaricQuantSampleInfo(
                @"C:\test.raw", "Control", 1, 1, 0, 1, "127-N", 127.0, false);
            Assert.That(sample.ToString(), Is.EqualTo("Plex1_127-N"));
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
        public void ToString_LargePlexId_ReturnsCorrectFormat()
        {
            var sample = new IsobaricQuantSampleInfo(
                @"C:\test.raw", "Control", 1, 1, 0, 999, "126", 126.0, false);
            Assert.That(sample.ToString(), Is.EqualTo("Plex999_126"));
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
        public void CompareTo_IdenticalPlexIdAndChannelLabel_ReturnsZero()
        {
            Assert.That(_sample1.CompareTo(_sampleIdenticalToSample1), Is.EqualTo(0));
        }

        [Test]
        public void CompareTo_DifferentPlexId_AscendingOrder()
        {
            var plex1 = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            var plex2 = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 2, "126", 126.0, false);
            var plex10 = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 10, "126", 126.0, false);

            Assert.That(plex1.CompareTo(plex2), Is.LessThan(0));
            Assert.That(plex2.CompareTo(plex1), Is.GreaterThan(0));
            Assert.That(plex1.CompareTo(plex10), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_SamePlexId_DifferentChannelLabel_AscendingOrder()
        {
            var channel126 = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            var channel127N = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "127N", 127.0, false);
            var channel131C = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "131C", 131.0, false);

            Assert.That(channel126.CompareTo(channel127N), Is.LessThan(0));
            Assert.That(channel127N.CompareTo(channel126), Is.GreaterThan(0));
            Assert.That(channel126.CompareTo(channel131C), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_PlexIdTakesPriorityOverChannelLabel()
        {
            var plex1Channel131 = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "131C", 131.0, false);
            var plex2Channel126 = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 2, "126", 126.0, false);

            // Plex 1 should come before Plex 2, even though 131C > 126 alphabetically
            Assert.That(plex1Channel131.CompareTo(plex2Channel126), Is.LessThan(0));
        }

        [Test]
        public void CompareTo_NonIsobaricSampleInfo_ReturnsPositive()
        {
            var spectraFileInfo = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            Assert.That(_sample1.CompareTo(spectraFileInfo), Is.GreaterThan(0));
        }

        [Test]
        public void CompareTo_NegativePlexId_OrderedCorrectly()
        {
            var plexNeg = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, -1, "126", 126.0, false);
            var plexZero = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 0, "126", 126.0, false);
            var plexPos = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);

            Assert.That(plexNeg.CompareTo(plexZero), Is.LessThan(0));
            Assert.That(plexZero.CompareTo(plexPos), Is.LessThan(0));
            Assert.That(plexNeg.CompareTo(plexPos), Is.LessThan(0));
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
        public void CompareTo_Sorting_ProducesCorrectOrder()
        {
            var samples = new List<IsobaricQuantSampleInfo>
            {
                new(@"C:\test.raw", "Control", 1, 1, 0, 2, "126", 126.0, false),
                new(@"C:\test.raw", "Control", 1, 1, 0, 1, "131C", 131.0, false),
                new(@"C:\test.raw", "Control", 1, 1, 0, 1, "127N", 127.0, false),
                new(@"C:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false),
                new(@"C:\test.raw", "Control", 1, 1, 0, 2, "127N", 127.0, false),
            };

            samples.Sort((a, b) => a.CompareTo(b));

            Assert.That(samples[0].PlexId, Is.EqualTo(1));
            Assert.That(samples[0].ChannelLabel, Is.EqualTo("126"));
            Assert.That(samples[1].PlexId, Is.EqualTo(1));
            Assert.That(samples[1].ChannelLabel, Is.EqualTo("127N"));
            Assert.That(samples[2].PlexId, Is.EqualTo(1));
            Assert.That(samples[2].ChannelLabel, Is.EqualTo("131C"));
            Assert.That(samples[3].PlexId, Is.EqualTo(2));
            Assert.That(samples[3].ChannelLabel, Is.EqualTo("126"));
            Assert.That(samples[4].PlexId, Is.EqualTo(2));
            Assert.That(samples[4].ChannelLabel, Is.EqualTo("127N"));
        }

        [Test]
        public void CompareTo_IsTransitive()
        {
            var a = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            var b = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, "127N", 127.0, false);
            var c = new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 2, "126", 126.0, false);

            // a < b (same plex, 126 < 127N)
            // b < c (plex 1 < plex 2)
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
        public void EdgeCase_WhitespaceChannelLabel()
        {
            var sample = new IsobaricQuantSampleInfo(
                @"C:\test.raw", "Control", 1, 1, 0, 1, "   ", 126.0, false);
            Assert.That(sample.ChannelLabel, Is.EqualTo("   "));
            Assert.That(sample.ToString(), Is.EqualTo("Plex1_   "));
        }

        [Test]
        public void EdgeCase_VeryLargePlexId()
        {
            var sample = new IsobaricQuantSampleInfo(
                @"C:\test.raw", "Control", 1, 1, 0, int.MaxValue, "126", 126.0, false);
            Assert.That(sample.PlexId, Is.EqualTo(int.MaxValue));
            Assert.That(sample.ToString(), Does.Contain(int.MaxValue.ToString()));
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
            Assert.That(sample.ToString(), Is.EqualTo("Plex1_αβγ"));
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