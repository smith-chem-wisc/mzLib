using NUnit.Framework;
using Omics.SequenceConversion;

namespace Test.Omics.SequenceConversion
{
    /// <summary>
    /// Tests for serializing canonical sequences into various format strings.
    /// </summary>
    [TestFixture]
    public class SerializationTests
    {
        private MzLibSequenceSerializer _mzLibSerializer;
        private MassShiftSequenceSerializer _massShiftSerializer;
        private ChronologerSequenceSerializer _chronologerSerializer;
        private MzLibSequenceParser _mzLibParser;
        private MassShiftSequenceParser _massShiftParser;

        [SetUp]
        public void Setup()
        {
            _mzLibSerializer = new MzLibSequenceSerializer();
            _massShiftSerializer = new MassShiftSequenceSerializer(new(4));
            _chronologerSerializer = new ChronologerSequenceSerializer();
            _mzLibParser = new MzLibSequenceParser();
            _massShiftParser = new MassShiftSequenceParser();
        }

        [Test]
        [TestCaseSource(typeof(GroundTruthTestData), nameof(GroundTruthTestData.CoreTestCases))]
        public void MzLibSerializer_CoreTestCases_SerializesCorrectly(GroundTruthTestData.TestCase testCase)
        {
            // Arrange - parse to get canonical form
            var canonical = _mzLibParser.Parse(testCase.MzLibFormat);
            Assert.That(canonical, Is.Not.Null);

            // Act
            var result = _mzLibSerializer.Serialize(canonical.Value);

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result, Is.EqualTo(testCase.MzLibFormat));
        }

        [Test]
        [TestCaseSource(typeof(GroundTruthTestData), nameof(GroundTruthTestData.CoreTestCases))]
        public void ChronologerSerializer_CoreTestCases_SerializesCorrectly(GroundTruthTestData.TestCase testCase)
        {
            // Arrange
            var canonical = _mzLibParser.Parse(testCase.MzLibFormat);
            Assert.That(canonical, Is.Not.Null);

            // Act
            var result = _chronologerSerializer.Serialize(canonical.Value);

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result, Is.EqualTo(testCase.ChronologerFormat));
        }

        [Test]
        public void MzLibSerializer_EmptySequence_ReturnsNull()
        {
            // Arrange
            var canonical = CanonicalSequence.Empty;

            // Act
            var result = _mzLibSerializer.Serialize(canonical, null, SequenceConversionHandlingMode.ReturnNull);

            // Assert
            Assert.That(result, Is.Null);
        }

        #region MassShift Serializer Tests

        [Test]
        [TestCaseSource(typeof(GroundTruthTestData), nameof(GroundTruthTestData.CoreTestCases))]
        public void MassShiftSerializer_CoreTestCases_SerializesCorrectly(GroundTruthTestData.TestCase testCase)
        {
            // Arrange - parse from MassShift format to get canonical form
            var canonical = _massShiftParser.Parse(testCase.MassShiftFormat);
            Assert.That(canonical, Is.Not.Null);

            // Act - serialize back to MassShift format
            var result = _massShiftSerializer.Serialize(canonical.Value);

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result, Is.EqualTo(testCase.MassShiftFormat));
        }

        [Test]
        [TestCaseSource(typeof(GroundTruthTestData), nameof(GroundTruthTestData.EdgeCases))]
        public void MassShiftSerializer_EdgeCases_SerializesCorrectly(GroundTruthTestData.TestCase testCase)
        {
            // Arrange
            var canonical = _massShiftParser.Parse(testCase.MassShiftFormat);
            Assert.That(canonical, Is.Not.Null);

            // Act
            var result = _massShiftSerializer.Serialize(canonical.Value);

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result, Is.EqualTo(testCase.MassShiftFormat));
        }

        [Test]
        [TestCaseSource(typeof(GroundTruthTestData), nameof(GroundTruthTestData.CoreTestCases))]
        public void MzLibToMassShift_CoreTestCases_ConvertsCorrectly(GroundTruthTestData.TestCase testCase)
        {
            // Arrange - parse from mzLib format
            var canonical = _mzLibParser.Parse(testCase.MzLibFormat);
            Assert.That(canonical, Is.Not.Null);

            // Act - serialize to MassShift format
            var result = _massShiftSerializer.Serialize(canonical.Value);

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result, Is.EqualTo(testCase.MassShiftFormat));
        }

        [Test]
        public void MassShiftSerializer_EmptySequence_ReturnsNull()
        {
            // Arrange
            var canonical = CanonicalSequence.Empty;

            // Act
            var result = _massShiftSerializer.Serialize(canonical, null, SequenceConversionHandlingMode.ReturnNull);

            // Assert
            Assert.That(result, Is.Null);
        }

        #endregion
    }
}
