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
        private ChronologerSequenceSerializer _chronologerSerializer;
        private MzLibSequenceParser _mzLibParser;
        private MassShiftSequenceParser _massShiftParser;

        [SetUp]
        public void Setup()
        {
            _mzLibSerializer = new MzLibSequenceSerializer();
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
        [TestCaseSource(typeof(GroundTruthTestData), nameof(GroundTruthTestData.CoreTestCases))]
        public void MassShiftToMzLib_CoreTestCases_ConvertsCorrectly(GroundTruthTestData.TestCase testCase)
        {
            // Arrange
            var canonical = _massShiftParser.Parse(testCase.MassShiftFormat);
            Assert.That(canonical, Is.Not.Null);

            // Act
            var result = _mzLibSerializer.Serialize(canonical.Value);

            // Assert
            Assert.That(result, Is.EqualTo(testCase.MzLibFormat));
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
    }
}
