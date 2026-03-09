using System.Linq;
using NUnit.Framework;
using Omics.SequenceConversion;

namespace Test.Omics.SequenceConversion
{
    /// <summary>
    /// Tests for round-trip conversions: format -> canonical -> format.
    /// </summary>
    [TestFixture]
    public class RoundTripTests
    {
        private MzLibSequenceParser _mzLibParser;
        private MzLibSequenceSerializer _mzLibSerializer;

        [SetUp]
        public void Setup()
        {
            _mzLibParser = new MzLibSequenceParser();
            _mzLibSerializer = new MzLibSequenceSerializer();
        }

        [Test]
        [TestCaseSource(typeof(GroundTruthTestData), nameof(GroundTruthTestData.CoreTestCases))]
        public void MzLibRoundTrip_CoreTestCases_PreservesSequence(GroundTruthTestData.TestCase testCase)
        {
            // Arrange
            var originalSequence = testCase.MzLibFormat;

            // Act - parse then serialize
            var parsed = _mzLibParser.Parse(originalSequence);
            Assert.That(parsed, Is.Not.Null);
            
            var serialized = _mzLibSerializer.Serialize(parsed.Value);

            // Assert
            Assert.That(serialized, Is.EqualTo(originalSequence));
        }

        [Test]
        public void MzLibRoundTrip_MultipleIterations_RemainsStable()
        {
            // Arrange
            var original = "PEPTM[Oxidation on M]IDE";

            // Act & Assert - perform multiple round trips
            var current = original;
            for (int i = 0; i < 5; i++)
            {
                var parsed = _mzLibParser.Parse(current);
                Assert.That(parsed, Is.Not.Null);
                
                var serialized = _mzLibSerializer.Serialize(parsed.Value);
                Assert.That(serialized, Is.EqualTo(original));
                
                current = serialized;
            }
        }
    }
}
