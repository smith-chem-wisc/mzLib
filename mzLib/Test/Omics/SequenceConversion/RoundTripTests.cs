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
        private MassShiftSequenceParser _massShiftParser;
        private MzLibSequenceSerializer _mzLibSerializer;
        private ChronologerSequenceSerializer _chronologerSerializer;

        [SetUp]
        public void Setup()
        {
            _mzLibParser = new MzLibSequenceParser();
            _massShiftParser = new MassShiftSequenceParser();
            _mzLibSerializer = new MzLibSequenceSerializer();
            _chronologerSerializer = new ChronologerSequenceSerializer();
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
        [TestCaseSource(typeof(GroundTruthTestData), nameof(GroundTruthTestData.CoreTestCases))]
        public void MassShiftToMzLibAndBack_CoreTestCases_PreservesSemantics(GroundTruthTestData.TestCase testCase)
        {
            // Arrange
            var originalMassShift = testCase.MassShiftFormat;

            // Act - MassShift to Canonical
            var canonical1 = _massShiftParser.Parse(originalMassShift);
            Assert.That(canonical1, Is.Not.Null);

            // Canonical to MzLib
            var mzLib = _mzLibSerializer.Serialize(canonical1.Value);
            Assert.That(mzLib, Is.Not.Null);
            Assert.That(mzLib, Is.EqualTo(testCase.MzLibFormat));

            // MzLib back to Canonical
            var canonical2 = _mzLibParser.Parse(mzLib);
            Assert.That(canonical2, Is.Not.Null);

            // Assert - canonical forms should be semantically equivalent
            Assert.That(canonical2.Value.BaseSequence, Is.EqualTo(canonical1.Value.BaseSequence));
            Assert.That(canonical2.Value.ResidueModifications.Count(), Is.EqualTo(canonical1.Value.ResidueModifications.Count()));
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
