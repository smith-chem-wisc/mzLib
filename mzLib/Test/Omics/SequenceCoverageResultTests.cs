using NUnit.Framework;
using Omics.BioPolymerGroup;
using System.Diagnostics.CodeAnalysis;

namespace Test.Omics
{
    /// <summary>
    /// Tests for BioPolymerGroup.SequenceCoverageResult nested class.
    /// Integration tests are in BioPolymerGroupSequenceCoverageTests.cs.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class SequenceCoverageResultTests
    {
        /// <summary>
        /// Verifies constructor initializes all lists as empty, non-null collections.
        /// Critical: Prevents null reference exceptions when CalculateSequenceCoverage populates data
        /// or when ToString accesses the lists before calculation.
        /// </summary>
        [Test]
        public void Constructor_InitializesAllListsAsEmptyNonNull()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            Assert.Multiple(() =>
            {
                Assert.That(result.SequenceCoverageFraction, Is.Not.Null.And.Empty);
                Assert.That(result.SequenceCoverageDisplayList, Is.Not.Null.And.Empty);
                Assert.That(result.SequenceCoverageDisplayListWithMods, Is.Not.Null.And.Empty);
                Assert.That(result.FragmentSequenceCoverageDisplayList, Is.Not.Null.And.Empty);
                Assert.That(result.ModsInfo, Is.Not.Null.And.Empty);
            });
        }
    }
}