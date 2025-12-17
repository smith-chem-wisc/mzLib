using NUnit.Framework;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using Omics;

namespace Test.Omics
{
    /// <summary>
    /// Minimal tests for ISpectralMatch interface behavior.
    /// Note: Most tests use TestSpectralMatch helper class defined in test files that need it.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class ISpectralMatchTests
    {
        /// <summary>
        /// Verifies ISpectralMatch.CompareTo orders by Score descending (higher scores first).
        /// Critical: Determines PSM ranking for FDR calculation and best-match selection.
        /// </summary>
        [Test]
        public void CompareTo_HigherScoreComesFirst()
        {
            var highScore = new TestSpectralMatch("file", "PEPTIDE", "PEPTIDE", score: 100, scanNumber: 1);
            var lowScore = new TestSpectralMatch("file", "PEPTIDE", "PEPTIDE", score: 50, scanNumber: 1);

            // Higher score should compare as greater (comes first when sorted descending)
            Assert.That(highScore.CompareTo(lowScore), Is.GreaterThan(0));
            Assert.That(lowScore.CompareTo(highScore), Is.LessThan(0));
        }

        /// <summary>
        /// Verifies GetIdentifiedBioPolymersWithSetMods returns empty enumeration (not null) when none provided.
        /// Critical: Prevents null reference exceptions in protein grouping code.
        /// </summary>
        [Test]
        public void GetIdentifiedBioPolymersWithSetMods_ReturnsEmptyNotNull()
        {
            var match = new TestSpectralMatch("file", "PEPTIDE", "PEPTIDE", score: 1, scanNumber: 1);

            var identified = match.GetIdentifiedBioPolymersWithSetMods();

            Assert.That(identified, Is.Not.Null);
            Assert.That(identified, Is.Empty);
        }

        #region Test Helper

        private class TestSpectralMatch : ISpectralMatch
        {
            public string FullFilePath { get; }
            public string FullSequence { get; }
            public string BaseSequence { get; }
            public double Score { get; }
            public int OneBasedScanNumber { get; }

            public TestSpectralMatch(string filePath, string fullSequence, string baseSequence, double score, int scanNumber)
            {
                FullFilePath = filePath ?? string.Empty;
                FullSequence = fullSequence ?? string.Empty;
                BaseSequence = baseSequence ?? string.Empty;
                Score = score;
                OneBasedScanNumber = scanNumber;
            }

            public int CompareTo(ISpectralMatch? other)
            {
                if (other is null) return 1;
                return Score.CompareTo(other.Score);
            }

            public IEnumerable<IBioPolymerWithSetMods> GetIdentifiedBioPolymersWithSetMods()
                => Enumerable.Empty<IBioPolymerWithSetMods>();
        }

        #endregion
    }
}