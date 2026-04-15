using NUnit.Framework;
using Chromatography.RetentionTimePrediction;
using Chromatography.RetentionTimePrediction.Chronologer;
using Chromatography.RetentionTimePrediction.SSRCalc;
using Chromatography.RetentionTimePrediction.CZE;
using Proteomics.ProteolyticDigestion;
using Omics.Modifications;
using System.Collections.Generic;
using System.Linq;
using System;
using Omics.SequenceConversion;

namespace Test.RetentionTimePrediction
{
    /// <summary>
    /// Tests for batch retention time equivalent prediction via PredictRetentionTimeEquivalents().
    /// Covers correctness, consistency with single-item predictions, mixed valid/invalid inputs,
    /// edge cases, and all three predictor types.
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class RetentionTimeEquivalentBatchTests
    {
        private List<IRetentionTimePredictor> _predictors;

        [SetUp]
        public void Setup()
        {
            _predictors = new List<IRetentionTimePredictor>
            {
                new SSRCalc3RetentionTimePredictor(),
                new ChronologerRetentionTimePredictor(SequenceConversionHandlingMode.RemoveIncompatibleElements),
                new CZERetentionTimePredictor(SequenceConversionHandlingMode.UsePrimarySequence, 1.0, 1.0)
            };
        }

        [TearDown]
        public void TearDown()
        {
            foreach (var predictor in _predictors)
                (predictor as IDisposable)?.Dispose();
        }

        #region Basic Batch Functionality

        [Test]
        public void PredictRetentionTimeEquivalents_ValidBatch_ReturnsSameCountAsInput()
        {
            var peptides = new[]
            {
                new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>()),
                new PeptideWithSetModifications("KALEIDIC", new Dictionary<string, Modification>()),
                new PeptideWithSetModifications("ANALYSIS", new Dictionary<string, Modification>())
            };

            foreach (var predictor in _predictors)
            {
                var results = predictor.PredictRetentionTimeEquivalents(peptides);

                Assert.That(results.Count, Is.EqualTo(peptides.Length),
                    $"{predictor.PredictorName} returned wrong count");
            }
        }

        [Test]
        public void PredictRetentionTimeEquivalents_ValidBatch_AllSucceed()
        {
            var peptides = new[]
            {
                new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>()),
                new PeptideWithSetModifications("PEPTIDER", new Dictionary<string, Modification>()),
                new PeptideWithSetModifications("ANALYSIS", new Dictionary<string, Modification>())
            };

            foreach (var predictor in _predictors)
            {
                var results = predictor.PredictRetentionTimeEquivalents(peptides);

                Assert.That(results.All(r => r.PredictedValue.HasValue), Is.True,
                    $"{predictor.PredictorName} had unexpected failures: " +
                    string.Join(", ", results.Where(r => !r.PredictedValue.HasValue).Select(r => r.FailureReason)));
            }
        }

        [Test]
        public void PredictRetentionTimeEquivalents_ContainsAllInputPeptides()
        {
            var sequences = new[] { "PEPTIDE", "PEPTIDER", "ANALYSIS", "TESTING", "PEPTIDER" };
            var peptides = sequences
                .Select(s => new PeptideWithSetModifications(s, new Dictionary<string, Modification>()))
                .ToList();

            foreach (var predictor in _predictors)
            {
                var results = predictor.PredictRetentionTimeEquivalents(peptides);

                // Results may be unordered and may collapse duplicates; only require that every
                // distinct input sequence is represented at least once in the output. A sorted-list
                // comparison would let a swapped-out duplicate pass silently, so compare distinct
                // sets instead.
                var expectedDistinct = new HashSet<string>(sequences);
                var resultDistinct = new HashSet<string>(results.Select(r => r.Peptide.BaseSequence));

                Assert.That(resultDistinct, Is.SupersetOf(expectedDistinct),
                    $"{predictor.PredictorName} did not return every distinct input peptide. " +
                    $"Missing: {string.Join(", ", expectedDistinct.Except(resultDistinct))}");
            }
        }

        #endregion

        #region Consistency With Single Predictions

        [Test]
        public void PredictRetentionTimeEquivalents_MatchesSinglePredictions()
        {
            var peptides = new[]
            {
                new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>()),
                new PeptideWithSetModifications("PEPTIDER", new Dictionary<string, Modification>()),
                new PeptideWithSetModifications("ANALYSIS", new Dictionary<string, Modification>())
            };

            foreach (var predictor in _predictors)
            {
                var batchResults = predictor.PredictRetentionTimeEquivalents(peptides)
                    .ToDictionary(r => r.Peptide.BaseSequence);

                foreach (var peptide in peptides)
                {
                    var singleResult = predictor.PredictRetentionTimeEquivalent(peptide, out var singleReason);
                    var batchResult = batchResults[peptide.BaseSequence];

                    Assert.That(batchResult.PredictedValue, Is.EqualTo(singleResult),
                        $"{predictor.PredictorName} batch result differed from single for {peptide.BaseSequence}");
                    Assert.That(batchResult.FailureReason, Is.EqualTo(singleReason),
                        $"{predictor.PredictorName} batch failure reason differed from single for {peptide.BaseSequence}");
                }
            }
        }

        [Test]
        public void PredictRetentionTimeEquivalents_CalledTwice_ReturnsSameValues()
        {
            var peptides = new[]
            {
                new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>()),
                new PeptideWithSetModifications("KALEIDIC", new Dictionary<string, Modification>())
            };

            foreach (var predictor in _predictors)
            {
                var results1 = predictor.PredictRetentionTimeEquivalents(peptides)
                    .ToDictionary(r => r.Peptide.BaseSequence);
                var results2 = predictor.PredictRetentionTimeEquivalents(peptides)
                    .ToDictionary(r => r.Peptide.BaseSequence);

                foreach (var peptide in peptides)
                {
                    Assert.That(results1[peptide.BaseSequence].PredictedValue,
                        Is.EqualTo(results2[peptide.BaseSequence].PredictedValue),
                        $"{predictor.PredictorName} gave inconsistent batch results for {peptide.BaseSequence}");
                }
            }
        }

        #endregion

        #region Mixed Valid and Invalid Inputs

        [Test]
        public void PredictRetentionTimeEquivalents_MixedBatch_CorrectSuccessAndFailureCounts()
        {
            var validPeptide = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());
            var shortPeptide = new PeptideWithSetModifications("PEP", new Dictionary<string, Modification>());
            var validPeptide2 = new PeptideWithSetModifications("PEPTIDER", new Dictionary<string, Modification>());
            var emptyPeptide = new PeptideWithSetModifications("", new Dictionary<string, Modification>());
            var validPeptide3 = new PeptideWithSetModifications("ANALYSIS", new Dictionary<string, Modification>());

            var peptides = new[] { validPeptide, shortPeptide, validPeptide2, emptyPeptide, validPeptide3 };

            foreach (var predictor in _predictors)
            {
                var results = predictor.PredictRetentionTimeEquivalents(peptides);

                Assert.That(results.Count, Is.EqualTo(5),
                    $"{predictor.PredictorName} returned wrong count for mixed batch");

                // Look up by reference equality on the original peptide objects rather than by
                // BaseSequence. Keying on BaseSequence is fragile: an empty or null key would
                // blow up ToDictionary, and any two peptides sharing a sequence would collide.
                var resultForValid = results.First(r => ReferenceEquals(r.Peptide, validPeptide));
                var resultForShort = results.First(r => ReferenceEquals(r.Peptide, shortPeptide));
                var resultForValid2 = results.First(r => ReferenceEquals(r.Peptide, validPeptide2));
                var resultForEmpty = results.First(r => ReferenceEquals(r.Peptide, emptyPeptide));
                var resultForValid3 = results.First(r => ReferenceEquals(r.Peptide, validPeptide3));

                Assert.That(resultForValid.PredictedValue.HasValue, Is.True, $"{predictor.PredictorName} PEPTIDE should succeed");
                Assert.That(resultForShort.PredictedValue.HasValue, Is.False, $"{predictor.PredictorName} PEP should fail (too short)");
                Assert.That(resultForShort.FailureReason, Is.EqualTo(RetentionTimeFailureReason.SequenceTooShort));
                Assert.That(resultForValid2.PredictedValue.HasValue, Is.True, $"{predictor.PredictorName} PEPTIDER should succeed");
                Assert.That(resultForEmpty.PredictedValue.HasValue, Is.False, $"{predictor.PredictorName} empty should fail");
                Assert.That(resultForEmpty.FailureReason, Is.EqualTo(RetentionTimeFailureReason.EmptySequence));
                Assert.That(resultForValid3.PredictedValue.HasValue, Is.True, $"{predictor.PredictorName} ANALYSIS should succeed");
            }
        }

        [Test]
        public void PredictRetentionTimeEquivalents_SequenceTooLong_ChronologerRejects()
        {
            // Only Chronologer has a MaxSequenceLength (50); SSRCalc3 and CZE do not
            var longPeptide = new PeptideWithSetModifications(new string('A', 60), new Dictionary<string, Modification>());

            using var chronologer = new ChronologerRetentionTimePredictor();
            var results = chronologer.PredictRetentionTimeEquivalents(new[] { longPeptide });

            Assert.That(results[0].PredictedValue.HasValue, Is.False);
            Assert.That(results[0].FailureReason, Is.EqualTo(RetentionTimeFailureReason.SequenceTooLong));
        }

        [Test]
        public void PredictRetentionTimeEquivalents_AllInvalid_AllFailWithReasons()
        {
            var peptides = new[]
            {
                new PeptideWithSetModifications("PEP", new Dictionary<string, Modification>()),  // too short for all
                new PeptideWithSetModifications("", new Dictionary<string, Modification>())       // empty for all
            };

            foreach (var predictor in _predictors)
            {
                var results = predictor.PredictRetentionTimeEquivalents(peptides);

                Assert.That(results.All(r => !r.PredictedValue.HasValue), Is.True,
                    $"{predictor.PredictorName} should have all failures");
                Assert.That(results.All(r => r.FailureReason.HasValue), Is.True,
                    $"{predictor.PredictorName} should provide failure reasons for all");
            }
        }

        #endregion

        #region Edge Cases

        [Test]
        public void PredictRetentionTimeEquivalents_EmptyList_ReturnsEmptyList()
        {
            var peptides = new List<PeptideWithSetModifications>();

            foreach (var predictor in _predictors)
            {
                var results = predictor.PredictRetentionTimeEquivalents(peptides);

                Assert.That(results.Count, Is.EqualTo(0),
                    $"{predictor.PredictorName} should return empty list for empty input");
            }
        }

        [Test]
        public void PredictRetentionTimeEquivalents_SinglePeptide_ReturnsSingleResult()
        {
            var peptides = new[]
            {
                new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>())
            };

            foreach (var predictor in _predictors)
            {
                var results = predictor.PredictRetentionTimeEquivalents(peptides);

                Assert.That(results.Count, Is.EqualTo(1),
                    $"{predictor.PredictorName} should return single result");
                Assert.That(results[0].PredictedValue.HasValue, Is.True,
                    $"{predictor.PredictorName} should succeed for valid single-item batch");
            }
        }

        [Test]
        public void PredictRetentionTimeEquivalents_DuplicatePeptides_ReturnsSameValueForEach()
        {
            var peptide = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());
            var peptides = Enumerable.Repeat(peptide, 5).ToList();

            foreach (var predictor in _predictors)
            {
                var results = predictor.PredictRetentionTimeEquivalents(peptides);

                Assert.That(results.Count, Is.EqualTo(5));
                var distinctValues = results.Select(r => r.PredictedValue).Distinct().Count();
                Assert.That(distinctValues, Is.EqualTo(1),
                    $"{predictor.PredictorName} should return same value for duplicate peptides");
            }
        }

        [Test]
        public void PredictRetentionTimeEquivalents_LargeBatch_CompletesSuccessfully()
        {
            var aminoAcids = "ACDEFGHIKLMNPQRSTVWY";
            var peptides = Enumerable.Range(0, 500)
                .Select(i =>
                {
                    var chars = Enumerable.Range(0, 7)
                        .Select(j => aminoAcids[(i * 7 + j) % aminoAcids.Length])
                        .ToArray();
                    return new PeptideWithSetModifications(
                        new string(chars), new Dictionary<string, Modification>());
                })
                .ToList();

            foreach (var predictor in _predictors)
            {
                var results = predictor.PredictRetentionTimeEquivalents(peptides);

                Assert.That(results.Count, Is.EqualTo(500),
                    $"{predictor.PredictorName} returned wrong count for large batch");

                var sampleIndices = new[] { 0, 1, 49, 50, 99, 100, 249, 250, 499 };
                foreach (var idx in sampleIndices)
                {
                    var singleValue = predictor.PredictRetentionTimeEquivalent(peptides[idx], out var singleReason);
                    var batchResult = results[idx];

                    Assert.That(batchResult.Peptide.BaseSequence, Is.EqualTo(peptides[idx].BaseSequence),
                        $"{predictor.PredictorName} peptide mismatch at index {idx}");
                    Assert.That(batchResult.PredictedValue, Is.EqualTo(singleValue),
                        $"{predictor.PredictorName} value mismatch at index {idx} for {peptides[idx].BaseSequence}");
                    Assert.That(batchResult.FailureReason, Is.EqualTo(singleReason),
                        $"{predictor.PredictorName} failure reason mismatch at index {idx}");
                }
            }
        }

        [Test]
        public void PredictRetentionTimeEquivalents_MultiThreaded_MatchesSingleThreaded()
        {
            var peptides = Enumerable.Range(0, 500)
                .Select(i => new PeptideWithSetModifications(
                    $"PEPTIDE{i % 20 + 1}", new Dictionary<string, Modification>()))
                .ToList();

            foreach (var predictor in _predictors)
            {
                var singleThreaded = predictor.PredictRetentionTimeEquivalents(peptides, maxThreads: 1)
                    .OrderBy(r => r.Peptide.BaseSequence)
                    .ThenBy(r => r.PredictedValue)
                    .ToList();

                var multiThreaded = predictor.PredictRetentionTimeEquivalents(peptides, maxThreads: 4)
                    .OrderBy(r => r.Peptide.BaseSequence)
                    .ThenBy(r => r.PredictedValue)
                    .ToList();

                Assert.That(multiThreaded.Count, Is.EqualTo(singleThreaded.Count),
                    $"{predictor.PredictorName} multi-threaded count mismatch");

                for (int i = 0; i < singleThreaded.Count; i++)
                {
                    Assert.That(multiThreaded[i].PredictedValue, Is.EqualTo(singleThreaded[i].PredictedValue),
                        $"{predictor.PredictorName} value mismatch at index {i} for {singleThreaded[i].Peptide.BaseSequence}");
                    Assert.That(multiThreaded[i].FailureReason, Is.EqualTo(singleThreaded[i].FailureReason),
                        $"{predictor.PredictorName} failure reason mismatch at index {i}");
                }
            }
        }

        [Test]
        public void PredictRetentionTimeEquivalents_EmptyList_MultiThreaded_ReturnsEmpty()
        {
            var peptides = new List<PeptideWithSetModifications>();

            foreach (var predictor in _predictors)
            {
                var results = predictor.PredictRetentionTimeEquivalents(peptides, maxThreads: 4);

                Assert.That(results.Count, Is.EqualTo(0),
                    $"{predictor.PredictorName} should return empty list for empty input with multiple threads");
            }
        }

        [Test]
        public void StreamRetentionTimeEquivalents_LargeBatch_ReturnsAllResults()
        {
            var peptides = Enumerable.Range(0, 500)
                .Select(i => new PeptideWithSetModifications(
                    $"PEPTIDE{i % 20 + 1}", new Dictionary<string, Modification>()))
                .ToList();

            foreach (var predictor in _predictors)
            {
                var streamedResults = predictor.StreamRetentionTimeEquivalents(peptides, maxThreads: 4).ToList();

                Assert.That(streamedResults.Count, Is.EqualTo(500),
                    $"{predictor.PredictorName} streamed wrong count");

                var batchResults = predictor.PredictRetentionTimeEquivalents(peptides, maxThreads: 1)
                    .OrderBy(r => r.Peptide.BaseSequence)
                    .ThenBy(r => r.PredictedValue)
                    .ToList();
                var sortedStreamed = streamedResults
                    .OrderBy(r => r.Peptide.BaseSequence)
                    .ThenBy(r => r.PredictedValue)
                    .ToList();

                for (int i = 0; i < batchResults.Count; i++)
                {
                    Assert.That(sortedStreamed[i].PredictedValue, Is.EqualTo(batchResults[i].PredictedValue),
                        $"{predictor.PredictorName} streamed value mismatch for {batchResults[i].Peptide.BaseSequence}");
                }
            }
        }

        [Test]
        public void StreamRetentionTimeEquivalents_EmptyList_MultiThreaded_YieldsNothing()
        {
            var peptides = new List<PeptideWithSetModifications>();

            foreach (var predictor in _predictors)
            {
                var results = predictor.StreamRetentionTimeEquivalents(peptides, maxThreads: 4).ToList();

                Assert.That(results.Count, Is.EqualTo(0),
                    $"{predictor.PredictorName} should yield nothing for empty input");
            }
        }

        [Test]
        public void StreamRetentionTimeEquivalents_SingleThread_MatchesBatchResults()
        {
            var peptides = new[]
            {
                new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>()),
                new PeptideWithSetModifications("PEPTIDER", new Dictionary<string, Modification>()),
                new PeptideWithSetModifications("ANALYSIS", new Dictionary<string, Modification>())
            };

            foreach (var predictor in _predictors)
            {
                var batchResults = predictor.PredictRetentionTimeEquivalents(peptides)
                    .OrderBy(r => r.Peptide.BaseSequence).ToList();
                var streamedResults = predictor.StreamRetentionTimeEquivalents(peptides)
                    .OrderBy(r => r.Peptide.BaseSequence).ToList();

                Assert.That(streamedResults.Count, Is.EqualTo(batchResults.Count),
                    $"{predictor.PredictorName} stream count differs from batch");

                for (int i = 0; i < batchResults.Count; i++)
                {
                    Assert.That(streamedResults[i].PredictedValue, Is.EqualTo(batchResults[i].PredictedValue),
                        $"{predictor.PredictorName} stream/batch value mismatch for {batchResults[i].Peptide.BaseSequence}");
                }
            }
        }

        [Test]
        public void PredictRetentionTimeEquivalents_ThreadCountExceedsItems_NoDeadlock()
        {
            var peptides = new[]
            {
                new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>()),
                new PeptideWithSetModifications("PEPTIDER", new Dictionary<string, Modification>())
            };

            foreach (var predictor in _predictors)
            {
                var results = predictor.PredictRetentionTimeEquivalents(peptides, maxThreads: 8);

                Assert.That(results.Count, Is.EqualTo(2),
                    $"{predictor.PredictorName} should handle maxThreads > item count gracefully");
                Assert.That(results.All(r => r.PredictedValue.HasValue), Is.True,
                    $"{predictor.PredictorName} all valid peptides should succeed");
            }
        }

        #endregion

        #region Factory and IDisposable Tests

        [Test]
        public void Factory_Create_SSRCalc3_ReturnsFunctionalPredictor()
        {
            using var predictor = RetentionTimePredictorFactory.Create(PredictorType.SSRCalc3);
            var peptides = new[] { new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>()) };

            var results = predictor.PredictRetentionTimeEquivalents(peptides);

            Assert.That(results.Count, Is.EqualTo(1));
            Assert.That(results[0].PredictedValue.HasValue, Is.True);
        }

        [Test]
        public void Factory_Create_Chronologer_ReturnsFunctionalPredictor()
        {
            using var predictor = RetentionTimePredictorFactory.Create(PredictorType.Chronologer);
            var peptides = new[] { new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>()) };

            var results = predictor.PredictRetentionTimeEquivalents(peptides);

            Assert.That(results.Count, Is.EqualTo(1));
            Assert.That(results[0].PredictedValue.HasValue, Is.True);
        }

        [Test]
        public void Factory_Create_CZE_ReturnsFunctionalPredictor()
        {
            using var predictor = RetentionTimePredictorFactory.Create(PredictorType.CZE);
            var peptides = new[] { new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>()) };

            var results = predictor.PredictRetentionTimeEquivalents(peptides);

            Assert.That(results.Count, Is.EqualTo(1));
            Assert.That(results[0].PredictedValue.HasValue, Is.True);
        }

        [Test]
        public void Factory_Create_InvalidType_ThrowsArgumentOutOfRangeException()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                RetentionTimePredictorFactory.Create((PredictorType)999));
        }

        [Test]
        public void Dispose_AfterBatchPrediction_DoesNotThrow()
        {
            var predictor = new ChronologerRetentionTimePredictor();
            var peptides = new[] { new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>()) };

            predictor.PredictRetentionTimeEquivalents(peptides);

            Assert.DoesNotThrow(() => predictor.Dispose());
        }

        #endregion
    }
}