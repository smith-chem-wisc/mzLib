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

                // Parallel production means results may be unordered; assert all sequences are present
                var resultSequences = results.Select(r => r.Peptide.BaseSequence).OrderBy(s => s).ToList();
                var expectedSequences = sequences.OrderBy(s => s).ToList();
                Assert.That(resultSequences, Is.EqualTo(expectedSequences),
                    $"{predictor.PredictorName} did not return all input peptides");
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
            var peptides = new[]
            {
                new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>()),  // valid
                new PeptideWithSetModifications("PEP", new Dictionary<string, Modification>()),       // too short for all
                new PeptideWithSetModifications("PEPTIDER", new Dictionary<string, Modification>()), // valid
                new PeptideWithSetModifications("", new Dictionary<string, Modification>()),          // empty for all
                new PeptideWithSetModifications("ANALYSIS", new Dictionary<string, Modification>())  // valid
            };

            foreach (var predictor in _predictors)
            {
                var results = predictor.PredictRetentionTimeEquivalents(peptides);

                Assert.That(results.Count, Is.EqualTo(5),
                    $"{predictor.PredictorName} returned wrong count for mixed batch");

                var bySequence = results.ToDictionary(r => r.Peptide.BaseSequence);

                Assert.That(bySequence["PEPTIDE"].PredictedValue.HasValue, Is.True, $"{predictor.PredictorName} PEPTIDE should succeed");
                Assert.That(bySequence["PEP"].PredictedValue.HasValue, Is.False, $"{predictor.PredictorName} PEP should fail (too short)");
                Assert.That(bySequence["PEP"].FailureReason, Is.EqualTo(RetentionTimeFailureReason.SequenceTooShort));
                Assert.That(bySequence["PEPTIDER"].PredictedValue.HasValue, Is.True, $"{predictor.PredictorName} PEPTIDER should succeed");
                Assert.That(bySequence[""].PredictedValue.HasValue, Is.False, $"{predictor.PredictorName} empty should fail");
                Assert.That(bySequence[""].FailureReason, Is.EqualTo(RetentionTimeFailureReason.EmptySequence));
                Assert.That(bySequence["ANALYSIS"].PredictedValue.HasValue, Is.True, $"{predictor.PredictorName} ANALYSIS should succeed");
            }
        }

        [Test]
        public void PredictRetentionTimeEquivalents_SequenceTooLong_ChronologerRejects()
        {
            // Only Chronologer has a MaxSequenceLength (50); SSRCalc3 and CZE do not
            var longPeptide = new PeptideWithSetModifications(new string('A', 60), new Dictionary<string, Modification>());

            var chronologer = new ChronologerRetentionTimePredictor();
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
            var peptides = Enumerable.Range(0, 500)
                .Select(i => new PeptideWithSetModifications(
                    $"PEPTIDE{i % 20 + 1}", new Dictionary<string, Modification>()))
                .ToList();

            foreach (var predictor in _predictors)
            {
                IReadOnlyList<(double? PredictedValue, IRetentionPredictable Peptide, RetentionTimeFailureReason? FailureReason)> results = null;
                Assert.DoesNotThrow(() => results = predictor.PredictRetentionTimeEquivalents(peptides),
                    $"{predictor.PredictorName} threw on large batch");
                Assert.That(results.Count, Is.EqualTo(500),
                    $"{predictor.PredictorName} returned wrong count for large batch");
            }
        }

        #endregion

        #region Factory and IDisposable Tests

        [Test]
        public void Factory_Create_SSRCalc3_ReturnsFunctionalPredictor()
        {
            var predictor = IRetentionTimePredictor.Create(IRetentionTimePredictor.PredictorType.SSRCalc3);
            var peptides = new[] { new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>()) };

            var results = predictor.PredictRetentionTimeEquivalents(peptides);

            Assert.That(results.Count, Is.EqualTo(1));
            Assert.That(results[0].PredictedValue.HasValue, Is.True);
        }

        [Test]
        public void Factory_Create_Chronologer_ReturnsFunctionalPredictor()
        {
            var predictor = IRetentionTimePredictor.Create(IRetentionTimePredictor.PredictorType.Chronologer);
            var peptides = new[] { new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>()) };

            var results = predictor.PredictRetentionTimeEquivalents(peptides);

            Assert.That(results.Count, Is.EqualTo(1));
            Assert.That(results[0].PredictedValue.HasValue, Is.True);
        }

        [Test]
        public void Factory_Create_CZE_ReturnsFunctionalPredictor()
        {
            var predictor = IRetentionTimePredictor.Create(IRetentionTimePredictor.PredictorType.CZE);
            var peptides = new[] { new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>()) };

            var results = predictor.PredictRetentionTimeEquivalents(peptides);

            Assert.That(results.Count, Is.EqualTo(1));
            Assert.That(results[0].PredictedValue.HasValue, Is.True);
        }

        [Test]
        public void Factory_Create_InvalidType_ThrowsArgumentOutOfRangeException()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                IRetentionTimePredictor.Create((IRetentionTimePredictor.PredictorType)999));
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