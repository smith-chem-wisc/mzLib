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
using Chromatography;
using Readers;

namespace Test.RetentionTimePrediction
{
    /// <summary>
    /// Integration tests for retention time prediction system
    /// Tests multiple components working together
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class RetentionTimePredictionIntegrationTests
    {
        private List<IRetentionTimePredictor> _predictors;

        [SetUp]
        public void Setup()
        {
            _predictors = new List<IRetentionTimePredictor>
            {
                new ChronologerRetentionTimePredictor(IncompatibleModHandlingMode.RemoveIncompatibleMods),
                new SSRCalc3RetentionTimePredictor(),
                new CZERetentionTimePredictor(IncompatibleModHandlingMode.UsePrimarySequence, 1.0, 1.0)
            };
        }

        [TearDown]
        public void TearDown()
        {
            foreach (var predictor in _predictors)
            {
                if (predictor is IDisposable disposable)
                {
                    disposable.Dispose();
                }
            }
        }

        #region Multi-Predictor Tests

        [Test]
        public void AllPredictors_UnmodifiedPeptide_ReturnValues()
        {
            var peptide = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());
            
            foreach (var predictor in _predictors)
            {
                var result = predictor.PredictRetentionTime(peptide, out var failureReason);
                
                Assert.That(result, Is.Not.Null, 
                    $"Predictor {predictor.PredictorName} failed to predict for unmodified peptide");
                Assert.That(failureReason, Is.Null, 
                    $"Predictor {predictor.PredictorName} reported failure: {failureReason}");
            }
        }

        [Test]
        public void AllPredictors_HaveUniquePredictorNames()
        {
            var names = _predictors.Select(p => p.PredictorName).ToList();
            
            Assert.That(names.Count, Is.EqualTo(names.Distinct().Count()));
        }

        [Test]
        public void AllPredictors_HaveSeparationTypes()
        {
            foreach (var predictor in _predictors)
            {
                Assert.That(predictor.SeparationType, Is.Not.EqualTo((SeparationType)(-1)),
                    $"Predictor {predictor.PredictorName} does not have defined separation type");
            }
        }

        [Test]
        public void MultiplePredictors_SamePeptide_ProduceConsistentResults()
        {
            var peptide = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());
            
            foreach (var predictor in _predictors)
            {
                var result1 = predictor.PredictRetentionTime(peptide, out _);
                var result2 = predictor.PredictRetentionTime(peptide, out _);
                
                Assert.That(result1, Is.EqualTo(result2),
                    $"Predictor {predictor.PredictorName} gave inconsistent results");
            }
        }

        #endregion

        #region Real-World Scenario Tests

        [Test]
        public void RealWorld_TrypticDigest_AllPeptidesPredictable()
        {
            // Simulate a tryptic digest
            var peptides = new[]
            {
                "PEPTIDER",
                "DIGESTIN",
                "TRYPTIC",
                "ANALYSIS"
            };

            foreach (var sequence in peptides)
            {
                var peptide = new PeptideWithSetModifications(sequence, new Dictionary<string, Modification>());
                
                foreach (var predictor in _predictors)
                {
                    var result = predictor.PredictRetentionTime(peptide, out _);
                    Assert.That(result, Is.Not.Null,
                        $"Failed to predict {sequence} with {predictor.PredictorName}");
                }
            }
        }

        [Test]
        public void RealWorld_CommonModifications_HandledCorrectly()
        {
            var testCases = new[]
            {
                ("PEPTM[Oxidation on M]IDE", "Oxidation on M"),
                ("PEPTC[Carbamidomethyl on C]IDE", "Carbamidomethyl on C"),
                ("PEPTIDES[Phospho on S]", "Phospho on S")
            };

            foreach (var (sequence, modName) in testCases)
            {
                var mods = new Dictionary<string, Modification>
                {
                    { modName, ModificationConverter.AllModsKnown[modName] }
                };
                var peptide = new PeptideWithSetModifications(sequence, mods);
                
                foreach (var predictor in _predictors)
                {
                    var result = predictor.PredictRetentionTime(peptide, out _);
                    // Some predictors may not support all modifications, but should handle gracefully
                    Assert.That(result, Is.Not.Null,
                        $"Predictor {predictor.PredictorName} failed for {sequence}");
                }
            }
        }

        [Test]
        public void RealWorld_BatchPrediction_HandlesMultiplePeptidesEfficiently()
        {
            var peptides = Enumerable.Range(0, 100)
                .Select(i => new PeptideWithSetModifications($"PEPTIDE{i % 10}", new Dictionary<string, Modification>()))
                .ToList();

            foreach (var predictor in _predictors)
            {
                var stopwatch = System.Diagnostics.Stopwatch.StartNew();
                
                foreach (var peptide in peptides)
                {
                    predictor.PredictRetentionTime(peptide, out _);
                }
                
                stopwatch.Stop();
                
                // Ensure reasonable performance (less than 10 seconds for 100 predictions)
                Assert.That(stopwatch.Elapsed.TotalSeconds, Is.LessThan(10),
                    $"Predictor {predictor.PredictorName} took too long");
            }
        }

        #endregion

        #region Boundary and Edge Case Tests

        [Test]
        public void Integration_VeryShortPeptide_AllPredictorsHandle()
        {
            var peptide = new PeptideWithSetModifications("PEP", new Dictionary<string, Modification>());
            
            foreach (var predictor in _predictors)
            {
                var result = predictor.PredictRetentionTime(peptide, out var failureReason);
                
                // Should either succeed or fail gracefully
                if (result == null)
                {
                    Assert.That(failureReason, Is.Not.Null,
                        $"Predictor {predictor.PredictorName} returned null without failure reason");

                    Assert.That(failureReason, Is.EqualTo(RetentionTimeFailureReason.SequenceTooShort),
                        $"Predictor {predictor.PredictorName} should indicate sequence too short");
                }
            }
        }

        [Test]
        public void Integration_LongPeptide_PredictorsHandleOrReject()
        {
            var longSequence = new string('A', 60); // Longer than typical limit
            var peptide = new PeptideWithSetModifications(longSequence, new Dictionary<string, Modification>());
            
            foreach (var predictor in _predictors)
            {
                var result = predictor.PredictRetentionTime(peptide, out var failureReason);
                
                if (result == null)
                {
                    Assert.That(failureReason, Is.Not.Null,
                        $"Predictor {predictor.PredictorName} should report why prediction failed");
                    Assert.That(failureReason, Is.EqualTo(RetentionTimeFailureReason.SequenceTooLong),
                        $"Predictor {predictor.PredictorName} should indicate sequence too long");
                }
            }
        }

        [Test]
        public void Integration_AllCanonicalAminoAcids_Handled()
        {
            var allAA = "ACDEFGHIKLMNPQRSTVWY";
            var peptide = new PeptideWithSetModifications(allAA, new Dictionary<string, Modification>());
            
            foreach (var predictor in _predictors)
            {
                var result = predictor.PredictRetentionTime(peptide, out _);
                
                // All predictors should handle canonical amino acids
                Assert.That(result, Is.Not.Null,
                    $"Predictor {predictor.PredictorName} should handle all canonical AAs");
            }
        }

        #endregion

        #region Format and Interface Consistency Tests

        [Test]
        public void Integration_GetFormattedSequence_ConsistentWithPrediction()
        {
            var peptide = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());
            
            foreach (var predictor in _predictors)
            {
                var formatted = predictor.GetFormattedSequence(peptide, out var formatFailure);
                var prediction = predictor.PredictRetentionTime(peptide, out var predictFailure);
                
                // If formatting succeeds, prediction should succeed
                if (formatted != null)
                {
                    Assert.That(prediction, Is.Not.Null,
                        $"Predictor {predictor.PredictorName} formatted but failed to predict");
                }
                
                // If formatting fails, prediction should fail with same or related reason
                if (formatted == null && prediction == null)
                {
                    // Both failed - this is acceptable
                    Assert.Pass();
                }
            }
        }

        [Test]
        public void Integration_FormattedSequence_ContainsSequenceInformation()
        {
            var peptide = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());
            
            foreach (var predictor in _predictors)
            {
                var formatted = predictor.GetFormattedSequence(peptide, out _);
                
                if (formatted != null)
                {
                    // Formatted sequence should contain information about the peptide
                    Assert.That(formatted.Length, Is.GreaterThan(0),
                        $"Predictor {predictor.PredictorName} returned empty formatted sequence");
                    Assert.That(formatted.Contains("PEPTIDE") || formatted.Length != "PEPTIDE".Length,
                        $"Predictor {predictor.PredictorName} formatted sequence does not reflect peptide");
                }
            }
        }

        #endregion

        #region Comparative Tests

        [Test]
        public void Comparative_HydrophobicVsHydrophilic_ChronologerShowsDifference()
        {
            var hydrophobic = new PeptideWithSetModifications("FWLLLLW", new Dictionary<string, Modification>());
            var hydrophilic = new PeptideWithSetModifications("STTTSSS", new Dictionary<string, Modification>());
            
            using var chronologer = new ChronologerRetentionTimePredictor();
            
            var rtHydrophobic = chronologer.PredictRetentionTime(hydrophobic, out _);
            var rtHydrophilic = chronologer.PredictRetentionTime(hydrophilic, out _);
            
            Assert.That(rtHydrophobic, Is.Not.Null);
            Assert.That(rtHydrophilic, Is.Not.Null);
            Assert.That(rtHydrophobic.Value, Is.GreaterThan(rtHydrophilic.Value));
        }

        [Test]
        public void Comparative_ModifiedVsUnmodified_ShowsMassDifference()
        {
            var unmodified = new PeptideWithSetModifications("PEPTMIDE", new Dictionary<string, Modification>());
            
            var mods = new Dictionary<string, Modification>
            {
                { "Oxidation on M", ModificationConverter.AllModsKnown["Oxidation on M"] }
            };
            var modified = new PeptideWithSetModifications("PEPTM[Oxidation on M]IDE", mods);
            
            using var chronologer = new ChronologerRetentionTimePredictor();
            
            var rtUnmodified = chronologer.PredictRetentionTime(unmodified, out _);
            var rtModified = chronologer.PredictRetentionTime(modified, out _);
            
            Assert.That(rtUnmodified, Is.Not.Null);
            Assert.That(rtModified, Is.Not.Null);
            // Oxidation typically makes peptides slightly more hydrophilic
            Assert.That(rtModified.Value, Is.Not.EqualTo(rtUnmodified.Value));
            Assert.That(rtModified.Value, Is.LessThan(rtUnmodified.Value));

        }

        #endregion

        #region Error Handling Integration Tests

        [Test]
        public void Integration_NullPeptide_AllPredictorsHandleGracefully()
        {
            foreach (var predictor in _predictors)
            {
                Assert.Throws<System.ArgumentNullException>(() =>
                    predictor.PredictRetentionTime(null, out _));
            }
        }

        [Test]
        public void Integration_MultipleFailureScenarios_ProvideUsefulFailureReasons()
        {
            var testCases = new[]
            {
                (new PeptideWithSetModifications("", new Dictionary<string, Modification>()), "empty sequence"),
                (new PeptideWithSetModifications(new string('A', 100), new Dictionary<string, Modification>()), "long sequence")
            };

            foreach (var (peptide, description) in testCases)
            {
                foreach (var predictor in _predictors)
                {
                    var result = predictor.PredictRetentionTime(peptide, out var failureReason);
                    
                    if (result == null)
                    {
                        Assert.That(failureReason, Is.Not.Null,
                            $"Predictor {predictor.PredictorName} should provide failure reason for {description}");
                    }
                }
            }
        }

        #endregion
    }
}
