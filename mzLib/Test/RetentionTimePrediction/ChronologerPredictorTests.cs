using System;
using NUnit.Framework;
using Chromatography.RetentionTimePrediction;
using Chromatography.RetentionTimePrediction.Chronologer;
using Proteomics.ProteolyticDigestion;
using Omics.Modifications;
using System.Collections.Generic;
using Chromatography;
using Readers;
using System.Collections.Concurrent;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;

namespace Test.RetentionTimePrediction
{
    /// <summary>
    /// Tests for Chronologer retention time predictor
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class ChronologerPredictorTests
    {
        private ChronologerRetentionTimePredictor _predictor;

        [SetUp]
        public void Setup()
        {
            _predictor = new ChronologerRetentionTimePredictor(IncompatibleModHandlingMode.RemoveIncompatibleMods);
        }

        [TearDown]
        public void TearDown()
        {
            _predictor?.Dispose();
        }

        #region Basic Functionality Tests

        [Test]
        public void Constructor_DefaultMode_CreatesPredictor()
        {
            using var predictor = new ChronologerRetentionTimePredictor();
            
            Assert.That(predictor, Is.Not.Null);
            Assert.That(predictor.PredictorName, Is.EqualTo("Chronologer"));
            Assert.That(predictor.SeparationType, Is.EqualTo(SeparationType.HPLC));
        }

        [Test]
        public void Constructor_WithMode_CreatesPredictor()
        {
            using var predictor = new ChronologerRetentionTimePredictor(IncompatibleModHandlingMode.ThrowException);
            
            Assert.That(predictor, Is.Not.Null);
        }

        [Test]
        public void PredictorName_ReturnsCorrectName()
        {
            Assert.That(_predictor.PredictorName, Is.EqualTo("Chronologer"));
        }

        [Test]
        public void SeparationType_ReturnsHPLC()
        {
            Assert.That(_predictor.SeparationType, Is.EqualTo(SeparationType.HPLC));
        }

        #endregion

        #region Sequence Validation Tests

        [Test]
        public void PredictRetentionTime_ValidUnmodifiedPeptide_ReturnsValue()
        {
            var peptide = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());
            
            var result = _predictor.PredictRetentionTime(peptide, out var failureReason);
            
            Assert.That(result, Is.Not.Null);
            Assert.That(result.Value, Is.GreaterThan(0));
            Assert.That(failureReason, Is.Null);
        }

        [Test]
        public void PredictRetentionTime_EmptySequence_ReturnsNull()
        {
            var peptide = new PeptideWithSetModifications("", new Dictionary<string, Modification>());
            
            var result = _predictor.PredictRetentionTime(peptide, out var failureReason);
            
            Assert.That(result, Is.Null);
            Assert.That(failureReason, Is.Not.Null);
        }

        [Test]
        public void PredictRetentionTime_SequenceTooLong_ReturnsNull()
        {
            // Create a sequence longer than 50 amino acids
            var longSequence = new string('A', 51);
            var peptide = new PeptideWithSetModifications(longSequence, new Dictionary<string, Modification>());
            
            var result = _predictor.PredictRetentionTime(peptide, out var failureReason);
            
            Assert.That(result, Is.Null);
            Assert.That(failureReason, Is.EqualTo(RetentionTimeFailureReason.SequenceTooLong));
        }

        [Test]
        public void PredictRetentionTime_SequenceWithSelenocysteine_ReturnsNull()
        {
            var peptide = new PeptideWithSetModifications("PEPTUDE", new Dictionary<string, Modification>());
            
            var result = _predictor.PredictRetentionTime(peptide, out var failureReason);
            
            Assert.That(result, Is.Null);
            Assert.That(failureReason, Is.EqualTo(RetentionTimeFailureReason.InvalidAminoAcid));
        }

        [Test]
        public void PredictRetentionTime_MaxLengthPeptide_ReturnsValue()
        {
            // Test at the boundary: exactly 50 amino acids
            var maxLengthSequence = new string('A', 50);
            var peptide = new PeptideWithSetModifications(maxLengthSequence, new Dictionary<string, Modification>());
            
            var result = _predictor.PredictRetentionTime(peptide, out var failureReason);
            
            Assert.That(result, Is.Not.Null);
            Assert.That(failureReason, Is.Null);
        }

        #endregion

        #region Supported Modification Tests

        [Test]
        public void PredictRetentionTime_OxidationOnMethionine_ReturnsValue()
        {
            var mods = new Dictionary<string, Modification>
            {
                { "Oxidation on M", ModificationConverter.AllModsKnown["Oxidation on M"] }
            };
            var peptide = new PeptideWithSetModifications("PEPTM[Oxidation on M]IDE", mods);
            
            var result = _predictor.PredictRetentionTime(peptide, out var failureReason);
            
            Assert.That(result, Is.Not.Null);
            Assert.That(failureReason, Is.Null);
        }

        [Test]
        public void PredictRetentionTime_CarbamidomethylOnCysteine_ReturnsValue()
        {
            var mods = new Dictionary<string, Modification>
            {
                { "Carbamidomethyl on C", ModificationConverter.AllModsKnown["Carbamidomethyl on C"] }
            };
            var peptide = new PeptideWithSetModifications("PEPTC[Carbamidomethyl on C]IDE", mods);
            
            var result = _predictor.PredictRetentionTime(peptide, out var failureReason);
            
            Assert.That(result, Is.Not.Null);
            Assert.That(failureReason, Is.Null);
        }

        [Test]
        public void PredictRetentionTime_PhosphorylationOnSerine_ReturnsValue()
        {
            var mods = new Dictionary<string, Modification>
            {
                { "Phosphorylation on S", ModificationConverter.AllModsKnown["Phosphorylation on S"] }
            };
            var peptide = new PeptideWithSetModifications("PEPTIDES[Phosphorylation on S]", mods);
            
            var result = _predictor.PredictRetentionTime(peptide, out var failureReason);
            
            Assert.That(result, Is.Not.Null);
            Assert.That(failureReason, Is.Null);
        }

        [Test]
        public void PredictRetentionTime_AcetylationOnLysine_ReturnsValue()
        {
            var mods = new Dictionary<string, Modification>
            {
                { "Acetylation on K", ModificationConverter.AllModsKnown["Acetylation on K"] }
            };
            var peptide = new PeptideWithSetModifications("PEPTIDEK[Acetylation on K]", mods);
            
            var result = _predictor.PredictRetentionTime(peptide, out var failureReason);
            
            Assert.That(result, Is.Not.Null);
            Assert.That(failureReason, Is.Null);
        }

        [Test]
        public void PredictRetentionTime_MultipleModifications_ReturnsValue()
        {
            var mods = new Dictionary<string, Modification>
            {
                { "Oxidation on M", ModificationConverter.AllModsKnown["Oxidation on M"] },
                { "Carbamidomethyl on C", ModificationConverter.AllModsKnown["Carbamidomethyl on C"] }
            };
            var peptide = new PeptideWithSetModifications("PEPTM[Oxidation on M]IDEC[Carbamidomethyl on C]", mods);
            
            var result = _predictor.PredictRetentionTime(peptide, out var failureReason);
            
            Assert.That(result, Is.Not.Null);
            Assert.That(failureReason, Is.Null);
        }

        #endregion

        #region Incompatible Modification Handling Tests

        [Test]
        public void PredictRetentionTime_IncompatibleMod_RemoveMode_ReturnsValue()
        {
            using var predictor = new ChronologerRetentionTimePredictor(IncompatibleModHandlingMode.RemoveIncompatibleMods);
            
            // Create a peptide with an unsupported modification
            var mods = new Dictionary<string, Modification>
            {
                { "HexNAc on N", ModificationConverter.AllModsKnown["HexNAc on N"] }
            };
            var peptide = new PeptideWithSetModifications("PEPTN[HexNAc on N]IDE", mods);
            
            var result = predictor.PredictRetentionTime(peptide, out var failureReason);
            
            Assert.That(result, Is.Not.Null);
            Assert.That(failureReason, Is.Null);
        }

        [Test]
        public void PredictRetentionTime_IncompatibleMod_UsePrimaryMode_ReturnsValue()
        {
            using var predictor = new ChronologerRetentionTimePredictor(IncompatibleModHandlingMode.UsePrimarySequence);
            
            var mods = new Dictionary<string, Modification>
            {
                { "HexNAc on N", ModificationConverter.AllModsKnown["HexNAc on N"] }
            };
            var peptide = new PeptideWithSetModifications("PEPTN[HexNAc on N]IDE", mods);
            
            var result = predictor.PredictRetentionTime(peptide, out var failureReason);
            
            Assert.That(result, Is.Not.Null);
            Assert.That(failureReason, Is.Null);
        }

        [Test]
        public void PredictRetentionTime_IncompatibleMod_ReturnNullMode_ReturnsNull()
        {
            using var predictor = new ChronologerRetentionTimePredictor(IncompatibleModHandlingMode.ReturnNull);
            
            var mods = new Dictionary<string, Modification>
            {
                { "HexNAc on N", ModificationConverter.AllModsKnown["HexNAc on N"] }
            };
            var peptide = new PeptideWithSetModifications("PEPTN[HexNAc on N]IDE", mods);
            
            var result = predictor.PredictRetentionTime(peptide, out var failureReason);
            
            Assert.That(result, Is.Null);
            Assert.That(failureReason, Is.EqualTo(RetentionTimeFailureReason.IncompatibleModifications));
        }

        [Test]
        public void PredictRetentionTime_IncompatibleMod_ThrowMode_ThrowsException()
        {
            using var predictor = new ChronologerRetentionTimePredictor(IncompatibleModHandlingMode.ThrowException);
            
            var mods = new Dictionary<string, Modification>
            {
                { "HexNAc on N", ModificationConverter.AllModsKnown["HexNAc on N"] }
            };
            var peptide = new PeptideWithSetModifications("PEPTN[HexNAc on N]IDE", mods);
            
            Assert.Throws<IncompatibleModificationException>(() =>
                predictor.PredictRetentionTime(peptide, out _));
        }

        #endregion

        #region Formatted Sequence Tests

        [Test]
        public void GetFormattedSequence_UnmodifiedPeptide_ReturnsFormattedSequence()
        {
            var peptide = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());
            
            var formatted = _predictor.GetFormattedSequence(peptide, out var failureReason);
            
            Assert.That(formatted, Is.Not.Null);
            Assert.That(formatted, Does.StartWith("-")); // Free N-terminus
            Assert.That(formatted, Does.EndWith("_")); // Free C-terminus
            Assert.That(formatted, Is.EqualTo("-PEPTIDE_")); // Check full formatted sequence
            Assert.That(failureReason, Is.Null);
        }

        [Test]
        public void GetFormattedSequence_OxidizedMethionine_ReplacesWithCode()
        {
            var mods = new Dictionary<string, Modification>
            {
                { "Oxidation on M", ModificationConverter.AllModsKnown["Oxidation on M"] }
            };
            var peptide = new PeptideWithSetModifications("PEPTM[Oxidation on M]IDE", mods);
            
            var formatted = _predictor.GetFormattedSequence(peptide, out var failureReason);
            
            Assert.That(formatted, Is.Not.Null);
            Assert.That(formatted, Does.Contain("m")); // Lowercase 'm' for oxidized M
            Assert.That(formatted, Does.Not.Contain("[")); // No brackets in formatted sequence
            Assert.That(formatted, Does.Not.Contain("]"));
            Assert.That(formatted, Is.EqualTo("-PEPTmIDE_")); // Check full formatted sequence
            Assert.That(failureReason, Is.Null);
        }

        [Test]
        public void GetFormattedSequence_CarbamidomethylCysteine_ReplacesWithCode()
        {
            var mods = new Dictionary<string, Modification>
            {
                { "Carbamidomethyl on C", ModificationConverter.AllModsKnown["Carbamidomethyl on C"] }
            };
            var peptide = new PeptideWithSetModifications("PEPTC[Carbamidomethyl on C]IDE", mods);
            
            var formatted = _predictor.GetFormattedSequence(peptide, out var failureReason);
            
            Assert.That(formatted, Is.Not.Null);
            Assert.That(formatted, Does.Contain("c")); // Lowercase 'c' for carbamidomethyl C
            Assert.That(formatted, Does.Not.Contain("[")); // No brackets in formatted sequence
            Assert.That(formatted, Does.Not.Contain("]"));
            Assert.That(formatted, Is.EqualTo("-PEPTcIDE_")); // Check full formatted sequence
            Assert.That(failureReason, Is.Null);
        }

        [Test]
        public void GetFormattedSequence_NTerminalAcetylation_AddsCorrectToken()
        {
            var mods = new Dictionary<string, Modification>
            {
                { "Acetylation on X", ModificationConverter.AllModsKnown["Acetylation on X"] }
            };
            var peptide = new PeptideWithSetModifications("[Acetylation on X]-PEPTIDE", mods);
            
            var formatted = _predictor.GetFormattedSequence(peptide, out var failureReason);
            
            Assert.That(formatted, Is.Not.Null);
            Assert.That(formatted, Does.StartWith("^")); // Acetylated N-terminus token
            Assert.That(formatted, Does.Not.Contain("[")); // No brackets in formatted sequence
            Assert.That(formatted, Does.Not.Contain("]"));
            Assert.That(formatted, Is.EqualTo("^PEPTIDE_")); 
            Assert.That(failureReason, Is.Null);
        }

        #endregion

        #region Prediction Consistency Tests

        [Test]
        public void PredictRetentionTime_SamePeptideTwice_ReturnsSameValue()
        {
            var peptide = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());
            
            var result1 = _predictor.PredictRetentionTime(peptide, out _);
            var result2 = _predictor.PredictRetentionTime(peptide, out _);
            
            Assert.That(result1, Is.EqualTo(result2));
        }

        [Test]
        public void PredictRetentionTime_DifferentSequences_ReturnsDifferentValues()
        {
            var peptide1 = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());
            var peptide2 = new PeptideWithSetModifications("PRETEIN", new Dictionary<string, Modification>());
            
            var result1 = _predictor.PredictRetentionTime(peptide1, out var failureReason1);
            var result2 = _predictor.PredictRetentionTime(peptide2, out var failureReason2);

            Assert.That(result1, Is.Not.Null);
            Assert.That(failureReason1, Is.Null);
            Assert.That(result2, Is.Not.Null);
            Assert.That(failureReason2, Is.Null);
            Assert.That(result1, Is.Not.EqualTo(result2));
        }

        [Test]
        public void PredictRetentionTime_HydrophobicPeptide_ReturnsHigherValue()
        {
            // More hydrophobic peptide (contains F, W, L)
            var hydrophobic = new PeptideWithSetModifications("FWLLLLW", new Dictionary<string, Modification>());
            // Less hydrophobic peptide (contains more polar residues)
            var hydrophilic = new PeptideWithSetModifications("STTTSSS", new Dictionary<string, Modification>());
            
            var resultHydrophobic = _predictor.PredictRetentionTime(hydrophobic, out _);
            var resultHydrophilic = _predictor.PredictRetentionTime(hydrophilic, out _);
            
            Assert.That(resultHydrophobic, Is.Not.Null);
            Assert.That(resultHydrophilic, Is.Not.Null);
            Assert.That(resultHydrophobic.Value, Is.GreaterThan(resultHydrophilic.Value));
        }

        #endregion

        #region Edge Cases

        [Test]
        public void PredictRetentionTime_AllAlanine_ReturnsValue()
        {
            var peptide = new PeptideWithSetModifications("AAAAAAAAAA", new Dictionary<string, Modification>());
            
            var result = _predictor.PredictRetentionTime(peptide, out var failureReason);
            
            Assert.That(result, Is.Not.Null);
            Assert.That(failureReason, Is.Null);
        }

        [Test]
        public void PredictRetentionTime_AllDifferentAminoAcids_ReturnsValue()
        {
            var peptide = new PeptideWithSetModifications("ACDEFGHIKLMNPQRSTVWY", new Dictionary<string, Modification>());
            
            var result = _predictor.PredictRetentionTime(peptide, out var failureReason);
            
            Assert.That(result, Is.Not.Null);
            Assert.That(failureReason, Is.Null);
        }

        [Test]
        public void PredictRetentionTime_ShortPeptide_ReturnsFailure()
        {
            var peptide = new PeptideWithSetModifications("PEP", new Dictionary<string, Modification>());
            
            var result = _predictor.PredictRetentionTime(peptide, out var failureReason);
            
            Assert.That(result, Is.Null);
            Assert.That(failureReason, Is.Not.Null);
            Assert.That(failureReason, Is.EqualTo(RetentionTimeFailureReason.SequenceTooShort));
        }

        #endregion

        #region Disposal Tests

        [Test]
        public void Dispose_CalledMultipleTimes_DoesNotThrow()
        {
            using var predictor = new ChronologerRetentionTimePredictor();
            
            Assert.DoesNotThrow(() => predictor.Dispose());
            Assert.DoesNotThrow(() => predictor.Dispose()); // Second dispose should not throw
        }

        [Test]
        public void PredictRetentionTime_AfterDispose_ThrowsObjectDisposedException()
        {
            var predictor = new ChronologerRetentionTimePredictor();
            var peptide = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());
            
            predictor.Dispose();
            
            Assert.That(predictor.PredictRetentionTime(peptide, out var failureReason), Is.Null);

            var method = predictor.GetType().GetMethod("PredictCore" , System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);

            Exception e = null;
            try
            {
                method!.Invoke(predictor, new object[] { peptide, null });
            }
            catch (Exception ex)
            {
                e = ex;
            }

            if (e == null)
            {
                Assert.Fail("Expected an exception when calling PredictCore after disposal, but no exception was thrown.");
            }
            else
            {
                Assert.That(e.InnerException, Is.TypeOf<ObjectDisposedException>());
                Assert.That(e.InnerException.Message, Does.Contain("Chronologer"));
            }
        }

        #endregion

        #region Thread Safety Tests

        [Test]
        public void PredictRetentionTime_ConcurrentCalls_WithLocking_ReturnsConsistentResults()
        {
            using var predictor = new ChronologerRetentionTimePredictor();
            var peptide = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());
            
            var results = new ConcurrentBag<double?>();
            var tasks = Enumerable.Range(0, 10).Select(_ => Task.Run(() =>
            {
                var result = predictor.PredictRetentionTime(peptide, out var _);
                results.Add(result);
            })).ToArray();

            Task.WaitAll(tasks);

            Assert.That(results, Has.Count.EqualTo(10));
            Assert.That(results.All(r => r.HasValue), Is.True);
            Assert.That(results.Distinct().Count(), Is.EqualTo(1)); // All results should be identical
        }

        [Test]
        public void PredictRetentionTime_ConcurrentDifferentPeptides_WithLocking_Succeeds()
        {
            using var predictor = new ChronologerRetentionTimePredictor();
            var peptides = new[]
            {
                new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>()),
                new PeptideWithSetModifications("PEEPPTIDE", new Dictionary<string, Modification>()),
                new PeptideWithSetModifications("TIDEPEPTIDE", new Dictionary<string, Modification>())
            };

            var results = new ConcurrentBag<(string sequence, double? rt)>();
            var tasks = Enumerable.Range(0, 30).Select(i => Task.Run(() =>
            {
                var peptide = peptides[i % 3];
                var result = predictor.PredictRetentionTime(peptide, out _);
                results.Add((peptide.BaseSequence, result));
            })).ToArray();

            Task.WaitAll(tasks);

            Assert.That(results, Has.Count.EqualTo(30));
            Assert.That(results.All(r => r.rt.HasValue), Is.True);
            
            // Verify same peptide always gives same result
            var groupedResults = results.GroupBy(r => r.sequence);
            foreach (var group in groupedResults)
            {
                var distinctRts = group.Select(r => r.rt).Distinct().Count();
                Assert.That(distinctRts, Is.EqualTo(1), $"Peptide {group.Key} had inconsistent results");
            }
        }

        #endregion
    }
}
