using NUnit.Framework;
using Chromatography.RetentionTimePrediction;
using Chromatography.RetentionTimePrediction.Chronologer;
using Chromatography.RetentionTimePrediction.Util;
using Proteomics.ProteolyticDigestion;
using Omics.Modifications;
using Proteomics;
using System.Collections.Generic;
using Readers;

namespace Test.RetentionTimePrediction
{
    /// <summary>
    /// Tests for sequence formatting and modification handling across different predictors
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class SequenceFormattingTests
    {
        #region Chronologer Formatting Tests

        [Test]
        public void ChronologerFormatting_UnmodifiedPeptide_AddsTermini()
        {
            using var predictor = new ChronologerRetentionTimePredictor();
            var peptide = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());
            
            var formatted = predictor.GetFormattedSequence(peptide, out var failureReason);
            
            Assert.That(formatted, Is.Not.Null);
            Assert.That(formatted, Does.StartWith("-")); // Free N-terminus
            Assert.That(formatted, Does.EndWith("_")); // Free C-terminus
            Assert.That(formatted, Does.Contain("PEPTIDE"));
            Assert.That(formatted, Is.EqualTo("-PEPTIDE_"));
            Assert.That(failureReason, Is.Null);
        }

        [Test]
        public void ChronologerFormatting_OxidizedMethionine_ConvertsTolowercase()
        {
            using var predictor = new ChronologerRetentionTimePredictor();
            var mods = new Dictionary<string, Modification>
            {
                { "Oxidation on M", ModificationConverter.AllModsKnown["Oxidation on M"] }
            };
            var peptide = new PeptideWithSetModifications("PEPTM[Oxidation on M]IDE", mods);
            
            var formatted = predictor.GetFormattedSequence(peptide, out var failureReason);
            
            Assert.That(formatted, Is.Not.Null);
            Assert.That(formatted, Does.Contain("m")); // Lowercase m for oxidation
            Assert.That(formatted, Does.Not.Contain("M")); // No uppercase M
            Assert.That(formatted, Does.Not.Contain("[")); // No brackets
            Assert.That(formatted, Is.EqualTo("-PEPTmIDE_"));
            Assert.That(failureReason, Is.Null);
        }

        [Test]
        public void ChronologerFormatting_CarbamidomethylCysteine_ConvertsToCode()
        {
            using var predictor = new ChronologerRetentionTimePredictor();
            var mods = new Dictionary<string, Modification>
            {
                { "Carbamidomethyl on C", ModificationConverter.AllModsKnown["Carbamidomethyl on C"] }
            };
            var peptide = new PeptideWithSetModifications("PEPTC[Carbamidomethyl on C]IDE", mods);
            
            var formatted = predictor.GetFormattedSequence(peptide, out var failureReason);
            
            Assert.That(formatted, Is.Not.Null);
            Assert.That(formatted, Does.Contain("c")); // Lowercase c for carbamidomethyl
            Assert.That(formatted, Does.Not.Contain("C")); // No uppercase C
            Assert.That(formatted, Does.Not.Contain("[")); // No brackets
            Assert.That(formatted, Is.EqualTo("-PEPTcIDE_"));
            Assert.That(failureReason, Is.Null);
        }

        [Test]
        public void ChronologerFormatting_PhosphorylationOnSerine_ConvertsToCode()
        {
            using var predictor = new ChronologerRetentionTimePredictor();
            var mods = new Dictionary<string, Modification>
            {
                { "Phospho on S", ModificationConverter.AllModsKnown["Phospho on S"] }
            };
            var peptide = new PeptideWithSetModifications("PEPTIDES[Phospho on S]", mods);
            
            var formatted = predictor.GetFormattedSequence(peptide, out var failureReason);
            
            Assert.That(formatted, Is.Not.Null);
            Assert.That(formatted, Does.Contain("s")); // Lowercase s for phosphorylation
            Assert.That(formatted, Is.EqualTo("-PEPTIDEs_"));
            Assert.That(failureReason, Is.Null);
        }

        [Test]
        public void ChronologerFormatting_PhosphorylationOnThreonine_ConvertsToCode()
        {
            using var predictor = new ChronologerRetentionTimePredictor();
            var mods = new Dictionary<string, Modification>
            {
                { "Phospho on T", ModificationConverter.AllModsKnown["Phospho on T"] }
            };
            var peptide = new PeptideWithSetModifications("PEPTIDET[Phospho on T]", mods);
            
            var formatted = predictor.GetFormattedSequence(peptide, out var failureReason);
            
            Assert.That(formatted, Is.Not.Null);
            Assert.That(formatted, Does.Contain("t")); // Lowercase t for phosphorylation
            Assert.That(formatted, Is.EqualTo("-PEPTIDEt_"));
            Assert.That(failureReason, Is.Null);
        }

        [Test]
        public void ChronologerFormatting_PhosphorylationOnTyrosine_ConvertsToCode()
        {
            using var predictor = new ChronologerRetentionTimePredictor();
            var mods = new Dictionary<string, Modification>
            {
                { "Phospho on Y", ModificationConverter.AllModsKnown["Phospho on Y"] }
            };
            var peptide = new PeptideWithSetModifications("PEPTIDEY[Phospho on Y]", mods);
            
            var formatted = predictor.GetFormattedSequence(peptide, out var failureReason);
            
            Assert.That(formatted, Is.Not.Null);
            Assert.That(formatted, Does.Contain("y")); // Lowercase y for phosphorylation
            Assert.That(formatted, Is.EqualTo("-PEPTIDEy_"));
            Assert.That(failureReason, Is.Null);
        }

        [Test]
        public void ChronologerFormatting_AcetylationOnLysine_ConvertsToCode()
        {
            using var predictor = new ChronologerRetentionTimePredictor();
            var mods = new Dictionary<string, Modification>
            {
                { "Acetylation on K", ModificationConverter.AllModsKnown["Acetylation on K"] }
            };
            var peptide = new PeptideWithSetModifications("PEPTIDEK[Acetylation on K]", mods);
            
            var formatted = predictor.GetFormattedSequence(peptide, out var failureReason);
            
            Assert.That(formatted, Is.Not.Null);
            Assert.That(formatted, Does.Contain("a")); // 'a' for acetylation
            Assert.That(formatted, Is.EqualTo("-PEPTIDEa_"));
            Assert.That(failureReason, Is.Null);
        }

        [Test]
        public void ChronologerFormatting_MultipleModifications_AllConverted()
        {
            using var predictor = new ChronologerRetentionTimePredictor();
            var mods = new Dictionary<string, Modification>
            {
                { "Oxidation on M", ModificationConverter.AllModsKnown["Oxidation on M"] },
                { "Carbamidomethyl on C", ModificationConverter.AllModsKnown["Carbamidomethyl on C"] },
                { "Phospho on S", ModificationConverter.AllModsKnown["Phospho on S"] }
            };
            var peptide = new PeptideWithSetModifications("M[Oxidation on M]C[Carbamidomethyl on C]S[Phospho on S]", mods);
            
            var formatted = predictor.GetFormattedSequence(peptide, out var failureReason);
            
            Assert.That(formatted, Is.Not.Null);
            Assert.That(formatted, Does.Contain("m")); // oxidized M
            Assert.That(formatted, Does.Contain("c")); // carbamidomethyl C
            Assert.That(formatted, Does.Contain("s")); // phospho S
            Assert.That(formatted, Does.Not.Contain("[")); // No brackets remain
            Assert.That(formatted, Is.EqualTo("-mcs_"));
            Assert.That(failureReason, Is.Null);
        }

        [Test]
        public void ChronologerFormatting_UnsupportedMod_RemoveMode_RemovesIt()
        {
            using var predictor = new ChronologerRetentionTimePredictor(IncompatibleModHandlingMode.RemoveIncompatibleMods);
            var mods = new Dictionary<string, Modification>
            {
                { "HexNAc on N", ModificationConverter.AllModsKnown["HexNAc on N"] }
            };
            var peptide = new PeptideWithSetModifications("PEPTN[HexNAc on N]IDE", mods);
            
            var formatted = predictor.GetFormattedSequence(peptide, out var failureReason);
            
            Assert.That(formatted, Is.Not.Null);
            Assert.That(formatted, Does.Not.Contain("[")); // Brackets removed
            Assert.That(formatted, Does.Contain("N")); // Base amino acid still there
            Assert.That(formatted, Is.EqualTo("-PEPTNIDE_"));
            Assert.That(failureReason, Is.Null);
        }

        [Test]
        public void ChronologerFormatting_UnsupportedMod_UsePrimaryMode_ReturnsBaseSequence()
        {
            using var predictor = new ChronologerRetentionTimePredictor(IncompatibleModHandlingMode.UsePrimarySequence);
            var mods = new Dictionary<string, Modification>
            {
                { "HexNAc on N", ModificationConverter.AllModsKnown["HexNAc on N"] }
            };
            var peptide = new PeptideWithSetModifications("PEPTN[HexNAc on N]IDE", mods);
            
            var formatted = predictor.GetFormattedSequence(peptide, out var failureReason);
            
            Assert.That(formatted, Is.Not.Null);
            Assert.That(formatted, Is.EqualTo("-PEPTNIDE_")); // Just base sequence
            Assert.That(failureReason, Is.Null);
        }

        [Test]
        public void ChronologerFormatting_UnsupportedMod_ReturnNullMode_ReturnsNull()
        {
            using var predictor = new ChronologerRetentionTimePredictor(IncompatibleModHandlingMode.ReturnNull);
            var mods = new Dictionary<string, Modification>
            {
                { "HexNAc on N", ModificationConverter.AllModsKnown["HexNAc on N"] }
            };
            var peptide = new PeptideWithSetModifications("PEPTN[HexNAc on N]IDE", mods);
            
            var formatted = predictor.GetFormattedSequence(peptide, out var failureReason);
            
            Assert.That(formatted, Is.Null);
            Assert.That(failureReason, Is.EqualTo(RetentionTimeFailureReason.IncompatibleModifications));
        }

        [Test]
        public void ChronologerFormatting_UnsupportedMod_ThrowMode_ThrowsException()
        {
            using var predictor = new ChronologerRetentionTimePredictor(IncompatibleModHandlingMode.ThrowException);
            var mods = new Dictionary<string, Modification>
            {
                { "HexNAc on N", ModificationConverter.AllModsKnown["HexNAc on N"] }
            };
            var peptide = new PeptideWithSetModifications("PEPTN[HexNAc on N]IDE", mods);
            
            Assert.Throws<IncompatibleModificationException>(() =>
                predictor.GetFormattedSequence(peptide, out _));
        }

        #endregion

        #region Mass Shift Formatting Tests

        [Test]
        public void MassShiftFormatting_PositiveMassShift_FormatsCorrectly()
        {
            var peptide = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());
            var massShiftSequence = peptide.FullSequenceWithMassShifts;
            
            // For unmodified peptide, should just be the sequence
            Assert.That(massShiftSequence, Is.EqualTo("PEPTIDE"));
        }

        [Test]
        public void MassShiftFormatting_OxidationOnM_IncludesPositiveMass()
        {
            var mods = new Dictionary<string, Modification>
            {
                { "Oxidation on M", ModificationConverter.AllModsKnown["Oxidation on M"] }
            };
            var peptide = new PeptideWithSetModifications("PEPTM[Oxidation on M]IDE", mods);
            var massShiftSequence = peptide.FullSequenceWithMassShifts;
            
            Assert.That(massShiftSequence, Does.Contain("[+"));
            Assert.That(massShiftSequence, Does.Contain("15.99")); // Oxidation mass
            Assert.That(massShiftSequence, Does.Contain("M"));
        }

        [Test]
        public void MassShiftFormatting_NegativeMassShift_FormatsWithMinus()
        {
            var mods = new Dictionary<string, Modification>
            {
                { "Glu to PyroGlu on Q", ModificationConverter.AllModsKnown["Glu to PyroGlu on Q"] }
            };
            var peptide = new PeptideWithSetModifications("Q[Glu to PyroGlu on Q]PEPTIDE", mods);
            var massShiftSequence = peptide.FullSequenceWithMassShifts;
            
            Assert.That(massShiftSequence, Does.Contain("[-"));
            Assert.That(massShiftSequence, Does.Contain("Q"));
        }

        [Test]
        public void MassShiftFormatting_PreservesSequenceOrder()
        {
            var mods = new Dictionary<string, Modification>
            {
                { "Oxidation on M", ModificationConverter.AllModsKnown["Oxidation on M"] }
            };
            var peptide = new PeptideWithSetModifications("PEPTM[Oxidation on M]IDE", mods);
            var massShiftSequence = peptide.FullSequenceWithMassShifts;
            
            // Remove mass shift annotations to get base sequence
            var withoutMods = System.Text.RegularExpressions.Regex.Replace(massShiftSequence, @"\[.*?\]", "");
            Assert.That(withoutMods, Is.EqualTo(peptide.BaseSequence));
        }

        #endregion

        #region Edge Cases

        [Test]
        public void Formatting_OnlyModifications_NoBaseSequence_HandlesGracefully()
        {
            // This should not be possible with normal construction, but test edge case
            using var predictor = new ChronologerRetentionTimePredictor();
            var peptide = new PeptideWithSetModifications("A", new Dictionary<string, Modification>());
            
            var formatted = predictor.GetFormattedSequence(peptide, out var failureReason);
            
            // Single amino acid peptides should still be formatted
            Assert.That(formatted, Is.Not.Null);
        }

        [Test]
        public void Formatting_AllModifiedResidues_AllConverted()
        {
            using var predictor = new ChronologerRetentionTimePredictor();
            var mods = new Dictionary<string, Modification>
            {
                { "Oxidation on M", ModificationConverter.AllModsKnown["Oxidation on M"] },
                { "Carbamidomethyl on C", ModificationConverter.AllModsKnown["Carbamidomethyl on C"] }
            };
            var peptide = new PeptideWithSetModifications("M[Oxidation on M]C[Carbamidomethyl on C]M[Oxidation on M]", mods);
            
            var formatted = predictor.GetFormattedSequence(peptide, out var failureReason);
            
            Assert.That(formatted, Is.Not.Null);
            // All modifications should be converted
            Assert.That(formatted, Does.Not.Contain("["));
            Assert.That(formatted, Does.Not.Contain("]"));
            // Should contain modified residue codes
            Assert.That(formatted, Does.Contain("m")); // oxidized M
            Assert.That(formatted, Does.Contain("c")); // carbamidomethyl C
            Assert.That(formatted, Is.EqualTo("-mcm_"));
            Assert.That(failureReason, Is.Null);
        }

        #endregion

        #region Consistency Tests

        [Test]
        public void Formatting_SameSequenceTwice_ReturnsSameResult()
        {
            using var predictor = new ChronologerRetentionTimePredictor();
            var peptide = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());
            
            var formatted1 = predictor.GetFormattedSequence(peptide, out _);
            var formatted2 = predictor.GetFormattedSequence(peptide, out _);
            
            Assert.That(formatted1, Is.EqualTo(formatted2));
        }

        [Test]
        public void Formatting_DifferentModHandlingModes_GiveDifferentResults()
        {
            var mods = new Dictionary<string, Modification>
            {
                { "HexNAc on N", ModificationConverter.AllModsKnown["HexNAc on N"] },
                {"Phospho on T", ModificationConverter.AllModsKnown["Phospho on T"] }
            };
            var peptide = new PeptideWithSetModifications("PEPT[Phospho on T]N[HexNAc on N]IDE", mods);
            
            using var removeMode = new ChronologerRetentionTimePredictor(IncompatibleModHandlingMode.RemoveIncompatibleMods);
            using var primaryMode = new ChronologerRetentionTimePredictor(IncompatibleModHandlingMode.UsePrimarySequence);
            
            var removeFormatted = removeMode.GetFormattedSequence(peptide, out _);
            var primaryFormatted = primaryMode.GetFormattedSequence(peptide, out _);
            
            Assert.That(removeFormatted, Is.Not.EqualTo(primaryFormatted));
        }

        #endregion
    }
}
