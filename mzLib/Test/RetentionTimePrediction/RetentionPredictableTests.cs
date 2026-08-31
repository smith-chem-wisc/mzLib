using NUnit.Framework;
using Chromatography.RetentionTimePrediction;
using Proteomics.ProteolyticDigestion;
using Omics.Modifications;
using System.Collections.Generic;
using Readers;

namespace Test.RetentionTimePrediction
{
    /// <summary>
    /// Tests for IRetentionPredictable interface implementation in PeptideWithSetModifications
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class RetentionPredictableTests
    {
        #region Basic Interface Implementation Tests

        [Test]
        public void BaseSequence_UnmodifiedPeptide_ReturnsCorrectSequence()
        {
            IRetentionPredictable peptide = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());
            
            Assert.That(peptide.BaseSequence, Is.EqualTo("PEPTIDE"));
        }

        [Test]
        public void BaseSequence_ModifiedPeptide_ReturnsSequenceWithoutModifications()
        {
            var mods = new Dictionary<string, Modification>
            {
                { "Oxidation on M", ModificationConverter.AllModsKnown["Oxidation on M"] }
            };
            IRetentionPredictable peptide = new PeptideWithSetModifications("PEPTM[Oxidation on M]IDE", mods);
            
            Assert.That(peptide.BaseSequence, Is.EqualTo("PEPTMIDE"));
        }

        [Test]
        public void FullSequence_UnmodifiedPeptide_ReturnsBaseSequence()
        {
            IRetentionPredictable peptide = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());
            
            Assert.That(peptide.FullSequence, Is.EqualTo("PEPTIDE"));
        }

        [Test]
        public void FullSequence_ModifiedPeptide_IncludesModifications()
        {
            var mods = new Dictionary<string, Modification>
            {
                { "Oxidation on M", ModificationConverter.AllModsKnown["Oxidation on M"] }
            };
            IRetentionPredictable peptide = new PeptideWithSetModifications("PEPTM[Oxidation on M]IDE", mods);
            
            Assert.That(peptide.FullSequence, Does.Contain("["));
            Assert.That(peptide.FullSequence, Does.Contain("]"));
            Assert.That(peptide.FullSequence, Does.Contain("PEPTM"));
        }

        [Test]
        public void MonoisotopicMass_ValidPeptide_ReturnsPositiveValue()
        {
            IRetentionPredictable peptide = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());
            
            Assert.That(peptide.MonoisotopicMass, Is.GreaterThan(0));
        }

        [Test]
        public void MonoisotopicMass_ModifiedPeptide_IncludesModificationMass()
        {
            var mods = new Dictionary<string, Modification>
            {
                { "Oxidation on M", ModificationConverter.AllModsKnown["Oxidation on M"] }
            };
            
            IRetentionPredictable unmodified = new PeptideWithSetModifications("PEPTMIDE", new Dictionary<string, Modification>());
            IRetentionPredictable modified = new PeptideWithSetModifications("PEPTM[Oxidation on M]IDE", mods);
            
            Assert.That(modified.MonoisotopicMass, Is.GreaterThan(unmodified.MonoisotopicMass));
        }

        #endregion

        #region FullSequenceWithMassShifts Tests

        [Test]
        public void FullSequenceWithMassShifts_UnmodifiedPeptide_ReturnsBaseSequence()
        {
            IRetentionPredictable peptide = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());
            
            Assert.That(peptide.FullSequenceWithMassShifts, Is.EqualTo("PEPTIDE"));
        }

        [Test]
        public void FullSequenceWithMassShifts_OxidizedMethionine_IncludesMassShift()
        {
            var mods = new Dictionary<string, Modification>
            {
                { "Oxidation on M", ModificationConverter.AllModsKnown["Oxidation on M"] }
            };
            IRetentionPredictable peptide = new PeptideWithSetModifications("PEPTM[Oxidation on M]IDE", mods);
            
            var withMassShifts = peptide.FullSequenceWithMassShifts;
            
            Assert.That(withMassShifts, Does.Contain("M[+15.99"));
        }

        [Test]
        public void FullSequenceWithMassShifts_Carbamidomethyl_IncludesMassShift()
        {
            var mods = new Dictionary<string, Modification>
            {
                { "Carbamidomethyl on C", ModificationConverter.AllModsKnown["Carbamidomethyl on C"] }
            };
            IRetentionPredictable peptide = new PeptideWithSetModifications("PEPTC[Carbamidomethyl on C]IDE", mods);
            
            var withMassShifts = peptide.FullSequenceWithMassShifts;
            
            Assert.That(withMassShifts, Does.Contain("PEPTC[+57.02"));
            Assert.That(withMassShifts, Does.Contain("]IDE"));
        }

        [Test]
        public void FullSequenceWithMassShifts_MultipleModifications_IncludesAllMassShifts()
        {
            var mods = new Dictionary<string, Modification>
            {
                { "Oxidation on M", ModificationConverter.AllModsKnown["Oxidation on M"] },
                { "Carbamidomethyl on C", ModificationConverter.AllModsKnown["Carbamidomethyl on C"] }
            };
            IRetentionPredictable peptide = new PeptideWithSetModifications("PEPTM[Oxidation on M]IDEC[Carbamidomethyl on C]", mods);
            
            var withMassShifts = peptide.FullSequenceWithMassShifts;
            
            Assert.That(withMassShifts, Does.Contain("PEPTM[+15.99"));
            Assert.That(withMassShifts, Does.Contain("]IDEC[+57.02"));
            Assert.That(withMassShifts, Does.EndWith("]"));
        }

        [Test]
        public void FullSequenceWithMassShifts_NegativeMassShift_FormatCorrectly()
        {
            var mods = new Dictionary<string, Modification>
            {
                { "Glu to PyroGlu on Q", ModificationConverter.AllModsKnown["Glu to PyroGlu on Q"] }
            };
            IRetentionPredictable peptide = new PeptideWithSetModifications("Q[Glu to PyroGlu on Q]PEPTIDE", mods);
            
            var withMassShifts = peptide.FullSequenceWithMassShifts;
            
            Assert.That(withMassShifts, Does.Contain("Q[-17.02"));
            Assert.That(withMassShifts, Does.Contain("]PEPTIDE"));
        }

        [Test]
        public void FullSequenceWithMassShifts_NTerminalModification_IncludesMassShift()
        {
            var mods = new Dictionary<string, Modification>
            {
                { "Acetylation on X", ModificationConverter.AllModsKnown["Acetylation on X"] }
            };
            IRetentionPredictable peptide = new PeptideWithSetModifications("[Acetylation on X]PEPTIDE", mods);
            
            var withMassShifts = peptide.FullSequenceWithMassShifts;
            
            Assert.That(withMassShifts, Does.Contain("[+42.01"));
            Assert.That(withMassShifts, Does.Contain("]PEPTIDE"));
        }

        [Test]
        public void FullSequenceWithMassShifts_CTerminalModification_IncludesMassShift()
        {
            var mods = new Dictionary<string, Modification>
            {
                { "Amidation on X", ModificationConverter.AllModsKnown["Amidation on X"] }
            };
            IRetentionPredictable peptide = new PeptideWithSetModifications("PEPTIDE-[Amidation on X]", mods);
            
            var withMassShifts = peptide.FullSequenceWithMassShifts;
            
            // Should include the mass shift
            Assert.That(withMassShifts, Does.Contain("["));
        }

        #endregion

        #region Consistency Tests

        [Test]
        public void BaseSequence_Length_MatchesMonoisotopicMassCalculation()
        {
            IRetentionPredictable peptide = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());
            
            Assert.That(peptide.BaseSequence.Length, Is.EqualTo(7));
            Assert.That(peptide.MonoisotopicMass, Is.GreaterThan(700)); // Approximate mass check
        }

        [Test]
        public void FullSequence_ContainsBaseSequence()
        {
            var mods = new Dictionary<string, Modification>
            {
                { "Oxidation on M", ModificationConverter.AllModsKnown["Oxidation on M"] }
            };
            IRetentionPredictable peptide = new PeptideWithSetModifications("PEPTM[Oxidation on M]IDE", mods);
            
            // Full sequence should contain all characters from base sequence
            foreach (char c in peptide.BaseSequence)
            {
                Assert.That(peptide.FullSequence, Does.Contain(c.ToString()));
            }
        }

        [Test]
        public void FullSequenceWithMassShifts_ContainsBaseSequenceCharacters()
        {
            var mods = new Dictionary<string, Modification>
            {
                { "Oxidation on M", ModificationConverter.AllModsKnown["Oxidation on M"] }
            };
            IRetentionPredictable peptide = new PeptideWithSetModifications("PEPTM[Oxidation on M]IDE", mods);
            
            // Remove brackets and their content to check base sequence
            var withoutBrackets = System.Text.RegularExpressions.Regex.Replace(
                peptide.FullSequenceWithMassShifts, @"\[.*?\]", "");
            
            Assert.That(withoutBrackets, Is.EqualTo(peptide.BaseSequence));
        }

        #endregion

        #region Edge Cases

        [Test]
        public void Interface_ShortPeptide_ImplementsCorrectly()
        {
            IRetentionPredictable peptide = new PeptideWithSetModifications("PEP", new Dictionary<string, Modification>());
            
            Assert.That(peptide.BaseSequence, Is.EqualTo("PEP"));
            Assert.That(peptide.FullSequence, Is.EqualTo("PEP"));
            Assert.That(peptide.FullSequenceWithMassShifts, Is.EqualTo("PEP"));
            Assert.That(peptide.MonoisotopicMass, Is.GreaterThan(0));
        }

        [Test]
        public void Interface_LongPeptide_ImplementsCorrectly()
        {
            var longSequence = "PEPTIDEPEPTIDEPEPTIDEPEPTIDE";
            IRetentionPredictable peptide = new PeptideWithSetModifications(longSequence, new Dictionary<string, Modification>());
            
            Assert.That(peptide.BaseSequence, Is.EqualTo(longSequence));
            Assert.That(peptide.FullSequence, Is.EqualTo(longSequence));
            Assert.That(peptide.MonoisotopicMass, Is.GreaterThan(0));
        }

        [Test]
        public void Interface_AllCanonicalAminoAcids_ImplementsCorrectly()
        {
            var allAA = "ACDEFGHIKLMNPQRSTVWY";
            IRetentionPredictable peptide = new PeptideWithSetModifications(allAA, new Dictionary<string, Modification>());
            
            Assert.That(peptide.BaseSequence, Is.EqualTo(allAA));
            Assert.That(peptide.MonoisotopicMass, Is.GreaterThan(0));
        }

        #endregion

        #region Mass Accuracy Tests

        [Test]
        public void MonoisotopicMass_Consistency_SameSequenceAlwaysSameMass()
        {
            IRetentionPredictable peptide1 = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());
            IRetentionPredictable peptide2 = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());
            
            Assert.That(peptide1.MonoisotopicMass, Is.EqualTo(peptide2.MonoisotopicMass));
        }

        [Test]
        public void MonoisotopicMass_ModificationAddsCorrectMass()
        {
            var unmodified = new PeptideWithSetModifications("PEPTMIDE", new Dictionary<string, Modification>());
            
            var mods = new Dictionary<string, Modification>
            {
                { "Oxidation on M", ModificationConverter.AllModsKnown["Oxidation on M"] }
            };
            var modified = new PeptideWithSetModifications("PEPTM[Oxidation on M]IDE", mods);
            
            var oxidationMass = ModificationConverter.AllModsKnown["Oxidation on M"].MonoisotopicMass;
            var expectedMass = unmodified.MonoisotopicMass + oxidationMass.Value;
            
            Assert.That(modified.MonoisotopicMass, Is.EqualTo(expectedMass).Within(0.01));
        }

        #endregion
    }
}
