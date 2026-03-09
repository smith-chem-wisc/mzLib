using System.Linq;
using NUnit.Framework;
using Omics.Modifications;
using Omics.SequenceConversion;

namespace Test.Omics.SequenceConversion
{
    /// <summary>
    /// Tests for parsing sequences from various formats into the canonical intermediate representation.
    /// </summary>
    [TestFixture]
    public class ParsingTests
    {
        private MzLibSequenceParser _mzLibParser;
        private MassShiftSequenceParser _massShiftParser;

        [SetUp]
        public void Setup()
        {
            _mzLibParser = new MzLibSequenceParser();
            _massShiftParser = new MassShiftSequenceParser();
        }

        [Test]
        [TestCaseSource(typeof(GroundTruthTestData), nameof(GroundTruthTestData.CoreTestCases))]
        public void MzLibParser_CoreTestCases_ParsesCorrectly(GroundTruthTestData.TestCase testCase)
        {
            // Act
            var result = _mzLibParser.Parse(testCase.MzLibFormat);

            // Assert
            Assert.That(result, Is.Not.Null, $"Failed to parse: {testCase.Description}");
            var canonical = result.Value;
            
            Assert.That(canonical.BaseSequence, Is.EqualTo(testCase.ExpectedBaseSequence));
            Assert.That(canonical.NTerminalModification.HasValue, Is.EqualTo(testCase.HasNTermMod));
            Assert.That(canonical.CTerminalModification.HasValue, Is.EqualTo(testCase.HasCTermMod));
            Assert.That(canonical.ResidueModifications.Count(), Is.EqualTo(testCase.ExpectedResidueMods));
        }

        [Test]
        [TestCaseSource(typeof(GroundTruthTestData), nameof(GroundTruthTestData.EdgeCases))]
        public void MzLibParser_EdgeCases_ParsesCorrectly(GroundTruthTestData.TestCase testCase)
        {
            // Act
            var result = _mzLibParser.Parse(testCase.MzLibFormat);

            // Assert
            Assert.That(result, Is.Not.Null, $"Failed to parse: {testCase.Description}");
            Assert.That(result.Value.BaseSequence, Is.EqualTo(testCase.ExpectedBaseSequence));
        }

        [Test]
        [TestCaseSource(typeof(GroundTruthTestData), nameof(GroundTruthTestData.CoreTestCases))]
        public void MassShiftParser_CoreTestCases_ParsesCorrectly(GroundTruthTestData.TestCase testCase)
        {
            // Act
            var result = _massShiftParser.Parse(testCase.MassShiftFormat);

            // Assert
            Assert.That(result, Is.Not.Null, $"Failed to parse: {testCase.Description}");
            var canonical = result.Value;
            
            Assert.That(canonical.BaseSequence, Is.EqualTo(testCase.ExpectedBaseSequence));
            Assert.That(canonical.ResidueModifications.Count(), Is.EqualTo(testCase.ExpectedResidueMods));
        }

        [Test]
        public void MzLibParser_EmptyString_ReturnsNull()
        {
            var result = _mzLibParser.Parse("", null, SequenceConversionHandlingMode.ReturnNull);
            Assert.That(result, Is.Null);
        }

        [Test]
        public void MzLibParser_ResidueModifications_HaveCorrectPositions()
        {
            // Arrange
            var sequence = "PEPS[Common Biological:Phosphorylation on S]TM[Common Variable:Oxidation on M]IDE";
            var phospho = Mods.GetModification("Phosphorylation on S");
            var oxidation = Mods.GetModification("Oxidation on M");

            // Act
            var result = _mzLibParser.Parse(sequence);
            
            // Assert
            Assert.That(result, Is.Not.Null);
            var canonical = result.Value;
            
            // Check position 3 (S)
            var sMod = canonical.GetModificationAt(3);
            Assert.That(sMod, Is.Not.Null);
            Assert.That(sMod.Value.IsResolved, Is.False);
            Assert.That(sMod.Value.ResidueIndex, Is.EqualTo(3));
            Assert.That(sMod.Value.TargetResidue, Is.EqualTo('S'));

            // Check position 5 (M)
            var mMod = canonical.GetModificationAt(5);
            Assert.That(mMod, Is.Not.Null);
            Assert.That(mMod.Value.IsResolved, Is.False);
            Assert.That(mMod.Value.ResidueIndex, Is.EqualTo(5));
            Assert.That(mMod.Value.TargetResidue, Is.EqualTo('M'));
        }

        [Test]
        public void MzLibParser_ResidueModifications_ResolvedProperly()
        {
            // Arrange
            var lookup = new MzLibModificationLookup();
            var sequence = "PEPS[Common Biological:Phosphorylation on S]TM[Common Variable:Oxidation on M]IDE";
            var phospho = Mods.GetModification("Phosphorylation on S");
            var oxidation = Mods.GetModification("Oxidation on M", ModificationNamingConvention.MetaMorpheus_Protein);

            // Act
            var result = _mzLibParser.Parse(sequence);

            // Assert
            Assert.That(result, Is.Not.Null);
            var canonical = result.Value;

            // Check position 3 (S)
            var sMod = canonical.GetModificationAt(3);
            Assert.That(sMod, Is.Not.Null);
            Assert.That(sMod.Value.ResidueIndex, Is.EqualTo(3));
            Assert.That(sMod.Value.TargetResidue, Is.EqualTo('S'));
            Assert.That(sMod.Value.IsResolved, Is.False);
            sMod = lookup.TryResolve(sMod.Value);

            Assert.That(sMod.Value.MzLibModification, Is.EqualTo(phospho));
            Assert.That(sMod.Value.ChemicalFormula, Is.EqualTo(phospho.ChemicalFormula));
            Assert.That(sMod.Value.HasMass);
            Assert.That(sMod.Value.EffectiveMass, Is.EqualTo(phospho.MonoisotopicMass));
            Assert.That(sMod.Value.UnimodId, Is.EqualTo(21));

            // Check position 5 (M)
            var mMod = canonical.GetModificationAt(5);
            Assert.That(mMod, Is.Not.Null);
            Assert.That(mMod.Value.ResidueIndex, Is.EqualTo(5));
            Assert.That(mMod.Value.TargetResidue, Is.EqualTo('M'));
            Assert.That(mMod.Value.IsResolved, Is.False);
            mMod = lookup.TryResolve(mMod.Value);


            Assert.That(mMod.Value.MzLibModification, Is.EqualTo(oxidation));
            Assert.That(mMod.Value.ChemicalFormula, Is.EqualTo(oxidation.ChemicalFormula));
            Assert.That(mMod.Value.HasMass);
            Assert.That(mMod.Value.EffectiveMass, Is.EqualTo(oxidation.MonoisotopicMass));
            Assert.That(mMod.Value.UnimodId, Is.EqualTo(35));
        }
    }
}
