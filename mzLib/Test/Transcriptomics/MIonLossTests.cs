using NUnit.Framework;
using Omics.Fragmentation;
using Chemistry;

namespace Test.Transcriptomics
{
    [TestFixture]
    public class MIonLossTests
    {
        [Test]
        public void PhosphoLoss_StaticInstance_IsInitializedCorrectly()
        {
            // Arrange
            var expectedName = "Phosphate Loss";
            var expectedAnnotation = "-P";
            var expectedFormula = ChemicalFormula.ParseFormula("H-1P-1O-3");

            // Act
            var phosphoLoss = MIonLoss.PhosphoLoss;

            // Assert
            Assert.That(phosphoLoss, Is.Not.Null);
            Assert.That(phosphoLoss.Annotation, Is.EqualTo(expectedAnnotation));
            Assert.That(phosphoLoss.Name, Is.EqualTo(expectedName));
            Assert.That(phosphoLoss.ThisChemicalFormula, Is.EqualTo(expectedFormula));
            Assert.That(phosphoLoss.MonoisotopicMass, Is.EqualTo(expectedFormula.MonoisotopicMass));
        }

        [Test]
        public void WaterLoss_StaticInstance_IsInitializedCorrectly()
        {
            // Arrange
            var expectedName = "Water Loss";
            var expectedAnnotation = "-H2O";
            var expectedFormula = ChemicalFormula.ParseFormula("H-2O-1");

            // Act
            var waterLoss = MIonLoss.WaterLoss;

            // Assert
            Assert.That(waterLoss, Is.Not.Null);
            Assert.That(waterLoss.Name, Is.EqualTo(expectedName));
            Assert.That(waterLoss.Annotation, Is.EqualTo(expectedAnnotation));
            Assert.That(waterLoss.ThisChemicalFormula, Is.EqualTo(expectedFormula));
            Assert.That(waterLoss.MonoisotopicMass, Is.EqualTo(expectedFormula.MonoisotopicMass));
        }

        [Test]
        public void Constructor_InitializesPropertiesCorrectly()
        {
            // Arrange
            var name = "TestLoss";
            var annotation = "-Test";
            var formula = ChemicalFormula.ParseFormula("C2H4O2");

            // Act
            var loss = new MIonLoss(name, annotation, formula);

            // Assert
            Assert.That(loss.Name, Is.EqualTo(name));
            Assert.That(loss.ThisChemicalFormula, Is.EqualTo(formula));
            Assert.That(loss.MonoisotopicMass, Is.EqualTo(formula.MonoisotopicMass));
        }

        [Test]
        public void ThisChemicalFormula_Setter_UpdatesFormulaAndMass()
        {
            // Arrange
            var loss = new MIonLoss("TestLoss", "TestLoss", ChemicalFormula.ParseFormula("H2O"));
            var newFormula = ChemicalFormula.ParseFormula("CO2");

            // Act
            loss.ThisChemicalFormula = newFormula;

            // Assert
            Assert.That(loss.ThisChemicalFormula, Is.EqualTo(newFormula));
            Assert.That(loss.MonoisotopicMass, Is.EqualTo(newFormula.MonoisotopicMass));
        }
    }
}
