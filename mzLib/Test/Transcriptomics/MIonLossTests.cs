using NUnit.Framework;
using Omics.Fragmentation;
using Chemistry;
using Transcriptomics.Digestion;
using Transcriptomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Readers;

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
            var expectedFormula = ChemicalFormula.ParseFormula("H1P1O3");

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
            var expectedFormula = ChemicalFormula.ParseFormula("H2O1");

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

        [Test]
        public void MatchedIonHasExpectedAnnotation()
        {
            var oligo = new OligoWithSetMods("GUACUG", []);
            var fragmentParams = new RnaFragmentationParams();
            fragmentParams.MIonLosses.Add(MIonLoss.WaterLoss);
            fragmentParams.MIonLosses.Add(MIonLoss.PhosphoLoss);

            List<Product> allProducts = new List<Product>();
            oligo.Fragment(MassSpectrometry.DissociationType.CID, FragmentationTerminus.Both, allProducts, fragmentParams);

            var mIons = allProducts.FindAll(p => p.ProductType == ProductType.M);
            Assert.That(mIons.Count, Is.EqualTo(3));

            List<MatchedFragmentIon> matchedIons = new();
            List<string> expectedAnnotation = new();
            foreach (var mIon in mIons)
            {
                matchedIons.Add(new MatchedFragmentIon(mIon, 1.0, 1.0, -2));
                expectedAnnotation.Add(mIon.Annotation + "-2");
            }

            for (int i = 0; i < matchedIons.Count; i++)
            {
                Assert.That(matchedIons[i].Annotation, Is.EqualTo(expectedAnnotation[i]));
            }
        }

        [Test]
        public static void ReadsInCorrectly()
        {
            string osmPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "OsmWithCustomMIons.osmtsv");
            var osm = SpectrumMatchTsvReader.ReadOsmTsv(osmPath, out var warnings);


            Assert.That(warnings.Count, Is.EqualTo(0));
            Assert.That(osm.Count, Is.EqualTo(1));
            var oligo = osm[0];
            // Extract all M ions (including custom)
            var mIons = oligo.MatchedIons
                .Where(ion => ion.NeutralTheoreticalProduct.ProductType == ProductType.M)
                .ToList();

            // Expected annotations as they would appear after parsing
            var expectedAnnotations = new[]
            {
                "M-5",         // Legacy M ion
                "M-5",         // new M
                "M-P-5",       // custom phosphate loss
                "M-A-5",       // custom adenosine loss
                "M-A-P-H20-5"  // custom multi-loss
            };
            Assert.That(mIons.Count, Is.EqualTo(expectedAnnotations.Length), "Unexpected number of M ions found.");

            // Check that all expected M ions are present
            foreach (var expected in expectedAnnotations)
            {
                Assert.That(mIons.Any(m => m.Annotation.Contains(expected)), $"Missing expected M ion: {expected}");
            }

            // Optionally, check charge, m/z, and intensity for each
            foreach (var mIon in mIons)
            {
                Assert.That(mIon.Charge, Is.EqualTo(-5)); // Should be nonzero
                Assert.That(mIon.Mz, Is.GreaterThan(0));
                Assert.That(mIon.Intensity, Is.GreaterThanOrEqualTo(0));
            }
        }
    }
}
