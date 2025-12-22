using MassSpectrometry;
using MassSpectrometry.ExperimentalDesign;
using NUnit.Framework;
using Omics;
using Quantification;
using Quantification.Interfaces;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;

namespace Test.Quantification
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class CollapseStrategyTests
    {
        #region Fraction Collapse Tests

        /// <summary>
        /// Verifies that SampleCollapseStrategy correctly collapses two fractions into one
        /// using Median aggregation.
        /// </summary>
        [Test]
        public void SampleCollapseStrategy_CollapseFractions_Median()
        {
            // Arrange: PSMs observed in both fractions
            var psm1 = new BaseSpectralMatch("file.raw", 100, 50.0, "PEPTIDEK", "PEPTIDEK");
            var psm2 = new BaseSpectralMatch("file.raw", 200, 60.0, "ANOTHERPEPTIDE", "ANOTHERPEPTIDE");
            var spectralMatches = new List<ISpectralMatch> { psm1, psm2 };

            // Two fractions from the same biological sample
            var fraction1 = new SpectraFileInfo("Sample_F1.raw", "Control", biorep: 1, techrep: 1, fraction: 1);
            var fraction2 = new SpectraFileInfo("Sample_F2.raw", "Control", biorep: 1, techrep: 1, fraction: 2);
            var samples = new List<ISampleInfo> { fraction1, fraction2 };

            var matrix = new SpectralMatchMatrix(spectralMatches, samples);
            matrix.SetRow(psm1, new double[] { 1000.0, 2000.0 }); // Median = 1500
            matrix.SetRow(psm2, new double[] { 500.0, 1500.0 });  // Median = 1000

            var strategy = new SampleCollapseStrategy(CollapseDimension.Fraction, AggregationType.Median);

            // Act
            var collapsed = strategy.CollapseSamples(matrix);

            // Assert
            Assert.That(collapsed.ColumnKeys.Count, Is.EqualTo(1));
            Assert.That(collapsed.GetRow(psm1)[0], Is.EqualTo(1500.0));
            Assert.That(collapsed.GetRow(psm2)[0], Is.EqualTo(1000.0));
            Assert.That(collapsed.ColumnKeys[0].Fraction, Is.EqualTo(0)); // Collapsed indicator
        }

        #endregion

        #region Technical Replicate Collapse Tests

        /// <summary>
        /// Verifies that SampleCollapseStrategy correctly collapses technical replicates.
        /// </summary>
        [Test]
        public void SampleCollapseStrategy_CollapseTechnicalReplicates_Median()
        {
            // Arrange
            var psm1 = new BaseSpectralMatch("file.raw", 100, 50.0, "PEPTIDEK", "PEPTIDEK");
            var spectralMatches = new List<ISpectralMatch> { psm1 };

            // Two technical replicates from the same biological sample (same fraction)
            var techRep1 = new SpectraFileInfo("Sample_T1.raw", "Control", biorep: 1, techrep: 1, fraction: 1);
            var techRep2 = new SpectraFileInfo("Sample_T2.raw", "Control", biorep: 1, techrep: 2, fraction: 1);
            var samples = new List<ISampleInfo> { techRep1, techRep2 };

            var matrix = new SpectralMatchMatrix(spectralMatches, samples);
            matrix.SetRow(psm1, new double[] { 1000.0, 3000.0 }); // Median = 2000

            var strategy = new SampleCollapseStrategy(CollapseDimension.TechnicalReplicate, AggregationType.Median);

            // Act
            var collapsed = strategy.CollapseSamples(matrix);

            // Assert
            Assert.That(collapsed.ColumnKeys.Count, Is.EqualTo(1));
            Assert.That(collapsed.GetRow(psm1)[0], Is.EqualTo(2000.0));
            Assert.That(collapsed.ColumnKeys[0].TechnicalReplicate, Is.EqualTo(0)); // Collapsed indicator
            Assert.That(collapsed.ColumnKeys[0].Fraction, Is.EqualTo(1)); // Preserved
        }

        /// <summary>
        /// Verifies that technical replicates from different bio reps are NOT collapsed together.
        /// </summary>
        [Test]
        public void SampleCollapseStrategy_CollapseTechnicalReplicates_DifferentBioReps_NotCollapsed()
        {
            // Arrange
            var psm1 = new BaseSpectralMatch("file.raw", 100, 50.0, "PEPTIDEK", "PEPTIDEK");
            var spectralMatches = new List<ISpectralMatch> { psm1 };

            // Two tech reps from DIFFERENT biological replicates
            var bioRep1_tech1 = new SpectraFileInfo("B1_T1.raw", "Control", biorep: 1, techrep: 1, fraction: 1);
            var bioRep2_tech1 = new SpectraFileInfo("B2_T1.raw", "Control", biorep: 2, techrep: 1, fraction: 1);
            var samples = new List<ISampleInfo> { bioRep1_tech1, bioRep2_tech1 };

            var matrix = new SpectralMatchMatrix(spectralMatches, samples);
            matrix.SetRow(psm1, new double[] { 1000.0, 2000.0 });

            var strategy = new SampleCollapseStrategy(CollapseDimension.TechnicalReplicate, AggregationType.Median);

            // Act
            var collapsed = strategy.CollapseSamples(matrix);

            // Assert: Should still have 2 columns (different bio reps not collapsed)
            Assert.That(collapsed.ColumnKeys.Count, Is.EqualTo(2));
            Assert.That(collapsed.GetRow(psm1)[0], Is.EqualTo(1000.0));
            Assert.That(collapsed.GetRow(psm1)[1], Is.EqualTo(2000.0));
        }

        #endregion

        #region Biological Replicate Collapse Tests

        /// <summary>
        /// Verifies that SampleCollapseStrategy correctly collapses biological replicates.
        /// </summary>
        [Test]
        public void SampleCollapseStrategy_CollapseBiologicalReplicates_Average()
        {
            // Arrange
            var psm1 = new BaseSpectralMatch("file.raw", 100, 50.0, "PEPTIDEK", "PEPTIDEK");
            var spectralMatches = new List<ISpectralMatch> { psm1 };

            // Three biological replicates from the same condition
            var bioRep1 = new SpectraFileInfo("B1.raw", "Treatment", biorep: 1, techrep: 1, fraction: 1);
            var bioRep2 = new SpectraFileInfo("B2.raw", "Treatment", biorep: 2, techrep: 1, fraction: 1);
            var bioRep3 = new SpectraFileInfo("B3.raw", "Treatment", biorep: 3, techrep: 1, fraction: 1);
            var samples = new List<ISampleInfo> { bioRep1, bioRep2, bioRep3 };

            var matrix = new SpectralMatchMatrix(spectralMatches, samples);
            matrix.SetRow(psm1, new double[] { 1000.0, 2000.0, 3000.0 }); // Average = 2000

            var strategy = new SampleCollapseStrategy(CollapseDimension.BiologicalReplicate, AggregationType.Average);

            // Act
            var collapsed = strategy.CollapseSamples(matrix);

            // Assert
            Assert.That(collapsed.ColumnKeys.Count, Is.EqualTo(1));
            Assert.That(collapsed.GetRow(psm1)[0], Is.EqualTo(2000.0));
            Assert.That(collapsed.ColumnKeys[0].BiologicalReplicate, Is.EqualTo(0)); // Collapsed indicator
            Assert.That(collapsed.ColumnKeys[0].Condition, Is.EqualTo("Treatment")); // Preserved
        }

        #endregion

        #region Aggregation Type Tests

        [Test]
        public void SampleCollapseStrategy_Sum_CorrectlyAggregates()
        {
            // Arrange
            var psm1 = new BaseSpectralMatch("file.raw", 100, 50.0, "PEPTIDEK", "PEPTIDEK");
            var fraction1 = new SpectraFileInfo("F1.raw", "Control", biorep: 1, techrep: 1, fraction: 1);
            var fraction2 = new SpectraFileInfo("F2.raw", "Control", biorep: 1, techrep: 1, fraction: 2);

            var matrix = new SpectralMatchMatrix(
                new List<ISpectralMatch> { psm1 },
                new List<ISampleInfo> { fraction1, fraction2 });
            matrix.SetRow(psm1, new double[] { 1000.0, 3000.0 }); // Sum = 4000

            var strategy = new SampleCollapseStrategy(CollapseDimension.Fraction, AggregationType.Sum);

            // Act
            var collapsed = strategy.CollapseSamples(matrix);

            // Assert
            Assert.That(collapsed.GetRow(psm1)[0], Is.EqualTo(4000.0));
        }

        #endregion

        #region Name Property Tests

        [Test]
        public void SampleCollapseStrategy_Name_ReflectsDimensionAndAggregation()
        {
            Assert.That(
                new SampleCollapseStrategy(CollapseDimension.Fraction, AggregationType.Median).Name,
                Is.EqualTo("Collapse_Fraction_Median"));

            Assert.That(
                new SampleCollapseStrategy(CollapseDimension.TechnicalReplicate, AggregationType.Average).Name,
                Is.EqualTo("Collapse_TechnicalReplicate_Average"));

            Assert.That(
                new SampleCollapseStrategy(CollapseDimension.BiologicalReplicate, AggregationType.Sum).Name,
                Is.EqualTo("Collapse_BiologicalReplicate_Sum"));
        }

        #endregion
    }
}