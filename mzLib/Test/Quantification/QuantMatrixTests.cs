using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Quantification;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Test.Quantification
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class QuantMatrixTests
    {
        /// <summary>
        /// Verifies that SpectralMatchMatrix can be constructed with real BaseSpectralMatch objects
        /// and that data can be set and retrieved correctly via SetRow/GetRow.
        /// </summary>
        [Test]
        public void SpectralMatchMatrix_PopulateWithBaseSpectralMatch_CorrectlyStoresAndRetrievesData()
        {
            // Arrange: Create real BaseSpectralMatch instances
            var psm1 = new BaseSpectralMatch(
                fullFilePath: "Sample1.raw",
                oneBasedScanNumber: 100,
                score: 50.5,
                fullSequence: "PEPTIDEK",
                baseSequence: "PEPTIDEK");

            var psm2 = new BaseSpectralMatch(
                fullFilePath: "Sample1.raw",
                oneBasedScanNumber: 200,
                score: 75.3,
                fullSequence: "ANOTHERPEPTIDE",
                baseSequence: "ANOTHERPEPTIDE");

            var psm3 = new BaseSpectralMatch(
                fullFilePath: "Sample2.raw",
                oneBasedScanNumber: 150,
                score: 60.0,
                fullSequence: "THIRDPEPTIDE",
                baseSequence: "THIRDPEPTIDE");

            var spectralMatches = new List<ISpectralMatch> { psm1, psm2, psm3 };

            // Create sample info (columns)
            var sample1 = new SpectraFileInfo("Sample1.raw", "Control", biorep: 1, techrep: 1, fraction: 1);
            var sample2 = new SpectraFileInfo("Sample2.raw", "Treatment", biorep: 1, techrep: 1, fraction: 1);
            var samples = new List<ISampleInfo> { sample1, sample2 };

            // Act: Create matrix and populate with intensity values
            var matrix = new SpectralMatchMatrix(spectralMatches, samples);

            // Set intensity values for each PSM across samples
            // Row 0 (psm1): intensities in Sample1=1000, Sample2=500
            matrix.SetRow(psm1, new double[] { 1000.0, 500.0 });
            // Row 1 (psm2): intensities in Sample1=2000, Sample2=1500
            matrix.SetRow(psm2, new double[] { 2000.0, 1500.0 });
            // Row 2 (psm3): intensities in Sample1=0 (not detected), Sample2=3000
            matrix.SetRow(psm3, new double[] { 0.0, 3000.0 });

            // Assert: Retrieve and verify data
            double[] row0 = matrix.GetRow(psm1);
            Assert.That(row0[0], Is.EqualTo(1000.0));
            Assert.That(row0[1], Is.EqualTo(500.0));

            double[] row1 = matrix.GetRow(psm2);
            Assert.That(row1[0], Is.EqualTo(2000.0));
            Assert.That(row1[1], Is.EqualTo(1500.0));

            double[] row2 = matrix.GetRow(psm3);
            Assert.That(row2[0], Is.EqualTo(0.0));
            Assert.That(row2[1], Is.EqualTo(3000.0));

            // Verify by index access as well
            double[] rowByIndex = matrix.GetRow(0);
            Assert.That(rowByIndex, Is.EqualTo(row0));
        }

    }
}
