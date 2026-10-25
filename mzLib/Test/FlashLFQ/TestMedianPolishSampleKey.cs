using FlashLFQ;
using MassSpectrometry;
using NUnit.Framework;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using Peptide = FlashLFQ.Peptide;

namespace Test.FlashLFQ
{
    /// <summary>
    /// Median polish counts its samples, one per (condition, biological replicate). It used to count them by the
    /// string Condition + BiologicalReplicate, so condition "A1" biorep 0 and condition "A" biorep 10 were one sample
    /// ("A10"): the peptide matrix got too few columns and filling it threw IndexOutOfRangeException.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class TestMedianPolishSampleKey
    {
        /// <summary>Protein intensities for three peptides over two samples, one file each.</summary>
        private static double[] ProteinIntensities(SpectraFileInfo first, SpectraFileInfo second)
        {
            var files = new List<SpectraFileInfo> { first, second };
            var results = new FlashLfqResults(files, new List<Identification>());
            var protein = new ProteinGroup("accession", "gene", "organism");
            results.ProteinGroups.Add(protein.ProteinGroupName, protein);

            double[][] intensities = { new[] { 1000.0, 2100 }, new[] { 3000.0, 5900 }, new[] { 500.0, 1050 } };
            for (int row = 0; row < intensities.Length; row++)
            {
                var peptide = new Peptide("PEPTIDE" + row, "PEPTIDE" + row, true, new HashSet<ProteinGroup> { protein });
                results.PeptideModifiedSequences.Add(peptide.Sequence, peptide);
                for (int col = 0; col < files.Count; col++)
                {
                    peptide.SetIntensity(files[col], intensities[row][col]);
                    peptide.SetDetectionType(files[col], DetectionType.MSMS);
                }
            }

            results.CalculateProteinResultsMedianPolish(useSharedPeptides: false);
            return files.Select(protein.GetIntensity).ToArray();
        }

        [Test]
        public static void ConditionA1Biorep0AndConditionABiorep10AreTwoSamples()
        {
            var colliding = ProteinIntensities(new SpectraFileInfo("a1", "A", 10, 0, 0), new SpectraFileInfo("b", "A1", 0, 0, 0));
            var distinct = ProteinIntensities(new SpectraFileInfo("a1", "A", 10, 0, 0), new SpectraFileInfo("b", "B", 0, 0, 0));

            Assert.That(colliding, Is.EqualTo(distinct), "the names do not change the quantification");
            Assert.That(colliding, Has.All.GreaterThan(0));
        }
    }
}
