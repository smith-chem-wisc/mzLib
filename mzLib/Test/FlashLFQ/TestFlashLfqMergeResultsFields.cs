using Chemistry;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Test.FileReadingTests;

namespace Test.FlashLFQ
{
    /// <summary>
    /// A peptide row in the QuantifiedPeptides output carries one intensity, one retention time and one
    /// detection type per file. MergeResultsWith has to carry all three across for a peptide both runs saw,
    /// not just intensity and detection type - the retention time it drops is a measurement belonging to
    /// the merged-in file, not something transferred between files.
    /// </summary>
    public class TestFlashLfqMergeResultsFields
    {
        private const string Sequence = "PEPTIDE";
        private const double RetentionTime = 1.02;

        private static string WriteMzml(string name)
        {
            ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(Sequence).GetChemicalFormula();
            IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);

            var scan = new MsDataScan(
                massSpectrum: new MzSpectrum(dist.Masses.Select(v => v.ToMz(1)).ToArray(),
                    dist.Intensities.Select(v => v * 1e6).ToArray(), false),
                oneBasedScanNumber: 1, msnOrder: 1, isCentroid: true, polarity: Polarity.Positive,
                retentionTime: RetentionTime, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 1e6, injectionTime: 1.0, noiseData: null,
                nativeId: "scan=1");

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, name + ".mzML");
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(new[] { scan }), path, false);
            return path;
        }

        private static FlashLfqResults QuantifyOneFile(string name, string proteinGroupName, out SpectraFileInfo file)
        {
            file = new SpectraFileInfo(WriteMzml(name), "a", 0, 0, 0);

            var identification = new Identification(file, Sequence, Sequence,
                new Proteomics.AminoAcidPolymer.Peptide(Sequence).MonoisotopicMass,
                RetentionTime, 1, new List<ProteinGroup> { new ProteinGroup(proteinGroupName, "gene", "org") });

            return new FlashLfqEngine(new List<Identification> { identification }, maxThreads: 1).Run();
        }

        [Test]
        public static void MergeCarriesRetentionTimeAndProteinGroupsForASharedPeptide()
        {
            var resultsA = QuantifyOneFile("merge_fields_a", "ProteinA", out SpectraFileInfo fileA);
            var resultsB = QuantifyOneFile("merge_fields_b", "ProteinB", out SpectraFileInfo fileB);

            double retentionTimeBeforeMerge = resultsB.PeptideModifiedSequences[Sequence].GetRetentionTime(fileB);
            Assert.That(retentionTimeBeforeMerge, Is.GreaterThan(0),
                "the peptide must be quantified in the second run, or this test proves nothing");

            resultsA.MergeResultsWith(resultsB);

            Peptide merged = resultsA.PeptideModifiedSequences[Sequence];

            Assert.That(merged.GetRetentionTime(fileA), Is.GreaterThan(0));
            Assert.That(merged.GetRetentionTime(fileB), Is.EqualTo(retentionTimeBeforeMerge),
                "the retention time measured in the merged-in file must survive the merge");
            Assert.That(merged.GetIntensity(fileB), Is.GreaterThan(0));
            Assert.That(merged.ProteinGroups.Select(p => p.ProteinGroupName),
                Is.EquivalentTo(new[] { "ProteinA", "ProteinB" }),
                "protein groups found only in the merged-in run must not be dropped");
        }
    }
}
