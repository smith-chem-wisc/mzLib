using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.BioPolymerGroup;
using Omics.SpectralMatch;
using Proteomics.ProteolyticDigestion;
using Quantification;
using Quantification.Strategies;
using Test.Omics;
using ChromatographicPeak = FlashLFQ.ChromatographicPeak;
using DetectionType = FlashLFQ.DetectionType;
using FlashLfqEngine = FlashLFQ.FlashLfqEngine;
using FlashLfqIdentification = FlashLFQ.Identification;
using FlashLfqProteinGroup = FlashLFQ.ProteinGroup;
using FlashLfqResults = FlashLFQ.FlashLfqResults;
using MbrChromatographicPeak = FlashLFQ.MbrChromatographicPeak;
using Protein = Proteomics.Protein;
using SpectrumMatchTsvReader = Readers.SpectrumMatchTsvReader;

namespace Test.Quantification
{
    /// <summary>
    /// Feeds FlashLFQ's chromatographic peaks through the Quantification engine and checks that the
    /// peptide and protein intensities it arrives at are the ones FlashLFQ computed internally.
    ///
    /// The engine is configured to do what FlashLFQ does, step for step:
    ///  - each qualifying peak is one spectral match, with the peak intensity as its only channel;
    ///  - peaks roll up to peptides by taking the maximum, as FlashLfqResults.CalculatePeptideResults does;
    ///  - no normalization or collapsing, since each file is its own sample here and FlashLFQ's
    ///    normalization (off in this run) already acts on the peaks themselves;
    ///  - peptides roll up to proteins with MedianPolishRollUp, using unique peptides only.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class FlashLfqParityTests
    {
        private static string DataDirectory => Path.Combine(
            TestContext.CurrentContext.TestDirectory, "FlashLFQ", "TestData");

        private class TestExperimentalDesign : IExperimentalDesign
        {
            public Dictionary<string, ISampleInfo[]> FileNameSampleInfoDictionary { get; }

            public TestExperimentalDesign(Dictionary<string, ISampleInfo[]> dict)
            {
                FileNameSampleInfoDictionary = dict;
            }
        }

        private static List<SpectraFileInfo> FileInfos() => new()
        {
            new SpectraFileInfo(Path.Combine(DataDirectory, "20100614_Velos1_TaGe_SA_K562_3.mzML"), "a", 0, 0, 0),
            new SpectraFileInfo(Path.Combine(DataDirectory, "20100614_Velos1_TaGe_SA_K562_4.mzML"), "b", 0, 0, 0),
        };

        /// <summary>
        /// One FlashLFQ protein group object per accession, shared by every identification that names it.
        /// FlashLFQ compares protein groups by reference, so a fresh object per PSM would split one
        /// protein's peptides across several groups.
        /// </summary>
        private static List<FlashLfqIdentification> ReadIdentifications(List<SpectraFileInfo> files)
        {
            var proteinGroups = new Dictionary<string, FlashLfqProteinGroup>();
            var identifications = new List<FlashLfqIdentification>();

            foreach (var psm in SpectrumMatchTsvReader.ReadPsmTsv(Path.Combine(DataDirectory, "AllPSMs.psmtsv"), out _))
            {
                SpectraFileInfo file = files.FirstOrDefault(f =>
                    Path.GetFileNameWithoutExtension(f.FullFilePathWithExtension) == psm.FileNameWithoutExtension);
                if (file == null) continue;

                var proteins = psm.ProteinAccession.Split('|')
                    .Select(accession =>
                    {
                        if (!proteinGroups.TryGetValue(accession, out var pg))
                        {
                            pg = new FlashLfqProteinGroup(accession, "", "");
                            proteinGroups[accession] = pg;
                        }
                        return pg;
                    })
                    .Distinct()
                    .ToList();

                identifications.Add(new FlashLfqIdentification(file, psm.BaseSeq, psm.FullSequence,
                    (double)psm.MonoisotopicMass, (double)psm.RetentionTime, psm.PrecursorCharge, proteins,
                    decoy: psm.DecoyContamTarget == "D"));
            }

            return identifications;
        }

        /// <summary>
        /// The peaks FlashLFQ itself counts toward a peptide's intensity -- the same predicate as
        /// FlashLfqResults.CalculatePeptideResults: unambiguous by full sequence, not a decoy, and for an
        /// MBR peak, below the q-value threshold and not a random retention time decoy.
        /// </summary>
        private static bool QualifiesForPeptideQuant(ChromatographicPeak peak, FlashLfqResults results) =>
            peak.NumIdentificationsByFullSeq == 1
            && !peak.Identifications.First().IsDecoy
            && (peak.DetectionType != DetectionType.MBR
                || (peak is MbrChromatographicPeak mbr && mbr.MbrQValue < results.MbrQValueThreshold && !mbr.RandomRt))
            && results.PeptideModifiedSequences.ContainsKey(peak.Identifications.First().ModifiedSequence);

        [Test]
        [NonParallelizable] // full FlashLFQ runs; keep them off a loaded thread pool
        [TestCase(false, TestName = "FlashLfqParity_WithoutMbr")]
        [TestCase(true, TestName = "FlashLfqParity_WithMbr")]
        public void QuantificationEngineReproducesFlashLfq(bool matchBetweenRuns)
        {
            // 1) FlashLFQ's own peptide and protein results
            List<SpectraFileInfo> files = FileInfos();
            FlashLfqResults flashLfq = new FlashLfqEngine(ReadIdentifications(files),
                matchBetweenRuns: matchBetweenRuns, maxThreads: 1, silent: true).Run();

            // 2) The same peptides and protein groups, as the Quantification engine's types
            var peptides = flashLfq.PeptideModifiedSequences.Values
                .ToDictionary(p => p, p => (IBioPolymerWithSetMods)new PeptideWithSetModifications(p.Sequence));

            var proteinGroups = new Dictionary<FlashLfqProteinGroup, IBioPolymerGroup>();
            foreach (FlashLfqProteinGroup flashLfqGroup in flashLfq.ProteinGroups.Values)
            {
                var all = new HashSet<IBioPolymerWithSetMods>();
                var unique = new HashSet<IBioPolymerWithSetMods>();
                foreach (var kvp in peptides.Where(kvp => kvp.Key.UseForProteinQuant && kvp.Key.ProteinGroups.Contains(flashLfqGroup)))
                {
                    all.Add(kvp.Value);
                    if (kvp.Key.ProteinGroups.Count == 1)
                    {
                        unique.Add(kvp.Value);
                    }
                }

                proteinGroups[flashLfqGroup] = new BioPolymerGroup(
                    new HashSet<IBioPolymer> { new Protein("", flashLfqGroup.ProteinGroupName) }, all, unique);
            }

            // 3) One spectral match per qualifying peak, carrying the peak intensity as its only channel
            var peptideBySequence = peptides.ToDictionary(kvp => kvp.Key.Sequence, kvp => kvp.Value);
            var spectralMatches = new List<ISpectralMatch>();
            int scanNumber = 1;
            foreach (var (file, peaks) in flashLfq.Peaks)
            {
                foreach (ChromatographicPeak peak in peaks.Where(p => QualifiesForPeptideQuant(p, flashLfq)))
                {
                    FlashLfqIdentification id = peak.Identifications.First();
                    spectralMatches.Add(new MockSpectralMatch(file.FullFilePathWithExtension, id.ModifiedSequence,
                        id.BaseSequence, 0, scanNumber++, new[] { peptideBySequence[id.ModifiedSequence] })
                    {
                        Intensities = new[] { peak.Intensity }
                    });
                }
            }

            var design = new TestExperimentalDesign(files.ToDictionary(
                f => Path.GetFileName(f.FullFilePathWithExtension), f => new ISampleInfo[] { f }));

            var parameters = new QuantificationParameters
            {
                SpectralMatchNormalizationStrategy = new NoNormalization(),
                SpectralMatchToPeptideRollUpStrategy = new MaxRollUp(),
                PeptideNormalizationStrategy = new NoNormalization(),
                CollapseStrategy = new NoCollapse(),
                CollapseAggregationStrategy = new SumAggregation(),
                PeptideToProteinRollUpStrategy = new MedianPolishRollUp(),
                ProteinNormalizationStrategy = new NoNormalization(),
                UseSharedPeptidesForProteinQuant = false,
                WriteRawInformation = false,
                WritePeptideInformation = false,
                WriteProteinInformation = false,
            };

            QuantificationResults quant = new QuantificationEngine(parameters, design, spectralMatches,
                peptides.Values.ToList(), proteinGroups.Values.ToList()).Run();

            Assert.That(quant.Success, Is.True, quant.Summary);

            // 4) Peptides: every peptide in every file
            var peptideMismatches = new List<string>();
            int quantifiedPeptideValues = 0;
            foreach (var (flashLfqPeptide, peptide) in peptides)
            {
                quant.PeptideIntensities.TryGetValue(peptide, out var row);
                foreach (SpectraFileInfo file in files)
                {
                    double expected = flashLfqPeptide.GetIntensity(file);
                    double actual = row != null && row.TryGetValue(file, out double v) ? v : 0;

                    if (expected > 0) quantifiedPeptideValues++;
                    if (expected != actual)
                    {
                        peptideMismatches.Add($"{flashLfqPeptide.Sequence} in {file.FilenameWithoutExtension}: " +
                            $"FlashLFQ {expected}, Quantification {actual}");
                    }
                }
            }

            // 5) Proteins: every protein group in every sample. FlashLFQ reports NaN for a sample it saw
            // but could not quantify; the matrix has one marker for that, 0, so NaN is expected as 0.
            var proteinMismatches = new List<string>();
            int quantifiedProteinValues = 0;
            foreach (var (flashLfqGroup, group) in proteinGroups)
            {
                quant.ProteinIntensities.TryGetValue(group, out var row);
                foreach (SpectraFileInfo file in files)
                {
                    double expected = flashLfqGroup.GetIntensity(file);
                    if (double.IsNaN(expected)) expected = 0;
                    double actual = row != null && row.TryGetValue(file, out double v) ? v : 0;

                    if (expected > 0) quantifiedProteinValues++;
                    if (expected != actual)
                    {
                        proteinMismatches.Add($"{flashLfqGroup.ProteinGroupName} in {file.FilenameWithoutExtension}: " +
                            $"FlashLFQ {expected}, Quantification {actual}");
                    }
                }
            }

            // Guard against agreeing about nothing
            Assert.That(quantifiedPeptideValues, Is.GreaterThan(100));
            Assert.That(quantifiedProteinValues, Is.GreaterThan(50));

            Assert.That(peptideMismatches, Is.Empty,
                $"{peptideMismatches.Count} of {peptides.Count * files.Count} peptide intensities differ:\n  " +
                string.Join("\n  ", peptideMismatches.Take(20)));
            Assert.That(proteinMismatches, Is.Empty,
                $"{proteinMismatches.Count} of {proteinGroups.Count * files.Count} protein intensities differ:\n  " +
                string.Join("\n  ", proteinMismatches.Take(20)));
        }
    }
}
