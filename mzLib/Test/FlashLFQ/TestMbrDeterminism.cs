using FlashLFQ;
using MassSpectrometry;
using NUnit.Framework;
using Readers;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test.FlashLFQ
{
    /// <summary>
    /// Match-between-runs has to give the same answer twice. Issue #1111 reported peptide intensities and
    /// protein results changing between runs on identical inputs, with the chromatographic peaks stable
    /// throughout - a peptide whose only peak was a real MBR transfer reported its intensity in one run and
    /// zero in the next.
    ///
    /// The path was: training rows for the PEP model were gathered by appending each parallel chunk as it
    /// finished, so their order depended on thread scheduling. FastTree runs single-threaded from a fixed
    /// seed, but row order changes its binning and split choices, so the model differed, so every PEP
    /// differed. PEP orders the target-decoy walk that assigns MbrQValue, and the peptide roll-up keeps an
    /// MBR peak only while its q-value is below the threshold.
    /// </summary>
    [TestFixture]
    public static class TestMbrDeterminism
    {
        private static string DataDirectory => Path.Combine(
            TestContext.CurrentContext.TestDirectory, "FlashLFQ", "TestData");

        private static List<Identification> ReadIdentifications(List<SpectraFileInfo> files)
        {
            var identifications = new List<Identification>();
            foreach (var psm in SpectrumMatchTsvReader.ReadPsmTsv(Path.Combine(DataDirectory, "AllPSMs.psmtsv"), out _))
            {
                SpectraFileInfo file = files.FirstOrDefault(f =>
                    Path.GetFileNameWithoutExtension(f.FullFilePathWithExtension) == psm.FileNameWithoutExtension);
                if (file == null) continue;

                var proteins = psm.ProteinAccession.Split('|').Select(a => new ProteinGroup(a, "", "")).ToList();
                identifications.Add(new Identification(file, psm.BaseSeq, psm.FullSequence,
                    (double)psm.MonoisotopicMass, (double)psm.RetentionTime, psm.PrecursorCharge, proteins,
                    decoy: psm.DecoyContamTarget == "D"));
            }
            return identifications;
        }

        /// <summary>
        /// Per-MBR-peak posterior error probability, keyed so the same transfer is comparable across runs.
        /// This is the sensitive measure: the peptide intensities only move when a changed PEP happens to
        /// push a q-value across the reporting threshold, which made an intensity-only test miss the defect
        /// four times in five.
        /// </summary>
        private static Dictionary<string, double?> QuantifyAndTakePeakPeps(out Dictionary<string, double> intensitiesOut)
        {
            var files = FileInfos();
            var results = new FlashLfqEngine(ReadIdentifications(files), matchBetweenRuns: true, maxThreads: -1).Run();

            var peps = new Dictionary<string, double?>();
            foreach (var fileAndPeaks in results.Peaks)
            {
                foreach (MbrChromatographicPeak peak in fileAndPeaks.Value
                    .Where(p => p.DetectionType == DetectionType.MBR)
                    .Cast<MbrChromatographicPeak>())
                {
                    Identification donor = peak.Identifications.First();
                    peps[$"{fileAndPeaks.Key.FilenameWithoutExtension}|{donor.ModifiedSequence}|{donor.PrecursorChargeState}"] = peak.MbrPep;
                }
            }

            intensitiesOut = new Dictionary<string, double>();
            foreach (Peptide peptide in results.PeptideModifiedSequences.Values)
            {
                foreach (SpectraFileInfo file in files)
                {
                    intensitiesOut[peptide.Sequence + "|" + file.FilenameWithoutExtension] = peptide.GetIntensity(file);
                }
            }

            return peps;
        }

        private static List<SpectraFileInfo> FileInfos() => new()
        {
            new SpectraFileInfo(Path.Combine(DataDirectory, "20100614_Velos1_TaGe_SA_K562_3.mzML"), "a", 0, 0, 0),
            new SpectraFileInfo(Path.Combine(DataDirectory, "20100614_Velos1_TaGe_SA_K562_4.mzML"), "a", 1, 0, 0),
        };

        /// <summary>
        /// Two runs over the same inputs, with threading left at its default, have to agree on every MBR
        /// peak's PEP and on every peptide intensity. Before the fix the PEPs differed on 5 to 142 of the
        /// 145 MBR peaks on every attempt, and six peptides swung between a real intensity and zero.
        /// </summary>
        [Test]
        [NonParallelizable] // two full quantifications; keep them off a loaded thread pool
        public static void RepeatedQuantificationGivesTheSameAnswer()
        {
            Dictionary<string, double?> firstPeps = QuantifyAndTakePeakPeps(out var firstIntensities);
            Dictionary<string, double?> secondPeps = QuantifyAndTakePeakPeps(out var secondIntensities);

            Assert.That(secondPeps.Count, Is.EqualTo(firstPeps.Count), "the same MBR peaks must be found");

            var pepsThatMoved = firstPeps
                .Where(kv => !secondPeps.TryGetValue(kv.Key, out double? v) || v != kv.Value)
                .Select(kv => $"{kv.Key}: {kv.Value} then {(secondPeps.TryGetValue(kv.Key, out double? v) ? v : null)}")
                .ToList();

            Assert.That(pepsThatMoved, Is.Empty,
                $"{pepsThatMoved.Count} of {firstPeps.Count} MBR peaks got a different posterior error " +
                "probability on the second run, which reorders the target-decoy walk behind MbrQValue:\n  " +
                string.Join("\n  ", pepsThatMoved.Take(10)));

            var intensitiesThatMoved = firstIntensities
                .Where(kv => !secondIntensities.TryGetValue(kv.Key, out double v) || v != kv.Value)
                .Select(kv => $"{kv.Key}: {kv.Value} then {(secondIntensities.TryGetValue(kv.Key, out double v) ? v.ToString() : "absent")}")
                .ToList();

            Assert.That(intensitiesThatMoved, Is.Empty,
                "peptide intensities changed between two runs on identical input:\n  " +
                string.Join("\n  ", intensitiesThatMoved));
        }
    }
}
