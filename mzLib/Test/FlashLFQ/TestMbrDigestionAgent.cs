using Chemistry;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Readers;
using Test.FileReadingTests;
using ChromatographicPeak = FlashLFQ.ChromatographicPeak;
using DetectionType = FlashLFQ.DetectionType;

namespace Test.FlashLFQ
{
    /// <summary>
    /// Match-between-runs may only transfer an identification into a file that could contain the analyte.
    /// Digestion decides that, so a transfer must not cross digestion agents - but it must still cross
    /// conditions and replicates, which is what match-between-runs exists to do.
    ///
    /// Every file here is written with the same MS1 peaks, so a refused transfer is a decision rather than
    /// an absence of signal: the peak the acceptor would have matched is present in all three files.
    /// </summary>
    public class TestMbrDigestionAgent
    {
        private const string TransferredSequence = "PEPTIDEVV";

        private static readonly List<string> PepSequences = new()
        {
            "PEPTIDE", "PEPTIDEV", TransferredSequence, "PEPTIDEVVV", "PEPTIDEVVVV", "PEPTIDEVVVVA", "PEPTIDEVVVVAA"
        };

        private static readonly double[] RetentionTimes = { 1.01, 1.02, 1.03, 1.035, 1.04, 1.045, 1.05 };

        private static string WriteMzml(string name)
        {
            var scans = new MsDataScan[PepSequences.Count];

            for (int p = 0; p < PepSequences.Count; p++)
            {
                ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(PepSequences[p]).GetChemicalFormula();
                IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);
                double[] mz = dist.Masses.Select(v => v.ToMz(1)).ToArray();
                double[] intensities = dist.Intensities.Select(v => v * 1e6).ToArray();

                scans[p] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: p + 1,
                    msnOrder: 1, isCentroid: true, polarity: Polarity.Positive, retentionTime: RetentionTimes[p],
                    scanWindowRange: new MzRange(400, 1600), scanFilter: "f", mzAnalyzer: MZAnalyzerType.Orbitrap,
                    totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (p + 1));
            }

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, name + ".mzML");
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(scans), path, false);
            return path;
        }

        private static Identification Id(SpectraFileInfo file, string sequence, string digestionAgent, double psmScore = 0)
        {
            int index = PepSequences.IndexOf(sequence);

            return new Identification(file, sequence, sequence,
                new Proteomics.AminoAcidPolymer.Peptide(sequence).MonoisotopicMass,
                RetentionTimes[index] + 0.001, 1,
                new List<ProteinGroup> { new ProteinGroup("MyProtein", "gene", "org") },
                psmScore: psmScore,
                digestionAgentName: digestionAgent);
        }

        /// <summary>
        /// Builds three files. The first two share a digestion agent and differ by biological replicate; the
        /// third uses a different agent. Only the first file identifies <see cref="TransferredSequence"/>.
        /// </summary>
        private static FlashLfqResults RunWithAgents(string prefix, string firstAgent, string secondAgent, string thirdAgent)
        {
            var file1 = new SpectraFileInfo(WriteMzml(prefix + "_1"), "a", 0, 0, 0);
            var file2 = new SpectraFileInfo(WriteMzml(prefix + "_2"), "a", 1, 0, 0);
            var file3 = new SpectraFileInfo(WriteMzml(prefix + "_3"), "a", 2, 0, 0);

            var identifications = new List<Identification>();

            // every file identifies these, giving the retention time calibration something to work with
            foreach (string shared in PepSequences.Where(s => s != TransferredSequence))
            {
                identifications.Add(Id(file1, shared, firstAgent));
                identifications.Add(Id(file2, shared, secondAgent));
                identifications.Add(Id(file3, shared, thirdAgent));
            }

            // only the first file identifies this one, so files 2 and 3 are candidates for a transfer
            identifications.Add(Id(file1, TransferredSequence, firstAgent));

            return new FlashLfqEngine(identifications, matchBetweenRuns: true, maxThreads: 1).Run();
        }

        private static bool WasTransferredInto(FlashLfqResults results, string fileNameWithoutExtension) => results.Peaks
            .Where(kvp => kvp.Key.FilenameWithoutExtension == fileNameWithoutExtension)
            .SelectMany(kvp => kvp.Value)
            .Any(peak => peak.DetectionType == DetectionType.MBR
                && peak.Identifications.Any(id => id.ModifiedSequence == TransferredSequence));

        [Test]
        public static void TransferCrossesReplicatesButNotDigestionAgents()
        {
            var results = RunWithAgents("mbr_agents", "trypsin", "trypsin", "Glu-C");

            Assert.That(WasTransferredInto(results, "mbr_agents_2"), Is.True,
                "a transfer between files sharing a digestion agent must still happen, or this test proves nothing");
            Assert.That(WasTransferredInto(results, "mbr_agents_3"), Is.False,
                "Glu-C cannot produce a tryptic peptide, so the identification must not be transferred into its file");
        }

        /// <summary>
        /// A donor is chosen per (sequence, digestion agent). Choosing one donor per sequence globally lets a
        /// file from the wrong agent win the slot, and the transfer that should have happened between the two
        /// files sharing an agent is then silently lost rather than merely misdirected.
        /// </summary>
        [Test]
        public static void DonorFromAnotherAgentDoesNotDisplaceTheSameAgentDonor()
        {
            var file1 = new SpectraFileInfo(WriteMzml("mbr_donor_1"), "a", 0, 0, 0);
            var file2 = new SpectraFileInfo(WriteMzml("mbr_donor_2"), "a", 1, 0, 0);
            var file3 = new SpectraFileInfo(WriteMzml("mbr_donor_3"), "a", 2, 0, 0);

            var identifications = new List<Identification>();

            foreach (string shared in PepSequences.Where(s => s != TransferredSequence))
            {
                identifications.Add(Id(file1, shared, "trypsin", psmScore: 1));
                identifications.Add(Id(file2, shared, "trypsin", psmScore: 1));
                identifications.Add(Id(file3, shared, "Glu-C", psmScore: 1));
            }

            // identified under both agents, but scored far higher under Glu-C, so it wins a global donor
            // contest. file2 shares trypsin with file1 and is the file that needs the transfer.
            identifications.Add(Id(file1, TransferredSequence, "trypsin", psmScore: 1));
            identifications.Add(Id(file3, TransferredSequence, "Glu-C", psmScore: 100));

            var results = new FlashLfqEngine(identifications, matchBetweenRuns: true, maxThreads: 1).Run();

            Assert.That(WasTransferredInto(results, "mbr_donor_2"), Is.True,
                "the trypsin donor in file 1 must still serve file 2, even though the Glu-C peak scores higher");
        }

        /// <summary>
        /// Callers that supply no digestion agent must see the behaviour they saw before the agent was added.
        /// </summary>
        [Test]
        public static void UnknownDigestionAgentLeavesTransfersUnrestricted()
        {
            var results = RunWithAgents("mbr_no_agents", null, null, null);

            Assert.That(WasTransferredInto(results, "mbr_no_agents_2"), Is.True);
            Assert.That(WasTransferredInto(results, "mbr_no_agents_3"), Is.True);
        }
    }
}
