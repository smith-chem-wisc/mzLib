using NUnit.Framework;
using Readers;
using System.Collections.Generic;
using System.Linq;
using FlashLFQ;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System.IO;
using FlashLFQ.PEP;
using System;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Test.FileReadingTests;
using UsefulProteomicsDatabases;


namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class TestPipEcho
    {
        [Test]
        [TestCase(3)]
        [TestCase(5)]
        public static void TestDonorGroupEqualizer(int numGroups)
        {
            SpectraFileInfo fakeFile = new SpectraFileInfo("fakeFile", "A", 1, 1, 1);
            Identification id = new Identification(fakeFile, "KPVGAAK", "KPVGAAK", 669.4173, 1.9398, 2, new List<ProteinGroup> { new ProteinGroup("P16403", "H12", "HUMAN") });

            MbrChromatographicPeak targetPeak = new MbrChromatographicPeak(id, fakeFile, 1, randomRt: false);
            MbrChromatographicPeak decoyPeak = new MbrChromatographicPeak(id, fakeFile, 1, randomRt: true);
            targetPeak.MbrScore = 100;

            Random random = new Random(42);
            List<DonorGroup> donorGroups = new List<DonorGroup>();
            for (int i = 0; i < 10000; i++)
            {
                int numberTargets = random.Next(0, 10);
                int numberDecoys = random.Next(0, 10);
                donorGroups.Add(new DonorGroup(id, Enumerable.Repeat(targetPeak, numberTargets).ToList(), Enumerable.Repeat(decoyPeak, numberDecoys).ToList()));
            }
            
            donorGroups = PepAnalysisEngine.OrderDonorGroups(donorGroups);
            var donorIndices = PepAnalysisEngine.GetDonorGroupIndices(donorGroups, numGroups: numGroups, scoreCutoff: 50);

            Assert.That(donorIndices.Count, Is.EqualTo(numGroups));
            List<int> targetPeakCounts = new();
            List<int> decoyPeakCounts = new();
            for (int i = 0; i < numGroups; i++)
            {
                int targetSum = 0;
                int decoySum = 0;
                foreach (int idx in donorIndices[i])
                {
                    targetSum += donorGroups[idx].TargetAcceptors.Count;
                    decoySum += donorGroups[idx].DecoyAcceptors.Count;
                }
                targetPeakCounts.Add(targetSum);
                decoyPeakCounts.Add(decoySum);
            } 

            // Assert that each group has an approximately equal number of target peaks
            Assert.That(targetPeakCounts.Max() - targetPeakCounts.Min(), Is.LessThanOrEqualTo(numGroups-1));
            // Assert that each group has an approximately equal number of decoy peaks
            Assert.That(decoyPeakCounts.Max() - decoyPeakCounts.Min(), Is.LessThanOrEqualTo(numGroups - 1));
        }

        [Test]
        public static void TestMbrScorer()
        {
            SpectraFileInfo fakeFile = new SpectraFileInfo("fakeFile", "A", 1, 1, 1);
            SpectraFileInfo fakeDonorFile = new SpectraFileInfo("fakeFile", "A", 1, 1, 1);

            double idMass = 669.4173;
            Identification id = new Identification(fakeFile, "KPVGAAK", "KPVGAAK", 669.4173, 1.9398, 2, new List<ProteinGroup> { new ProteinGroup("P16403", "H12", "HUMAN") });
            Identification id2 = new Identification(fakeFile, "KPVGK", "KPVGK", 669.4173, 1.9398, 2, new List<ProteinGroup> { new ProteinGroup("P16403", "H12", "HUMAN") });
            Identification donorId = new Identification(fakeFile, "KPVGK", "KPVGK", 669.4173, 1.9398, 2, new List<ProteinGroup> { new ProteinGroup("P16403", "H12", "HUMAN") });
            id.PeakfindingMass = idMass;
            id2.PeakfindingMass = idMass;
            donorId.PeakfindingMass = idMass;

            var peak1 = new ChromatographicPeak(id, fakeFile);
            var peak2 = new ChromatographicPeak(id, fakeFile);
            var peak3 = new ChromatographicPeak(id2, fakeFile);
            var peak4 = new ChromatographicPeak(id, fakeFile);
            var donorPeak = new ChromatographicPeak(donorId, fakeDonorFile);
            var acceptorPeak = new MbrChromatographicPeak(donorId, fakeFile, 1, false);

            IndexedMassSpectralPeak imsPeak = new IndexedMassSpectralPeak((idMass + 0.001).ToMz(1), 1.1, 1, 25);
            IndexedMassSpectralPeak imsPeak2 = new IndexedMassSpectralPeak((idMass - 0.001).ToMz(1), 1, 2, 26);
            var iso1 = new FlashLFQ.IsotopicEnvelope(imsPeak, 1, 1, 0.98);
            var iso2 = new FlashLFQ.IsotopicEnvelope(imsPeak2, 1, 1, 0.9);

            peak1.IsotopicEnvelopes.Add(iso1);
            peak1.IsotopicEnvelopes.Add(iso2);
            peak1.CalculateIntensityForThisFeature(false);

            peak4.IsotopicEnvelopes.Add(iso2);
            peak4.CalculateIntensityForThisFeature(false);

            donorPeak.IsotopicEnvelopes.Add(iso2);
            donorPeak.CalculateIntensityForThisFeature(false);

            acceptorPeak.IsotopicEnvelopes.Add(iso1);
            acceptorPeak.CalculateIntensityForThisFeature(false);


            var peakList = new List<ChromatographicPeak> { peak1, peak4 };
            var peakDict = peakList.ToDictionary(keySelector: p => p.Apex.IndexedPeak, elementSelector: p => p);

            // Builds a scorer. Ppm Error and Intensity distributions both have mean and std-dev of 1
            MbrScorer scorer = new MbrScorer(peakDict, peakList, new MathNet.Numerics.Distributions.Normal(1, 1), new MathNet.Numerics.Distributions.Normal(1,1));
            
            scorer.AddRtPredErrorDistribution(fakeDonorFile, new List<double> { 0.5, 0.6, 0.5, 0.6, 0.5, 0.6, 0.5 }, 2);

            acceptorPeak.MbrScore = scorer.ScoreMbr(acceptorPeak, donorPeak, predictedRt: 25.1);

            Assert.That(acceptorPeak.MbrScore, Is.EqualTo(58.7).Within(0.1));
            Assert.That(acceptorPeak.PpmScore, Is.EqualTo(0.62).Within(0.01));
            Assert.That(acceptorPeak.IntensityScore, Is.EqualTo(0.32).Within(0.01));
            Assert.That(acceptorPeak.RtScore, Is.EqualTo(0.96).Within(0.01));
            Assert.That(acceptorPeak.ScanCountScore, Is.EqualTo(0.5).Within(0.01));
            Assert.That(acceptorPeak.IsotopicDistributionScore, Is.EqualTo(0.74).Within(0.01));
        }

        [Test]
        public static void TestSpectraFileInfoString()
        {
            SpectraFileInfo fakeFile = new SpectraFileInfo(@"C:\Users\xyz\data\fakeFile.raw", "A", 1, 1, 1); 
            Assert.AreEqual("fakeFile.raw", fakeFile.ToString());
        }

        [Test]
        public static void TestChromatographicPeakEquals()
        {
            SpectraFileInfo fakeFile = new SpectraFileInfo("fakeFile", "A", 1, 1, 1);
            Identification id = new Identification(fakeFile, "KPVGAAK", "KPVGAAK", 669.4173, 1.9398, 2, new List<ProteinGroup> { new ProteinGroup("P16403", "H12", "HUMAN") });
            Identification id2 = new Identification(fakeFile, "KPVGK", "KPVGK", 669.4173, 1.9398, 2, new List<ProteinGroup> { new ProteinGroup("P16403", "H12", "HUMAN") });

            var peak1 = new MbrChromatographicPeak(id,  fakeFile, 1, randomRt: false);
            var peak2 = new MbrChromatographicPeak(id,  fakeFile, 1, randomRt: false);
            var peak3 = new MbrChromatographicPeak(id2, fakeFile, 1, randomRt: false);
            var peak4 = new MbrChromatographicPeak(id, fakeFile, 1, randomRt: false);


            IndexedMassSpectralPeak imsPeak = new IndexedMassSpectralPeak(1, 1, 1, 25);
            IndexedMassSpectralPeak imsPeak2 = new IndexedMassSpectralPeak(1, 1, 1, 50);
            var iso1 = new FlashLFQ.IsotopicEnvelope(imsPeak, 1, 1, 1);
            var iso2 = new FlashLFQ.IsotopicEnvelope(imsPeak2, 1, 1, 1);

            peak1.IsotopicEnvelopes.Add(iso1);
            peak1.CalculateIntensityForThisFeature(false);

            peak2.IsotopicEnvelopes.Add(iso1);
            peak2.CalculateIntensityForThisFeature(false);

            peak3.IsotopicEnvelopes.Add(iso1);
            peak3.CalculateIntensityForThisFeature(false);

            peak4.IsotopicEnvelopes.Add(iso2);
            peak4.CalculateIntensityForThisFeature(false);

            Assert.That(peak1.Equals(peak2));
            Assert.That(!peak1.Equals(peak3));
            Assert.That(!peak1.Equals(peak4));

        }

        /// <summary>
        /// This test MatchBetweenRuns by creating two fake mzML files and a list of fake IDs. 
        /// There are multiple sets of IDs, where most are shared between the two runs but one+ is/are missing
        /// MBR is tested by ensuring that IDs are transferred between runs
        /// </summary>
        [Test]
        public static void TestFlashLfqMatchBetweenRunsNearestNeighborDonors()
        {
            List<string> filesToWrite = new List<string> { "mzml_1", "mzml_2", "mzml_3" };
            List<string> pepSequences = new List<string>
                {
                "PEPTIDE",
                "PEPTIDEV",
                "PEPTIDEVV",
                "TARGETPEP",
                "PEPTIDEVVV",
                "PEPTIDEVVVV",
                "PEPTIDEVVVVA",
                "PEPTIDEVVVVAA"
            };
            double intensity = 1e6;

            double[] file1Rt = new double[] { 1.01, 1.02, 1.03, 1.033, 1.035, 1.04, 1.045, 1.05 };
            double[] file2Rt = new double[] { 1.00, 1.025, 1.03, 1.031, 1.035, 1.04, 1.055, 1.07 };

            // generate mzml files (5 peptides each)
            for (int f = 0; f < filesToWrite.Count; f++)
            {
                // 1 MS1 scan per peptide
                MsDataScan[] scans = new MsDataScan[8];

                for (int p = 0; p < pepSequences.Count; p++)
                {
                    ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(pepSequences[p]).GetChemicalFormula();
                    IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);
                    double[] mz = dist.Masses.Select(v => v.ToMz(1)).ToArray();
                    double[] intensities = dist.Intensities.Select(v => v * intensity).ToArray();
                    if(f == 2)
                    {
                        // Make file 3 the most intense
                        intensities = intensities.Select(v => v * 5).ToArray();
                    }
                    double rt;
                    if (f == 1)
                    {
                        rt = file2Rt[p];
                    }
                    else
                    {
                        rt = file1Rt[p];
                    }

                    // add the scan
                    scans[p] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: p + 1, msnOrder: 1, isCentroid: true,
                        polarity: Polarity.Positive, retentionTime: rt, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                        mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (p + 1));
                }

                // write the .mzML
                Readers.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(scans),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, filesToWrite[f] + ".mzML"), false);
            }

            // set up spectra file info
            SpectraFileInfo file1 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, filesToWrite[0] + ".mzML"), "a", 0, 0, 0);
            SpectraFileInfo file2 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, filesToWrite[1] + ".mzML"), "a", 1, 0, 0);
            SpectraFileInfo file3 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, filesToWrite[2] + ".mzML"), "a", 2, 0, 0);

            // create some PSMs
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            Identification id1 = new Identification(file1, "PEPTIDE", "PEPTIDE",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDE").MonoisotopicMass, file1Rt[0] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id2 = new Identification(file1, "PEPTIDEV", "PEPTIDEV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEV").MonoisotopicMass, file1Rt[1] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id3 = new Identification(file1, "PEPTIDEVV", "PEPTIDEVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVV").MonoisotopicMass, file1Rt[2] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id4 = new Identification(file1, "PEPTIDEVVV", "PEPTIDEVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVVV").MonoisotopicMass, file1Rt[4] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id5 = new Identification(file1, "PEPTIDEVVVV", "PEPTIDEVVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVVVV").MonoisotopicMass, file1Rt[5] + 0.001, 1, new List<ProteinGroup> { pg });

            Identification id6 = new Identification(file2, "PEPTIDE", "PEPTIDE",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDE").MonoisotopicMass, file2Rt[0] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id7 = new Identification(file2, "PEPTIDEV", "PEPTIDEV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEV").MonoisotopicMass, file2Rt[1] + 0.001, 1, new List<ProteinGroup> { pg });
            // missing ID 8 - MBR feature - "PEPTIDEVV"

            Identification id9 = new Identification(file2, "PEPTIDEVVV", "PEPTIDEVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVVV").MonoisotopicMass, file2Rt[4] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id10 = new Identification(file2, "PEPTIDEVVVV", "PEPTIDEVVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVVVV").MonoisotopicMass, file2Rt[5] + 0.001, 1, new List<ProteinGroup> { pg });


            Identification id11 = new Identification(file3, "PEPTIDEV", "PEPTIDEV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEV").MonoisotopicMass, file1Rt[1] + 0.001, 1, new List<ProteinGroup> { pg }); // same as peak 2
            Identification id12 = new Identification(file3, "PEPTIDEVV", "PEPTIDEVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVV").MonoisotopicMass, file1Rt[2] + 0.001, 1, new List<ProteinGroup> { pg }); // same as peak 3, but higher intensity
            Identification id13 = new Identification(file3, "PEPTIDEVVV", "PEPTIDEVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVVV").MonoisotopicMass, file1Rt[4] + 0.001, 1, new List<ProteinGroup> { pg }); // same as peak 4


            // create the FlashLFQ engine
            FlashLfqEngine neighborsEngine = new FlashLfqEngine(new List<Identification> { id1, id2, id3, id4, id5, id6, id7, id9, id10, id11, id12, id13 }, 
                matchBetweenRuns: true, donorCriterion: DonorCriterion.Neighbors);

            //run the engine
            var results = neighborsEngine.Run();

            Assert.That(results.Peaks[file2].Count == 5);
            Assert.That(results.Peaks[file2].Count(p => p.DetectionType == DetectionType.MBR) == 1);

            var peak = results.Peaks[file2].First(p => p.DetectionType == DetectionType.MBR);
            var otherFilePeak = results.Peaks[file1].Where(p => p.Identifications.First().BaseSequence ==
                peak.Identifications.First().BaseSequence).First();


            Assert.That(peak.Intensity > 0);
            Assert.That(peak.Intensity == otherFilePeak.Intensity);
            Assert.That(peak.Identifications.First().FileInfo == file1); // assure that the ID came from file 1, ie, the donor with the most neighboring peaks

            // create the FlashLFQ engine
            FlashLfqEngine intensityEngine = new FlashLfqEngine(new List<Identification> { id1, id2, id3, id4, id5, id6, id7, id9, id10, id11, id12, id13 },
                matchBetweenRuns: true, donorCriterion: DonorCriterion.Intensity);

            //run the engine
            results = intensityEngine.Run();

            Assert.That(results.Peaks[file2].Count == 5);
            Assert.That(results.Peaks[file2].Count(p => p.DetectionType == DetectionType.MBR) == 1);

            peak = results.Peaks[file2].First(p => p.DetectionType == DetectionType.MBR);
            otherFilePeak = results.Peaks[file3].Where(p => p.Identifications.First().BaseSequence ==
                peak.Identifications.First().BaseSequence).First();


            Assert.That(peak.Intensity > 0);
            Assert.That(peak.Intensity, Is.EqualTo(otherFilePeak.Intensity/5).Within(1)); // file 3 is five times more intense than file 2
            Assert.That(peak.Identifications.First().FileInfo == file3); // assure that the ID came from file 3, ie, the most intense donor peaks

        }

    }
}
