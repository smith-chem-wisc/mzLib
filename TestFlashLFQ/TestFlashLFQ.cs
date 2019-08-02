using Chemistry;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;
using ChromatographicPeak = FlashLFQ.ChromatographicPeak;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test
{
    [TestFixture]
    internal class TestFlashLFQ
    {
        private static Stopwatch Stopwatch { get; set; }

        [SetUp]
        public static void Setuppp()
        {
            Stopwatch = new Stopwatch();
            Stopwatch.Start();
        }

        [TearDown]
        public static void TearDown()
        {
            Console.WriteLine($"Analysis time: {Stopwatch.Elapsed.Hours}h {Stopwatch.Elapsed.Minutes}m {Stopwatch.Elapsed.Seconds}s");
        }

        [Test]
        public static void TestFlashLfq()
        {
            // get the raw file paths
            SpectraFileInfo raw = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-raw.raw"), "a", 0, 0, 0);
            SpectraFileInfo mzml = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-mzml.mzml"), "a", 0, 1, 0);

            // create some PSMs
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            Identification id1 = new Identification(raw, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            Identification id2 = new Identification(raw, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.05811, 2, new List<ProteinGroup> { pg });
            Identification id3 = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            Identification id4 = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.05811, 2, new List<ProteinGroup> { pg });

            // create the FlashLFQ engine
            FlashLfqEngine engine = new FlashLfqEngine(new List<Identification> { id1, id2, id3, id4 }, normalize: true);

            // run the engine
            var results = engine.Run();

            // check raw results
            Assert.That(results.Peaks[raw].Count == 1);
            Assert.That(results.Peaks[raw].First().Intensity > 0);
            Assert.That(!results.Peaks[raw].First().IsMbrPeak);
            Assert.That(results.PeptideModifiedSequences["EGFQVADGPLYR"].GetIntensity(raw) > 0);
            Assert.That(results.ProteinGroups["MyProtein"].GetIntensity(raw) > 0);

            // check mzml results
            Assert.That(results.Peaks[mzml].Count == 1);
            Assert.That(results.Peaks[mzml].First().Intensity > 0);
            Assert.That(!results.Peaks[mzml].First().IsMbrPeak);
            Assert.That(results.PeptideModifiedSequences["EGFQVADGPLYR"].GetIntensity(mzml) > 0);
            Assert.That(results.ProteinGroups["MyProtein"].GetIntensity(mzml) > 0);

            // check that condition normalization worked
            int int1 = (int)System.Math.Round(results.Peaks[mzml].First().Intensity, 0);
            int int2 = (int)System.Math.Round(results.Peaks[raw].First().Intensity, 0);
            Assert.That(int1 == int2);

            // test peak output
            results.WriteResults(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"peaks.tsv"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"modSeq.tsv"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"protein.tsv"),
                null);
        }

        [Test]
        public static void TestEvelopQuantification()
        {
            Loaders.LoadElements();

            double monoIsotopicMass = 1350.65681;
            double massOfAveragine = 111.1254;
            double numberOfAveragines = monoIsotopicMass / massOfAveragine;

            double averageC = 4.9384 * numberOfAveragines;
            double averageH = 7.7583 * numberOfAveragines;
            double averageO = 1.4773 * numberOfAveragines;
            double averageN = 1.3577 * numberOfAveragines;
            double averageS = 0.0417 * numberOfAveragines;

            ChemicalFormula myFormula = ChemicalFormula.ParseFormula(
                "C" + (int)Math.Round(averageC) +
                "H" + (int)Math.Round(averageH) +
                "O" + (int)Math.Round(averageO) +
                "N" + (int)Math.Round(averageN) +
                "S" + (int)Math.Round(averageS));


            // get the raw file paths
            SpectraFileInfo mzml = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-mzml.mzml"), "a", 0, 0, 0);

            // create some PSMs
            var pg = new ProteinGroup("MyProtein", "gene", "org");

            Identification id3 = new Identification(mzml, "", "1", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg }, myFormula);
            Identification id4 = new Identification(mzml, "", "2", 1350.65681, 94.05811, 2, new List<ProteinGroup> { pg }, myFormula);

            // create the FlashLFQ engine
            FlashLfqEngine engine = new FlashLfqEngine(new List<Identification> { id3, id4 }, normalize: true);

            // run the engine
            var results = engine.Run();

            Assert.IsTrue(results.Peaks.First().Value.First().Intensity > 0);
        }

        [Test]
        public static void TestFlashLfqNormalization()
        {
            // ********************************* check biorep normalization *********************************
            // get the raw file paths
            SpectraFileInfo raw = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-raw.raw"), "a", 0, 0, 0);
            SpectraFileInfo mzml = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-mzml.mzml"), "a", 1, 0, 0);

            // create some PSMs
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            Identification id1 = new Identification(raw, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            Identification id2 = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });

            // create the FlashLFQ engine
            var results = new FlashLfqEngine(new List<Identification> { id1, id2 }, normalize: true).Run();

            // check that biorep normalization worked
            int int1 = (int)System.Math.Round(results.Peaks[mzml].First().Intensity, 0);
            int int2 = (int)System.Math.Round(results.Peaks[raw].First().Intensity, 0);
            Assert.That(int1 > 0);
            Assert.That(int1 == int2);

            // ********************************* check condition normalization *********************************
            raw = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-raw.raw"), "a", 0, 0, 0);
            mzml = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-mzml.mzml"), "b", 0, 0, 0);

            id1 = new Identification(raw, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            id2 = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });

            results = new FlashLfqEngine(new List<Identification> { id1, id2 }, normalize: true).Run();

            int int3 = (int)System.Math.Round(results.Peaks[mzml].First().Intensity, 0);
            int int4 = (int)System.Math.Round(results.Peaks[raw].First().Intensity, 0);
            Assert.That(int3 > 0);
            Assert.That(int3 == int4);

            // ********************************* check techrep normalization *********************************
            raw = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-raw.raw"), "a", 0, 0, 0);
            mzml = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-mzml.mzml"), "a", 0, 1, 0);

            id1 = new Identification(raw, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            id2 = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });

            results = new FlashLfqEngine(new List<Identification> { id1, id2 }, normalize: true).Run();

            int int5 = (int)System.Math.Round(results.Peaks[mzml].First().Intensity, 0);
            int int6 = (int)System.Math.Round(results.Peaks[raw].First().Intensity, 0);
            Assert.That(int5 > 0);
            Assert.That(int5 == int6);

            Assert.That(int1 == int3);
            Assert.That(int1 == int5);

            // ********************************* check fraction normalization *********************************
            raw = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-raw.raw"), "a", 0, 0, 0);
            var raw2 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-raw.raw"), "a", 0, 0, 1);
            mzml = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-mzml.mzml"), "a", 1, 0, 0);
            var mzml2 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-mzml.mzml"), "a", 1, 0, 1);

            id1 = new Identification(raw, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            id2 = new Identification(raw2, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            var id3 = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            var id4 = new Identification(mzml2, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });

            results = new FlashLfqEngine(new List<Identification> { id1, id2, id3, id4 }, normalize: true).Run();

            int int7 = (int)System.Math.Round(results.PeptideModifiedSequences["EGFQVADGPLYR"].GetIntensity(raw) + results.PeptideModifiedSequences["EGFQVADGPLYR"].GetIntensity(raw2));
            int int8 = (int)System.Math.Round(results.PeptideModifiedSequences["EGFQVADGPLYR"].GetIntensity(mzml) + results.PeptideModifiedSequences["EGFQVADGPLYR"].GetIntensity(mzml2));
            Assert.That(int7 > 0);
            Assert.That(int7 == int8);
        }

        [Test]
        public static void TestFlashLfqMergeResults()
        {
            SpectraFileInfo rawA = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-raw.raw"), "a", 0, 0, 0);
            SpectraFileInfo mzmlA = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-mzml.mzml"), "a", 0, 1, 0);

            // create some PSMs
            var pgA = new ProteinGroup("MyProtein", "gene", "org");
            Identification id1A = new Identification(rawA, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pgA });
            Identification id2A = new Identification(rawA, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.05811, 2, new List<ProteinGroup> { pgA });
            Identification id3A = new Identification(mzmlA, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pgA });
            Identification id4A = new Identification(mzmlA, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.05811, 2, new List<ProteinGroup> { pgA });

            // create the FlashLFQ engine
            FlashLfqEngine engineA = new FlashLfqEngine(new List<Identification> { id1A, id2A, id3A, id4A });

            // run the engine
            var resultsA = engineA.Run();

            SpectraFileInfo rawB = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-raw.raw"), "b", 0, 0, 0);
            SpectraFileInfo mzmlB = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-mzml.mzml"), "b", 0, 1, 0);

            // create some PSMs
            var pgB = new ProteinGroup("MyProtein", "gene", "org");
            Identification id1 = new Identification(rawB, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pgB });
            Identification id2 = new Identification(rawB, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.05811, 2, new List<ProteinGroup> { pgB });
            Identification id3 = new Identification(mzmlB, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pgB });
            Identification id4 = new Identification(mzmlB, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.05811, 2, new List<ProteinGroup> { pgB });

            // create the FlashLFQ engine
            FlashLfqEngine engineB = new FlashLfqEngine(new List<Identification> { id1, id2, id3, id4 });

            // run the engine
            var resultsB = engineB.Run();

            resultsA.MergeResultsWith(resultsB);
            Assert.AreEqual(4, resultsA.Peaks.Count);
            Assert.AreEqual(1, resultsA.PeptideModifiedSequences.Count);
            Assert.AreEqual(1, resultsA.ProteinGroups.Count);
            Assert.AreEqual(4, resultsA.SpectraFiles.Count);
        }

        [Test]
        public static void TestFlashLfqMatchBetweenRuns()
        {
            List<string> filesToWrite = new List<string> { "mzml_1", "mzml_2" };
            List<string> pepSequences = new List<string> { "PEPTIDE", "PEPTIDEV", "PEPTIDEVV", "PEPTIDEVVV", "PEPTIDEVVVV" };
            double intensity = 1e6;

            double[] file1Rt = new double[] { 1.01, 1.02, 1.03, 1.04, 1.05 };
            double[] file2Rt = new double[] { 1.00, 1.025, 1.04, 1.055, 1.070 };

            Loaders.LoadElements();

            // generate mzml files (5 peptides each)
            for (int f = 0; f < filesToWrite.Count; f++)
            {
                // 1 MS1 scan per peptide
                MsDataScan[] scans = new MsDataScan[5];

                for (int p = 0; p < pepSequences.Count; p++)
                {
                    ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(pepSequences[p]).GetChemicalFormula();
                    IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);
                    double[] mz = dist.Masses.Select(v => v.ToMz(1)).ToArray();
                    double[] intensities = dist.Intensities.Select(v => v * intensity).ToArray();
                    double rt;
                    if (f == 0)
                    {
                        rt = file1Rt[p];
                    }
                    else
                    {
                        rt = file2Rt[p];
                    }

                    // add the scan
                    scans[p] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: p + 1, msnOrder: 1, isCentroid: true,
                        polarity: Polarity.Positive, retentionTime: rt, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                        mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (p + 1));
                }

                // write the .mzML
                IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(scans),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, filesToWrite[f] + ".mzML"), false);
            }

            // set up spectra file info
            SpectraFileInfo file1 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, filesToWrite[0] + ".mzML"), "a", 0, 0, 0);
            SpectraFileInfo file2 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, filesToWrite[1] + ".mzML"), "a", 1, 0, 0);

            // create some PSMs
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            Identification id1 = new Identification(file1, "PEPTIDE", "PEPTIDE",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDE").MonoisotopicMass, file1Rt[0] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id2 = new Identification(file1, "PEPTIDEV", "PEPTIDEV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEV").MonoisotopicMass, file1Rt[1] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id3 = new Identification(file1, "PEPTIDEVV", "PEPTIDEVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVV").MonoisotopicMass, file1Rt[2] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id4 = new Identification(file1, "PEPTIDEVVV", "PEPTIDEVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVVV").MonoisotopicMass, file1Rt[3] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id5 = new Identification(file1, "PEPTIDEVVVV", "PEPTIDEVVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVVVV").MonoisotopicMass, file1Rt[4] + 0.001, 1, new List<ProteinGroup> { pg });

            Identification id6 = new Identification(file2, "PEPTIDE", "PEPTIDE",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDE").MonoisotopicMass, file2Rt[0] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id7 = new Identification(file2, "PEPTIDEV", "PEPTIDEV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEV").MonoisotopicMass, file2Rt[1] + 0.001, 1, new List<ProteinGroup> { pg });
            // missing ID 8 - MBR feature
            Identification id9 = new Identification(file2, "PEPTIDEVVV", "PEPTIDEVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVVV").MonoisotopicMass, file2Rt[3] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id10 = new Identification(file2, "PEPTIDEVVVV", "PEPTIDEVVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVVVV").MonoisotopicMass, file2Rt[4] + 0.001, 1, new List<ProteinGroup> { pg });

            // create the FlashLFQ engine
            FlashLfqEngine engine = new FlashLfqEngine(new List<Identification> { id1, id2, id3, id4, id5, id6, id7, id9, id10 }, matchBetweenRuns: true);

            // run the engine
            var results = engine.Run();

            Assert.That(results.Peaks[file2].Count == 5);
            Assert.That(results.Peaks[file2].Where(p => p.IsMbrPeak).Count() == 1);

            var peak = results.Peaks[file2].Where(p => p.IsMbrPeak).First();
            var otherFilePeak = results.Peaks[file1].Where(p => p.Identifications.First().BaseSequence ==
                peak.Identifications.First().BaseSequence).First();

            Assert.That(peak.Intensity > 0);
            Assert.That(peak.Intensity == otherFilePeak.Intensity);

            Assert.That(results.Peaks[file1].Count == 5);
            Assert.That(results.Peaks[file1].Where(p => p.IsMbrPeak).Count() == 0);

            Assert.That(results.ProteinGroups["MyProtein"].GetIntensity(file1) > 0);
            Assert.That(results.ProteinGroups["MyProtein"].GetIntensity(file2) > 0);
        }

        [Test]
        public static void TestPeakSplittingLeft()
        {
            string fileToWrite = "myMzml.mzML";
            string peptide = "PEPTIDE";
            double intensity = 1e6;

            Loaders.LoadElements();

            // generate mzml file

            // 1 MS1 scan per peptide
            MsDataScan[] scans = new MsDataScan[10];
            double[] intensityMultipliers = { 1, 3, 1, 1, 3, 5, 10, 5, 3, 1 };

            for (int s = 0; s < scans.Length; s++)
            {
                ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(peptide).GetChemicalFormula();
                IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);
                double[] mz = dist.Masses.Select(v => v.ToMz(1)).ToArray();
                double[] intensities = dist.Intensities.Select(v => v * intensity * intensityMultipliers[s]).ToArray();

                // add the scan
                scans[s] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: s + 1, msnOrder: 1, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: 1.0 + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }

            // write the .mzML
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(scans),
                Path.Combine(TestContext.CurrentContext.TestDirectory, fileToWrite), false);

            // set up spectra file info
            SpectraFileInfo file1 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, fileToWrite), "", 0, 0, 0);

            // create some PSMs
            var pg = new ProteinGroup("MyProtein", "gene", "org");

            Identification id1 = new Identification(file1, peptide, peptide,
                new Proteomics.AminoAcidPolymer.Peptide(peptide).MonoisotopicMass, 1.7 + 0.001, 1, new List<ProteinGroup> { pg });

            // create the FlashLFQ engine
            FlashLfqEngine engine = new FlashLfqEngine(new List<Identification> { id1 });

            // run the engine
            var results = engine.Run();
            ChromatographicPeak peak = results.Peaks.First().Value.First();

            Assert.That(peak.Apex.IndexedPeak.RetentionTime == 1.6);
            Assert.That(peak.SplitRT == 1.3);
            Assert.That(!peak.IsotopicEnvelopes.Any(p => p.IndexedPeak.RetentionTime < 1.3));
            Assert.That(peak.IsotopicEnvelopes.Count == 6);
        }

        [Test]
        public static void TestPeakSplittingRight()
        {
            string fileToWrite = "myMzml.mzML";
            string peptide = "PEPTIDE";
            double intensity = 1e6;

            Loaders.LoadElements();

            // generate mzml file

            // 1 MS1 scan per peptide
            MsDataScan[] scans = new MsDataScan[10];
            double[] intensityMultipliers = { 1, 3, 5, 10, 5, 3, 1, 1, 3, 1 };

            for (int s = 0; s < scans.Length; s++)
            {
                ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(peptide).GetChemicalFormula();
                IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);
                double[] mz = dist.Masses.Select(v => v.ToMz(1)).ToArray();
                double[] intensities = dist.Intensities.Select(v => v * intensity * intensityMultipliers[s]).ToArray();

                // add the scan
                scans[s] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: s + 1, msnOrder: 1, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: 1.0 + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }

            // write the .mzML
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(scans),
                Path.Combine(TestContext.CurrentContext.TestDirectory, fileToWrite), false);

            // set up spectra file info
            SpectraFileInfo file1 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, fileToWrite), "", 0, 0, 0);

            // create some PSMs
            var pg = new ProteinGroup("MyProtein", "gene", "org");

            Identification id1 = new Identification(file1, peptide, peptide,
                new Proteomics.AminoAcidPolymer.Peptide(peptide).MonoisotopicMass, 1.3 + 0.001, 1, new List<ProteinGroup> { pg });

            // create the FlashLFQ engine
            FlashLfqEngine engine = new FlashLfqEngine(new List<Identification> { id1 });

            // run the engine
            var results = engine.Run();
            ChromatographicPeak peak = results.Peaks.First().Value.First();

            Assert.That(peak.Apex.IndexedPeak.RetentionTime == 1.3);
            Assert.That(peak.SplitRT == 1.6);
            Assert.That(!peak.IsotopicEnvelopes.Any(p => p.IndexedPeak.RetentionTime > 1.6));
            Assert.That(peak.IsotopicEnvelopes.Count == 6);
        }

        [Test]
        public static void TestPeakSplittingRightWithEmptyScan()
        {
            string fileToWrite = "myMzml.mzML";
            string peptide = "PEPTIDE";
            double intensity = 1e6;

            Loaders.LoadElements();

            // generate mzml file

            // 1 MS1 scan per peptide
            MsDataScan[] scans = new MsDataScan[10];
            double[] intensityMultipliers = { 1, 3, 5, 10, 5, 3, 1, 1, 3, 1 };

            for (int s = 0; s < scans.Length; s++)
            {
                ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(peptide).GetChemicalFormula();
                IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);
                double[] mz = dist.Masses.Select(v => v.ToMz(1)).ToArray();
                double[] intensities = dist.Intensities.Select(v => v * intensity * intensityMultipliers[s]).ToArray();

                if (s == 7)
                {
                    mz = new[] { 401.0 };
                    intensities = new[] { 1000.0 };
                }

                // add the scan
                scans[s] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: s + 1, msnOrder: 1, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: 1.0 + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }

            // write the .mzML
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(scans),
                Path.Combine(TestContext.CurrentContext.TestDirectory, fileToWrite), false);

            // set up spectra file info
            SpectraFileInfo file1 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, fileToWrite), "", 0, 0, 0);

            // create some PSMs
            var pg = new ProteinGroup("MyProtein", "gene", "org");

            Identification id1 = new Identification(file1, peptide, peptide,
                new Proteomics.AminoAcidPolymer.Peptide(peptide).MonoisotopicMass, 1.3 + 0.001, 1, new List<ProteinGroup> { pg });

            // create the FlashLFQ engine
            FlashLfqEngine engine = new FlashLfqEngine(new List<Identification> { id1 });

            // run the engine
            var results = engine.Run();
            ChromatographicPeak peak = results.Peaks.First().Value.First();

            Assert.That(peak.Apex.IndexedPeak.RetentionTime == 1.3);
            Assert.That(peak.SplitRT == 1.6);
            Assert.That(!peak.IsotopicEnvelopes.Any(p => p.IndexedPeak.RetentionTime > 1.6));
            Assert.That(peak.IsotopicEnvelopes.Count == 6);
        }

        [Test]
        public static void TestPeakSplittingLeftWithEmptyScan()
        {
            string fileToWrite = "myMzml.mzML";
            string peptide = "PEPTIDE";
            double intensity = 1e6;

            Loaders.LoadElements();

            // generate mzml file

            // 1 MS1 scan per peptide
            MsDataScan[] scans = new MsDataScan[10];
            double[] intensityMultipliers = { 1, 3, 1, 1, 3, 5, 10, 5, 3, 1 };

            for (int s = 0; s < scans.Length; s++)
            {
                ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(peptide).GetChemicalFormula();
                IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);
                double[] mz = dist.Masses.Select(v => v.ToMz(1)).ToArray();
                double[] intensities = dist.Intensities.Select(v => v * intensity * intensityMultipliers[s]).ToArray();

                if (s == 2)
                {
                    mz = new[] { 401.0 };
                    intensities = new[] { 1000.0 };
                }

                // add the scan
                scans[s] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: s + 1, msnOrder: 1, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: 1.0 + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }

            // write the .mzML
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(scans),
                Path.Combine(TestContext.CurrentContext.TestDirectory, fileToWrite), false);

            // set up spectra file info
            SpectraFileInfo file1 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, fileToWrite), "", 0, 0, 0);

            // create some PSMs
            var pg = new ProteinGroup("MyProtein", "gene", "org");

            Identification id1 = new Identification(file1, peptide, peptide,
                new Proteomics.AminoAcidPolymer.Peptide(peptide).MonoisotopicMass, 1.3 + 0.001, 1, new List<ProteinGroup> { pg });

            // create the FlashLFQ engine
            FlashLfqEngine engine = new FlashLfqEngine(new List<Identification> { id1 });

            // run the engine
            var results = engine.Run();
            ChromatographicPeak peak = results.Peaks.First().Value.First();

            Assert.That(peak.Apex.IndexedPeak.RetentionTime == 1.6);
            Assert.That(peak.SplitRT == 1.3);
            Assert.That(!peak.IsotopicEnvelopes.Any(p => p.IndexedPeak.RetentionTime < 1.3));
            Assert.That(peak.IsotopicEnvelopes.Count == 6);
        }

        [Test]
        public static void TestToString()
        {
            // many of these are just to check that the ToString methods don't cause crashes
            var indexedPeak = new IndexedMassSpectralPeak(1.0, 2.0, 4, 5.0);
            Assert.That(indexedPeak.ToString().Equals("1.000; 4"));

            var spectraFile = new SpectraFileInfo("myFullPath", "", 0, 0, 0);
            string spectraString = spectraFile.ToString();

            var proteinGroup = new ProteinGroup("Accession", "Gene", "Organism");
            string pgString = proteinGroup.ToString(new List<SpectraFileInfo> { spectraFile });

            var identification = new Identification(
                spectraFile, "PEPTIDE", "PEPTIDE", 1.0, 2.0, 3,
                new List<ProteinGroup> { proteinGroup });
            string idString = identification.ToString();

            var chromPeak = new ChromatographicPeak(identification, false, spectraFile);
            string chromPeakString = chromPeak.ToString();
            chromPeak.CalculateIntensityForThisFeature(true);
            string peakAfterCalculatingIntensity = chromPeak.ToString();
        }

        [Test]
        public static void TestNotFound()
        {
            FlashLFQ.Peptide p = new FlashLFQ.Peptide("Seq", true, new HashSet<ProteinGroup>());
            var notFound = p.GetDetectionType(new SpectraFileInfo("", "", 0, 0, 0));
            Assert.That(notFound == DetectionType.NotDetected);
        }
        
        [Test]
        public static void TestAmbiguous()
        {
            // get the raw file paths
            SpectraFileInfo mzml = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-mzml.mzml"), "a", 0, 1, 0);

            // create some PSMs
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            Identification id3 = new Identification(mzml, "EGFQVADGPLRY", "EGFQVADGPLRY", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            Identification id4 = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.05811, 2, new List<ProteinGroup> { pg });

            // create the FlashLFQ engine
            FlashLfqEngine engine = new FlashLfqEngine(new List<Identification> { id3, id4 });

            // run the engine
            var results = engine.Run();

            Assert.That(results.Peaks[mzml].Count == 1);
            Assert.That(results.Peaks[mzml].First().Intensity > 0);
            Assert.That(!results.Peaks[mzml].First().IsMbrPeak);
            Assert.That(results.Peaks[mzml].First().NumIdentificationsByFullSeq == 2);
            Assert.That(results.PeptideModifiedSequences["EGFQVADGPLYR"].GetIntensity(mzml) == 0);
            Assert.That(results.PeptideModifiedSequences["EGFQVADGPLRY"].GetIntensity(mzml) == 0);
            Assert.That(results.ProteinGroups["MyProtein"].GetIntensity(mzml) == 0);

            // test peak output
            results.WriteResults(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"peaks.tsv"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"modSeq.tsv"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"protein.tsv"),
                null);
        }
        
        [Test]
        public static void TestMatchBetweenRunsWithNoIdsInCommon()
        {
            List<string> filesToWrite = new List<string> { "mzml_1", "mzml_2" };
            List<string> pepSequences = new List<string> { "PEPTIDE", "PEPTIDEV", "PEPTIDEVV", "PEPTIDEVVV", "PEPTIDEVVVV" };
            double intensity = 1e6;

            double[] file1Rt = new double[] { 1.01, 1.02, 1.03, 1.04, 1.05 };
            double[] file2Rt = new double[] { 1.015, 1.030, 1.036, 1.050, 1.065 };

            Loaders.LoadElements();

            // generate mzml files (5 peptides each)
            for (int f = 0; f < filesToWrite.Count; f++)
            {
                // 1 MS1 scan per peptide
                MsDataScan[] scans = new MsDataScan[5];

                for (int p = 0; p < pepSequences.Count; p++)
                {
                    ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(pepSequences[p]).GetChemicalFormula();
                    IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);
                    double[] mz = dist.Masses.Select(v => v.ToMz(1)).ToArray();
                    double[] intensities = dist.Intensities.Select(v => v * intensity).ToArray();
                    double rt;
                    if (f == 0)
                    {
                        rt = file1Rt[p];
                    }
                    else
                    {
                        rt = file2Rt[p];
                    }

                    // add the scan
                    scans[p] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: p + 1, msnOrder: 1, isCentroid: true,
                        polarity: Polarity.Positive, retentionTime: rt, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                        mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (p + 1));
                }

                // write the .mzML
                IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(scans),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, filesToWrite[f] + ".mzML"), false);
            }

            // set up spectra file info
            SpectraFileInfo file1 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, filesToWrite[0] + ".mzML"), "a", 0, 0, 0);
            SpectraFileInfo file2 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, filesToWrite[1] + ".mzML"), "a", 1, 0, 0);

            // create some PSMs
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            var myMbrProteinGroup = new ProteinGroup("MyMbrProtein", "MbrGene", "org");

            Identification id1 = new Identification(file1, "PEPTIDE", "PEPTIDE",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDE").MonoisotopicMass, file1Rt[0] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id2 = new Identification(file1, "PEPTIDEV", "PEPTIDEV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEV").MonoisotopicMass, file1Rt[1] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id3 = new Identification(file1, "PEPTIDEVV", "PEPTIDEVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVV").MonoisotopicMass, file1Rt[2] + 0.001, 1, new List<ProteinGroup> { myMbrProteinGroup });
            Identification id4 = new Identification(file1, "PEPTIDEVVV", "PEPTIDEVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVVV").MonoisotopicMass, file1Rt[3] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id5 = new Identification(file1, "PEPTIDEVVVV", "PEPTIDEVVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVVVV").MonoisotopicMass, file1Rt[4] + 0.001, 1, new List<ProteinGroup> { pg });

            Identification id6 = new Identification(file2, "PEPTIED", "PEPTIED",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIED").MonoisotopicMass, file2Rt[0] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id7 = new Identification(file2, "PEPTIEDV", "PEPTIEDV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIEDV").MonoisotopicMass, file2Rt[1] + 0.001, 1, new List<ProteinGroup> { pg });
            // missing ID 8 - MBR feature
            Identification id9 = new Identification(file2, "PEPTIEDVVV", "PEPTIEDVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIEDVVV").MonoisotopicMass, file2Rt[3] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id10 = new Identification(file2, "PEPTIEDVVVV", "PEPTIEDVVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIEDVVVV").MonoisotopicMass, file2Rt[4] + 0.001, 1, new List<ProteinGroup> { pg });

            FlashLfqEngine engine = new FlashLfqEngine(new List<Identification> { id1, id2, id3, id4, id5, id6, id7, id9, id10 }, matchBetweenRuns: true);
            var results = engine.Run();

            // no assertions - just don't crash
        }

        [Test]
        public static void TestFlashLfqDoesNotRemovePeptides()
        {
            Loaders.LoadElements();

            Residue x = new Residue("a", 'a', "a", Chemistry.ChemicalFormula.ParseFormula("C{13}6H12N{15}2O"), ModificationSites.All); //+8 lysine
            Residue lightLysine = Residue.GetResidue('K');

            Residue.AddNewResiduesToDictionary(new List<Residue> { new Residue("heavyLysine", 'a', "a", x.ThisChemicalFormula, ModificationSites.All) });

            SpectraFileInfo fileInfo = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, @"SilacTest.mzML"), "", 0, 0, 0);
            FlashLfqEngine engine = new FlashLfqEngine(
                new List<Identification>
                {
                    new Identification(fileInfo,"RDILSSNNQHGILPLSWNIPELVNMGQWK","RDILSSNNQHGILPLSWNIPELVNM[Common Variable:Oxidation on M]GQWK",3374.7193792,98.814005,3,new List<FlashLFQ.ProteinGroup>{new FlashLFQ.ProteinGroup("P01027","C3","Mus") },null, true),
                    new Identification(fileInfo,"RDILSSNNQHGILPLSWNIPELVNMGQWa","RDILSSNNQHGILPLSWNIPELVNM[Common Variable:Oxidation on M]GQWa",3382.733578,98.814005,3,new List<FlashLFQ.ProteinGroup>{new FlashLFQ.ProteinGroup("P01027+8.014","C3","Mus") },null, true),
                    new Identification(fileInfo,"RDILSSNNQHGILPLSWNIPELVNMGQWK","RDILSSNNQHGILPLSWNIPELVNM[Common Variable:Oxidation on M]GQWK",3374.7193792,98.7193782,4,new List<FlashLFQ.ProteinGroup>{new FlashLFQ.ProteinGroup("P01027","C3","Mus") },null, true),
                    new Identification(fileInfo,"RDILSSNNQHGILPLSWNIPELVNMGQWa","RDILSSNNQHGILPLSWNIPELVNM[Common Variable:Oxidation on M]GQWa",3382.733578,98.7193782,4,new List<FlashLFQ.ProteinGroup>{new FlashLFQ.ProteinGroup("P01027+8.014","C3","Mus") },null, true),
                },
                ppmTolerance: 5,
                silent: true,
                maxThreads: 7
                );
            var results = engine.Run();
            Assert.IsTrue(results.PeptideModifiedSequences.Count == 2);
        }

        [Test]
        public static void TestBayesianProteinQuantification()
        {
            // this mostly just tests that the Bayesian quant algorithm produces "reasonable"
            // estimates of the posterior error probability (PEP) and mean differences.

            // the idea here is to create a protein, assign some peptides with intensities to it,
            // and compute the probability that the protein is changing between the conditions "a" and "b".
            // the intensities in condition "b" are about double that of condition "a", and the results
            // of the Bayesian estimation should reflect that.
            ProteinGroup pg = new ProteinGroup("Accession", "Gene", "Organism");
            var p = new FlashLFQ.Peptide("PEPTIDE", true, new HashSet<ProteinGroup> { pg });

            var files = new List<SpectraFileInfo>
            {
                new SpectraFileInfo("a1", "a", 0, 0, 0),
                new SpectraFileInfo("a2", "a", 1, 0, 0),
                new SpectraFileInfo("a3", "a", 2, 0, 0),
                new SpectraFileInfo("b1", "b", 0, 0, 0),
                new SpectraFileInfo("b2", "b", 1, 0, 0),
                new SpectraFileInfo("b3", "b", 2, 0, 0)
            };

            var res = new FlashLfqResults(files, new List<Identification>
            {
                new Identification(null, "SEQUENCE", "SEQUENCE", 0, 0, 0, new List<ProteinGroup>{ new ProteinGroup("Accession", "Gene", "Organism") })
            });

            FlashLFQ.Peptide peptide = res.PeptideModifiedSequences.First().Value;
            ProteinGroup proteinGroup = res.ProteinGroups.First().Value;

            peptide.SetIntensity(files[0], 900);
            peptide.SetIntensity(files[1], 1000);
            peptide.SetIntensity(files[2], 1100);

            peptide.SetIntensity(files[3], 1950);
            peptide.SetIntensity(files[4], 2000);
            peptide.SetIntensity(files[5], 2050);

            var engine = new ProteinQuantificationEngine(res, maxThreads: 1, baseCondition: "a", randomSeed: 0);
            engine.Run();

            var quantResult = proteinGroup.ConditionToQuantificationResults["b"];

            Assert.That(Math.Round(quantResult.NullHypothesisCutoff.Value, 3) == 0.366);
            Assert.That(Math.Round(quantResult.PosteriorErrorProbability, 3) == 0.082);
            Assert.That(Math.Round(quantResult.FoldChangePointEstimate, 3) == 0.996);
            Assert.That(quantResult.PeptideFoldChangeMeasurements.Count == 1);
            Assert.That(quantResult.PeptideFoldChangeMeasurements.SelectMany(v => v.foldChanges).Count() == 3);

            string filepath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"bayesianProteinQuant.tsv");
            res.WriteResults(null, null, null, filepath);

            var textResults = File.ReadAllLines(filepath);
            Assert.That(textResults.Length == 2);
            var line = textResults[1].Split(new char[] { '\t' });
            Assert.That(Math.Round(double.Parse(line[11]), 3) == 0.082);
            File.Delete(filepath);

            // try with defined fold-change cutoff
            proteinGroup.ConditionToQuantificationResults.Clear();
            engine = new ProteinQuantificationEngine(res, maxThreads: 1, baseCondition: "a", randomSeed: 1, foldChangeCutoff: 0.8);
            engine.Run();

            quantResult = proteinGroup.ConditionToQuantificationResults["b"];

            Assert.That(Math.Round(quantResult.PosteriorErrorProbability, 3) == 0.179);
            Assert.That(Math.Round(quantResult.FoldChangePointEstimate, 3) == 1.014);
            Assert.That(quantResult.PeptideFoldChangeMeasurements.Count == 1);
            Assert.That(quantResult.PeptideFoldChangeMeasurements.SelectMany(v => v.foldChanges).Count() == 3);

            // try with some missing values
            peptide.SetIntensity(files[1], 0);
            peptide.SetIntensity(files[5], 0);

            proteinGroup.ConditionToQuantificationResults.Clear();
            engine = new ProteinQuantificationEngine(res, maxThreads: 1, baseCondition: "a", randomSeed: 2, foldChangeCutoff: 0.5);
            engine.Run();

            quantResult = proteinGroup.ConditionToQuantificationResults["b"];

            Assert.That(Math.Round(quantResult.PosteriorErrorProbability, 3) == 0.478);
            Assert.That(Math.Round(quantResult.FoldChangePointEstimate, 3) == 1.357);
            Assert.That(quantResult.PeptideFoldChangeMeasurements.Count == 1);
            Assert.That(quantResult.PeptideFoldChangeMeasurements.SelectMany(v => v.foldChanges).Count() == 2);

            // try with paired samples
            peptide.SetIntensity(files[0], 100);
            peptide.SetIntensity(files[1], 1000);
            peptide.SetIntensity(files[2], 10000);

            peptide.SetIntensity(files[3], 210);
            peptide.SetIntensity(files[4], 2200);
            peptide.SetIntensity(files[5], 21500);

            proteinGroup.ConditionToQuantificationResults.Clear();
            engine = new ProteinQuantificationEngine(res, maxThreads: 1, baseCondition: "a", randomSeed: 3, pairedSamples: true);
            engine.Run();

            quantResult = proteinGroup.ConditionToQuantificationResults["b"];

            Assert.That(Math.Round(quantResult.PosteriorErrorProbability, 3) == 0.077);
            Assert.That(Math.Round(quantResult.FoldChangePointEstimate, 3) == 1.110);
            Assert.That(quantResult.PeptideFoldChangeMeasurements.Count == 1);
            Assert.That(quantResult.PeptideFoldChangeMeasurements.SelectMany(v => v.foldChanges).Count() == 3);
        }

        [Test]
        public static void RealDataMbrTest()
        {
            string psmFile = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"PSMsForMbrTest.psmtsv");

            SpectraFileInfo f1r1 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"f1r1_sliced_mbr.raw"), "a", 0, 0, 0);
            SpectraFileInfo f1r2 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"f1r2_sliced_mbr.raw"), "a", 1, 0, 0);

            List<Identification> ids = new List<Identification>();
            foreach (string line in File.ReadAllLines(psmFile))
            {
                var split = line.Split(new char[] { '\t' });

                if (split.Contains("File Name") || string.IsNullOrWhiteSpace(line))
                {
                    continue;
                }

                SpectraFileInfo file = null;

                if (split[0].Contains("f1r1"))
                {
                    file = f1r1;
                }
                else if (split[0].Contains("f1r2"))
                {
                    file = f1r2;
                }

                string baseSequence = split[12];
                string fullSequence = split[13];
                double monoMass = double.Parse(split[21]);
                double rt = double.Parse(split[2]);
                int z = (int)double.Parse(split[6]);
                
                Identification id = new Identification(file, baseSequence, fullSequence, monoMass, rt, z, new List<ProteinGroup>());
                ids.Add(id);
            }

            var engine = new FlashLfqEngine(ids, matchBetweenRuns: true, maxThreads: 1);
            var results = engine.Run();
            
            var f1r1MbrResults = results.PeptideModifiedSequences
                .Where(p => p.Value.GetDetectionType(f1r1) == DetectionType.MSMS && p.Value.GetDetectionType(f1r2) == DetectionType.MBR).ToList();

            Assert.That(f1r1MbrResults.Count >= 82);
            var res1 = f1r1MbrResults.First(p => p.Key == "HVAVTNMNEHSSR").Value;
            double intensity1 = Math.Log(res1.GetIntensity(f1r1), 2);
            double intensity2 = Math.Log(res1.GetIntensity(f1r2), 2);
            Assert.That(Math.Abs(intensity1 - intensity2) < 0.5);

            var f1r2MbrResults = results
                .PeptideModifiedSequences
                .Where(p => p.Value.GetDetectionType(f1r1) == DetectionType.MBR && p.Value.GetDetectionType(f1r2) == DetectionType.MSMS).ToList();

            Assert.That(f1r2MbrResults.Count >= 142);
            res1 = f1r2MbrResults.First(p => p.Key == "YGIDKDKVDK").Value;
            intensity1 = Math.Log(res1.GetIntensity(f1r1), 2);
            intensity2 = Math.Log(res1.GetIntensity(f1r2), 2);
            Assert.That(Math.Abs(intensity1 - intensity2) < 0.5);
        }
    }
}