using Chemistry;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    internal class Test
    {
        [Test]
        public static void TestFlashLfqAdvancedProteinQuant()
        {
            List<string> filesToWrite = new List<string> { "mzml_1", "mzml_2" };
            List<string> pepSequences = new List<string> { "PEPTIDE", "MYPEPTIDE", "VVVVVPEPTIDE" };
            double[,] amounts = new double[2, 3] { { 1000000, 1000000, 1000000 },
                                                   { 2000000, 2000000, 900000 } };
            Loaders.LoadElements(Path.Combine(TestContext.CurrentContext.TestDirectory, @"elements.dat"));

            // generate mzml files (3 peptides each)
            for(int f = 0; f < filesToWrite.Count; f++)
            {
                // 1 MS1 scan per peptide
                MsDataScan[] scans = new MsDataScan[3];
                
                for(int p = 0; p < pepSequences.Count; p++)
                {
                    ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(pepSequences[p]).GetChemicalFormula();
                    IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);
                    double[] mz = dist.Masses.Select(v => v.ToMz(1)).ToArray();
                    double[] intensities = dist.Intensities.Select(v => v * amounts[f, p]).ToArray();

                    // add the scan
                    scans[p] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: p + 1, msnOrder: 1, isCentroid: true, 
                        polarity: Polarity.Positive, retentionTime: 1.0 + (p / 10.0), scanWindowRange: new MzRange(400, 1600), scanFilter: "f", 
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
            Identification id1 = new Identification(file1, "PEPTIDE", "PEPTIDE", 799.35996, 1.01, 1, new List<ProteinGroup> { pg });
            Identification id2 = new Identification(file1, "MYPEPTIDE", "MYPEPTIDE", 1093.46377, 1.11, 1, new List<ProteinGroup> { pg });
            Identification id3 = new Identification(file1, "VVVVVPEPTIDE", "VVVVVPEPTIDE", 1294.70203, 1.21, 1, new List<ProteinGroup> { pg });

            Identification id4 = new Identification(file2, "PEPTIDE", "PEPTIDE", 799.35996, 1.01, 1, new List<ProteinGroup> { pg });
            Identification id5 = new Identification(file2, "MYPEPTIDE", "MYPEPTIDE", 1093.46377, 1.11, 1, new List<ProteinGroup> { pg });
            Identification id6 = new Identification(file2, "VVVVVPEPTIDE", "VVVVVPEPTIDE", 1294.70203, 1.21, 1, new List<ProteinGroup> { pg });

            // create the FlashLFQ engine
            FlashLFQEngine engine = new FlashLFQEngine(new List<Identification> { id1, id2, id3, id4, id5, id6 }, normalize: false, advancedProteinQuant: true);

            // run the engine
            var results = engine.Run();

            // third peptide should be low-weighted
            // protein should be ~sum of first two peptide intensities (a little lower, because some smaller isotope peaks get skipped)
            double file1ProteinIntensity = results.proteinGroups["MyProtein"].GetIntensity(file1);
            Assert.That(file1ProteinIntensity < 2e6);
            Assert.That(file1ProteinIntensity > 1e6);
            
            double file2ProteinIntensity = results.proteinGroups["MyProtein"].GetIntensity(file2);
            Assert.That(file2ProteinIntensity < 4e6);
            Assert.That(file2ProteinIntensity > 3e6);
        }

        [Test]
        public static void TestFlashLFQ()
        {
            // get the raw file paths
            SpectraFileInfo raw = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-raw.raw"), "a", 0, 0, 0);
            SpectraFileInfo mzml = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-mzml.mzml"), "a", 0, 1, 0);

            // create some PSMs
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            Identification id1 = new Identification(raw, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            Identification id2 = new Identification(raw, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.05811, 2, new List<ProteinGroup> { pg });
            Identification id3 = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            Identification id4 = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.05811, 2, new List<ProteinGroup> { pg });

            // create the FlashLFQ engine
            FlashLFQEngine engine = new FlashLFQEngine(new List<Identification> { id1, id2, id3, id4 }, normalize: true);

            // run the engine
            var results = engine.Run();

            // check raw results
            Assert.That(results.peaks[raw].Count == 1);
            Assert.That(results.peaks[raw].First().Intensity > 0);
            Assert.That(!results.peaks[raw].First().IsMbrFeature);
            Assert.That(results.peptideBaseSequences["EGFQVADGPLYR"].GetIntensity(raw) > 0);
            Assert.That(results.peptideModifiedSequences["EGFQVADGPLYR"].GetIntensity(raw) > 0);
            Assert.That(results.proteinGroups["MyProtein"].GetIntensity(raw) > 0);

            // check mzml results
            Assert.That(results.peaks[mzml].Count == 1);
            Assert.That(results.peaks[mzml].First().Intensity > 0);
            Assert.That(!results.peaks[mzml].First().IsMbrFeature);
            Assert.That(results.peptideBaseSequences["EGFQVADGPLYR"].GetIntensity(mzml) > 0);
            Assert.That(results.peptideModifiedSequences["EGFQVADGPLYR"].GetIntensity(mzml) > 0);
            Assert.That(results.proteinGroups["MyProtein"].GetIntensity(mzml) > 0);

            // check that condition normalization worked
            int int1 = (int)System.Math.Round(results.peaks[mzml].First().Intensity, 0);
            int int2 = (int)System.Math.Round(results.peaks[raw].First().Intensity, 0);
            Assert.That(int1 == int2);

            // test peak output
            results.WriteResults(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"peaks.tsv"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"modSeq.tsv"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"baseSeq.tsv"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"protein.tsv"));
        }

        [Test]
        public static void TestFlashLFQNormalization()
        {
            // ********************************* check biorep normalization *********************************
            // get the raw file paths
            SpectraFileInfo raw = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-raw.raw"), "a", 0, 0, 0);
            SpectraFileInfo mzml = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-mzml.mzml"), "a", 1, 0, 0);

            // create some PSMs
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            Identification id1 = new Identification(raw, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            Identification id2 = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });

            // create the FlashLFQ engine
            var results = new FlashLFQEngine(new List<Identification> { id1, id2 }, normalize: true).Run();

            // check that biorep normalization worked
            int int1 = (int)System.Math.Round(results.peaks[mzml].First().Intensity, 0);
            int int2 = (int)System.Math.Round(results.peaks[raw].First().Intensity, 0);
            Assert.That(int1 > 0);
            Assert.That(int1 == int2);

            // ********************************* check condition normalization *********************************
            raw = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-raw.raw"), "a", 0, 0, 0);
            mzml = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-mzml.mzml"), "b", 0, 0, 0);

            id1 = new Identification(raw, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            id2 = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });

            results = new FlashLFQEngine(new List<Identification> { id1, id2 }, normalize: true).Run();

            int int3 = (int)System.Math.Round(results.peaks[mzml].First().Intensity, 0);
            int int4 = (int)System.Math.Round(results.peaks[raw].First().Intensity, 0);
            Assert.That(int3 > 0);
            Assert.That(int3 == int4);

            // ********************************* check techrep normalization *********************************
            raw = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-raw.raw"), "a", 0, 0, 0);
            mzml = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-mzml.mzml"), "a", 0, 1, 0);

            id1 = new Identification(raw, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            id2 = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });

            results = new FlashLFQEngine(new List<Identification> { id1, id2 }, normalize: true).Run();

            int int5 = (int)System.Math.Round(results.peaks[mzml].First().Intensity, 0);
            int int6 = (int)System.Math.Round(results.peaks[raw].First().Intensity, 0);
            Assert.That(int5 > 0);
            Assert.That(int5 == int6);

            Assert.That(int1 == int3);
            Assert.That(int1 == int5);


            // ********************************* check fraction normalization *********************************
            raw = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-raw.raw"), "a", 0, 0, 0);
            var raw2 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-raw.raw"), "a", 0, 0, 1);
            mzml = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-mzml.mzml"), "a", 1, 0, 0);
            var mzml2 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-mzml.mzml"), "a", 1, 0, 1);

            id1 = new Identification(raw, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            id2 = new Identification(raw2, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            var id3 = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            var id4 = new Identification(mzml2, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });

            results = new FlashLFQEngine(new List<Identification> { id1, id2, id3, id4 }, normalize: true).Run();

            int int7 = (int)System.Math.Round(results.peptideBaseSequences["EGFQVADGPLYR"].GetIntensity(raw) + results.peptideBaseSequences["EGFQVADGPLYR"].GetIntensity(raw2));
            int int8 = (int)System.Math.Round(results.peptideBaseSequences["EGFQVADGPLYR"].GetIntensity(mzml) + results.peptideBaseSequences["EGFQVADGPLYR"].GetIntensity(mzml2));
            Assert.That(int7 > 0);
            Assert.That(int7 == int8);
        }

        [Test]
        public static void MergeFlashLFQResults()
        {
            SpectraFileInfo rawA = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-raw.raw"), "a", 0, 0, 0);
            SpectraFileInfo mzmlA = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-mzml.mzml"), "a", 0, 1, 0);

            // create some PSMs
            var pgA = new ProteinGroup("MyProtein", "gene", "org");
            Identification id1A = new Identification(rawA, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pgA });
            Identification id2A = new Identification(rawA, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.05811, 2, new List<ProteinGroup> { pgA });
            Identification id3A = new Identification(mzmlA, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pgA });
            Identification id4A = new Identification(mzmlA, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.05811, 2, new List<ProteinGroup> { pgA });

            // create the FlashLFQ engine
            FlashLFQEngine engineA = new FlashLFQEngine(new List<Identification> { id1A, id2A, id3A, id4A });

            // run the engine
            var resultsA = engineA.Run();

            SpectraFileInfo rawB = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-raw.raw"), "b", 0, 0, 0);
            SpectraFileInfo mzmlB = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-mzml.mzml"), "b", 0, 1, 0);

            // create some PSMs
            var pgB = new ProteinGroup("MyProtein", "gene", "org");
            Identification id1 = new Identification(rawB, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pgB });
            Identification id2 = new Identification(rawB, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.05811, 2, new List<ProteinGroup> { pgB });
            Identification id3 = new Identification(mzmlB, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pgB });
            Identification id4 = new Identification(mzmlB, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.05811, 2, new List<ProteinGroup> { pgB });

            // create the FlashLFQ engine
            FlashLFQEngine engineB = new FlashLFQEngine(new List<Identification> { id1, id2, id3, id4 });

            // run the engine
            var resultsB = engineB.Run();

            resultsA.MergeResultsWith(resultsB);
            Assert.AreEqual(4, resultsA.peaks.Count);
            Assert.AreEqual(1, resultsA.peptideBaseSequences.Count);
            Assert.AreEqual(1, resultsA.peptideModifiedSequences.Count);
            Assert.AreEqual(1, resultsA.proteinGroups.Count);
            Assert.AreEqual(4, resultsA.spectraFiles.Count);
        }
    }
}