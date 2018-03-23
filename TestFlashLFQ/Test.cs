using FlashLFQ;
using IO.MzML;
using IO.Thermo;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    internal class Test
    {
        #region Public Methods

        [Test]
        public static void TestFlashLFQ()
        {
            // get the raw file paths
            RawFileInfo raw = new RawFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-raw.raw"));
            RawFileInfo mzml = new RawFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-mzml.mzml"));

            // create some PSMs
            Identification id1 = new Identification(raw, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<string> { "MyProtein" });
            Identification id2 = new Identification(raw, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.05811, 2, new List<string> { "MyProtein" });
            Identification id3 = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<string> { "MyProtein" });
            Identification id4 = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.05811, 2, new List<string> { "MyProtein" });

            // create the FlashLFQ engine
            FlashLFQEngine engine = new FlashLFQEngine(new List<Identification> { id1, id2, id3, id4 });

            // run the engine
            var results = engine.Run();

            // check raw results
            Assert.That(results.peaks[raw].Count == 1);
            Assert.That(results.peaks[raw].First().intensity > 0);
            Assert.That(!results.peaks[raw].First().isMbrFeature);
            Assert.That(results.peptideBaseSequences["EGFQVADGPLYR"].intensities[raw] > 0);
            Assert.That(results.peptideModifiedSequences["EGFQVADGPLYR"].intensities[raw] > 0);
            Assert.That(results.proteinGroups["MyProtein"].intensities[raw] > 0);

            // check mzml results
            Assert.That(results.peaks[mzml].Count == 1);
            Assert.That(results.peaks[mzml].First().intensity > 0);
            Assert.That(!results.peaks[mzml].First().isMbrFeature);
            Assert.That(results.peptideBaseSequences["EGFQVADGPLYR"].intensities[mzml] > 0);
            Assert.That(results.peptideModifiedSequences["EGFQVADGPLYR"].intensities[mzml] > 0);
            Assert.That(results.proteinGroups["MyProtein"].intensities[mzml] > 0);

            // test peak output
            List<string> output = new List<string>() { FlashLFQ.ChromatographicPeak.TabSeparatedHeader };
            foreach (var peak in results.peaks.SelectMany(p => p.Value))
                output.Add(peak.ToString());
            Assert.That(output.Count == 3);

            // test peptide base sequence output
            output = new List<string>() { Peptide.TabSeparatedHeader };
            foreach (var pep in results.peptideBaseSequences)
                output.Add(pep.Value.ToString());
            Assert.That(output.Count == 2);

            // test peptide mod sequence output
            output = new List<string>() { Peptide.TabSeparatedHeader };
            foreach (var pep in results.peptideModifiedSequences)
                output.Add(pep.Value.ToString());
            Assert.That(output.Count == 2);

            // test protein output
            output = new List<string>() { ProteinGroup.TabSeparatedHeader };
            foreach (var protein in results.proteinGroups)
                output.Add(protein.Value.ToString());
            Assert.That(output.Count == 2);
        }

        [Test]
        public static void TestFlashLFQWithPassedFile()
        {
            // read periodic table - needed to open the raw files
            PeriodicTableLoader.Load(Path.Combine(TestContext.CurrentContext.TestDirectory, @"elements.dat"));

            // get the raw files
            string rawPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-raw.raw");
            string mzmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-mzml.mzml");

            var rawFile = ThermoDynamicData.InitiateDynamicConnection(rawPath);
            var mzmlFile = Mzml.LoadAllStaticData(mzmlPath);

            RawFileInfo raw = new RawFileInfo(rawPath, rawFile);
            RawFileInfo mzml = new RawFileInfo(mzmlPath, mzmlFile);

            // create some PSMs
            Identification id1 = new Identification(raw, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<string> { "MyProtein" });
            Identification id2 = new Identification(raw, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.05811, 2, new List<string> { "MyProtein" });
            Identification id3 = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<string> { "MyProtein" });
            Identification id4 = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.05811, 2, new List<string> { "MyProtein" });

            // create the FlashLFQ engine
            FlashLFQEngine engine = new FlashLFQEngine(new List<Identification> { id1, id2, id3, id4 });

            // run the engine
            var results = engine.Run();

            // check raw results
            Assert.That(results.peaks[raw].Count == 1);
            Assert.That(results.peaks[raw].First().intensity > 0);
            Assert.That(!results.peaks[raw].First().isMbrFeature);
            Assert.That(results.peptideBaseSequences["EGFQVADGPLYR"].intensities[raw] > 0);
            Assert.That(results.peptideModifiedSequences["EGFQVADGPLYR"].intensities[raw] > 0);
            Assert.That(results.proteinGroups["MyProtein"].intensities[raw] > 0);

            // check mzml results
            Assert.That(results.peaks[mzml].Count == 1);
            Assert.That(results.peaks[mzml].First().intensity > 0);
            Assert.That(!results.peaks[mzml].First().isMbrFeature);
            Assert.That(results.peptideBaseSequences["EGFQVADGPLYR"].intensities[mzml] > 0);
            Assert.That(results.peptideModifiedSequences["EGFQVADGPLYR"].intensities[mzml] > 0);
            Assert.That(results.proteinGroups["MyProtein"].intensities[mzml] > 0);

            // test peak output
            List<string> output = new List<string>() { FlashLFQ.ChromatographicPeak.TabSeparatedHeader };
            foreach (var peak in results.peaks.SelectMany(p => p.Value))
                output.Add(peak.ToString());
            Assert.That(output.Count == 3);

            // test peptide base sequence output
            output = new List<string>() { Peptide.TabSeparatedHeader };
            foreach (var pep in results.peptideBaseSequences)
                output.Add(pep.Value.ToString());
            Assert.That(output.Count == 2);

            // test peptide mod sequence output
            output = new List<string>() { Peptide.TabSeparatedHeader };
            foreach (var pep in results.peptideModifiedSequences)
                output.Add(pep.Value.ToString());
            Assert.That(output.Count == 2);

            // test protein output
            output = new List<string>() { ProteinGroup.TabSeparatedHeader };
            foreach (var protein in results.proteinGroups)
                output.Add(protein.Value.ToString());
            Assert.That(output.Count == 2);
        }

        [Test]
        public static void NotARealTest()
        {
            // get the raw file paths
            RawFileInfo raw1 = new RawFileInfo(@"C:\Data\ionstarSample\B02_06_161103_A1_HCD_OT_4ul.raw");
            RawFileInfo raw2 = new RawFileInfo(@"C:\Data\ionstarSample\B02_07_161103_A2_HCD_OT_4ul.raw");

            List<RawFileInfo> rawfiles = new List<RawFileInfo> { raw1, raw2 };
            // create some PSMs
            var myFile = File.ReadAllLines(@"C:\Data\ionstarSample\msms.txt");
            List<Identification> ids = new List<Identification>();
            int lineNum = 0;
            int fileNameCol = 0;
            int baseSequCol = 0;
            int fullSequCol = 0;
            int monoMassCol = 0;
            int msmsRetnCol = 0;
            int chargeStCol = 0;
            int protNameCol = 0;

            foreach(var line in myFile)
            {
                var param = line.Split('\t');

                if (lineNum == 0)
                {
                    fileNameCol = Array.IndexOf(param, "Raw file");
                    baseSequCol = Array.IndexOf(param, "Sequence");
                    fullSequCol = Array.IndexOf(param, "Modified sequence");
                    monoMassCol = Array.IndexOf(param, "Mass");
                    msmsRetnCol = Array.IndexOf(param, "Retention time");
                    chargeStCol = Array.IndexOf(param, "Charge");
                    protNameCol = Array.IndexOf(param, "Proteins");
                }
                else
                {
                    string fileName = param[fileNameCol];
                    string BaseSequence = param[baseSequCol];
                    string ModSequence = param[fullSequCol];
                    double monoisotopicMass = double.Parse(param[monoMassCol]);
                    double ms2RetentionTime = double.Parse(param[msmsRetnCol]);
                    int chargeState = (int)double.Parse(param[chargeStCol]);

                    List<string> proteinGroups = new List<string>();
                    var g = param[protNameCol].Split(new char[] { '\t' }, StringSplitOptions.RemoveEmptyEntries);
                    if (g.Any())
                        proteinGroups.Add(g.First().Trim());

                    // construct id
                    var fileNameNoExt = Path.GetFileNameWithoutExtension(fileName);
                    var rawFileInfoToUse = rawfiles.Where(p => p.filenameWithoutExtension.Equals(fileNameNoExt)).FirstOrDefault();
                    if (rawFileInfoToUse == null)
                    {
                        // skip PSMs for files with no spectrum data input
                        lineNum++;
                        continue;
                    }

                    var ident = new Identification(rawFileInfoToUse, BaseSequence, ModSequence, monoisotopicMass, ms2RetentionTime, chargeState, proteinGroups);
                    ids.Add(ident);
                }

                lineNum++;
            }

            // create the FlashLFQ engine
            FlashLFQEngine engine = new FlashLFQEngine(ids, 5, 5, true);

            // run the engine
            var results = engine.Run();

            string myPeakOutput = @"C:\Users\rmillikin\Desktop\myPeaks.tsv";
            string myPeptideOutput = @"C:\Users\rmillikin\Desktop\myPeptides.tsv";
            string myPeptideModOutput = @"C:\Users\rmillikin\Desktop\myPeptideMods.tsv";
            string myProteinOutput = @"C:\Users\rmillikin\Desktop\myProteins.tsv";

            // test peak output
            List<string> output = new List<string>() { FlashLFQ.ChromatographicPeak.TabSeparatedHeader };
            foreach (var peak in results.peaks.SelectMany(p => p.Value))
                output.Add(peak.ToString());
            File.WriteAllLines(myPeakOutput, output);

            // test peptide base sequence output
            output = new List<string>() { Peptide.TabSeparatedHeader };
            foreach (var pep in results.peptideBaseSequences)
                output.Add(pep.Value.ToString());
            File.WriteAllLines(myPeptideOutput, output);

            // test peptide mod sequence output
            output = new List<string>() { Peptide.TabSeparatedHeader };
            foreach (var pep in results.peptideModifiedSequences)
                output.Add(pep.Value.ToString());
            File.WriteAllLines(myPeptideModOutput, output);

            // test protein output
            output = new List<string>() { ProteinGroup.TabSeparatedHeader };
            foreach (var protein in results.proteinGroups)
                output.Add(protein.Value.ToString());
            File.WriteAllLines(myProteinOutput, output);
        }

        [Test]
        public static void TestFlashLFQNormalization()
        {
            double[] peptide1 = new double[] { 960, 768, 800 };
            double[] peptide2 = new double[] { 2093, 1674, 1700 };

            double averagePep1 = (peptide1[0] + peptide1[1] + peptide1[2]) / 3.0;
            double averagePep2 = (peptide2[0] + peptide2[1] + peptide2[2]) / 3.0;

            double[] diffFromAverage1 = new double[] { (peptide1[0] - averagePep1) / averagePep1, (peptide1[1] - averagePep1) / averagePep1 };
            double[] diffFromAverage2 = new double[] { (peptide2[0] - averagePep2) / averagePep2, (peptide2[1] - averagePep2) / averagePep2 };
            


            double[] normalizedPep1 = new double[] { peptide1[0] * (1 - diffFromAverage1[0]), peptide1[1] * (1 - diffFromAverage1[1]) };
            double[] normalizedPep2 = new double[] { peptide2[0] * (1 - diffFromAverage2[0]), peptide2[1] * (1 - diffFromAverage2[1]) };

            double error = Math.Abs(normalizedPep1[0] - normalizedPep1[1]) + Math.Abs(normalizedPep2[0] - normalizedPep2[1]);
        }

        #endregion Public Methods
    }
}