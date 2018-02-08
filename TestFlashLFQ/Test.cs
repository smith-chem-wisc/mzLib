using FlashLFQ;
using IO.MzML;
using IO.Thermo;
using MassSpectrometry;
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

        #endregion Public Methods
    }
}