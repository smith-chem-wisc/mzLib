using FlashLFQ;
using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test
{
    [TestFixture]
    internal class Test
    {
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
            FlashLFQEngine engine = new FlashLFQEngine(new List<Identification> { id1, id2, id3, id4 });

            // run the engine
            var results = engine.Run();

            // check raw results
            Assert.That(results.peaks[raw].Count == 1);
            Assert.That(results.peaks[raw].First().Intensity > 0);
            Assert.That(!results.peaks[raw].First().IsMbrFeature);
            Assert.That(results.peptideBaseSequences["EGFQVADGPLYR"].GetIntensity(raw) > 0);
            Assert.That(results.peptideModifiedSequences["EGFQVADGPLYR"].GetIntensity(raw) > 0);
            Assert.That(results.proteinGroups["MyProtein"].intensities[raw] > 0);

            // check mzml results
            Assert.That(results.peaks[mzml].Count == 1);
            Assert.That(results.peaks[mzml].First().Intensity > 0);
            Assert.That(!results.peaks[mzml].First().IsMbrFeature);
            Assert.That(results.peptideBaseSequences["EGFQVADGPLYR"].GetIntensity(mzml) > 0);
            Assert.That(results.peptideModifiedSequences["EGFQVADGPLYR"].GetIntensity(mzml) > 0);
            Assert.That(results.proteinGroups["MyProtein"].intensities[mzml] > 0);

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
            var results = new FlashLFQEngine(new List<Identification> { id1, id2 }).Run();

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

            results = new FlashLFQEngine(new List<Identification> { id1, id2 }).Run();

            int int3 = (int)System.Math.Round(results.peaks[mzml].First().Intensity, 0);
            int int4 = (int)System.Math.Round(results.peaks[raw].First().Intensity, 0);
            Assert.That(int3 > 0);
            Assert.That(int3 == int4);

            // ********************************* check techrep normalization *********************************
            raw = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-raw.raw"), "a", 0, 0, 0);
            mzml = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-mzml.mzml"), "a", 0, 1, 0);

            id1 = new Identification(raw, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            id2 = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });

            results = new FlashLFQEngine(new List<Identification> { id1, id2 }).Run();

            int int5 = (int)System.Math.Round(results.peaks[mzml].First().Intensity, 0);
            int int6 = (int)System.Math.Round(results.peaks[raw].First().Intensity, 0);
            Assert.That(int5 > 0);
            Assert.That(int5 == int6);

            Assert.That(int1 == int3);
            Assert.That(int1 == int5);
        }
    }
}