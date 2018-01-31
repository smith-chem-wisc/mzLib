using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;
using FlashLFQ;
using IO.Thermo;
using System.IO;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    class Test
    {
        [Test]
        public static void TestEverything()
        {
            Console.WriteLine("UNIT TEST - Entering unit test");
            string elements = Path.Combine(TestContext.CurrentContext.TestDirectory, "elements.dat");
            string files = TestContext.CurrentContext.TestDirectory;
            string ident = Path.Combine(TestContext.CurrentContext.TestDirectory, "aggregatePSMs_5ppmAroundZero.psmtsv");
            
            FlashLFQEngine engine = new FlashLFQEngine();
            Console.WriteLine("UNIT TEST - About to load elements");
            Loaders.LoadElements(elements);
            Console.WriteLine("UNIT TEST - Finished loading elements");

            Assert.That(engine.ParseArgs(new string[] {
                        "--idt " + ident,
                        "--rep " + files,
                        "--ppm 5",
                        "--sil false",
                        "--pau false",
                        "--mbr true" }
                    ));
            Console.WriteLine("UNIT TEST - Done making engine");
            engine.globalStopwatch.Start();
            Assert.That(engine.outputFolder != null);
            engine.SetParallelization(1);

            //Assert.That(engine.ReadPeriodicTable());

            Console.WriteLine("UNIT TEST - About to read TSV file");
            Assert.That(engine.ReadIdentificationsFromTSV());
            Console.WriteLine("UNIT TEST - Finished reading TSV");
            engine.ConstructIndexTemplateFromIdentifications();
            Console.WriteLine("UNIT TEST - Finished constructing bins");
            Assert.That(engine.observedMzsToUseForIndex.Count > 0);
            Assert.That(engine.baseSequenceToIsotopicDistribution.Count > 0);
            Console.WriteLine("UNIT TEST - Bins are OK");

            for (int i = 0; i < engine.filePaths.Length; i++)
            {
                Console.WriteLine("UNIT TEST - Quantifying file " + (i + 1));
                try
                {
                    Assert.That(engine.Quantify(null, engine.filePaths[i]));
                }
                catch (AssertionException)
                {
                    Console.WriteLine("UNIT TEST - Could not quantify file \"" + engine.filePaths[i] + "\"");
                }
            }

            //if (engine.mbr)
            //    engine.RetentionTimeCalibrationAndErrorCheckMatchedFeatures();

            Console.WriteLine("UNIT TEST - Quantifying proteins ");
            engine.QuantifyProteins();

            Console.WriteLine("UNIT TEST - Asserting results");
            Assert.That(engine.SumFeatures(engine.allFeaturesByFile.SelectMany(p => p), true).Any());
            Assert.That(engine.SumFeatures(engine.allFeaturesByFile.SelectMany(p => p), false).Any());

            Assert.That(engine.allFeaturesByFile[0].First().intensity > 0);
            Assert.That(engine.allFeaturesByFile[1].First().intensity > 0);

            Assert.That(engine.allFeaturesByFile[0].Count == 1);
            Assert.That(engine.allFeaturesByFile[1].Count == 1);

            Assert.That(!engine.allFeaturesByFile[0].First().isMbrFeature);
            Assert.That(!engine.allFeaturesByFile[1].First().isMbrFeature);
            Console.WriteLine("UNIT TEST - All passed");
        }

        [Test]
        public static void TestExternalNoPassedFile()
        {
            Console.WriteLine("UNIT TEST - Entering unit test");
            string[] filePaths = Directory.GetFiles(TestContext.CurrentContext.TestDirectory).Where(f => f.Substring(f.IndexOf('.')).ToUpper().Equals(".RAW") || f.Substring(f.IndexOf('.')).ToUpper().Equals(".MZML")).ToArray();
            string elements = Path.Combine(TestContext.CurrentContext.TestDirectory, "elements.dat");
            string ident = Path.Combine(TestContext.CurrentContext.TestDirectory, "aggregatePSMs_5ppmAroundZero.psmtsv");

            FlashLFQEngine engine = new FlashLFQEngine();
            Console.WriteLine("UNIT TEST - About to load elements");
            Loaders.LoadElements(elements);
            Console.WriteLine("UNIT TEST - Finished loading elements");

            engine.PassFilePaths(filePaths);
            Assert.That(engine.ParseArgs(new string[] {
                        "--ppm 5",
                        "--sil false",
                        "--pau false",
                        "--mbr true" }
                    ));
            Console.WriteLine("UNIT TEST - Done making engine");
            engine.globalStopwatch.Start();
            engine.SetParallelization(1);
            
            Console.WriteLine("UNIT TEST - Adding identifications");
            var ids = File.ReadAllLines(ident);
            int lineCount = 1;
            foreach(var line in ids)
            {
                if(lineCount != 1)
                {
                    var splitLine = line.Split('\t');
                    engine.AddIdentification(Path.GetFileNameWithoutExtension(splitLine[0]), splitLine[20], splitLine[21], double.Parse(splitLine[27]), double.Parse(splitLine[2]), (int) double.Parse(splitLine[6]), new List<string> { splitLine[14] });
                }
                lineCount++;
            }
            Console.WriteLine("UNIT TEST - Finished adding IDs");

            engine.ConstructIndexTemplateFromIdentifications();
            Console.WriteLine("UNIT TEST - Finished constructing bins");
            Assert.That(engine.observedMzsToUseForIndex.Count > 0);
            Assert.That(engine.baseSequenceToIsotopicDistribution.Count > 0);
            Console.WriteLine("UNIT TEST - Bins are OK");

            for (int i = 0; i < engine.filePaths.Length; i++)
            {
                Console.WriteLine("UNIT TEST - Quantifying file " + (i + 1));
                try
                {
                    Assert.That(engine.Quantify(null, engine.filePaths[i]));
                }
                catch (AssertionException)
                {
                    Console.WriteLine("UNIT TEST - Could not quantify file \"" + engine.filePaths[i] + "\"");
                }
            }

            //if (engine.mbr)
            //    engine.RetentionTimeCalibrationAndErrorCheckMatchedFeatures();

            Console.WriteLine("UNIT TEST - Quantifying proteins ");
            engine.QuantifyProteins();

            Console.WriteLine("UNIT TEST - Asserting results");
            Assert.That(engine.SumFeatures(engine.allFeaturesByFile.SelectMany(p => p), true).Any());
            Assert.That(engine.SumFeatures(engine.allFeaturesByFile.SelectMany(p => p), false).Any());

            Assert.That(engine.allFeaturesByFile[0].First().intensity > 0);
            Assert.That(engine.allFeaturesByFile[1].First().intensity > 0);

            Assert.That(!engine.allFeaturesByFile[0].First().isMbrFeature);
            Assert.That(!engine.allFeaturesByFile[1].First().isMbrFeature);
            Console.WriteLine("UNIT TEST - All passed");
        }

        [Test]
        public static void TestExternalPassedFile()
        {
            Console.WriteLine("UNIT TEST - Entering unit test");
            string[] filePaths = Directory.GetFiles(TestContext.CurrentContext.TestDirectory).Where(f => f.Substring(f.IndexOf('.')).ToUpper().Equals(".RAW")).ToArray();
            var thermoFile = ThermoStaticData.LoadAllStaticData(filePaths[0]);

            string elements = Path.Combine(TestContext.CurrentContext.TestDirectory, "elements.dat");
            string ident = Path.Combine(TestContext.CurrentContext.TestDirectory, "aggregatePSMs_5ppmAroundZero.psmtsv");

            FlashLFQEngine engine = new FlashLFQEngine();
            Console.WriteLine("UNIT TEST - About to load elements");
            Loaders.LoadElements(elements);
            Console.WriteLine("UNIT TEST - Finished loading elements");

            engine.PassFilePaths(filePaths);
            Assert.That(engine.ParseArgs(new string[] {
                        "--ppm 5",
                        "--sil false",
                        "--pau false",
                        "--mbr true" }
                    ));
            Console.WriteLine("UNIT TEST - Done making engine");
            engine.globalStopwatch.Start();
            engine.SetParallelization(1);

            Console.WriteLine("UNIT TEST - Adding identifications");
            var ids = File.ReadAllLines(ident);
            int lineCount = 1;
            foreach (var line in ids)
            {
                if (lineCount != 1)
                {
                    var splitLine = line.Split('\t');
                    engine.AddIdentification(Path.GetFileNameWithoutExtension(splitLine[0]), splitLine[20], splitLine[21], double.Parse(splitLine[27]), double.Parse(splitLine[2]), (int)double.Parse(splitLine[6]), new List<string> { splitLine[14] });
                }
                lineCount++;
            }
            Console.WriteLine("UNIT TEST - Finished adding IDs");

            engine.ConstructIndexTemplateFromIdentifications();
            Console.WriteLine("UNIT TEST - Finished constructing bins");
            Assert.That(engine.observedMzsToUseForIndex.Count > 0);
            Assert.That(engine.baseSequenceToIsotopicDistribution.Count > 0);
            Console.WriteLine("UNIT TEST - Bins are OK");

            for (int i = 0; i < engine.filePaths.Length; i++)
            {
                Console.WriteLine("UNIT TEST - Quantifying file " + (i + 1));
                try
                {
                    Assert.That(engine.Quantify(thermoFile, engine.filePaths[i]));
                }
                catch (AssertionException)
                {
                    Console.WriteLine("UNIT TEST - Could not quantify file \"" + engine.filePaths[i] + "\"");
                }
            }

            //if (engine.mbr)
            //    engine.RetentionTimeCalibrationAndErrorCheckMatchedFeatures();

            Console.WriteLine("UNIT TEST - Quantifying proteins ");
            engine.QuantifyProteins();

            Console.WriteLine("UNIT TEST - Asserting results");
            Assert.That(engine.SumFeatures(engine.allFeaturesByFile.SelectMany(p => p), true).Any());
            Assert.That(engine.SumFeatures(engine.allFeaturesByFile.SelectMany(p => p), false).Any());

            Assert.That(engine.allFeaturesByFile[0].First().intensity > 0);

            Assert.That(!engine.allFeaturesByFile[0].First().isMbrFeature);
            Console.WriteLine("UNIT TEST - All passed");
        }
    }
}
