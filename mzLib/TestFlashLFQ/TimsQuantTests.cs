using NUnit.Framework;
using Readers;
using System.Collections.Generic;
using System.Linq;
using FlashLFQ;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System.IO;
using FlashLFQ.PEP;
using System;
using System.Windows.Documents;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Test.FileReadingTests;
using UsefulProteomicsDatabases;
using NUnit.Framework.Legacy;
using FlashLFQ.PeakIndexing;


namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]

    public class TimsQuantTests
    {

        [Test]
        public static void LocalDataTinyTest()
        {
            string testDataDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData");
            string outputDirectory = @"D:\timsTOF_Data_Bruker\ddaPASEF_data\MM_15_20_Tolerance\FlashLFQ_Test";
            Directory.CreateDirectory(outputDirectory);

            string psmFile = @"D:\timsTOF_Data_Bruker\ddaPASEF_data\MM_15_20_Tolerance\Task1-SearchTask\AllPSMs.psmtsv";

            SpectraFileInfo f1r1 = new SpectraFileInfo(@"D:\timsTOF_Data_Bruker\ddaPASEF_data\200ngHeLaPASEF_1min.d", "one", 1, 1, 1);
            //SpectraFileInfo f1r2 = new SpectraFileInfo(@"D:\timsTOF_Data_Bruker\ddaPASEF_data\200ngHeLaPASEF_1min-Copy.d", "one", 1, 1, 1);
            //SpectraFileInfo f1r2 = new SpectraFileInfo(@"D:\timsTOF_Data_Bruker\ddaPASEF_data\50ng_K562_extreme_3min.d", "two", 1, 1, 1);
            //SpectraFileInfo f1r3 = new SpectraFileInfo(@"D:\timsTOF_Data_Bruker\ddaPASEF_data\20230505_TIMS05_PaSk_MA_HeLa_6min_ddaP_S1-F2_1_2352.d", "three", 1, 1, 1);

            //List<string> acceptableProteinGroupAccessions = new() { "Q7KZF4", "Q15149", "P52298" };

            List<Identification> ids = new List<Identification>();
            Dictionary<string, ProteinGroup> allProteinGroups = new Dictionary<string, ProteinGroup>();
            foreach (string line in File.ReadAllLines(psmFile))
            {
                var split = line.Split(new char[] { '\t' });

                //skip the header
                if (split.Contains("File Name") || string.IsNullOrWhiteSpace(line))
                {
                    continue;
                }

                SpectraFileInfo file = null;

                if (split[0].Contains("HeLaPASEF_1min"))
                {
                    file = f1r1;
                }
                else continue;

                if (split[23].Contains("|") || double.Parse(split[56]) > 0.01)
                {
                    continue;
                }

                string baseSequence = split[13];
                string fullSequence = split[14];
                double monoMass = double.Parse(split[23]);
                double rt = double.Parse(split[2]) / 60;
                int z = (int)double.Parse(split[6]);
                var proteinSubset = split[26].Split(new char[] { '|' });
                string organism = split[29];
                string gene = split[27];

                double score = double.Parse(split[10]);
                string targetContamDecoy = split[39];
                List<ProteinGroup> proteinGroups = new();


                foreach (var protein in proteinSubset)
                {
                    if (allProteinGroups.TryGetValue(protein, out var proteinGroup))
                    {
                        proteinGroups.Add(proteinGroup);
                    }
                    else
                    {
                        allProteinGroups.Add(protein, new ProteinGroup(protein, gene, organism));
                        proteinGroups.Add(allProteinGroups[protein]);
                    }
                }

                Identification id = new Identification(file, baseSequence, fullSequence, monoMass, rt, z, proteinGroups,
                    psmScore: score, decoy: targetContamDecoy.Contains("D"));
                ids.Add(id);
            }

            var engine = new FlashLfqEngine(ids,
                matchBetweenRuns: true,
                requireMsmsIdInCondition: false,
                useSharedPeptidesForProteinQuant: true,
                ppmTolerance: 15,
                isotopeTolerancePpm: 15,
                maxThreads: 5);
            var results = engine.Run();

            results.WriteResults(Path.Combine(outputDirectory, "peaks.tsv"), Path.Combine(outputDirectory, "peptides.tsv"), Path.Combine(outputDirectory, "proteins.tsv"), Path.Combine(outputDirectory, "bayesian.tsv"), true);

            var peaks = results.Peaks.Values.ToList();
            var peptides = results.PeptideModifiedSequences.Values.ToList();
            var proteins = results.ProteinGroups.Values.ToList();

            Assert.Pass();

        }

        [Test]
        public static void LocalDataFakeMBRTest()
        {
            string testDataDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData");
            string outputDirectory = @"D:\timsTOF_Data_Bruker\ddaPASEF_data\MM_15_20_Tolerance\FlashLFQ_Test";
            Directory.CreateDirectory(outputDirectory);

            string psmFile = @"D:\timsTOF_Data_Bruker\ddaPASEF_data\MM_15_20_Tolerance\Task1-SearchTask\AllPSMs.psmtsv";

            SpectraFileInfo f1r1 = new SpectraFileInfo(@"D:\timsTOF_Data_Bruker\ddaPASEF_data\200ngHeLaPASEF_1min.d", "one", 1, 1, 1);
            SpectraFileInfo f1r2 = new SpectraFileInfo(@"D:\timsTOF_Data_Bruker\ddaPASEF_data\200ngHeLaPASEF_1min-Copy.d", "one", 2, 1, 1);
            //SpectraFileInfo f1r2 = new SpectraFileInfo(@"D:\timsTOF_Data_Bruker\ddaPASEF_data\50ng_K562_extreme_3min.d", "two", 1, 1, 1);
            //SpectraFileInfo f1r3 = new SpectraFileInfo(@"D:\timsTOF_Data_Bruker\ddaPASEF_data\20230505_TIMS05_PaSk_MA_HeLa_6min_ddaP_S1-F2_1_2352.d", "three", 1, 1, 1);

            //List<string> acceptableProteinGroupAccessions = new() { "Q7KZF4", "Q15149", "P52298" };

            List<Identification> ids = new List<Identification>();
            Dictionary<string, ProteinGroup> allProteinGroups = new Dictionary<string, ProteinGroup>();
            int iterator = 0;
            foreach (string line in File.ReadAllLines(psmFile))
            {
                var split = line.Split(new char[] { '\t' });

                //skip the header
                if (split.Contains("File Name") || string.IsNullOrWhiteSpace(line))
                {
                    continue;
                }

                SpectraFileInfo file = null;

                if (split[0].Contains("HeLaPASEF_1min"))
                {
                    file = f1r1;
                }
                else continue;

                if (split[23].Contains("|") || double.Parse(split[56]) > 0.01)
                {
                    continue;
                }

                string baseSequence = split[13];
                string fullSequence = split[14];
                double monoMass = double.Parse(split[23]);
                double rt = double.Parse(split[2]) / 60;
                int z = (int)double.Parse(split[6]);
                var proteinSubset = split[26].Split(new char[] { '|' });
                string organism = split[29];
                string gene = split[27];

                double score = double.Parse(split[10]);
                string targetContamDecoy = split[39];
                List<ProteinGroup> proteinGroups = new();


                foreach (var protein in proteinSubset)
                {
                    if (allProteinGroups.TryGetValue(protein, out var proteinGroup))
                    {
                        proteinGroups.Add(proteinGroup);
                    }
                    else
                    {
                        allProteinGroups.Add(protein, new ProteinGroup(protein, gene, organism));
                        proteinGroups.Add(allProteinGroups[protein]);
                    }
                }

                Identification id = new Identification(file, baseSequence, fullSequence, monoMass, rt, z, proteinGroups,
                    psmScore: score, decoy: targetContamDecoy.Contains("D"));
                ids.Add(id);
                iterator++;
                if(iterator % 5 > 0)
                {
                    id = new Identification(f1r2, baseSequence, fullSequence, monoMass, rt, z, proteinGroups,
                        psmScore: score, decoy: targetContamDecoy.Contains("D"));
                    ids.Add(id);
                }
            }

            var engine = new FlashLfqEngine(ids,
                matchBetweenRuns: true,
                requireMsmsIdInCondition: false,
                useSharedPeptidesForProteinQuant: true,
                ppmTolerance: 15,
                isotopeTolerancePpm: 15,
                maxThreads: 20);
            var results = engine.Run();

            results.WriteResults(Path.Combine(outputDirectory, "peaks.tsv"), Path.Combine(outputDirectory, "peptides.tsv"), Path.Combine(outputDirectory, "proteins.tsv"), Path.Combine(outputDirectory, "bayesian.tsv"), true);

            var peaks = results.Peaks.Values.ToList();
            var peptides = results.PeptideModifiedSequences.Values.ToList();
            var proteins = results.ProteinGroups.Values.ToList();

            Assert.Pass();

        }

        [Test]
        public static void LocalDataSmallTest()
        {
            string testDataDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData");
            string outputDirectory = @"D:\timsTOF_Data_Bruker\ddaPASEF_data\MM_15_20_Tolerance\FlashLFQ_Test";
            Directory.CreateDirectory(outputDirectory);

            string psmFile = @"D:\timsTOF_Data_Bruker\ddaPASEF_data\MM_15_20_Tolerance\Task1-SearchTask\AllPSMs.psmtsv";

            //SpectraFileInfo f1r1 = new SpectraFileInfo(@"D:\timsTOF_Data_Bruker\ddaPASEF_data\200ngHeLaPASEF_1min.d", "one", 1, 1, 1);
            //SpectraFileInfo f1r2 = new SpectraFileInfo(@"D:\timsTOF_Data_Bruker\ddaPASEF_data\50ng_K562_extreme_3min.d", "two", 1, 1, 1);
            SpectraFileInfo f1r3 = new SpectraFileInfo(@"D:\timsTOF_Data_Bruker\ddaPASEF_data\20230505_TIMS05_PaSk_MA_HeLa_6min_ddaP_S1-F2_1_2352.d", "three", 1, 1, 1);

            //List<string> acceptableProteinGroupAccessions = new() { "Q7KZF4", "Q15149", "P52298" };

            List<Identification> ids = new List<Identification>();
            Dictionary<string, ProteinGroup> allProteinGroups = new Dictionary<string, ProteinGroup>();
            foreach (string line in File.ReadAllLines(psmFile))
            {
                var split = line.Split(new char[] { '\t' });

                //skip the header
                if (split.Contains("File Name") || string.IsNullOrWhiteSpace(line))
                {
                    continue;
                }

                SpectraFileInfo file = null;

                //if (split[0].Contains("HeLaPASEF_1min"))
                //{
                //    file = f1r1;
                //}
                //else if (split[0].Contains("_extreme_3min"))
                //{
                //    file = f1r2;
                //}
                //else 
                if (split[0].Contains("_6min_ddaP"))
                {
                    file = f1r3;
                }
                else
                {
                    continue;
                }

                if (split[23].Contains("|") || double.Parse(split[56]) > 0.01)
                {
                    continue;
                }

                string baseSequence = split[13];
                string fullSequence = split[14];
                double monoMass = double.Parse(split[23]);
                double rt = double.Parse(split[2]) / 60;
                int z = (int)double.Parse(split[6]);
                var proteinSubset = split[26].Split(new char[] { '|' });
                string organism = split[29];
                string gene = split[27];

                double score = double.Parse(split[10]);
                string targetContamDecoy = split[39];
                List<ProteinGroup> proteinGroups = new();


                foreach (var protein in proteinSubset)
                {
                    if (allProteinGroups.TryGetValue(protein, out var proteinGroup))
                    {
                        proteinGroups.Add(proteinGroup);
                    }
                    else
                    {
                        allProteinGroups.Add(protein, new ProteinGroup(protein, gene, organism));
                        proteinGroups.Add(allProteinGroups[protein]);
                    }
                }

                Identification id = new Identification(file, baseSequence, fullSequence, monoMass, rt, z, proteinGroups, 
                    psmScore: score, decoy: targetContamDecoy.Contains("D"));
                ids.Add(id);
            }

            var engine = new FlashLfqEngine(ids,
                matchBetweenRuns: false,
                requireMsmsIdInCondition: false,
                useSharedPeptidesForProteinQuant: true,
                maxThreads: -1);
            var results = engine.Run();

            results.WriteResults(Path.Combine(outputDirectory, "peaks.tsv"), Path.Combine(outputDirectory, "peptides.tsv"), Path.Combine(outputDirectory, "proteins.tsv"), Path.Combine(outputDirectory, "bayesian.tsv"), true);

            var peaks = results.Peaks.Values.ToList();
            var peptides = results.PeptideModifiedSequences.Values.ToList();
            var proteins = results.ProteinGroups.Values.ToList();

            Assert.Pass();

        }

        [Test]
        public static void LocalDataBigTest()
        {
            string testDataDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData");
            string outputDirectory = @"D:\PXD014777_timsTOF_spikeIn\FlashLFQ_LocalDataBigTest";
            Directory.CreateDirectory(outputDirectory);

            string psmFile = @"D:\PXD014777_timsTOF_spikeIn\MMSearch_2_18_25\Task1-SearchTask\AllPSMs.psmtsv";

            SpectraFileInfo f1r1 = new SpectraFileInfo(@"D:\PXD014777_timsTOF_spikeIn\20180809_120min_200ms_WEHI25_brute20k_timsON_100ng_HYE124A_Slot1-7_1_890.d", 
                "A", 1, 1, 1);
            SpectraFileInfo f1r2 = new SpectraFileInfo(@"D:\PXD014777_timsTOF_spikeIn\20180809_120min_200ms_WEHI25_brute20k_timsON_100ng_HYE124A_Slot1-7_1_891.d", "A", 2, 1, 1);
            SpectraFileInfo f1r3 = new SpectraFileInfo(@"D:\PXD014777_timsTOF_spikeIn\20180809_120min_200ms_WEHI25_brute20k_timsON_100ng_HYE124A_Slot1-7_1_892.d", "A", 3, 1, 1);

            SpectraFileInfo g1r1 = new SpectraFileInfo(@"D:\PXD014777_timsTOF_spikeIn\20180809_120min_200ms_WEHI25_brute20k_timsON_100ng_HYE124B_Slot1-8_1_893.d", "B", 1, 1, 1);
            SpectraFileInfo g1r2 = new SpectraFileInfo(@"D:\PXD014777_timsTOF_spikeIn\20180809_120min_200ms_WEHI25_brute20k_timsON_100ng_HYE124B_Slot1-8_1_894.d", "B", 2, 1, 1);
            SpectraFileInfo g1r3 = new SpectraFileInfo(@"D:\PXD014777_timsTOF_spikeIn\20180809_120min_200ms_WEHI25_brute20k_timsON_100ng_HYE124B_Slot1-8_1_895.d", "B", 3, 1, 1);

            //List<string> acceptableProteinGroupAccessions = new() { "Q7KZF4", "Q15149", "P52298" };

            List<Identification> ids = new List<Identification>();
            Dictionary<string, ProteinGroup> allProteinGroups = new Dictionary<string, ProteinGroup>();
            foreach (string line in File.ReadAllLines(psmFile))
            {
                var split = line.Split(new char[] { '\t' });

                //skip the header
                if (split.Contains("File Name") || string.IsNullOrWhiteSpace(line))
                {
                    continue;
                }

                SpectraFileInfo file = null;

                if (split[0].Contains("_890"))
                {
                    file = f1r1;
                }
                else if (split[0].Contains("_893"))
                {
                    file = g1r1;
                }
                else
                    continue;

                //else if (split[0].Contains("_891"))
                //{
                //    file = f1r2;
                //}
                //else if (split[0].Contains("_892"))
                //{
                //    file = f1r3;
                //}
                
                //else if (split[0].Contains("_894"))
                //{
                //    file = g1r2;
                //}
                //else if (split[0].Contains("_895"))
                //{
                //    file = g1r3;
                //}

                if (split[23].Contains("|") || double.Parse(split[56]) > 0.01)
                {
                    continue;
                }

                string baseSequence = split[13];
                string fullSequence = split[14];
                double monoMass = double.Parse(split[23]);
                double rt = double.Parse(split[2]) / 60;
                int z = (int)double.Parse(split[6]);
                var proteinSubset = split[26].Split(new char[] { '|' });
                string organism = split[29];
                string gene = split[27];

                double score = double.Parse(split[10]);
                string targetContamDecoy = split[39];
                List<ProteinGroup> proteinGroups = new();


                foreach (var protein in proteinSubset)
                {
                    if (allProteinGroups.TryGetValue(protein, out var proteinGroup))
                    {
                        proteinGroups.Add(proteinGroup);
                    }
                    else
                    {
                        allProteinGroups.Add(protein, new ProteinGroup(protein, gene, organism));
                        proteinGroups.Add(allProteinGroups[protein]);
                    }
                }

                Identification id = new Identification(file, baseSequence, fullSequence, monoMass, rt, z, proteinGroups,
                    psmScore: score, decoy: targetContamDecoy.Contains("D"));
                ids.Add(id);
            }

            var engine = new FlashLfqEngine(ids,
                matchBetweenRuns: true,
                ppmTolerance: 15,
                isotopeTolerancePpm: 15,
                matchBetweenRunsFdrThreshold: 0.01,
                requireMsmsIdInCondition: false,
                useSharedPeptidesForProteinQuant: true,
                maxThreads: 20);
            var results = engine.Run();

            results.WriteResults(Path.Combine(outputDirectory, "peaks.tsv"), Path.Combine(outputDirectory, "peptides.tsv"), Path.Combine(outputDirectory, "proteins.tsv"), Path.Combine(outputDirectory, "bayesian.tsv"), true);

            var peaks = results.Peaks.Values.ToList();
            var peptides = results.PeptideModifiedSequences.Values.ToList();
            var proteins = results.ProteinGroups.Values.ToList();

            Assert.Pass();

        }

        //[Test]
        //public static void TestPeakSplittingLeft()
        //{
        //    int[] intensityMultipliers = { 1, 3, 1, 1, 3, 5, 10, 5, 3, 1 };

        //    TraceableTimsTofPeak testPeak = new TraceableTimsTofPeak(1, 1);
        //    for (int i = 0; i < 10; i++)
        //    {
        //        testPeak.IonMobilityPeaks.Add(new IonMobilityPeak(1.2, i + 1, intensityMultipliers[i]));
        //    }

        //    var indexedPeaks = testPeak.GetIndexedPeaks();
        //    Assert.That(indexedPeaks.Count(), Is.EqualTo(2));
        //    Assert.That(indexedPeaks.First().Mz, Is.EqualTo(1.2).Within(0.01));
        //}

        //[Test]
        //public static void TestPeakSplittingDoubleSided()
        //{
        //    int[] intensityMultipliers = { 1, 3, 1, 1, 3, 5, 10, 5, 3, 1, 1, 5, 1 };

        //    TraceableTimsTofPeak testPeak = new TraceableTimsTofPeak(1, 1);
        //    for (int i = 0; i < 13; i++)
        //    {
        //        testPeak.IonMobilityPeaks.Add(new IonMobilityPeak(1.2, i + 1, intensityMultipliers[i]));
        //    }

        //    var indexedPeaks = testPeak.GetIndexedPeaks();
        //    Assert.That(indexedPeaks.Count(), Is.EqualTo(3));
        //    Assert.That(indexedPeaks.First().Mz, Is.EqualTo(1.2).Within(0.01));
        //}
    }
}
