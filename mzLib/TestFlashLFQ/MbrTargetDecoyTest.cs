using Chemistry;
using FlashLFQ;
using MassSpectrometry;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using NUnit.Framework;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Easy.Common.Extensions;
using Test.FileReadingTests;
using UsefulProteomicsDatabases;
using ChromatographicPeak = FlashLFQ.ChromatographicPeak;
using Stopwatch = System.Diagnostics.Stopwatch;
using Peptide = Proteomics.AminoAcidPolymer.Peptide;
using System.Windows.Shapes;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    internal class MbrTargetDecoyTest
    {
        [Test]
        [TestCase(0, ExpectedResult = 0)]
        [TestCase(1, ExpectedResult = -1)]
        [TestCase(2, ExpectedResult = 1)]
        [TestCase(3, ExpectedResult = -2)]
        [TestCase(5, ExpectedResult = -3)]
        [TestCase(6, ExpectedResult = 3)]
        public static int TestDecoySearchFlipFlop(int searchCount)
        {
            // Integer division take ceiling: https://stackoverflow.com/questions/17944/how-to-round-up-the-result-of-integer-division
            int result = (searchCount + 2 - 1) / 2;
            result = searchCount % 2 == 0 ? result : -1 * result;

            return result;
        }

        [Test]
        // This is gonna have a bunch of local file references, just a heads up. Dont make github try and build this one
        public static void TwoFileMbrTest()
        {
            //string psmFile = @"D:\SingleCellDataSets\Organoid\TwoFileSearch\Task1-SearchTask\subset_psms.psmtsv";
            string psmFile = @"D:\SingleCellDataSets\Organoid\TwoFileSearch\Task1-SearchTask\AllPSMs_1PercentFdr.psmtsv";

            SpectraFileInfo j5 = new SpectraFileInfo(@"D:\SingleCellDataSets\Organoid\raw_files\HFL1SC_Unhealthy_CH2_J5.raw", "a", 0, 0, 0);
            SpectraFileInfo j6 = new SpectraFileInfo(@"D:\SingleCellDataSets\Organoid\raw_files\HFL1SC_Unhealthy_CH2_J6.raw", "a", 1, 0, 0);

            List<Identification> ids = new List<Identification>();
            Dictionary<string, ProteinGroup> allProteinGroups = new Dictionary<string, ProteinGroup>();
            foreach (string line in File.ReadAllLines(psmFile))
            {
                var split = line.Split(new char[] { '\t' });

                if (split.Contains("File Name") || string.IsNullOrWhiteSpace(line))
                {
                    continue;
                }

                SpectraFileInfo file = null;

                if (split[0].Contains("J5"))
                {
                    file = j5;
                }
                else if (split[0].Contains("J6"))
                {
                    file = j6;
                }

                string baseSequence = split[12];
                string fullSequence = split[13];
                double monoMass = double.Parse(split[22]);
                double rt = double.Parse(split[2]);
                int z = (int)double.Parse(split[6]);
                var proteins = split[25].Split(new char[] { '|' });
                List<ProteinGroup> proteinGroups = new List<ProteinGroup>();
                foreach (var protein in proteins)
                {
                    if (allProteinGroups.TryGetValue(protein, out var proteinGroup))
                    {
                        proteinGroups.Add(proteinGroup);
                    }
                    else
                    {
                        allProteinGroups.Add(protein, new ProteinGroup(protein, "", "Homo Sapiens"));
                        proteinGroups.Add(allProteinGroups[protein]);
                    }
                }

                bool isDecoy = split[32] == "Y";
                double score = double.TryParse(split[9], out var s) ? s : 0;

                Identification id = new Identification(file, baseSequence, fullSequence, monoMass, rt, z, proteinGroups, decoy: isDecoy, psmScore: score);
                ids.Add(id);
            }

            
            var engine = new FlashLfqEngine(ids, matchBetweenRuns: true, requireMsmsIdInCondition: false, maxThreads: 5, donorCriterion: 'S');
            var results = engine.Run();

            var test = results.Peaks.Values.SelectMany(peakList => peakList).ToList();
            int place = 0;

            List<ChromatographicPeak> mbrPeaks = new();

            mbrPeaks.AddRange(test.Where(peak => peak.IsMbrPeak && peak.RandomRt && !peak.DecoyPeptide).ToList());
            mbrPeaks.AddRange(test.Where(peak => peak.IsMbrPeak && peak.DecoyPeptide && !peak.RandomRt).ToList());
            mbrPeaks.AddRange(test.Where(peak => peak.IsMbrPeak && peak.DecoyPeptide && peak.RandomRt).ToList());
            mbrPeaks.AddRange(test.Where(peak => peak.IsMbrPeak && !peak.DecoyPeptide & !peak.RandomRt).ToList());


            using (StreamWriter writer = new StreamWriter(@"D:\SingleCellDataSets\Organoid\TwoFileSearch\Task1-SearchTask\RealMBR\MbrResults_1PepTest.tsv"))
            {
                writer.WriteLine(ChromatographicPeak.TabSeparatedHeader);
                foreach (var peak in mbrPeaks)
                {
                    writer.WriteLine(peak);
                }
            }

            //using (StreamWriter writer = new StreamWriter(@"D:\SingleCellDataSets\Organoid\TwoFileSearch\Task1-SearchTask\RealMBR\AllDecoys_minRtDiff.tsv"))
            //{
            //    writer.WriteLine(ChromatographicPeak.TabSeparatedHeader);
            //    foreach (var peak in engine.DecoyPeaks)
            //    {
            //        writer.WriteLine(peak);
            //    }
            //}

            var f1r1MbrResults = results
                .PeptideModifiedSequences
                .Where(p => p.Value.GetDetectionType(j5) == DetectionType.MBR && p.Value.GetDetectionType(j6) == DetectionType.MSMS).ToList();

            Assert.That(f1r1MbrResults.Count >= 132);

            results.WriteResults(peaksOutputPath: @"C:\Users\Alex\Desktop\FlashTest\AllPeaks.tsv", null, null, null, true);

            var f1r2MbrResults = results.PeptideModifiedSequences
                .Where(p => p.Value.GetDetectionType(j5) == DetectionType.MSMS && p.Value.GetDetectionType(j6) == DetectionType.MBR).ToList();

            Assert.That(f1r2MbrResults.Count >= 77);

            List<(double, double)> peptideIntensities = new List<(double, double)>();

            foreach (var peptide in f1r1MbrResults)
            {
                double mbrIntensity = Math.Log(peptide.Value.GetIntensity(j5));
                double msmsIntensity = Math.Log(peptide.Value.GetIntensity(j6));
                peptideIntensities.Add((mbrIntensity, msmsIntensity));
            }

            double corr = Correlation.Pearson(peptideIntensities.Select(p => p.Item1), peptideIntensities.Select(p => p.Item2));
            Assert.Greater(corr, 0.8);

            peptideIntensities.Clear();
            foreach (var peptide in f1r2MbrResults)
            {
                double mbrIntensity = Math.Log(peptide.Value.GetIntensity(j6));
                double msmsIntensity = Math.Log(peptide.Value.GetIntensity(j5));
                peptideIntensities.Add((mbrIntensity, msmsIntensity));
            }

            corr = Correlation.Pearson(peptideIntensities.Select(p => p.Item1), peptideIntensities.Select(p => p.Item2));

            Assert.That(corr > 0.7);

            // the "requireMsmsIdInCondition" field requires that at least one MS/MS identification from a protein
            // has to be observed in a condition for match-between-runs
            j5.Condition = "b";
            engine = new FlashLfqEngine(ids, matchBetweenRuns: true, requireMsmsIdInCondition: true, maxThreads: 5);
            results = engine.Run();
            var proteinsObservedInF1 = ids.Where(p => p.FileInfo == j5).SelectMany(p => p.ProteinGroups).Distinct().ToList();
            var proteinsObservedInF2 = ids.Where(p => p.FileInfo == j6).SelectMany(p => p.ProteinGroups).Distinct().ToList();
            var proteinsObservedInF1ButNotF2 = proteinsObservedInF1.Except(proteinsObservedInF2).ToList();
            foreach (ProteinGroup protein in proteinsObservedInF1ButNotF2)
            {
                Assert.That(results.ProteinGroups[protein.ProteinGroupName].GetIntensity(j6) == 0);
            }
        }


        [Test]
        // This is gonna have a bunch of local file references, just a heads up. Dont make github try and build this one
        public static void ThreeFileMbrTest()
        {
            //string psmFile = @"D:\SingleCellDataSets\Organoid\TwoFileSearch\Task1-SearchTask\subset_psms.psmtsv";
            string psmFile = @"D:\SingleCellDataSets\Organoid\TwoFileSearch\Task1-SearchTask\AllPSMs_1PercentFdr.psmtsv";
            string psmFile2 = @"D:\SingleCellDataSets\Organoid\Search_MM_320\Task1-SearchTask\Individual File Results\HFL1SC_Unhealthy_CH2_J7-calib_PSMs.psmtsv";

            SpectraFileInfo j5 = new SpectraFileInfo(@"D:\SingleCellDataSets\Organoid\raw_files\HFL1SC_Unhealthy_CH2_J5.raw", "a", 0, 0, 0);
            SpectraFileInfo j6 = new SpectraFileInfo(@"D:\SingleCellDataSets\Organoid\raw_files\HFL1SC_Unhealthy_CH2_J6.raw", "a", 1, 0, 0);
            SpectraFileInfo j7 = new SpectraFileInfo(@"D:\SingleCellDataSets\Organoid\Calibration_MM_320\Task1-CalibrateTask\HFL1SC_Unhealthy_CH2_J7-calib.mzML", "a", 2, 0, 0);


            List<Identification> ids = new List<Identification>();
            Dictionary<string, ProteinGroup> allProteinGroups = new Dictionary<string, ProteinGroup>();
            foreach (string line in File.ReadAllLines(psmFile))
            {
                var split = line.Split(new char[] { '\t' });

                if (split.Contains("File Name") || string.IsNullOrWhiteSpace(line))
                {
                    continue;
                }

                SpectraFileInfo file = null;

                if (split[0].Contains("J5"))
                {
                    file = j5;
                }
                else if (split[0].Contains("J6"))
                {
                    file = j6;
                }

                string baseSequence = split[12];
                string fullSequence = split[13];
                double monoMass = double.Parse(split[22]);
                double rt = double.Parse(split[2]);
                int z = (int)double.Parse(split[6]);
                var proteins = split[25].Split(new char[] { '|' });
                List<ProteinGroup> proteinGroups = new List<ProteinGroup>();
                foreach (var protein in proteins)
                {
                    if (allProteinGroups.TryGetValue(protein, out var proteinGroup))
                    {
                        proteinGroups.Add(proteinGroup);
                    }
                    else
                    {
                        allProteinGroups.Add(protein, new ProteinGroup(protein, "", ""));
                        proteinGroups.Add(allProteinGroups[protein]);
                    }
                }

                bool isDecoy = split[32] == "Y";

                Identification id = new Identification(file, baseSequence, fullSequence, monoMass, rt, z, proteinGroups, decoy: isDecoy);
                ids.Add(id);
            }
            foreach(string line in File.ReadAllLines(psmFile2))
            {
                var split = line.Split(new char[] { '\t' });

                if (split.Contains("File Name") || string.IsNullOrWhiteSpace(line))
                {
                    continue;
                }

                SpectraFileInfo file = j7;

                double qval = Double.Parse(split[50]);
                if (qval > 0.01) continue;

                string baseSequence = split[12];
                string fullSequence = split[13];
                if(!double.TryParse(split[22], out var x))
                {
                    continue; // Occurs for ambiguous peptides
                }
                double monoMass = double.Parse(split[22]);
                double rt = double.Parse(split[2]);
                int z = (int)double.Parse(split[6]);
                var proteins = split[25].Split(new char[] { '|' });
                List<ProteinGroup> proteinGroups = new List<ProteinGroup>();
                foreach (var protein in proteins)
                {
                    if (allProteinGroups.TryGetValue(protein, out var proteinGroup))
                    {
                        proteinGroups.Add(proteinGroup);
                    }
                    else
                    {
                        allProteinGroups.Add(protein, new ProteinGroup(protein, "", ""));
                        proteinGroups.Add(allProteinGroups[protein]);
                    }
                }

                bool isDecoy = split[32] == "Y";
                double score = double.TryParse(split[9], out var s) ? s : 0;

                Identification id = new Identification(file, baseSequence, fullSequence, monoMass, rt, z, proteinGroups, decoy: isDecoy, psmScore: score);
                ids.Add(id);
            }


            var engine = new FlashLfqEngine(ids, matchBetweenRuns: true, requireMsmsIdInCondition: false, maxThreads: 5, donorCriterion: 'S');
            var results = engine.Run();

            var test = results.Peaks.Values.SelectMany(peakList => peakList).ToList();
            int place = 0;

            List<ChromatographicPeak> mbrPeaks = new();

            mbrPeaks.AddRange(test.Where(peak => peak.IsMbrPeak && peak.RandomRt && !peak.DecoyPeptide).ToList());
            mbrPeaks.AddRange(test.Where(peak => peak.IsMbrPeak && peak.DecoyPeptide && !peak.RandomRt).ToList());
            mbrPeaks.AddRange(test.Where(peak => peak.IsMbrPeak && peak.DecoyPeptide && peak.RandomRt).ToList());
            mbrPeaks.AddRange(test.Where(peak => peak.IsMbrPeak && !peak.DecoyPeptide & !peak.RandomRt).ToList());


            using (StreamWriter writer = new StreamWriter(@"D:\SingleCellDataSets\Organoid\TwoFileSearch\Task1-SearchTask\RealMBR\MbrResults_1PeakPerPepScore.tsv"))
            {
                writer.WriteLine(ChromatographicPeak.TabSeparatedHeader);
                foreach (var peak in mbrPeaks)
                {
                    writer.WriteLine(peak);
                }
            }

            //using (StreamWriter writer = new StreamWriter(@"D:\SingleCellDataSets\Organoid\TwoFileSearch\Task1-SearchTask\RealMBR\AllDecoys_minRtDiff.tsv"))
            //{
            //    writer.WriteLine(ChromatographicPeak.TabSeparatedHeader);
            //    foreach (var peak in engine.DecoyPeaks)
            //    {
            //        writer.WriteLine(peak);
            //    }
            //}

            var f1r1MbrResults = results
                .PeptideModifiedSequences
                .Where(p => p.Value.GetDetectionType(j5) == DetectionType.MBR && p.Value.GetDetectionType(j6) == DetectionType.MSMS).ToList();

            Assert.That(f1r1MbrResults.Count >= 132);

            var f1r2MbrResults = results.PeptideModifiedSequences
                .Where(p => p.Value.GetDetectionType(j5) == DetectionType.MSMS && p.Value.GetDetectionType(j6) == DetectionType.MBR).ToList();

            Assert.That(f1r2MbrResults.Count >= 77);

            List<(double, double)> peptideIntensities = new List<(double, double)>();

            foreach (var peptide in f1r1MbrResults)
            {
                double mbrIntensity = Math.Log(peptide.Value.GetIntensity(j5));
                double msmsIntensity = Math.Log(peptide.Value.GetIntensity(j6));
                peptideIntensities.Add((mbrIntensity, msmsIntensity));
            }

            double corr = Correlation.Pearson(peptideIntensities.Select(p => p.Item1), peptideIntensities.Select(p => p.Item2));
            Assert.Greater(corr, 0.8);

            peptideIntensities.Clear();
            foreach (var peptide in f1r2MbrResults)
            {
                double mbrIntensity = Math.Log(peptide.Value.GetIntensity(j6));
                double msmsIntensity = Math.Log(peptide.Value.GetIntensity(j5));
                peptideIntensities.Add((mbrIntensity, msmsIntensity));
            }

            corr = Correlation.Pearson(peptideIntensities.Select(p => p.Item1), peptideIntensities.Select(p => p.Item2));

            Assert.That(corr > 0.7);

            // the "requireMsmsIdInCondition" field requires that at least one MS/MS identification from a protein
            // has to be observed in a condition for match-between-runs
            j5.Condition = "b";
            engine = new FlashLfqEngine(ids, matchBetweenRuns: true, requireMsmsIdInCondition: true, maxThreads: 5);
            results = engine.Run();
            var proteinsObservedInF1 = ids.Where(p => p.FileInfo == j5).SelectMany(p => p.ProteinGroups).Distinct().ToList();
            var proteinsObservedInF2 = ids.Where(p => p.FileInfo == j6).SelectMany(p => p.ProteinGroups).Distinct().ToList();
            var proteinsObservedInF1ButNotF2 = proteinsObservedInF1.Except(proteinsObservedInF2).ToList();
            foreach (ProteinGroup protein in proteinsObservedInF1ButNotF2)
            {
                Assert.That(results.ProteinGroups[protein.ProteinGroupName].GetIntensity(j6) == 0);
            }
        }
    }
}
