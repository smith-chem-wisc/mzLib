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
        public static void DecoyPeakFindTrial()
        {
            string decoyPeptidePath = @"C:\Users\Alex\Source\Repos\chronologer\chronologer-main(noChange)\Arabidopsis_Tryptic_Peptides_Slice.tsv";
            string spectraFilePath = @"D:\SingleCellDataSets\Organoid\Calibration_MM_320\Task1-CalibrateTask\HFL1SC_Unhealthy_CH2_J5-calib.mzML";
            ProteinGroup pg = new ProteinGroup("xyz", "x", "z");
            List<ProteinGroup> pgs = new List<ProteinGroup> { pg };
            SpectraFileInfo j5 = new SpectraFileInfo(spectraFilePath, "A", 1, 1, 1);
            double rtRange = 4.0;

            List<Identification> decoys = new();
            Loaders.LoadElements();

            using (StreamReader reader = new StreamReader(decoyPeptidePath))
            {
                reader.ReadLine();
                while(!reader.EndOfStream)
                {
                    string[] lineSplit = reader.ReadLine().Split('\t');
                    if (double.Parse(lineSplit[4]) > 57.5) continue;

                    Peptide peptide = new Peptide(sequence: lineSplit[0]);
                    Identification decoyId = new Identification(j5, peptide.BaseSequence, peptide.BaseSequence,
                        peptide.MonoisotopicMass, ms2RetentionTimeInMinutes: double.Parse(lineSplit[4]),
                        chargeState: peptide.MonoisotopicMass > 1600 ? 3 : 2, pgs);
                    decoys.Add(decoyId);
                }
            }


            FlashLfqEngine engine = new FlashLfqEngine(
                decoys
                );
            
            engine.CalculateTheoreticalIsotopeDistributions();
            engine._ms1Scans = new Dictionary<SpectraFileInfo, Ms1ScanInfo[]>();
            engine._peakIndexingEngine.IndexMassSpectralPeaks(j5, true, engine._ms1Scans);


            List<ChromatographicPeak> decoyPeaks = new();
            List<double> massDifferences = new();
            Dictionary<IndexedMassSpectralPeak, ChromatographicPeak> apexPeakDict = new();
            int decoysConsidered = 0;
            Random rnd = new Random();
            foreach (Identification decoy in decoys)
            {

                int rndInt = rnd.Next(1, 13);
                // Eliminate ~ half of the decoys with mass greater than 2000 daltons
                // This is an ad-hoc way of matching target and decoy mass distribution
                if (decoy.PeakfindingMass > 1800 && rndInt % 2 == 0) continue;
                else if (decoy.PeakfindingMass > 1400 && rndInt % 3 == 0) continue;

                if (decoy.Ms2RetentionTimeInMinutes < 8 && rndInt < 11) continue;

                PpmTolerance tolerance = new PpmTolerance(10);
                var foundPeak = engine.FindDecoyPeak(
                    j5,
                    apexPeakDict,
                    tolerance,
                    (decoy.Ms2RetentionTimeInMinutes, rtRange, null, null),
                    decoy);

                if (foundPeak != null)
                {
                    decoyPeaks.Add(foundPeak);
                    massDifferences.Add(
                        Math.Abs(
                            decoy.PeakfindingMass - foundPeak.Apex.IndexedPeak.Mz.ToMass(foundPeak.Apex.ChargeState)
                            ));
                }

                decoysConsidered++;
                if (decoyPeaks.Count >= 750) break;
            }

            int placeholder = 0;

            double massDiffMean = massDifferences.Select(m => Math.Abs(m)).Average();
            double envelopeCountMean = decoyPeaks.Select(peak => peak.IsotopicEnvelopes.Count).Average();
            double intensityMean = decoyPeaks.Select(peak => peak.IsotopicEnvelopes.Select(e => e.Intensity).Average()).Average();

            placeholder = 1;

            // Repeat, but for targets
            string targetPeptidePath = @"C:\Users\Alex\Source\Repos\chronologer\chronologer-main(noChange)\Unhealthy_CH2_J5_MBR_Predicted.tsv";
            //For MBR Predicted file
            int fullSeqCol = 3;
            int massCol = 5;
            int rtColumn = 24;
            
            List<Identification> targetIDs = new();
            using (StreamReader reader = new StreamReader(targetPeptidePath))
            {
                reader.ReadLine();
                while (!reader.EndOfStream)
                {
                    string[] lineSplit = reader.ReadLine().Split('\t');
                    if (lineSplit[fullSeqCol].Contains('[')) continue;
                    if (double.Parse(lineSplit[rtColumn]) > 60) continue;
                    

                    Peptide peptide = new Peptide(sequence: lineSplit[fullSeqCol]);
                    Identification targetId = new Identification(j5, peptide.BaseSequence, peptide.BaseSequence,
                        peptide.MonoisotopicMass, ms2RetentionTimeInMinutes: double.Parse(lineSplit[rtColumn]),
                        chargeState: peptide.MonoisotopicMass > 1600 ? 3 : 2, pgs);
                    targetIDs.Add(targetId);
                }
            }


            engine = new FlashLfqEngine(
                targetIDs
                );

            engine.CalculateTheoreticalIsotopeDistributions();
            engine._ms1Scans = new Dictionary<SpectraFileInfo, Ms1ScanInfo[]>();
            engine._peakIndexingEngine.IndexMassSpectralPeaks(j5, true, engine._ms1Scans);


            List<ChromatographicPeak> targetPeaks = new();
            List<double> massDifferencesTarget = new();
            Dictionary<IndexedMassSpectralPeak, ChromatographicPeak> apexPeakDictTarget = new();
            int targetsConsidered = 0;
            foreach (Identification target in targetIDs)
            {
                PpmTolerance tolerance = new PpmTolerance(10);
                var foundPeak = engine.FindDecoyPeak(
                    j5,
                    apexPeakDictTarget,
                    tolerance,
                    (target.Ms2RetentionTimeInMinutes, rtRange, null, null),
                    target);

                if (foundPeak != null)
                {
                    targetPeaks.Add(foundPeak);
                    massDifferencesTarget.Add(
                        Math.Abs(
                            target.PeakfindingMass - foundPeak.Apex.IndexedPeak.Mz.ToMass(foundPeak.Apex.ChargeState)
                            ));
                }

                targetsConsidered++;
                if (targetPeaks.Count >= 750) break;
            }

            double massDiffMeanT = massDifferencesTarget.Select(m => Math.Abs(m)).Average();
            double envelopeCountMeanT = targetPeaks.Select(peak => peak.IsotopicEnvelopes.Count).Average();
            double intensityMeanT = targetPeaks.Select(peak => peak.IsotopicEnvelopes.Select(e => e.Intensity).Average()).Average();

            placeholder = 2;

            using(StreamWriter writer = new StreamWriter(@"C:\Users\Alex\Desktop\MBR_10_30\Take8_MbrTargetRT_4MinWindow.tsv"))
            {
                string[] header = new string[]
                    {
                        "Sequence",
                        "Target/Decoy",
                        "Theoretical Peak Finding Mass",
                        "Found Mass",
                        "Number of Scans",
                        "Apex Intensity",
                        "Retention Time",
                        "Predicted Retention Time"
                    };
                writer.WriteLine(string.Join('\t', header));

                foreach (var decoy in decoyPeaks)
                {
                    double peakFindingMass = decoy.Identifications.First().PeakfindingMass;
                    header = new string[]
                    {
                        decoy.Identifications.First().BaseSequence,
                        "D",
                        decoy.Identifications.First().PeakfindingMass.ToString(),
                        decoy.Apex.IndexedPeak.Mz.ToMass(decoy.Apex.ChargeState).ToString(),
                        decoy.IsotopicEnvelopes.Count.ToString(),
                        decoy.Apex.Intensity.ToString(),
                        decoy.Apex.IndexedPeak.RetentionTime.ToString(),
                        decoy.Identifications.First().Ms2RetentionTimeInMinutes.ToString()
                    };
                    writer.WriteLine(string.Join('\t', header));
                }

                foreach (var target in targetPeaks)
                {
                    double peakFindingMass = target.Identifications.First().PeakfindingMass;
                    header = new string[]
                    {
                        target.Identifications.First().BaseSequence,
                        "T",
                        target.Identifications.First().PeakfindingMass.ToString(),
                        target.Apex.IndexedPeak.Mz.ToMass(target.Apex.ChargeState).ToString(),
                        target.IsotopicEnvelopes.Count.ToString(),
                        target.Apex.Intensity.ToString(),
                        target.Apex.IndexedPeak.RetentionTime.ToString(),
                        target.Identifications.First().Ms2RetentionTimeInMinutes.ToString()
                    };
                    writer.WriteLine(string.Join('\t', header));
                }
            }

        }

    }
}
