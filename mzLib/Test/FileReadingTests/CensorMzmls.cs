using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using Proteomics.PSM;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using Proteomics;
using Readers;
using MassSpectrometry;

namespace Test.FileReadingTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class CensorMzmls
    {

        [Test]
        public static void CensorKellyTwoProteomeData()
        {
            string peptideFilePath = @"D:\Kelly_TwoProteomeData\NineFileSearch_MM105\Task1-SearchTask\AllPeptides.psmtsv";
            List<PsmFromTsv> parsedPeptides = SpectrumMatchTsvReader.ReadPsmTsv(peptideFilePath, out var warnings).Where(psm => psm.PEP_QValue < 0.002).ToList();
            HashSet<string> allPeptideSeqs = parsedPeptides.Select(pep => pep.FullSequence).ToHashSet();
            Assert.That(warnings.Count, Is.EqualTo(0));
            Assert.That(parsedPeptides.Count, Is.GreaterThan(100));

            string psmFilePath = @"D:\Kelly_TwoProteomeData\NineFileSearch_MM105\Task1-SearchTask\AllPSMs.psmtsv";
            List<PsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFilePath, out warnings).Where(psm => psm.PEP_QValue < 0.01).ToList();

            var psmHeaderDict = SpectrumMatchTsvReader.ParseHeader(File.ReadLines(psmFilePath).First());
            var psmHeaderString = string.Join('\t', psmHeaderDict.Select(kvp => kvp).OrderBy(kvp => kvp.Value).Select(kvp => kvp.Value));

            var psmsByFile = parsedPsms.GroupBy(psm => psm.FileNameWithoutExtension).ToDictionary(group => group.Key, group => group.ToList());
            var peptideSequencesByFile = parsedPeptides
                .GroupBy(psm => psm.FileNameWithoutExtension)
                .ToDictionary(
                    group => group.Key, 
                    group => group.Select(psm => psm.FullSequence).ToHashSet());

            string mzmlFolder = @"D:\Kelly_TwoProteomeData\TenFile_NewDataMM105\Task1-CalibrateTask";
            MzSpectrum blankSpectrum = new MzSpectrum(mz: new double[] { 150 }, intensities: new double[] { 10000 }, false);
            string outputFolder = @"D:\Kelly_TwoProteomeData\CensoredDataFiles";

            List<PsmFromTsv> censoredPsms = new();

            foreach (var kvp in psmsByFile)
            {
                if (kvp.Key.Contains("HeYe")) continue;

                var psmsToRemove = kvp.Value
                   // If this file contains the top peptide instance, we can't remove it
                   .Where(psm => !peptideSequencesByFile[kvp.Key].Contains(psm.FullSequence) && psm.OrganismName.Equals("Homo sapiens") && psm.DecoyContamTarget == "T")
                   .GroupBy(psm => psm.FullSequence) 
                   .Where(group => allPeptideSeqs.Contains(group.Key) && group.Count() == 1)
                   .OrderBy(group => Guid.NewGuid()) // Order randomly, that way we don't only select the top scoring psms
                   .SelectMany(group => group.Take(1))
                   .Take(500)
                   .ToList();
                censoredPsms.AddRange(psmsToRemove);

                // Open up the mzMl file specified by the psm FileName
                // Select all scans referenced by the psms
                // set the intensity of all peaks but 1 to zero
                // Save the new mzMl file
                string mzmlFullPath = Path.Combine(mzmlFolder, psmsToRemove.First().FileNameWithoutExtension + ".mzML");
                var file = MsDataFileReader.GetDataFile(mzmlFullPath);
                file.LoadAllStaticData();

                foreach(var psm in psmsToRemove)
                {
                    MsDataScan scan = file.GetOneBasedScan(psm.Ms2ScanNumber);
                    scan.MassSpectrum = blankSpectrum;
                }

                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(file,
                    Path.Combine(outputFolder, psmsToRemove.First().FileNameWithoutExtension + "-censored.mzML"),
                    writeIndexed: true);

            }

            using (StreamWriter writer = new StreamWriter(Path.Combine(outputFolder, "CensoredPsms.psmtsv")))
            {
                writer.WriteLine(PsmFromTsv.GetHeader());
                foreach(var psm in censoredPsms)
                {
                    writer.WriteLine(psm);
                }
            }

        }

        [Test]
        public static void CensorFraggerKelly()
        {
            string peptideFilePath = @"D:\Kelly_TwoProteomeData\IonQuant1Percent\combined_peptide.tsv";

            var peptideFile = new MsFraggerPeptideFile(peptideFilePath);
            peptideFile.LoadResults();
            List<MsFraggerPeptide> parsedPeptides = peptideFile.Results;

            int placeholder = 0;
            // Filter on peptideProphet probability, 0.99 or higher
        }

        [Test]
        public static void CensorGygiTwoProteomeData()
        {
            string peptideFilePath = @"D:\GygiTwoProteome_PXD014415\Yeast_Human_Arabad_Contam_search\AllPeptides.psmtsv";
            List<PsmFromTsv> parsedPeptides = SpectrumMatchTsvReader.ReadPsmTsv(peptideFilePath, out var warnings).Where(psm => psm.PEP_QValue < 0.002).ToList();
            HashSet<string> allPeptideSeqs = parsedPeptides.Select(pep => pep.FullSequence).ToHashSet();
            Assert.That(warnings.Count, Is.EqualTo(0));
            Assert.That(parsedPeptides.Count, Is.GreaterThan(100));

            string psmFilePath = @"D:\GygiTwoProteome_PXD014415\Yeast_Human_Arabad_Contam_search\AllPSMs.psmtsv";
            List<PsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFilePath, out warnings).Where(psm => psm.PEP_QValue < 0.01).ToList();

            var psmsByFile = parsedPsms.GroupBy(psm => psm.FileNameWithoutExtension).ToDictionary(group => group.Key, group => group.ToList());
            var peptideSequencesByFile = parsedPeptides
                .GroupBy(psm => psm.FileNameWithoutExtension)
                .ToDictionary(
                    group => group.Key,
                    group => group.Select(psm => psm.FullSequence).ToHashSet());

            var psmHeaderDict = SpectrumMatchTsvReader.ParseHeader(File.ReadLines(psmFilePath).First());
            var psmHeaderString = string.Join('\t', psmHeaderDict.Select(kvp => kvp).OrderBy(kvp => kvp.Value).Select(kvp => kvp.Value));

            string mzmlFolder = @"D:\GygiTwoProteome_PXD014415\MM105_Calibration_Search_MBR\Task1-CalibrateTask";
            MzSpectrum blankSpectrum = new MzSpectrum(mz: new double[] { 150 }, intensities: new double[] { 10000 }, false);
            string outputFolder = @"D:\GygiTwoProteome_PXD014415\CensoredDataFiles";
            List<PsmFromTsv> censoredPsms = new();


            foreach (var kvp in psmsByFile)
            {
                if (kvp.Key.Contains("yaeast")) continue;

                var psmsToRemove = kvp.Value
                   // If this file contains the top peptide instance, we can't remove it
                   .Where(psm => !peptideSequencesByFile[kvp.Key].Contains(psm.FullSequence) && psm.OrganismName.Equals("Homo sapiens") && psm.DecoyContamTarget == "T")
                   .GroupBy(psm => psm.FullSequence)
                   .Where(group => allPeptideSeqs.Contains(group.Key) && group.Count() == 1)
                   .OrderBy(group => Guid.NewGuid()) // Order randomly, ensure different things are selected for each file
                   .SelectMany(group => group.Take(1))
                   .Take(500)
                   .ToList();
                censoredPsms.AddRange(psmsToRemove);

                // Open up the mzMl file specified by the psm FileName
                // Select all scans referenced by the psms
                // set the intensity of all peaks but 1 to zero
                // Save the new mzMl file
                string mzmlFullPath = Path.Combine(mzmlFolder, psmsToRemove.First().FileNameWithoutExtension + ".mzML");
                var file = MsDataFileReader.GetDataFile(mzmlFullPath);
                file.LoadAllStaticData();

                foreach (var psm in psmsToRemove)
                {
                    MsDataScan scan = file.GetOneBasedScan(psm.Ms2ScanNumber);
                    scan.MassSpectrum = blankSpectrum;
                }

                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(file,
                    Path.Combine(outputFolder, psmsToRemove.First().FileNameWithoutExtension + "-censored.mzML"),
                    writeIndexed: true);

            }

            using (StreamWriter writer = new StreamWriter(Path.Combine(outputFolder, "CensoredPsms.psmtsv")))
            {
                writer.WriteLine(PsmFromTsv.GetHeader());
                foreach (var psm in censoredPsms)
                {
                    writer.WriteLine(psm);
                }
            }

        }

        [Test]
        public static void CensorInHouseHumanTwoProteomeData()
        {
            string peptideFilePath = @"D:\Human_Ecoli_TwoProteome_60minGradient\MM105_Human_Mixed\Task1-SearchTask\AllPeptides.psmtsv";
            List<PsmFromTsv> parsedPeptides = SpectrumMatchTsvReader.ReadPsmTsv(peptideFilePath, out var warnings).Where(psm => psm.PEP_QValue < 0.002).ToList();
            HashSet<string> allPeptideSeqs = parsedPeptides.Select(pep => pep.FullSequence).ToHashSet();
            Assert.That(warnings.Count, Is.EqualTo(0));
            Assert.That(parsedPeptides.Count, Is.GreaterThan(100));

            string psmFilePath = @"D:\Human_Ecoli_TwoProteome_60minGradient\MM105_Human_Mixed\Task1-SearchTask\AllPSMs.psmtsv";
            List<PsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFilePath, out warnings).Where(psm => psm.PEP_QValue < 0.01).ToList();

            var psmsByFile = parsedPsms.GroupBy(psm => psm.FileNameWithoutExtension).ToDictionary(group => group.Key, group => group.ToList());
            var peptideSequencesByFile = parsedPeptides
                .GroupBy(psm => psm.FileNameWithoutExtension)
                .ToDictionary(
                    group => group.Key,
                    group => group.Select(psm => psm.FullSequence).ToHashSet());

            var psmHeaderDict = SpectrumMatchTsvReader.ParseHeader(File.ReadLines(psmFilePath).First());
            var psmHeaderString = string.Join('\t', psmHeaderDict.Select(kvp => kvp).OrderBy(kvp => kvp.Value).Select(kvp => kvp.Value));

            string mzmlFolder = @"D:\Human_Ecoli_TwoProteome_60minGradient\CalibrateSearch_4_19_24\Human_Calibrated_Files";
            MzSpectrum blankSpectrum = new MzSpectrum(mz: new double[] { 150 }, intensities: new double[] { 10000 }, false);
            string outputFolder = @"D:\Human_Ecoli_TwoProteome_60minGradient\CensoredHumanData";
            List<PsmFromTsv> censoredPsms = new();


            foreach (var kvp in psmsByFile)
            {
                if (kvp.Key.Contains("yaeast")) continue;

                var psmsToRemove = kvp.Value
                   // If this file contains the top peptide instance, we can't remove it
                   .Where(psm => !peptideSequencesByFile[kvp.Key].Contains(psm.FullSequence) && psm.OrganismName.Equals("Homo sapiens") && psm.DecoyContamTarget == "T")
                   .GroupBy(psm => psm.FullSequence)
                   .Where(group => allPeptideSeqs.Contains(group.Key) && group.Count() == 1)
                   .OrderBy(group => Guid.NewGuid()) // Order randomly, ensure different things are selected for each file
                   .SelectMany(group => group.Take(1))
                   .Take(1000)
                   .ToList();
                censoredPsms.AddRange(psmsToRemove);

                // Open up the mzMl file specified by the psm FileName
                // Select all scans referenced by the psms
                // set the intensity of all peaks but 1 to zero
                // Save the new mzMl file
                string mzmlFullPath = Path.Combine(mzmlFolder, psmsToRemove.First().FileNameWithoutExtension + ".mzML");
                var file = MsDataFileReader.GetDataFile(mzmlFullPath);
                file.LoadAllStaticData();

                foreach (var psm in psmsToRemove)
                {
                    MsDataScan scan = file.GetOneBasedScan(psm.Ms2ScanNumber);
                    scan.MassSpectrum = blankSpectrum;
                }

                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(file,
                    Path.Combine(outputFolder, psmsToRemove.First().FileNameWithoutExtension + "-censored.mzML"),
                    writeIndexed: true);

            }

            using (StreamWriter writer = new StreamWriter(Path.Combine(outputFolder, "CensoredPsms.psmtsv")))
            {
                writer.WriteLine(PsmFromTsv.GetHeader());
                foreach (var psm in censoredPsms)
                {
                    writer.WriteLine(psm);
                }
            }

        }

    }
}
