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
using CsvHelper;
using static Nett.TomlObjectFactory;

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
        public static void TinyTest()
        {
            string filePath = @"D:\Human_Ecoli_TwoProteome_60minGradient\RawData\04-12-24_Human_C18_3mm_50msec_stnd-60min_1_centroid.mzML";

            var file = MsDataFileReader.GetDataFile(filePath);
            file.LoadAllStaticData();

            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(file, @"D:\Human_Ecoli_TwoProteome_60minGradient\RawData\04-12-24_Human_C18_3mm_50msec_stnd-60min_1_mzLib_convert.mzML", writeIndexed: true);
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


        [Test]
        public static void CensorFraggerKelly()
        {
            string peptideFilePath = @"D:\Kelly_TwoProteomeData\MsConvertMzMls\IonQuant_1Percent\combined_modified_peptide.tsv";

            var peptideFile = new MsFraggerPeptideFile(peptideFilePath);
            peptideFile.LoadResults();
            List<MsFraggerPeptide> parsedPeptides = peptideFile.Results;

            var testPept = parsedPeptides.First();

            MzSpectrum blankSpectrum = new MzSpectrum(mz: new double[] { 150 }, intensities: new double[] { 10000 }, false);
            int placeholder = 0;
            // Filter on peptideProphet probability, 0.99 or higher

            // For each peptide, select the run with the highest intensity. That run is the master donor and will not be censored
            // Construct a dictionary linking file/run to the list of master donor peptides
            Dictionary<string, List<string>> masterDonorPeptidesByFile = new();
            foreach (var peptide in parsedPeptides)
            {
                string masterDonor = peptide.IntensityByFile.MaxBy(kvp => kvp.Value).Key;
                if (!masterDonorPeptidesByFile.ContainsKey(masterDonor))
                {
                    masterDonorPeptidesByFile.Add(masterDonor, new List<string> { peptide.FullSequence });
                }
                else
                {
                    masterDonorPeptidesByFile[masterDonor].Add(peptide.FullSequence);
                }
            }

            placeholder += 1;
            List<string> acceptorSamples = new List<string> { "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9" };
            var expFile = new MsFraggerExperimentFile(@"D:\Kelly_TwoProteomeData\MsConvertMzMls\IonQuant_1Percent\experiment_annotation.tsv");
            List<MsFraggerExperiment> experiments = expFile.Results;

            List<MsFraggerPsm> censoredPsms = new();


            // Iterate through the psm.tsv file for each run
            // Select only human peptides with peptide prophet probability > 0.99
            // and group by modified peptide (create some field that defaults to base sequence if mod seq is NA)
            // Fragger only gives scan retention times (in seconds), so we'll need to map the RT to the scan number, which will suck
            // Randomly select 500 peptides to censor
            // Load in the raw file, iterate through each peptide group
            // for every psm in the group, censor the corresponding spectrum
            foreach (var exp in experiments.Where(exp => exp.FullFilePathWithExtension.Contains("02nguL")))
            {
                var psms = new MsFraggerPsmFile(Path.Combine(@"D:\Kelly_TwoProteomeData\MsConvertMzMls\IonQuant_1Percent\", exp.Sample, "psm.tsv"));

                var psmsToRemove = psms
                    .Where(pep => pep.Protein.Contains("HUMAN") && !pep.Protein.Contains("REV") && !pep.Protein.Contains("CON")
                        && pep.PeptideProphetProbability > 0.99
                        && !masterDonorPeptidesByFile[exp.Sample].Contains(pep.ModifiedSequence))
                    .GroupBy(pep => pep.ModifiedSequence)
                    .Where(group => group.Count() == 1)
                    .OrderBy(group => Guid.NewGuid())
                    .SelectMany(group => group.Take(1))
                    .Take(500)
                    .ToList();

                censoredPsms.AddRange(psmsToRemove);

                // Open up the mzMl file specified by the psm FileName
                // Select all scans referenced by the psms
                // set the intensity of all peaks but 1 to zero
                // Save the new mzMl file
                string rawFilePath = Path.Combine(exp.FullFilePathWithExtension);
                var file = MsDataFileReader.GetDataFile(rawFilePath);
                file.LoadAllStaticData();

                foreach (var psm in psmsToRemove)
                {
                    MsDataScan scan = file.GetOneBasedScan(psm.OneBasedScanNumber);
                    scan.MassSpectrum = blankSpectrum;
                }

                string fileName = psmsToRemove.First().FileNameWithoutExtension;

                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(file,
                                       Path.Combine(@"D:\Kelly_TwoProteomeData\MsConvertMzMls\CensoredFiles", fileName + "-censored.mzML"),
                                                          writeIndexed: true);
            }

            // Write out the list of censored PSMs 
            using (var csv =
                new CsvWriter(new StreamWriter(Path.Combine(@"D:\Kelly_TwoProteomeData\MsConvertMzMls\CensoredFiles", "CensoredPsms.tsv")), MsFraggerPsm.CsvConfiguration))
            {
                csv.WriteHeader<MsFraggerPsm>();
                foreach (var result in censoredPsms)
                {
                    csv.NextRecord();
                    csv.WriteRecord(result);
                }
            }

        }

        [Test]
        public static void CensorFraggerGygi()
        {
            string peptideFilePath = @"D:\GygiTwoProteome_PXD014415\MsConvertmzMLs\IonQuant_1Percent\combined_modified_peptide.tsv";

            var peptideFile = new MsFraggerPeptideFile(peptideFilePath);
            peptideFile.LoadResults();
            List<MsFraggerPeptide> parsedPeptides = peptideFile.Results;

            var testPept = parsedPeptides.First();

            MzSpectrum blankSpectrum = new MzSpectrum(mz: new double[] { 150 }, intensities: new double[] { 10000 }, false);
            int placeholder = 0;

            // For each peptide, select the run with the highest intensity. That run is the master donor and will not be censored
            // Construct a dictionary linking file/run to the list of master donor peptides
            Dictionary<string, List<string>> masterDonorPeptidesByFile = new();
            foreach (var peptide in parsedPeptides)
            {
                string masterDonor = peptide.IntensityByFile.MaxBy(kvp => kvp.Value).Key;
                if (!masterDonorPeptidesByFile.ContainsKey(masterDonor))
                {
                    masterDonorPeptidesByFile.Add(masterDonor, new List<string> { peptide.FullSequence });
                }
                else
                {
                    masterDonorPeptidesByFile[masterDonor].Add(peptide.FullSequence);
                }
            }

            placeholder += 1;
            List<string> acceptorSamples = new List<string>
            {
                "A_14", "A_15", "A_3", "A_4", "A_5",
                "A_6", "A_7", "A_8", "A_9", "A_10",
                "A_11", "A_12", "A_13", "A_14", "A_15",
                "A_16", "A_17", "A_18", "A_19", "A_20"
            };
            var expFile = new MsFraggerExperimentFile(@"D:\GygiTwoProteome_PXD014415\MsConvertmzMLs\IonQuant_1Percent\experiment_annotation.tsv");
            List<MsFraggerExperiment> experiments = expFile.Results;

            List<MsFraggerPsm> censoredPsms = new();


            // Iterate through the psm.tsv file for each run
            // Select only human peptides with peptide prophet probability > 0.99
            // and group by modified peptide (create some field that defaults to base sequence if mod seq is NA)
            // Fragger only gives scan retention times (in seconds), so we'll need to map the RT to the scan number, which will suck
            // Randomly select 500 peptides to censor
            // Load in the raw file, iterate through each peptide group
            // for every psm in the group, censor the corresponding spectrum
            foreach (var exp in experiments.Where(exp => exp.FullFilePathWithExtension.Contains("_human_90min")))
            {
                var psms = new MsFraggerPsmFile(Path.Combine(@"D:\GygiTwoProteome_PXD014415\MsConvertmzMLs\IonQuant_1Percent", exp.Sample, "psm.tsv"));

                var psmsToRemove = psms
                    .Where(pep => pep.Protein.Contains("HUMAN") && !pep.Protein.Contains("REV") && !pep.Protein.Contains("CON")
                        && pep.PeptideProphetProbability > 0.99
                        && !masterDonorPeptidesByFile[exp.Sample].Contains(pep.ModifiedSequence))
                    .GroupBy(pep => pep.ModifiedSequence)
                    .Where(group => group.Count() == 1)
                    .OrderBy(group => Guid.NewGuid())
                    .SelectMany(group => group.Take(1))
                    .Take(500)
                    .ToList();

                censoredPsms.AddRange(psmsToRemove);

                // Open up the mzMl file specified by the psm FileName
                // Select all scans referenced by the psms
                // set the intensity of all peaks but 1 to zero
                // Save the new mzMl file
                string rawFilePath = Path.Combine(exp.FullFilePathWithExtension);
                var file = MsDataFileReader.GetDataFile(rawFilePath);
                file.LoadAllStaticData();

                foreach (var psm in psmsToRemove)
                {
                    MsDataScan scan = file.GetOneBasedScan(psm.OneBasedScanNumber);
                    scan.MassSpectrum = blankSpectrum;
                }

                string fileName = psmsToRemove.First().FileNameWithoutExtension;

                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(file,
                                       Path.Combine(@"D:\GygiTwoProteome_PXD014415\MsConvertmzMLs\CensoredFiles", fileName + "-censored.mzML"),
                                                          writeIndexed: true);
            }

            // Write out the list of censored PSMs 
            using (var csv =
                new CsvWriter(new StreamWriter(Path.Combine(@"D:\GygiTwoProteome_PXD014415\MsConvertmzMLs\CensoredFiles", "CensoredPsms.tsv")), MsFraggerPsm.CsvConfiguration))
            {
                csv.WriteHeader<MsFraggerPsm>();
                foreach (var result in censoredPsms)
                {
                    csv.NextRecord();
                    csv.WriteRecord(result);
                }
            }

        }

        [Test]
        public static void CensorFraggerInHouseHuman()
            
        {
            string peptideFilePath = @"D:\Human_Ecoli_TwoProteome_60minGradient\IonQuant_1Percent_mzML\combined_modified_peptide.tsv";

            var peptideFile = new MsFraggerPeptideFile(peptideFilePath);
            peptideFile.LoadResults();
            List<MsFraggerPeptide> parsedPeptides = peptideFile.Results;

            var testPept = parsedPeptides.First();

            MzSpectrum blankSpectrum = new MzSpectrum(mz: new double[] { 150 }, intensities: new double[] { 10000 }, false);
            int placeholder = 0;
            // Filter on peptideProphet probability, 0.99 or higher

            // For each peptide, select the run with the highest intensity. That run is the master donor and will not be censored
            // Construct a dictionary linking file/run to the list of master donor peptides
            Dictionary<string, List<string>> masterDonorPeptidesByFile = new();
            foreach (var peptide in parsedPeptides)
            {
                string masterDonor = peptide.IntensityByFile.MaxBy(kvp => kvp.Value).Key;
                if (!masterDonorPeptidesByFile.ContainsKey(masterDonor))
                {
                    masterDonorPeptidesByFile.Add(masterDonor, new List<string> { peptide.FullSequence });
                }
                else
                {
                    masterDonorPeptidesByFile[masterDonor].Add(peptide.FullSequence);
                }
            }

            placeholder += 1;
            List<string> acceptorSamples = new List<string> { "a_1", "a_2", "a_3", "a_4", "a_5", "a_6", "a_7", "a_8", "a_9", "a_10" };
            var expFile = new MsFraggerExperimentFile(@"D:\Human_Ecoli_TwoProteome_60minGradient\IonQuant_1Percent_mzML\experiment_annotation.tsv");
            List<MsFraggerExperiment> experiments = expFile.Results;

            List<MsFraggerPsm> censoredPsms = new();


            // Iterate through the psm.tsv file for each run
            // Select only human peptides with peptide prophet probability > 0.99
            // and group by modified peptide (create some field that defaults to base sequence if mod seq is NA)
            // Fragger only gives scan retention times (in seconds), so we'll need to map the RT to the scan number, which will suck
            // Randomly select 500 peptides to censor
            // Load in the raw file, iterate through each peptide group
            // for every psm in the group, censor the corresponding spectrum
            foreach (var exp in experiments.Where(exp => exp.FullFilePathWithExtension.Contains("Human_C18")))
            {
                var psms = new MsFraggerPsmFile(Path.Combine(@"D:\Human_Ecoli_TwoProteome_60minGradient\IonQuant_1Percent_mzML", exp.Sample, "psm.tsv"));

                var psmsToRemove = psms
                    .Where(pep => pep.Protein.Contains("HUMAN") && !pep.Protein.Contains("REV") && !pep.Protein.Contains("CON")
                        && pep.PeptideProphetProbability > 0.99
                        && !masterDonorPeptidesByFile[exp.Sample].Contains(pep.ModifiedSequence))
                    .GroupBy(pep => pep.ModifiedSequence)
                    .Where(group => group.Count() == 1)
                    .OrderBy(group => Guid.NewGuid())
                    .SelectMany(group => group.Take(1))
                    .Take(500)
                    .ToList();

                censoredPsms.AddRange(psmsToRemove);

                // Open up the mzMl file specified by the psm FileName
                // Select all scans referenced by the psms
                // set the intensity of all peaks but 1 to zero
                // Save the new mzMl file
                string rawFilePath = Path.Combine(exp.FullFilePathWithExtension);
                var file = MsDataFileReader.GetDataFile(rawFilePath);
                file.LoadAllStaticData();

                foreach (var psm in psmsToRemove)
                {
                    MsDataScan scan = file.GetOneBasedScan(psm.OneBasedScanNumber);
                    scan.MassSpectrum = blankSpectrum;
                }

                string fileName = psmsToRemove.First().FileNameWithoutExtension;

                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(file,
                                       Path.Combine(@"D:\Human_Ecoli_TwoProteome_60minGradient\RawData\CensoredFiles", fileName + "-censored.mzML"),
                                                          writeIndexed: true);
            }

            // Write out the list of censored PSMs 
            using (var csv =
                new CsvWriter(new StreamWriter(Path.Combine(@"D:\Human_Ecoli_TwoProteome_60minGradient\RawData\CensoredFiles", "CensoredPsms.tsv")), MsFraggerPsm.CsvConfiguration))
            {
                csv.WriteHeader<MsFraggerPsm>();
                foreach (var result in censoredPsms)
                {
                    csv.NextRecord();
                    csv.WriteRecord(result);
                }
            }

        }

    }
}
