using CsvHelper;
using Easy.Common.Extensions;
using MassSpectrometry;
using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace Test.FileReadingTests
{
    public class TestMaxQuantResultFile
    {
        [TestFixture]
        [ExcludeFromCodeCoverage]
        public class MaxQuantEvidenceTests
        {
            [Test]
            public void TestMaxQuantEvidenceResultFile()
            {
                var path = @"D:\Kelly_TwoProteomeData\combined_MaxQuant_Kelly\txt\evidence.txt";
                var file = new MaxQuantEvidenceFile(path);
                file.LoadResults();
                var results = file.Results.ToList();

                Assert.That(results.IsNotNullOrEmpty(), "No results were loaded from the file");

                var mbrHits = results.Where(r => r.MatchScore != null).ToList();

                var humanProteins = ProteinDbLoader.LoadProteinFasta(@"D:\Kelly_TwoProteomeData\combined_MaxQuant_Kelly\txt\uniprot_HSapiens_80k_03_2024.fasta",
                    true, DecoyType.None, isContaminant: false, out var errors);

                HashSet<string> proteinAccesions = new HashSet<string>(humanProteins.Select(p => p.Accession));

                var humanMbrHits = mbrHits
                    .Where(r => r.ProteinList.All(p => proteinAccesions.Contains(p)))
                    .ToList();

                int placeholder = 0;

            }

            [Test]
            public void CensorKellyMaxQuant()
            {
                // Read in evidence
                // Keep only human peptides with Q_value above cutoff (0.01)
                // Group by peptide mod sequence
                // Keep only groups with more than one peptide
               
                var resultsPath = @"D:\Kelly_TwoProteomeData\combined_MaxQuant_Kelly\txt\evidence.txt";
                var file = new MaxQuantEvidenceFile(resultsPath);
                file.LoadResults();
                var results = file.Results.ToList();

                var humanProteins = ProteinDbLoader.LoadProteinFasta(@"D:\Kelly_TwoProteomeData\combined_MaxQuant_Kelly\txt\uniprot_HSapiens_80k_03_2024.fasta",
                    true, DecoyType.None, isContaminant: false, out var errors);

                HashSet<string> proteinAccesions = new HashSet<string>(humanProteins.Select(p => p.Accession));

                results = results.Where(r => r.ProteinList.All(p => 
                    proteinAccesions.Contains(p)) 
                    && r.PEP < 0.01
                    && !r.PotentialContaminant.Equals("+")
                    && !r.Reverse.Equals("+"))
                   .ToList();

                // Then, write a dictionary that pairs the files with the list of top-scoring peptides that originated from that file
                Dictionary<string, List<string>> masterDonorPeptidesByFile = new();
                foreach (var peptide in results.GroupBy(r => r.ModifiedSequence))
                {
                    var bestPep = peptide.MaxBy(kvp => kvp.Score);
                    if (!masterDonorPeptidesByFile.ContainsKey(bestPep.FileNameWithoutExtension))
                    {
                        masterDonorPeptidesByFile.Add(bestPep.FileNameWithoutExtension, new List<string> { bestPep.ModifiedSequence });
                    }
                    else
                    {
                        masterDonorPeptidesByFile[bestPep.FileNameWithoutExtension].Add(bestPep.ModifiedSequence);
                    }
                }


                // Finally, iterate through all peptides for each file,
                // take 500 randomly if they're not in the top-scoring list for that file

                string rawFileFolder = @"D:\Kelly_TwoProteomeData\NewDataApril2024\DDA_sc_hela";
                string censoredFileFolder = @"D:\Kelly_TwoProteomeData\CensoredDataFiles_MaxQuant";
                MzSpectrum blankSpectrum = new MzSpectrum(mz: new double[] { 150 }, intensities: new double[] { 10000 }, false);
                List<MaxQuantEvidence> censoredResults = new();

                foreach (var filePeptideGroup in results.GroupBy(r => r.FileNameWithoutExtension))
                {
                    if (filePeptideGroup.Key.Contains("HeYe"))
                    {
                        continue;
                    }

                    var filePeptides = filePeptideGroup
                        .Where(
                            r => r.PsmCount == 1 && 
                            !masterDonorPeptidesByFile[filePeptideGroup.Key].Contains(r.Modifications))
                        .OrderBy(r => Guid.NewGuid())
                        .Take(500)
                        .ToList();
           
                    censoredResults.AddRange(filePeptides);

                    // Open the raw file
                    // for each peptide in filePeptide, find the associated scan and replace it with a blank spectrum
                    
                    var rawFile = Path.Combine(rawFileFolder, filePeptideGroup.Key + ".raw");
                    var reader = MsDataFileReader.GetDataFile(rawFile);
                    reader.LoadAllStaticData();

                    foreach(var peptide in filePeptides)
                    {
                        var scan = reader.GetOneBasedScan((int)peptide.Ms2ScanNumber);
                        scan.MassSpectrum = blankSpectrum;
                    }

                    string fileName = Path.Combine(censoredFileFolder, filePeptideGroup.Key + "_censored.mzML");
                    MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(reader, fileName, writeIndexed: true);
                }


                // Write out the list of censored PSMs 
                using (var csv =
                    new CsvWriter(new StreamWriter(Path.Combine(@"D:\Kelly_TwoProteomeData\CensoredDataFiles_MaxQuant", "CensoredPsms.tsv")), MaxQuantEvidence.CsvConfiguration))
                {
                    csv.WriteHeader<MaxQuantEvidence>();
                    foreach (var result in censoredResults)
                    {
                        csv.NextRecord();
                        csv.WriteRecord(result);
                    }
                }

            }


            [Test]
            public void CensorGygiMaxQuant()
            {
                // Read in evidence
                // Keep only human peptides with Q_value above cutoff (0.01)
                // Group by peptide mod sequence
                // Keep only groups with more than one peptide

                var resultsPath = @"D:\GygiTwoProteome_PXD014415\combined_MaxQuant_Gygi\txt\evidence.txt";
                var file = new MaxQuantEvidenceFile(resultsPath);
                file.LoadResults();
                var results = file.Results.ToList();

                var humanProteins = ProteinDbLoader.LoadProteinFasta(@"D:\Kelly_TwoProteomeData\combined_MaxQuant_Kelly\txt\uniprot_HSapiens_80k_03_2024.fasta",
                    true, DecoyType.None, isContaminant: false, out var errors);

                HashSet<string> proteinAccesions = new HashSet<string>(humanProteins.Select(p => p.Accession));

                results = results.Where(r => r.ProteinList.All(p =>
                    proteinAccesions.Contains(p))
                    && r.PEP < 0.01
                    && !r.PotentialContaminant.Equals("+")
                    && !r.Reverse.Equals("+"))
                   .ToList();

                // Then, write a dictionary that pairs the files with the list of top-scoring peptides that originated from that file
                Dictionary<string, List<string>> masterDonorPeptidesByFile = new();
                foreach (var peptide in results.GroupBy(r => r.ModifiedSequence))
                {
                    var bestPep = peptide.MaxBy(kvp => kvp.Score);
                    if (!masterDonorPeptidesByFile.ContainsKey(bestPep.FileNameWithoutExtension))
                    {
                        masterDonorPeptidesByFile.Add(bestPep.FileNameWithoutExtension, new List<string> { bestPep.ModifiedSequence });
                    }
                    else
                    {
                        masterDonorPeptidesByFile[bestPep.FileNameWithoutExtension].Add(bestPep.ModifiedSequence);
                    }
                }


                // Finally, iterate through all peptides for each file,
                // take 500 randomly if they're not in the top-scoring list for that file

                string rawFileFolder = @"D:\GygiTwoProteome_PXD014415\Raw_Files";
                string censoredFileFolder = @"D:\GygiTwoProteome_PXD014415\CensoredDataFiles_MaxQuant";
                MzSpectrum blankSpectrum = new MzSpectrum(mz: new double[] { 150 }, intensities: new double[] { 10000 }, false);
                List<MaxQuantEvidence> censoredResults = new();

                foreach (var filePeptideGroup in results.GroupBy(r => r.FileNameWithoutExtension))
                {
                    if (filePeptideGroup.Key.Contains("yaeast"))
                    {
                        continue;
                    }

                    var filePeptides = filePeptideGroup
                        .Where(
                            r => r.PsmCount == 1 &&
                            !masterDonorPeptidesByFile[filePeptideGroup.Key].Contains(r.Modifications))
                        .OrderBy(r => Guid.NewGuid())
                        .Take(500)
                        .ToList();

                    censoredResults.AddRange(filePeptides);

                    // Open the raw file
                    // for each peptide in filePeptide, find the associated scan and replace it with a blank spectrum

                    var rawFile = Path.Combine(rawFileFolder, filePeptideGroup.Key + ".raw");
                    var reader = MsDataFileReader.GetDataFile(rawFile);
                    reader.LoadAllStaticData();

                    foreach (var peptide in filePeptides)
                    {
                        var scan = reader.GetOneBasedScan((int)peptide.Ms2ScanNumber);
                        scan.MassSpectrum = blankSpectrum;
                    }

                    string fileName = Path.Combine(censoredFileFolder, filePeptideGroup.Key + "_censored.mzML");
                    MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(reader, fileName, writeIndexed: true);
                }


                // Write out the list of censored PSMs 
                using (var csv =
                    new CsvWriter(new StreamWriter(Path.Combine(censoredFileFolder, "CensoredPsms.tsv")), MaxQuantEvidence.CsvConfiguration))
                {
                    csv.WriteHeader<MaxQuantEvidence>();
                    foreach (var result in censoredResults)
                    {
                        csv.NextRecord();
                        csv.WriteRecord(result);
                    }
                }

            }


            [Test]
            public void CensorInHouseMaxQuant()
            {
                // Read in evidence
                // Keep only human peptides with Q_value above cutoff (0.01)
                // Group by peptide mod sequence
                // Keep only groups with more than one peptide

                var resultsPath = @"D:\Human_Ecoli_TwoProteome_60minGradient\combined_MaxQuant_InHouse\txt\evidence.txt";
                var file = new MaxQuantEvidenceFile(resultsPath);
                file.LoadResults();
                var results = file.Results.ToList();

                var humanProteins = ProteinDbLoader.LoadProteinFasta(@"D:\Kelly_TwoProteomeData\combined_MaxQuant_Kelly\txt\uniprot_HSapiens_80k_03_2024.fasta",
                    true, DecoyType.None, isContaminant: false, out var errors);

                HashSet<string> proteinAccesions = new HashSet<string>(humanProteins.Select(p => p.Accession));

                results = results.Where(r => r.ProteinList.All(p =>
                    proteinAccesions.Contains(p))
                    && r.PEP < 0.01
                    && !r.PotentialContaminant.Equals("+")
                    && !r.Reverse.Equals("+"))
                   .ToList();

                // Then, write a dictionary that pairs the files with the list of top-scoring peptides that originated from that file
                Dictionary<string, List<string>> masterDonorPeptidesByFile = new();
                foreach (var peptide in results.GroupBy(r => r.ModifiedSequence))
                {
                    var bestPep = peptide.MaxBy(kvp => kvp.Score);
                    if (!masterDonorPeptidesByFile.ContainsKey(bestPep.FileNameWithoutExtension))
                    {
                        masterDonorPeptidesByFile.Add(bestPep.FileNameWithoutExtension, new List<string> { bestPep.ModifiedSequence });
                    }
                    else
                    {
                        masterDonorPeptidesByFile[bestPep.FileNameWithoutExtension].Add(bestPep.ModifiedSequence);
                    }
                }


                // Finally, iterate through all peptides for each file,
                // take 500 randomly if they're not in the top-scoring list for that file

                string rawFileFolder = @"D:\Human_Ecoli_TwoProteome_60minGradient\RawData";
                string censoredFileFolder = @"D:\Human_Ecoli_TwoProteome_60minGradient\CensoredHumanData_MaxQuant";
                MzSpectrum blankSpectrum = new MzSpectrum(mz: new double[] { 150 }, intensities: new double[] { 10000 }, false);
                List<MaxQuantEvidence> censoredResults = new();

                foreach (var filePeptideGroup in results.GroupBy(r => r.FileNameWithoutExtension))
                {
                    if (filePeptideGroup.Key.Contains("Ecoli_10to1"))
                    {
                        continue;
                    }

                    var filePeptides = filePeptideGroup
                        .Where(
                            r => r.PsmCount == 1 &&
                            !masterDonorPeptidesByFile[filePeptideGroup.Key].Contains(r.Modifications))
                        .OrderBy(r => Guid.NewGuid())
                        .Take(500)
                        .ToList();

                    censoredResults.AddRange(filePeptides);

                    // Open the raw file
                    // for each peptide in filePeptide, find the associated scan and replace it with a blank spectrum

                    var rawFile = Path.Combine(rawFileFolder, filePeptideGroup.Key + ".raw");
                    var reader = MsDataFileReader.GetDataFile(rawFile);
                    reader.LoadAllStaticData();

                    foreach (var peptide in filePeptides)
                    {
                        var scan = reader.GetOneBasedScan((int)peptide.Ms2ScanNumber);
                        scan.MassSpectrum = blankSpectrum;
                    }

                    string fileName = Path.Combine(censoredFileFolder, filePeptideGroup.Key + "_censored.mzML");
                    MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(reader, fileName, writeIndexed: true);
                }


                // Write out the list of censored PSMs 
                using (var csv =
                    new CsvWriter(new StreamWriter(Path.Combine(censoredFileFolder, "CensoredPsms.tsv")), MaxQuantEvidence.CsvConfiguration))
                {
                    csv.WriteHeader<MaxQuantEvidence>();
                    foreach (var result in censoredResults)
                    {
                        csv.NextRecord();
                        csv.WriteRecord(result);
                    }
                }

            }
        }   
    }
}
