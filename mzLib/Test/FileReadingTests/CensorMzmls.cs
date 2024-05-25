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
                   .Where(psm => !peptideSequencesByFile[kvp.Key].Contains(psm.FullSequence) && psm.OrganismName.Equals("Homo sapiens"))
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

            string mzmlFolder = @"D:\GygiTwoProteome_PXD014415\MM105_Calibration_Search_MBR\Task1-CalibrateTask";
            MzSpectrum blankSpectrum = new MzSpectrum(mz: new double[] { 150 }, intensities: new double[] { 10000 }, false);
            string outputFolder = @"D:\GygiTwoProteome_PXD014415\CensoredDataFiles";

            foreach (var kvp in psmsByFile)
            {
                if (kvp.Key.Contains("yaeast")) continue;

                var peptideGroups = kvp.Value
                   // If this file contains the top peptide instance, we can't remove it
                   .Where(psm => !peptideSequencesByFile[kvp.Key].Contains(psm.FullSequence) && psm.OrganismName.Equals("Homo sapiens"))
                   .GroupBy(psm => psm.FullSequence)
                   .OrderBy(group => group.Key) // Order alphabetically, that way we don't only select the top scoring psms
                   .Where(group => allPeptideSeqs.Contains(group.Key) && group.Count() == 1)
                   .OrderBy(group => Guid.NewGuid())
                   .Take(1000);

                // Open up the mzMl file specified by the psm FileName
                // Select all scans referenced by the psms
                // set the intensity of all peaks but 1 to zero
                // Save the new mzMl file
                string mzmlFullPath = Path.Combine(mzmlFolder, peptideGroups.First().First().FileNameWithoutExtension + ".mzML");
                var file = MsDataFileReader.GetDataFile(mzmlFullPath);
                file.LoadAllStaticData();

                foreach (var group in peptideGroups)
                {
                    foreach (var psm in group)
                    {
                        MsDataScan scan = file.GetOneBasedScan(psm.Ms2ScanNumber);
                        scan.MassSpectrum = blankSpectrum;
                    }
                }

                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(file,
                    Path.Combine(outputFolder, peptideGroups.First().First().FileNameWithoutExtension + "-censored.mzML"),
                    writeIndexed: true);

            }

        }

    }
}
