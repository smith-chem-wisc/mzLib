using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
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

namespace Test.FileReadingTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class TestPsmFromTsv
    {

        [Test]
        [TestCase("oglycoSinglePsms.psmtsv", 2)] // oglyco
        [TestCase("nglyco_f5.psmtsv", 5)] // nglyco
        [TestCase("VariantCrossTest.psmtsv", 15)] // variant crossing
        [TestCase("XL_Intralinks.tsv", 6)] // variant crossing
        [TestCase("XLink.psmtsv", 19)] // variant crossing
        public static void TestPsmReaderWithMultipleEntryPoints(string path, int expected)
        {
            string psmFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\SearchResults",
                path);
            List<PsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFilePath, out var warnings);
            Assert.That(warnings.Count, Is.EqualTo(0));
            Assert.That(parsedPsms.Count, Is.EqualTo(expected));

            PsmFromTsvFile psmFile = FileReader.ReadFile<PsmFromTsvFile>(psmFilePath);
            List<PsmFromTsv> psmFilePsms = psmFile.Results;
            Assert.That(psmFilePsms.Count, Is.EqualTo(parsedPsms.Count));
            for (int i = 0; i < parsedPsms.Count; i++)
            {
                Assert.That(parsedPsms[i].FullSequence, Is.EqualTo(psmFilePsms[i].FullSequence));
            }


            SpectrumMatchFromTsvFile spectralMatchFile = FileReader.ReadFile<SpectrumMatchFromTsvFile>(psmFilePath);
            List<SpectrumMatchFromTsv> spectralMatchFilePsms = spectralMatchFile.Results;
            Assert.That(psmFilePsms.Count, Is.EqualTo(parsedPsms.Count));
            for (int i = 0; i < parsedPsms.Count; i++)
            {
                PsmFromTsv casted = spectralMatchFilePsms[i] as PsmFromTsv;
                if (casted == null) Assert.Fail();

                Assert.That(parsedPsms[i].FullSequence, Is.EqualTo(spectralMatchFilePsms[i].FullSequence));
            }
        }

        [Test]
        public static void ReadInternalIonsPsms()
        {
            string psmFile = @"FileReadingTests\SearchResults\internalIons.psmtsv";
            List<PsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFile, out var warnings);
            Assert.AreEqual(10, parsedPsms.Count);
            var internalIons = parsedPsms[0].MatchedIons.Where(i => i.Annotation.Contains("yIb")).ToList();
            Assert.AreEqual(97, internalIons.Count);
        }

        [Test]
        public static void ReadOGlycoPsmsLocalizedGlycans()
        {
            string psmFile = @"FileReadingTests\SearchResults\oglyco.psmtsv";
            List<PsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFile, out var warnings);
            Assert.AreEqual(9, parsedPsms.Count);

            // read glycans if applicable
            List<Tuple<int, string, double>> localGlycans = null;
            if (parsedPsms[0].GlycanLocalizationLevel != null)
            {
                localGlycans = PsmFromTsv.ReadLocalizedGlycan(parsedPsms[0].LocalizedGlycan);
            }

            Assert.AreEqual(1, localGlycans.Count);

        }

        [Test]
        public static void ReadExcelEditedPsms()
        {
            string psmFile = @"FileReadingTests\SearchResults\ExcelEditedPeptide.psmtsv";
            List<PsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFile, out var warnings);
            Assert.AreEqual(1, parsedPsms.Count);
            IEnumerable<string> expectedIons = new string[] { "y3+1", "y4+1", "b4+1", "b5+1", "b6+1", "b8+1" };
            Assert.That(6 == parsedPsms[0].MatchedIons.Select(p => p.Annotation).Intersect(expectedIons).Count());
        }

        [Test]
        public static void TestPsmFromTsvDisambiguatingConstructor()
        {
            // initialize values
            string psmTsvPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\SearchResults\TDGPTMDSearchResults.psmtsv");
            List<string> warnings = new();
            List<PsmFromTsv> psms = SpectrumMatchTsvReader.ReadPsmTsv(psmTsvPath, out warnings);
            PsmFromTsv psm = psms.First();

            // non ambiguous construction should not be successful
            string fullSeq = psm.FullSequence;
            fullSeq = fullSeq.Substring(0, fullSeq.Length - 1);
            PsmFromTsv modifiedPsm = new(psm, fullSeq);
            Assert.That(modifiedPsm.FullSequence == fullSeq);

            // disambiguation construction
            var ambiguousPsms = psms.Where(p => p.FullSequence.Contains('|'));
            PsmFromTsv ambiguousPsm = ambiguousPsms.First();
            var fullSeqStrings = ambiguousPsm.FullSequence.Split('|');

            PsmFromTsv modifiedAmbiguousPsm = new(ambiguousPsm, fullSeqStrings[0]);
            List<string[]> test = new();
            foreach (var ambPsm in ambiguousPsms)
            {
                PsmFromTsv disambiguatedPSM = new(ambPsm, ambPsm.FullSequence.Split("|")[0]);
                Assert.That(disambiguatedPSM.StartAndEndResiduesInProtein == ambPsm.StartAndEndResiduesInProtein.Split("|")[0]);
                Assert.That(disambiguatedPSM.BaseSeq == ambPsm.BaseSeq.Split("|")[0]);
                Assert.That(disambiguatedPSM.EssentialSeq == ambPsm.EssentialSeq.Split("|")[0]);
                Assert.That(disambiguatedPSM.ProteinAccession == ambPsm.ProteinAccession.Split("|")[0]);
                Assert.That(disambiguatedPSM.PeptideMonoMass == ambPsm.PeptideMonoMass.Split("|")[0]);
                Assert.That(disambiguatedPSM.MassDiffDa == ambPsm.MassDiffDa.Split("|")[0]);
                Assert.That(disambiguatedPSM.MassDiffPpm == ambPsm.MassDiffPpm.Split("|")[0]);
                Assert.That(disambiguatedPSM.ProteinName == ambPsm.ProteinName.Split("|")[0]);
                Assert.That(disambiguatedPSM.GeneName == ambPsm.GeneName.Split("|")[0]);

                for (int i = 0; i < ambPsm.MatchedIons.Count; i++)
                {
                    Assert.That(disambiguatedPSM.MatchedIons[i] == ambPsm.MatchedIons[i]);
                }

                if (ambPsm.StartAndEndResiduesInProtein.Split("|").Count() > 1)
                {
                    for (int i = 1; i < ambPsm.StartAndEndResiduesInProtein.Split("|").Count(); i++)
                    {
                        disambiguatedPSM = new(ambPsm, ambPsm.FullSequence.Split("|")[i], i);
                        Assert.That(disambiguatedPSM.StartAndEndResiduesInProtein == ambPsm.StartAndEndResiduesInProtein.Split("|")[i]);
                        Assert.That(disambiguatedPSM.BaseSeq == ambPsm.BaseSeq.Split("|")[i]);
                        Assert.That(disambiguatedPSM.EssentialSeq == ambPsm.EssentialSeq.Split("|")[i]);
                        Assert.That(disambiguatedPSM.ProteinAccession == ambPsm.ProteinAccession.Split("|")[i]);
                        Assert.That(disambiguatedPSM.ProteinName == ambPsm.ProteinName.Split("|")[i]);
                        Assert.That(disambiguatedPSM.GeneName == ambPsm.GeneName.Split("|")[i]);

                        if (ambPsm.PeptideMonoMass.Split("|").Count() == 1)
                        {
                            Assert.That(disambiguatedPSM.PeptideMonoMass == ambPsm.PeptideMonoMass.Split("|")[0]);
                            Assert.That(disambiguatedPSM.MassDiffDa == ambPsm.MassDiffDa.Split("|")[0]);
                            Assert.That(disambiguatedPSM.MassDiffPpm == ambPsm.MassDiffPpm.Split("|")[0]);
                        }
                        else
                        {
                            Assert.That(disambiguatedPSM.PeptideMonoMass == ambPsm.PeptideMonoMass.Split("|")[i]);
                            Assert.That(disambiguatedPSM.MassDiffDa == ambPsm.MassDiffDa.Split("|")[i]);
                            Assert.That(disambiguatedPSM.MassDiffPpm == ambPsm.MassDiffPpm.Split("|")[i]);
                        }
                    }
                }
            }
        }

        [Test]
        public static void TestParseModification()
        {
            string psmTsvPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\SearchResults\TDGPTMDSearchResults.psmtsv");
            List<string> warnings = new();
            List<PsmFromTsv> psms = SpectrumMatchTsvReader.ReadPsmTsv(psmTsvPath, out warnings).Take(20).ToList();
            Assert.That(warnings.Count == 2);
            Assert.AreEqual("Could not read line: 2", warnings[0]);
            Assert.AreEqual("Warning: 1 PSMs were not read.", warnings[1]);

            // psm with single modificaiton
            PsmFromTsv singleMod = psms[0];
            var modDict = Omics.SpectrumMatch.SpectrumMatchFromTsv.ParseModifications(singleMod.FullSequence);
            Assert.That(modDict.Count == 1);
            Assert.That(modDict.ContainsKey(37));
            Assert.That(modDict.Values.First().Contains("Common Fixed:Carbamidomethyl on C"));

            // psm with two modifications
            PsmFromTsv twoMods = psms[15];
            modDict = Omics.SpectrumMatch.SpectrumMatchFromTsv.ParseModifications(twoMods.FullSequence);
            Assert.That(modDict.Count == 2);
            Assert.That(modDict.ContainsKey(0) && modDict.ContainsKey(104));
            Assert.That(modDict[0].Count == 1);
            Assert.That(modDict[0].Contains("UniProt:N-acetylserine on S"));
            Assert.That(modDict[104].Count == 1);
            Assert.That(modDict[104].Contains("UniProt:N5-methylglutamine on Q"));


            // psm with two mods on the same amino acid
            string fullSeq = "[Common Fixed:Carbamidomethyl on C]|[UniProt:N-acetylserine on S]KPRKIEEIKDFLLTARRKDAKSVKIKKNKDNVKFK";
            modDict = Omics.SpectrumMatch.SpectrumMatchFromTsv.ParseModifications(fullSeq);
            Assert.That(modDict.Count == 1);
            Assert.That(modDict.ContainsKey(0));
            Assert.That(modDict[0].Count == 2);
            Assert.That(modDict[0].Contains("Common Fixed:Carbamidomethyl on C"));
            Assert.That(modDict[0].Contains("UniProt:N-acetylserine on S"));
        }

        [Test]
        public static void TestRemoveSpecialCharacters()
        {
            // successful removal of the default character
            string toRemove = "ANDVHAO|CNVASDF|ABVCUAE";
            int length = toRemove.Length;
            Omics.SpectrumMatch.SpectrumMatchFromTsv.RemoveSpecialCharacters(ref toRemove);
            Assert.That(toRemove.Length == length - 2);
            Assert.That(toRemove.Equals("ANDVHAOCNVASDFABVCUAE"));

            // does not remove default character when prompted otherwise
            toRemove = "ANDVHAO|CNVASDF|ABVCUAE";
            Omics.SpectrumMatch.SpectrumMatchFromTsv.RemoveSpecialCharacters(ref toRemove, specialCharacter: @"\[");
            Assert.That(toRemove.Length == length);
            Assert.That(toRemove.Equals("ANDVHAO|CNVASDF|ABVCUAE"));

            // replaces default symbol when prompted
            Omics.SpectrumMatch.SpectrumMatchFromTsv.RemoveSpecialCharacters(ref toRemove, replacement: @"%");
            Assert.That(toRemove.Length == length);
            Assert.That(toRemove.Equals("ANDVHAO%CNVASDF%ABVCUAE"));

            // replaces inputted symbol with non-default symbol
            Omics.SpectrumMatch.SpectrumMatchFromTsv.RemoveSpecialCharacters(ref toRemove, replacement: @"=", specialCharacter: @"%");
            Assert.That(toRemove.Length == length);
            Assert.That(toRemove.Equals("ANDVHAO=CNVASDF=ABVCUAE"));
        }

        [Test]
        public static void TestToString()
        {
            string psmTsvPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\SearchResults\TDGPTMDSearchResults.psmtsv");
            List<string> warnings = new();
            List<PsmFromTsv> psms = SpectrumMatchTsvReader.ReadPsmTsv(psmTsvPath, out warnings).Take(3).ToList();
            Assert.That(warnings.Count == 2);
            Assert.AreEqual("Could not read line: 2", warnings[0]);
            Assert.AreEqual("Warning: 1 PSMs were not read.", warnings[1]);

            Assert.That(psms[0].FullSequence.Equals(psms[0].ToString()));
            Assert.That(psms[1].FullSequence.Equals(psms[1].ToString()));
            Assert.That(psms[2].FullSequence.Equals(psms[2].ToString()));
        }

        [Test]
        public static void TestParenthesesRemovalForSilac()
        {
            string baseSequence = "ASDF(+8.01)ASDF";
            string cleanedSequence = Omics.SpectrumMatch.SpectrumMatchFromTsv.RemoveParentheses(baseSequence);
            Assert.IsTrue(cleanedSequence.Equals("ASDFASDF"));
        }

        [Test]
        public static void TestSimpleToLibrarySpectrum()
        {
            string psmTsvPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\SearchResults\TDGPTMDSearchResults.psmtsv");
            List<string> warnings = new();
            List<PsmFromTsv> psms = SpectrumMatchTsvReader.ReadPsmTsv(psmTsvPath, out warnings).Take(3).ToList();
            Assert.That(warnings.Count == 2);
            Assert.AreEqual("Could not read line: 2", warnings[0]);
            Assert.AreEqual("Warning: 1 PSMs were not read.", warnings[1]);

            string librarySpectrum = psms[0].ToLibrarySpectrum().ToString();

            string expectedLibrarySpectrum = File.ReadAllText(Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\SearchResults\simple.msp"));

            //not a great way to test equality but we are experiencing a great deal of 10th digit rounding differences
            Assert.AreEqual(Regex.Matches(expectedLibrarySpectrum, "ppm").Count, Regex.Matches(librarySpectrum, "ppm").Count);


            //the code below tests the addition and correct output for neutral loss fragments
            Product p = new Product(ProductType.bWaterLoss, FragmentationTerminus.N, 1, 1, 1, 18);
            MatchedFragmentIon matchedIon = new(p, 1, 1, 1);
            psms[0].MatchedIons.Add(matchedIon);
            string librarySpectrumWithNeutralLoss = psms[0].ToLibrarySpectrum().ToString();

            Assert.That(librarySpectrumWithNeutralLoss.Contains("WaterLoss"));
        }

        [Test]
        [TestCase("FileReader - PsmFromTsv")]
        [TestCase("FileReader - SpectrumMatchFromTsv")]
        [TestCase("File Construction - PsmFromTsv")]
        [TestCase("File Construction - SpectrumMatchFromTsv")]
        public static void TestPsmFiles(string fileLoadingType)
        {
            string psmTsvPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\SearchResults\TDGPTMDSearchResults.psmtsv");
            List<PsmFromTsv> psms = SpectrumMatchTsvReader.ReadPsmTsv(psmTsvPath, out var warnings);
            Assert.That(warnings.Count == 2);

            IResultFile loadedFile = null;
            switch (fileLoadingType)
            {
                case "FileReader - PsmFromTsv":
                    var file = FileReader.ReadFile<PsmFromTsvFile>(psmTsvPath);
                    file.LoadResults();
                    Assert.That(file.Results.Count == psms.Count);
                    loadedFile = file;
                    break;

                case "FileReader - SpectrumMatchFromTsv":
                    var file2 = FileReader.ReadFile<SpectrumMatchFromTsvFile>(psmTsvPath);
                    file2.LoadResults();
                    Assert.That(file2.Results.Count == psms.Count);
                    loadedFile = file2;
                    break;

                case "File Construction - PsmFromTsv":
                    var file3 = new PsmFromTsvFile(psmTsvPath);
                    file3.LoadResults();
                    Assert.That(file3.Results.Count == psms.Count);
                    loadedFile = file3;
                    break;

                case "File Construction - SpectrumMatchFromTsv":
                    var file4 = new SpectrumMatchFromTsvFile(psmTsvPath);
                    file4.LoadResults();
                    Assert.That(file4.Results.Count == psms.Count);
                    loadedFile = file4;
                    break;

                default:
                    Assert.Fail();
                    break;
            }

            Assert.That(loadedFile.FileType == SupportedFileType.psmtsv);
            Assert.Throws<NotImplementedException>(() => { loadedFile.WriteResults(""); });
        }

        [Test]
        public static void ReadPsmFromTsvWithNewHeaderTerms()
        {
            string psmTsvPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\SearchResults\NewHeader.psmtsv");
            List<PsmFromTsv> psms = SpectrumMatchTsvReader.ReadPsmTsv(psmTsvPath, out var warnings);
            Assert.That(warnings.Count == 2);

            IResultFile loadedFile = null;
            var file = FileReader.ReadFile<PsmFromTsvFile>(psmTsvPath);
            file.LoadResults();
            Assert.That(file.Results.Count == psms.Count);
            loadedFile = file;
        }

        [Test]
        public static void ReadXLinkPsmFromTsv()
        {
            string psmTsvPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\SearchResults\XLink.psmtsv");
            List<PsmFromTsv> psms = SpectrumMatchTsvReader.ReadPsmTsv(psmTsvPath, out var warnings);

            IResultFile loadedFile = null;
            var file = FileReader.ReadFile<PsmFromTsvFile>(psmTsvPath);
            file.LoadResults();
            Assert.That(file.Results.Count == psms.Count);
            loadedFile = file;
        }
    }
}
