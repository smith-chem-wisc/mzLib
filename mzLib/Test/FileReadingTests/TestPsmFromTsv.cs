﻿using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using NUnit.Framework.Legacy;
using Omics.Fragmentation;
using Readers;
using MzLibUtil;
using FlashLFQ;
using NUnit.Framework.Interfaces;

namespace Test.FileReadingTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class TestPsmFromTsv
    {

        [Test]
        [TestCase("oglycoSinglePsms.psmtsv", 2)] // oglyco
        [TestCase("oGlycoAllPsms.psmtsv", 10)] // oglyco - AllPsms
        [TestCase("nglyco_f5.psmtsv", 5)] // nglyco
        [TestCase("nglyco_f5_NewVersion.psmtsv", 5)]
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

            var spectralMatchFile = FileReader.ReadFile<SpectrumMatchFromTsvFile>(psmFilePath);
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
        public static void ReadInternalIonsPsms_IntensitiesCorrect()
        {
            string psmFile = @"FileReadingTests\SearchResults\internalIons.psmtsv";
            List<PsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFile, out var warnings);
            Assert.AreEqual(10, parsedPsms.Count);

            var internalIons = parsedPsms[0].MatchedIons.Where(i => i.Annotation.Contains("yIb")).ToList();
            Assert.AreEqual(97, internalIons.Count);

            // Example: check a few known intensities from the file
            // yIb[2-6]+1:286.11465:1:177
            var ion = internalIons.FirstOrDefault(i => i.Annotation.StartsWith("yIb[2-6]+1"));
            Assert.IsNotNull(ion);
            Assert.AreEqual(177, ion.Intensity, 1e-6);

            // yIb[2-9]+1:457.17603:1:140
            ion = internalIons.FirstOrDefault(i => i.Annotation.StartsWith("yIb[2-9]+1"));
            Assert.IsNotNull(ion);
            Assert.AreEqual(140, ion.Intensity, 1e-6);

            // yIb[2-10]+1:514.20691:1:168
            ion = internalIons.FirstOrDefault(i => i.Annotation.StartsWith("yIb[2-10]+1"));
            Assert.IsNotNull(ion);
            Assert.AreEqual(168, ion.Intensity, 1e-6);
        }

        [Test]
        public static void ReadInternalIonsPsms_MassErrorCorrect()
        {
            string psmFile = @"FileReadingTests\SearchResults\internalIons.psmtsv";
            List<PsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFile, out var warnings);
            Assert.AreEqual(10, parsedPsms.Count);

            var internalIons = parsedPsms[0].MatchedIons.Where(i => i.Annotation.Contains("yIb")).ToList();
            Assert.AreEqual(97, internalIons.Count);

            // Example: check a few known mass errors from the file
            // yIb[2-6]+1:286.11465:1:177:-0.00006
            var ion = internalIons.FirstOrDefault(i => i.Annotation.StartsWith("yIb[2-6]+1"));
            Assert.IsNotNull(ion);
            Assert.AreEqual(0.00006, ion.MassErrorDa, 1e-6);

            // yIb[2-9]+1:457.17603:1:140:-0.00296
            ion = internalIons.FirstOrDefault(i => i.Annotation.StartsWith("yIb[2-9]+1"));
            Assert.IsNotNull(ion);
            Assert.AreEqual(-0.00296, ion.MassErrorDa, 1e-6);

            // yIb[2-10]+1:514.20691:1:168:0.00646
            ion = internalIons.FirstOrDefault(i => i.Annotation.StartsWith("yIb[2-10]+1"));
            Assert.IsNotNull(ion);
            Assert.AreEqual(0.00646, ion.MassErrorDa, 1e-6);

            // Check MassErrorPpm as well (values from file: 0.21, -6.49, 12.59, etc.)
            Assert.AreEqual(0.21, internalIons.First(i => i.Annotation.StartsWith("yIb[2-6]+1")).MassErrorPpm, 1e-2);
            Assert.AreEqual(-6.49, internalIons.First(i => i.Annotation.StartsWith("yIb[2-9]+1")).MassErrorPpm, 1e-2);
            Assert.AreEqual(12.59, internalIons.First(i => i.Annotation.StartsWith("yIb[2-10]+1")).MassErrorPpm, 1e-2);
        }

        [Test]
        public static void TestOneOverK0Reading()
        {
            string psmFile = @"FileReadingTests\SearchResults\OneOverK0Example.psmtsv";
            List<PsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFile, out var warnings);
            Assert.AreEqual(14, parsedPsms.Count);
            var oneOverK0 = parsedPsms[0].OneOverK0;
            Assert.That(oneOverK0, Is.EqualTo(1.1068).Within(0.0001));
            var clonedPsm = new PsmFromTsv(parsedPsms[0], parsedPsms[0].FullSequence, 0);
            Assert.That(clonedPsm.OneOverK0, Is.EqualTo(1.1068).Within(0.0001));
        }

        [Test]
        [TestCase(@"FileReadingTests\SearchResults\oglyco.psmtsv")]
        [TestCase(@"FileReadingTests\SearchResults\oglyco_NewVersion.psmtsv")]
        public static void ReadOGlycoPsms_LegacyVersion(string psmFile)
        {
            // In this test, we are testing that the glyco-specific fields are read correctly from older versions of the oglyco.psmtsv file.
            // There are few headers edited in the newer versions. In order to maintain backward compatibility, we need to ensure that the older files are still read correctly.
            List<GlycoPsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadGlycoPsmTsv(psmFile, out var warnings);
            Assert.AreEqual(9, parsedPsms.Count);
            Assert.That(parsedPsms.All(p=>p.FlankingResidues != null));
            Assert.That(parsedPsms.All(p => p.TotalGlycanSites != null));
            Assert.That(parsedPsms.All(p => p.NGlycanMotifCheck != null));
            Assert.That(parsedPsms.All(p => p.AllPotentialGlycanLocalization != null));
            Assert.That(parsedPsms.All(p => p.AllSiteSpecificLocalizationProbability != null));
        }

        [Test]
        [TestCase(@"FileReadingTests\SearchResults\oglyco.psmtsv")]
        [TestCase(@"FileReadingTests\SearchResults\oglyco_NewVersion.psmtsv")]
        [TestCase(@"FileReadingTests\SearchResults\oglycoWithWrongExtension.tsv")]
        public static void ReadOGlycoPsmsLocalizedGlycans_Glyco(string psmFile)
        {
            // Using direct glycopsm reader
            List<GlycoPsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadGlycoPsmTsv(psmFile, out var warnings);
            Assert.AreEqual(9, parsedPsms.Count);

            // read glycans if applicable
            List<Tuple<int, string, double>> localGlycans = null;
            if (parsedPsms[0].GlycanLocalizationLevel != null)
            {
                localGlycans = GlycoPsmFromTsv.ReadLocalizedGlycan(parsedPsms[0].LocalizedGlycanInPeptide);
            }

            Assert.AreEqual(1, localGlycans.Count);
        }

        [Test]
        [TestCase(@"FileReadingTests\SearchResults\oglyco.psmtsv")]
        [TestCase(@"FileReadingTests\SearchResults\oglyco_NewVersion.psmtsv")]
        [TestCase(@"FileReadingTests\SearchResults\oglycoWithWrongExtension.tsv")]
        public static void ReadOGlycoPsmsLocalizedGlycans_Generic(string psmFile)
        {
            // Using generic psm reader and casting
            var parsedPsms = SpectrumMatchTsvReader.ReadTsv(psmFile, out var warnings).Cast<GlycoPsmFromTsv>().ToList();
            Assert.AreEqual(9, parsedPsms.Count);

            // read glycans if applicable
            List<Tuple<int, string, double>> localGlycans = null;
            if (parsedPsms[0].GlycanLocalizationLevel != null)
            {
                localGlycans = GlycoPsmFromTsv.ReadLocalizedGlycan(parsedPsms[0].LocalizedGlycanInPeptide);
            }

            Assert.AreEqual(1, localGlycans.Count);
        }

        [Test]
        [TestCase(@"FileReadingTests\SearchResults\oglyco.psmtsv")]
        [TestCase(@"FileReadingTests\SearchResults\oglyco_NewVersion.psmtsv")]
        [TestCase(@"FileReadingTests\SearchResults\oglycoWithWrongExtension.tsv")]
        public static void ReadOGlycoPsmsLocalizedGlycans_AsPsm(string psmFile)
        {
            // Using generic psm reader and casting
            var parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFile, out var warnings).Cast<GlycoPsmFromTsv>().ToList();
            Assert.AreEqual(9, parsedPsms.Count);

            // read glycans if applicable
            List<Tuple<int, string, double>> localGlycans = null;
            if (parsedPsms[0].GlycanLocalizationLevel != null)
            {
                localGlycans = GlycoPsmFromTsv.ReadLocalizedGlycan(parsedPsms[0].LocalizedGlycanInPeptide);
            }

            Assert.AreEqual(1, localGlycans.Count);
        }

        [Test]
        [TestCase(@"FileReadingTests\SearchResults\oglyco.psmtsv")]
        [TestCase(@"FileReadingTests\SearchResults\oglyco_NewVersion.psmtsv")]
        public static void GlcyoPsmFromTsv_ModifiedPsmFunction(string path)
        {
            // Test that the GlycoPsmFromTsv modified PSM constructor works as intended
            // The modified PSM should have the same glycan information as the unmodified PSM
            var unModifiedPsm = SpectrumMatchTsvReader.ReadPsmTsv(path, out var warnings).Cast<GlycoPsmFromTsv>().First();
            string newFullSequence = "REVEDPQVAQLELGGGPGAGDLQT[O-Glycosylation:H1N1A2 on T]LALEVAQQKR";
            var modifiedPsm = new GlycoPsmFromTsv(unModifiedPsm, newFullSequence);
            Assert.AreEqual(newFullSequence, modifiedPsm.FullSequence); // make sure the full sequence is updated
            // make sure all glycan information is the same
            Assert.AreEqual(unModifiedPsm.BaseSeq, modifiedPsm.BaseSeq);
            Assert.AreEqual(unModifiedPsm.GlycanComposition, modifiedPsm.GlycanComposition);
            Assert.AreEqual(unModifiedPsm.GlycanMass, modifiedPsm.GlycanMass);
            Assert.AreEqual(unModifiedPsm.GlycanStructure, modifiedPsm.GlycanStructure);
            Assert.AreEqual(unModifiedPsm.GlycanLocalizationLevel, modifiedPsm.GlycanLocalizationLevel);
            Assert.AreEqual(unModifiedPsm.LocalizedGlycanInPeptide, modifiedPsm.LocalizedGlycanInPeptide);
            Assert.AreEqual(unModifiedPsm.LocalizedGlycanInProtein, modifiedPsm.LocalizedGlycanInProtein);
            Assert.AreEqual(unModifiedPsm.R138144, modifiedPsm.R138144);
            Assert.AreEqual(unModifiedPsm.LocalizedScores, modifiedPsm.LocalizedScores);
            Assert.AreEqual(unModifiedPsm.YionScore, modifiedPsm.YionScore);
            Assert.AreEqual(unModifiedPsm.DiagonosticIonScore, modifiedPsm.DiagonosticIonScore);
            Assert.AreEqual(unModifiedPsm.TotalGlycanSites, modifiedPsm.TotalGlycanSites);
        }

        [Test]
        public static void GlycoPsmReader_CreatesNormalPsmsForAllPSMFileWhenAppropriate()
        {
            int expectedGlyco = 9;
            int expectedNormal = 1;
            string path = @"FileReadingTests\SearchResults\oGlycoAllPsms.psmtsv";

            var parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(path, out var warnings);
            Assert.That(warnings.Count, Is.EqualTo(0));
            Assert.AreEqual(expectedGlyco + expectedNormal, parsedPsms.Count);

            var glyco = parsedPsms.OfType<GlycoPsmFromTsv>().ToList();
            var normal = parsedPsms.Except(glyco).ToList();
            Assert.AreEqual(expectedGlyco, glyco.Count);
            Assert.AreEqual(expectedNormal, normal.Count);

            var parsedSpectralMatches = SpectrumMatchTsvReader.ReadTsv(path, out warnings);
            Assert.AreEqual(expectedGlyco + expectedNormal, parsedSpectralMatches.Count);

            glyco = parsedSpectralMatches.OfType<GlycoPsmFromTsv>().ToList();
            var norm= parsedSpectralMatches.Except(glyco).ToList();
            Assert.AreEqual(expectedGlyco, glyco.Count);
            Assert.AreEqual(expectedNormal, norm.Count);
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
        public static void FailWhenPathIsBadAndContainsNoting()
        {
            string psmFile = @"DatabaseTests\bad4.fasta";
            var lines = File.ReadAllLines(psmFile).Length;

            // Throws the right type of exception
            Assert.Throws<MzLibException>(() =>
            {
                SpectrumMatchTsvReader.ReadPsmTsv(psmFile, out _);
            });

            // Wrapped up the parsing exception in a MzLibException
            bool caught = false;
            List<string> warnings = [];
            List<PsmFromTsv> parsedPsms = [];
            try
            {
                parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFile, out warnings);
            }
            catch (MzLibException e)
            {
                caught = true;
                Assert.That(parsedPsms.Count == 0);
                Assert.That(warnings.Count == lines);
                Assert.That(e is MzLibException);
                Assert.That(e.InnerException, Is.Not.Null);
                Assert.That(e.InnerException, Is.TypeOf<MzLibException>());
                Assert.That(e.InnerException.Message, Does.Contain("type not supported"));
            }

            Assert.That(caught);
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

            // psm with single modification
            PsmFromTsv singleMod = psms[0];
            var modDict = SpectrumMatchFromTsv.ParseModifications(singleMod.FullSequence);
            Assert.That(modDict.Count == 1);
            Assert.That(modDict.ContainsKey(37));
            Assert.That(modDict.Values.First().Contains("Common Fixed:Carbamidomethyl on C"));

            // psm with two modifications
            PsmFromTsv twoMods = psms[15];
            modDict = SpectrumMatchFromTsv.ParseModifications(twoMods.FullSequence);
            Assert.That(modDict.Count == 2);
            Assert.That(modDict.ContainsKey(0) && modDict.ContainsKey(104));
            Assert.AreEqual(modDict[0], "UniProt:N-acetylserine on S");
            Assert.AreEqual(modDict[104], "UniProt:N5-methylglutamine on Q");

            // test sequence with mods at all relevant positions. Mods include cation mod to test bracket selectivity when parsing mods. 
            string allPosSeq = "[mod type1: testmodW on N-term]S[mod type2: testmodX on S]TGTSQ[Common Artifact: Fe[II] on Q]ADDC[mod type3: testmodY on C]-[mod type4: testmodZ on C-Term]";
            // test if the locations of mods are correct
            modDict = SpectrumMatchFromTsv.ParseModifications(allPosSeq);
            Assert.That(modDict.Count == 5);
            Assert.AreEqual(modDict.Keys.Order(), new List<int>{ 0, 1, 6, 10, 11 });
            Assert.AreEqual(modDict[0], "mod type1: testmodW on N-term");
            Assert.AreEqual(modDict[1], "mod type2: testmodX on S");
            Assert.AreEqual(modDict[6], "Common Artifact: Fe[II] on Q");
            Assert.AreEqual(modDict[10], "mod type3: testmodY on C");
            Assert.AreEqual(modDict[11], "mod type4: testmodZ on C-Term");
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
            string cleanedSequence = SpectrumMatchFromTsv.RemoveParentheses(baseSequence);
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

            var testResult = file.First();
            Assert.That(testResult != null);
            Assert.That(!Double.IsNaN(testResult.PEP));
            Assert.That(!Double.IsNaN(testResult.PEP_QValue));
            Assert.That(!Double.IsNaN(testResult.RetentionTime));
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

        [Test]
        [TestCase("BottomUpExample.psmtsv")] // Bottom-Up
        [TestCase("TDGPTMDSearchResults.psmtsv")] // TopDown
        public static void TestProteinGroupInfos(string path)
        {
            string psmFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\SearchResults",
                path);
            List<PsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFilePath, out var warnings);

            SpectrumMatchFromTsv ambiguousResult = parsedPsms.First(psm => psm.Accession.Contains("|"));
            Assert.That(ambiguousResult.Accession.Contains("|"));
            int ambiguityCount = ambiguousResult.Accession.Count(c => c == '|') + 1;
            Assert.AreEqual(ambiguityCount, ambiguousResult.ProteinGroupInfos.Count);

            CollectionAssert.AreEquivalent(ambiguousResult.ProteinGroupInfos.Select(g => g.proteinAccessions).ToList(),
                ambiguousResult.Accession.Split('|').ToList());
            CollectionAssert.AreEquivalent(ambiguousResult.ProteinGroupInfos.Select(g => g.geneName).ToList(),
                ambiguousResult.GeneName.Split('|').ToList());
            Assert.AreEqual(ambiguousResult.ProteinGroupInfos.Select(g => g.organism).First(),
                ambiguousResult.OrganismName); // All these results are specific to one Organism, so no ambiguity at the organismLevel
        }

        [Test]
        [TestCase("BottomUpExample.psmtsv", "04-30-13_CAST_Frac5_4uL.raw", @"D:/Data/This/Is/A/Folder/Tree/04-30-13_CAST_Frac4_6uL.mzML")] // Bottom-Up
        [TestCase("TDGPTMDSearchResults.psmtsv", "TDGPTMDSearchSingleSpectra.raw", @"C:/Data/FXN3_tr1_032017-calib.mzML")] // TopDown
        public static void TestFileNameFilePathDictionary(string path, string filePathA, string filePathB)
        {
            string psmFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\SearchResults",
                path);
            IQuantifiableResultFile file = FileReader.ReadQuantifiableResultFile(psmFilePath);
            List<string> filePaths = new List<string> { filePathA, filePathB };

            var dictionary = file.FileNameToFilePath(filePaths);
            Assert.That(dictionary.Count == 2);
            Assert.That(file.GetQuantifiableResults().Any(p => p.FileName.Equals(dictionary.Keys.First())));
            Assert.That(file.GetQuantifiableResults().Any(p => p.FileName.Equals(dictionary.Keys.Last())));
        }

        [Test]
        public static void TestPsmFromTsvIdentifications()
        {
            string psmFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\SearchResults", "BottomUpExample.psmtsv");

            IQuantifiableResultFile file = FileReader.ReadQuantifiableResultFile(psmFilePath);
            SpectraFileInfo f1 = new SpectraFileInfo("04-30-13_CAST_Frac5_4uL.raw", "A", 0, 0, 0);
            SpectraFileInfo f2 = new SpectraFileInfo(@"D:/Data/This/Is/A/Folder/Tree/04-30-13_CAST_Frac4_6uL.mzML", "A", 0, 0, 0);
            var ids = file.MakeIdentifications(new() { f1, f2});

            CollectionAssert.AreEquivalent(ids.Select(i => i.ModifiedSequence), file.GetQuantifiableResults().Select(i => i.FullSequence));
        }
    }
}
