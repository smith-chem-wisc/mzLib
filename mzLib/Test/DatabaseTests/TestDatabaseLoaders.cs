// Copyright 2016 Stefan Solntsev
//
// This file (ChemicalFormula.cs) is part of Chemistry Library.
//
// Chemistry Library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Chemistry Library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with Chemistry Library. If not, see <http://www.gnu.org/licenses/>

using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Omics.Modifications;
using UsefulProteomicsDatabases;
using Stopwatch = System.Diagnostics.Stopwatch;
using NUnit.Framework.Legacy;
using Omics;
using Omics.BioPolymer;

namespace Test.DatabaseTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class TestDatabaseLoaders
    {
        private static Stopwatch Stopwatch { get; set; }

        [SetUp]
        public static void Setuppp()
        {
            Stopwatch = new Stopwatch();
            Stopwatch.Start();
        }

        [TearDown]
        public static void TearDown()
        {
            Console.WriteLine($"Analysis time: {Stopwatch.Elapsed.Hours}h {Stopwatch.Elapsed.Minutes}m {Stopwatch.Elapsed.Seconds}s");
        }

        [Test]
        public static void LoadIsoforms()
        {
            var protein = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "Isoform.fasta"), true, DecoyType.None,
                false, out var errors, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                ProteinDbLoader.UniprotOrganismRegex);
            Assert.AreEqual("Q13409", protein[0].Accession);
            Assert.AreEqual("Q13409-2", protein[1].Accession);
            Assert.AreEqual("Q13409-3", protein[2].Accession);
            Assert.AreEqual("Q13813", protein[3].Accession);
            Assert.AreEqual("Q13813-2", protein[4].Accession);
            Assert.AreEqual("Q13813-3", protein[5].Accession);
            Assert.AreEqual("Q14103", protein[6].Accession);
            Assert.AreEqual("Q14103-2", protein[7].Accession);
            Assert.AreEqual("Q14103-3", protein[8].Accession);
            Assert.AreEqual("Q14103-4", protein[9].Accession);
            ProteinDbWriter.WriteXmlDatabase(protein, Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "IsoformTest.xml"));
            var proteinXml = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "IsoformTest.xml"), true, DecoyType.None, null, false, null, out var unknownMod);
            Assert.AreEqual("Q13409", proteinXml[0].Accession);
            Assert.AreEqual("Q13409-2", proteinXml[1].Accession);
            Assert.AreEqual("Q13409-3", proteinXml[2].Accession);
            Assert.AreEqual("Q13813", proteinXml[3].Accession);
            Assert.AreEqual("Q13813-2", proteinXml[4].Accession);
            Assert.AreEqual("Q13813-3", proteinXml[5].Accession);
            Assert.AreEqual("Q14103", proteinXml[6].Accession);
            Assert.AreEqual("Q14103-2", proteinXml[7].Accession);
            Assert.AreEqual("Q14103-3", proteinXml[8].Accession);
            Assert.AreEqual("Q14103-4", proteinXml[9].Accession);
        }

        [Test]
        [TestCase("cRAP_databaseGPTMD.xml", DecoyType.None)]
        [TestCase("uniprot_aifm1.fasta", DecoyType.None)]
        [TestCase("cRAP_databaseGPTMD.xml", DecoyType.Reverse)]
        [TestCase("uniprot_aifm1.fasta", DecoyType.Reverse)]
        public void LoadingIsReproducible(string fileName, DecoyType decoyType)
        {
            // Load in proteins
            var dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", fileName);
            List<Protein> proteins1 = null;
            List<Protein> proteins2 = null;
            if(fileName.Contains(".xml"))
            {
                proteins1 = ProteinDbLoader.LoadProteinXML(dbPath, true, decoyType, null, false, null, out var unknownModifications);
                proteins2 = ProteinDbLoader.LoadProteinXML(dbPath, true, decoyType, null, false, null, out unknownModifications);
            }
            else if (fileName.Contains(".fasta"))
            {
                proteins1 = ProteinDbLoader.LoadProteinFasta(dbPath, true, decoyType, false, out var unknownModifications);
                proteins2 = ProteinDbLoader.LoadProteinFasta(dbPath, true, decoyType, false, out unknownModifications);
            }
            else
            {
                Assert.Fail("Unknown file type");
            }

            // check are equivalent lists of proteins
            Assert.AreEqual(proteins1.Count, proteins2.Count);
            // Because decoys are sorted before they are returned, the order should be identical
            Assert.AreEqual(proteins1, proteins2);
        }

        [Test]
        [TestCase("proteinEntryLipidMoietyBindingRegion.xml", DecoyType.Reverse)]
        public void LoadingLipidAsMod(string fileName, DecoyType decoyType)
        {
            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml"));
            Dictionary<string, int> formalChargesDictionary = Loaders.GetFormalChargesDictionary(psiModDeserialized);
            List<Modification> UniProtPtms = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt"), formalChargesDictionary).ToList();

            // Load in proteins
            var dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", fileName);
            List<Protein> proteins1 = ProteinDbLoader.LoadProteinXML(dbPath, true, decoyType, UniProtPtms, false, null, out var unknownModifications);
            List<Protein> proteins2 = ProteinDbLoader.LoadProteinXML(dbPath, true, decoyType, UniProtPtms, false, null, out unknownModifications);

            // check are equivalent lists of proteins
            Assert.AreEqual(proteins1.Count, proteins2.Count);
            // Because decoys are sorted before they are returned, the order should be identical
            Assert.AreEqual(proteins1, proteins2);
            var oneBasedPossibleLocalizedModifications = proteins1[0].OneBasedPossibleLocalizedModifications[36];
            var firstMod = oneBasedPossibleLocalizedModifications.First();
            Assert.AreEqual("LIPID", firstMod.FeatureType);
            Assert.AreEqual("Anywhere.", firstMod.LocationRestriction);
            Assert.AreEqual("S-palmitoyl cysteine on C", firstMod.IdWithMotif);
        }

        [Test]
        public static void LoadModWithNl()
        {
            var hah = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "cfInNL.txt"), out var errors).First() as Modification;
            int count = 0;
            foreach (KeyValuePair<MassSpectrometry.DissociationType, List<double>> item in hah.NeutralLosses)
            {
                foreach (double loos in item.Value)
                {
                    count++;
                }
            }
            Assert.AreEqual(2, count);
        }

        [Test]
        public void TestUpdateUnimod()
        {
            var unimodLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "unimod_tables.xml");
            Loaders.UpdateUnimod(unimodLocation);
            Loaders.UpdateUnimod(unimodLocation);
        }

        [Test]
        public void TestUpdatePsiMod()
        {
            var psimodLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "lal.xml");
            Loaders.UpdatePsiMod(psimodLocation);
            Loaders.UpdatePsiMod(psimodLocation);
        }
        [Test]
        public void TestUpdatePsiModObo()
        {
            string testDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "obo");
            Directory.CreateDirectory(testDirectory);
            var psiModOboLocation = Path.Combine(testDirectory, "psi-mod.obo");

            using (StringWriter sw = new())
            {
                Console.SetOut(sw);
                Loaders.UpdatePsiModObo(psiModOboLocation);

                string expected = "psi-mod.obo database did not exist, writing to disk\r\n";
                Assert.AreEqual(expected, sw.ToString());
                sw.Close();
            }

            using (StringWriter sw = new())
            {
                Console.SetOut(sw);
                Loaders.UpdatePsiModObo(psiModOboLocation);

                string expected = "psi-mod.obo database is up to date, doing nothing\r\n";
                Assert.AreEqual(expected, sw.ToString());
                sw.Close();
            }

            //create and empty obo that will be seen as different from the downloaded file and then be updated.
            File.WriteAllText(psiModOboLocation, "");

            using (StringWriter sw = new())
            {
                Console.SetOut(sw);
                Loaders.UpdatePsiModObo(psiModOboLocation);

                string expected = "psi-mod.obo database updated, saving old version as backup\r\n";
                Assert.AreEqual(expected, sw.ToString());
                sw.Close();
            }

            string[] files = Directory.GetFiles(testDirectory);
            foreach (string file in files)
            {
                File.SetAttributes(file, FileAttributes.Normal);
                File.Delete(file);
            }
            Directory.Delete(testDirectory, false);

            // Now you have to restore default output stream
            var standardOutput = new StreamWriter(Console.OpenStandardOutput())
            {
                AutoFlush = true
            };
            Console.SetOut(standardOutput);
        }
        [Test]
        public void TestUpdateElements()
        {
            Assert.IsTrue(PeriodicTable.ValidateAbundances(1e-15));
            Assert.IsTrue(PeriodicTable.ValidateAverageMasses(1e-2));
        }

        [Test]
        public void TestUpdateUniprot()
        {
            var uniprotLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist.txt");
            Loaders.UpdateUniprot(uniprotLocation);
            Loaders.UpdateUniprot(uniprotLocation);
        }

        [Test]
        public void FilesEqualHash()
        {
            var fake = Path.Combine(TestContext.CurrentContext.TestDirectory, "fake.txt");
            using (StreamWriter file = new StreamWriter(fake))
                file.WriteLine("fake");
            Loaders.UpdateUniprot(fake);
            fake = Path.Combine(TestContext.CurrentContext.TestDirectory, "fake1.txt");
            using (StreamWriter file = new StreamWriter(fake))
                file.WriteLine("fake");
            Loaders.UpdateUnimod(fake);
            fake = Path.Combine(TestContext.CurrentContext.TestDirectory, "fake2.txt");
            using (StreamWriter file = new StreamWriter(fake))
                file.WriteLine("fake");
            Loaders.UpdatePsiMod(fake);
            fake = Path.Combine(TestContext.CurrentContext.TestDirectory, "fake3.txt");
            using (StreamWriter file = new StreamWriter(fake))
                file.WriteLine("fake");
        }

        [Test]
        public void TestPsiModLoading()
        {
            string psiModPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "PSI-MOD.obo");

            var psiMods = Loaders.ReadPsiModFile(psiModPath);

            // N6,N6,N6-trimethyllysine
            var trimethylLysine = psiMods.First(b => b.Id.Equals("MOD:00083"));
            Assert.AreEqual("1+",
                trimethylLysine.ValuePairs
                    .First(b => b.Value.Contains("FormalCharge")).GetFormalChargeString());

            // Phosphoserine
            bool resultBool = psiMods.First(b => b.Id.Equals("MOD:00046"))
                .ValuePairs.Any(i => i.Value.Contains("FormalCharge"));
            Assert.IsFalse(resultBool);

            // ensure that there are negative numbers in the formal charges
            Dictionary<string, int> formalChargesDictionary = Loaders.GetFormalChargesDictionary(psiMods);
            bool anyNegativeValue = formalChargesDictionary.Values.Any(i => i < 0);
            Assert.IsTrue(anyNegativeValue);
        }

        [Test]
        public void FilesLoading() //delete mzLib\Test\bin\x64\Debug to update your local unimod list
        {
            string uniModPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "unimod_tables2.xml");
            string psiModPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml");
            string uniProtPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt");

            // UniModPTMs
            var unimodMods = Loaders.LoadUnimod(uniModPath).ToList();
            Assert.Greater(unimodMods.Count, 2700); // UniMod PTM list may be updated at some point, causing the unit test to fail

            List<Modification> myList = unimodMods.Where(m => m.OriginalId.Equals("HexNAc(2)")).ToList();

            Modification testMod = myList.First();
            int neutralLossCount = 0;
            if (testMod.NeutralLosses.Count != 0)
            {
                foreach (KeyValuePair<MassSpectrometry.DissociationType, List<double>> item in testMod.NeutralLosses)
                {
                    foreach (double loss in item.Value)
                    {
                        neutralLossCount++;
                    }
                }
            }

            Assert.AreEqual(2, neutralLossCount);

            // PsiMod PTMs
            var psiModDeserialized = Loaders.LoadPsiMod(psiModPath);

            // N6,N6,N6-trimethyllysine
            var trimethylLysine = psiModDeserialized.Items.OfType<UsefulProteomicsDatabases.Generated.oboTerm>().First(b => b.id.Equals("MOD:00083"));
            Assert.AreEqual("1+", trimethylLysine.xref_analog.First(b => b.dbname.Equals("FormalCharge")).name);

            // Phosphoserine
            Assert.IsFalse(psiModDeserialized.Items.OfType<UsefulProteomicsDatabases.Generated.oboTerm>().First(b => b.id.Equals("MOD:00046")).xref_analog.Any(b => b.dbname.Equals("FormalCharge")));

            Dictionary<string, int> formalChargesDictionary = Loaders.GetFormalChargesDictionary(psiModDeserialized);

            // UniProt PTMs
            var uniprotPtms = Loaders.LoadUniprot(uniProtPath, formalChargesDictionary).ToList();
            Assert.LessOrEqual(300, uniprotPtms.Count()); // UniProt PTM list may be updated at some point, causing the unit test to fail

            // write UniProt and UniMod PTMs to a file
            using (StreamWriter w = new StreamWriter(Path.Combine(TestContext.CurrentContext.TestDirectory, "test.txt")))
            {
                foreach (var nice in uniprotPtms)
                {
                    w.WriteLine(nice.ToString());
                    w.WriteLine("//");
                }
                foreach (var nice in unimodMods)
                {
                    w.WriteLine(nice.ToString());
                    w.WriteLine("//");
                }
            }

            // read in the file and make sure that it has the same number of PTMs
            var sampleModList = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "test.txt"), out var errors).ToList();

            Assert.AreEqual(uniprotPtms.Count + unimodMods.Count, sampleModList.Count());

            List<Modification> myOtherList = new List<Modification>();
            foreach (Modification mod in sampleModList)
            {
                if (mod.IdWithMotif != null && mod.IdWithMotif.Contains("Acetyl"))
                {
                    myOtherList.Add(mod);
                }
            }

            var thisMod = myOtherList.First();
            Assert.IsTrue(thisMod.MonoisotopicMass > 42);
            Assert.IsTrue(thisMod.MonoisotopicMass < 43);

            File.Delete(uniModPath);
            File.Delete(psiModPath);
            File.Delete(uniProtPath);
        }

        /// <summary>
        /// Tests loading an annotated PTM with a longer known motif (>1 character in the motif)
        /// </summary>
        [Test]
        public void SampleLoadModWithLongMotif()
        {
            ModificationMotif.TryGetMotif("msgRgk", out var motif);
            Modification testMod = new Modification(_originalId: "Asymmetric dimethylarginine", _modificationType: "Test", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 100.0);
            List<Modification> allKnownMods = new List<Modification> { testMod };

            Assert.That(testMod.ValidModification);
            Assert.That(testMod.Target.ToString().Equals("msgRgk"));

            Protein protein = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "modified_start.xml"), true, DecoyType.None, allKnownMods, false, new List<string>(), out var unk).First();

            Assert.That(protein.BaseSequence.StartsWith("MSGRGK"));
            Assert.That(protein.OneBasedPossibleLocalizedModifications.Count == 1);
            Assert.That(protein.OneBasedPossibleLocalizedModifications.First().Value.First() == testMod);
        }

        [Test]
        public void SampleModFileLoading()
        {
            PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFile.txt"), out var errors);
        }

        [Test]
        public void SampleModFileLoadingFail1()
        {
            var b = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFileFail1.txt"), out var errors);
            Assert.AreEqual(0, b.Count());
        }

        [Test]
        public void SampleModFileLoadingFail2()
        {
            var b = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFileFail2.txt"), out var errors);
            Assert.AreEqual(0, b.Count());
        }

        [Test]
        public void SampleModFileLoadingFail3()
        {
            Assert.That(() => PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFileFail3.txt"), out var errors).ToList(),
                                            Throws.TypeOf<MzLibException>()
                                            .With.Property("Message")
                                            .EqualTo("Input string for chemical formula was in an incorrect format: $%&$%"));
        }

        [Test]
        public void SampleModFileLoadingFail4()
        {
            Assert.That(() => PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "m.txt"), out var errors).ToList(),
                                            Throws.TypeOf<MzLibException>()
                                            .With.Property("Message")
                                            .EqualTo("0 or 238.229666 is not a valid monoisotopic mass"));
        }

        [Test]
        public void SampleModFileLoadingFail5()
        {
            var b = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFileFail5.txt"), out var errors);
            Assert.AreEqual(0, b.Count());
        }

        [Test]
        public void SampleModFileLoadingFail6()
        {
            var b = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFileFail5.txt"), out var errors);
            Assert.AreEqual(0, b.Count());
        }

        [Test]
        public void CompactFormReading()
        {
            Assert.AreEqual(2, PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFileDouble.txt"), out var errors).Count());
        }

        [Test]
        public void CompactFormReading2()
        {
            Assert.AreEqual(2, PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFileDouble2.txt"), out var errors).Count());
        }

        [Test]
        public void Modification_read_write_into_proteinDb()
        {
            var sampleModList = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "z.txt"), out var errors).ToList();
            Assert.AreEqual(1, sampleModList.OfType<Modification>().Count());
            Protein protein = new Protein("MCSSSSSSSSSS", "accession", "organism", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>> { { 2, sampleModList.OfType<Modification>().ToList() } }, null, "name", "full_name", false, false, new List<DatabaseReference>(), new List<SequenceVariation>(), disulfideBonds: new List<DisulfideBond>());
            Assert.AreEqual(1, protein.OneBasedPossibleLocalizedModifications[2].OfType<Modification>().Count());
            ProteinDbWriter.WriteXmlDatabase(new List<Protein> { protein }, Path.Combine(TestContext.CurrentContext.TestDirectory, "test_modifications_with_proteins.xml"));
            List<Protein> new_proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "test_modifications_with_proteins.xml"),
                true, DecoyType.None, new List<Modification>(), false, new List<string>(), out Dictionary<string, Modification> um);
            Assert.AreEqual(1, new_proteins.Count);
            Assert.AreEqual(1, new_proteins[0].OneBasedPossibleLocalizedModifications.Count);
            Assert.AreEqual(1, new_proteins[0].OneBasedPossibleLocalizedModifications.SelectMany(kv => kv.Value).Count());
            Assert.AreEqual("Type", new_proteins[0].OneBasedPossibleLocalizedModifications.SelectMany(kv => kv.Value).OfType<Modification>().First().ModificationType);
            Assert.AreEqual("Palmitoylation on C", new_proteins[0].OneBasedPossibleLocalizedModifications[2][0].IdWithMotif);
            Assert.AreEqual(1, new_proteins[0].OneBasedPossibleLocalizedModifications[2].OfType<Modification>().Count());

            // Check that Modifications were saved after last load
            Assert.AreEqual(1, ProteinDbLoader.GetPtmListFromProteinXml(Path.Combine(TestContext.CurrentContext.TestDirectory, @"test_modifications_with_proteins.xml")).Count);
            Assert.True(ProteinDbLoader.GetPtmListFromProteinXml(Path.Combine(TestContext.CurrentContext.TestDirectory, @"test_modifications_with_proteins.xml"))[0] == new_proteins[0].OneBasedPossibleLocalizedModifications.SelectMany(kv => kv.Value).First());

            //But that we can still read modifications from other protein XMLs that exist
            Assert.AreEqual(0, ProteinDbLoader.GetPtmListFromProteinXml(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "xml.xml")).Count);

            // Check that Modifications were saved after last load
            var b = ProteinDbLoader.GetPtmListFromProteinXml(Path.Combine(TestContext.CurrentContext.TestDirectory, @"test_modifications_with_proteins.xml"));
            Assert.AreEqual(1, b.Count);

            var c = ProteinDbLoader.GetPtmListFromProteinXml(Path.Combine(TestContext.CurrentContext.TestDirectory, @"test_modifications_with_proteins.xml"))[0];
            var d = new_proteins[0].OneBasedPossibleLocalizedModifications.SelectMany(kv => kv.Value).First();

            Assert.IsTrue(c.Equals(d));

            //But that we can still read modifications from other protein XMLs that exist
            Assert.AreEqual(0, ProteinDbLoader.GetPtmListFromProteinXml(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "xml.xml")).Count);
        }

        [Test]
        public void MultiMod_ProteinDbWriter()
        {
            var sampleModList = PtmListLoader
                .ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "z.txt"),
                    out var errors).ToList();
            var currentMod = sampleModList.First();
            // create slightly different modifications
            var newMod = new Modification(_originalId: "1" + currentMod.OriginalId, _target: currentMod.Target,
                _modificationType: currentMod.ModificationType,
                _accession: currentMod.Accession, _locationRestriction: currentMod.LocationRestriction,
                _featureType: currentMod.FeatureType,
                _chemicalFormula: currentMod.ChemicalFormula);
            var newMod2 = new Modification(_originalId: "2" + currentMod.OriginalId, _target: currentMod.Target,
                _modificationType: currentMod.ModificationType,
                _accession: currentMod.Accession, _locationRestriction: currentMod.LocationRestriction,
                _featureType: currentMod.FeatureType,
                _chemicalFormula: currentMod.ChemicalFormula);
            var newMod3 = new Modification(_originalId: "3" + currentMod.OriginalId, _target: currentMod.Target,
                _modificationType: currentMod.ModificationType,
                _accession: currentMod.Accession, _locationRestriction: currentMod.LocationRestriction,
                _featureType: currentMod.FeatureType,
                _chemicalFormula: currentMod.ChemicalFormula);
            var newMod4 = new Modification(_originalId: "4" + currentMod.OriginalId, _target: currentMod.Target,
                _modificationType: currentMod.ModificationType,
                _accession: currentMod.Accession, _locationRestriction: currentMod.LocationRestriction,
                _featureType: currentMod.FeatureType,
                _chemicalFormula: currentMod.ChemicalFormula);
            var newMod5 = new Modification(_originalId: "5" + currentMod.OriginalId, _target: currentMod.Target,
                _modificationType: currentMod.ModificationType,
                _accession: currentMod.Accession, _locationRestriction: currentMod.LocationRestriction,
                _featureType: currentMod.FeatureType,
                _chemicalFormula: currentMod.ChemicalFormula);
            sampleModList.AddRange(new List<Modification>() { newMod, newMod2, newMod3, newMod4, newMod5 });
            Assert.AreEqual(6, sampleModList.OfType<Modification>().Count());
            // Create a protein with all possible modifications
            Protein protein = new Protein(
                "MCMCMCSSSSSSSS",
                "accession",
                "organism",
                new List<Tuple<string, string>>(),
                new Dictionary<int, List<Modification>>
                {
                    { 2, sampleModList.OfType<Modification>().ToList() },
                    { 4, sampleModList.OfType<Modification>().ToList() },
                    { 6, sampleModList.OfType<Modification>().ToList() },
                },
                null,
                "name",
                "full_name",
                false,
                false,
                new List<DatabaseReference>(),
                new List<SequenceVariation>(),
                disulfideBonds: new List<DisulfideBond>());

            Assert.AreEqual(6, protein.OneBasedPossibleLocalizedModifications[2].OfType<Modification>().Count());
            Assert.AreEqual(18, protein.OneBasedPossibleLocalizedModifications.SelectMany(kvp => kvp.Value).Count());
            ProteinDbWriter.WriteXmlDatabase(new List<Protein> { protein },
                Path.Combine(TestContext.CurrentContext.TestDirectory, "test_modifications_with_proteins.xml"));
            List<Protein> newProteins = ProteinDbLoader.LoadProteinXML(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "test_modifications_with_proteins.xml"),
                true, DecoyType.None, new List<Modification>(), false, new List<string>(),
                out Dictionary<string, Modification> um);

            // Create a second protein with the same modifications, but listed in a different order.
            sampleModList.Reverse();
            Protein modShuffledProtein = new Protein(
                "MCMCMCSSSSSSSS",
                "accession",
                "organism",
                new List<Tuple<string, string>>(),
                new Dictionary<int, List<Modification>>
                {
                    { 2, sampleModList.OfType<Modification>().ToList() },
                    { 4, sampleModList.OfType<Modification>().ToList() },
                    { 6, sampleModList.OfType<Modification>().ToList() },
                },
                null,
                "name",
                "full_name",
                false,
                false,
                new List<DatabaseReference>(),
                new List<SequenceVariation>(),
                disulfideBonds: new List<DisulfideBond>());
            string shuffledProteinFileName = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "test_shuffled_modifications_with_proteins.xml");
            ProteinDbWriter.WriteXmlDatabase(new List<Protein> { modShuffledProtein }, shuffledProteinFileName);
            List<Protein> newShuffledProteins = ProteinDbLoader.LoadProteinXML(shuffledProteinFileName,
                true, DecoyType.None, new List<Modification>(), false, new List<string>(), out um);

            // We've read in proteins from both databases. Assert that they are equal
            Assert.AreEqual(newShuffledProteins.First().Accession, newProteins.First().Accession);
            Assert.AreEqual(newShuffledProteins.First(), newProteins.First());

            // Now, ensure that the modification dictionaries for each are equivalent (contain the same mods) and equal (contain the same mods in the same order)
            for(int i = 1; i<4; i++)
            {
                int oneBasedResidue = i * 2;

                Assert.That(newShuffledProteins.First().OneBasedPossibleLocalizedModifications[oneBasedResidue],
                    Is.EquivalentTo(newProteins.First().OneBasedPossibleLocalizedModifications[oneBasedResidue]));

                Assert.That(newShuffledProteins.First().OneBasedPossibleLocalizedModifications[oneBasedResidue],
                    Is.EqualTo(newProteins.First().OneBasedPossibleLocalizedModifications[oneBasedResidue]));
            }
        }

        [Test]
        public static void Test_MetaMorpheusStyleProteinDatabaseWriteAndREad()
        {
            string proteinDbFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestProteinSplitAcrossFiles.xml");

            ModificationMotif.TryGetMotif("D", out ModificationMotif motif);
            Modification mod = new Modification(_originalId: "mod1", _modificationType: "mt", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 10);

            IDictionary<int, List<Modification>> oneBasedModification = new Dictionary<int, List<Modification>>
            {
                { 3, new List<Modification>{ mod } }
            };

            Protein prot1 = new Protein("MEDEEK", "prot1", oneBasedModifications: oneBasedModification);
            List<Protein> proteinList = new List<Protein> { prot1 };
            ProteinDbWriter.WriteXmlDatabase( proteinList, proteinDbFilePath);

            var lines = File.ReadAllLines(proteinDbFilePath);
            List<Protein> newProteinList = ProteinDbLoader.LoadProteinXML(proteinDbFilePath, true, DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out var um, -1);
        }

        [Test]
        public void DoNotWriteSameModTwiceAndDoNotWriteInHeaderSinceDifferent()
        {
            var sampleModList = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "z.txt"), out var errors).ToList();
            Protein protein = new Protein("MCSSSSSSSSSS", "accession", "organism", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>> { { 2, sampleModList.OfType<Modification>().ToList() } }, null, "name", "full_name", false, false, new List<DatabaseReference>(), new List<SequenceVariation>(), disulfideBonds: new List<DisulfideBond>());
            Assert.AreEqual(1, protein.OneBasedPossibleLocalizedModifications[2].OfType<Modification>().Count());

            Dictionary<string, HashSet<Tuple<int, Modification>>> dictWithThisMod = new Dictionary<string, HashSet<Tuple<int, Modification>>>();

            HashSet<Tuple<int, Modification>> value = new HashSet<Tuple<int, Modification>>();

            var modReadFromFile = sampleModList.First() as Modification;
            ModificationMotif.TryGetMotif("C", out ModificationMotif motif);
            Modification newMod = new Modification(_originalId: "Palmitoylation of C", _modificationType: "Type", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: modReadFromFile.ChemicalFormula, _monoisotopicMass: modReadFromFile.MonoisotopicMass, _featureType: "MOD_RES", _fileOrigin: "E:\\GitClones\\mzLib\\Test\\bin\\x64\\Debug\\DatabaseTests\\z.txt");

            Assert.IsTrue(newMod.Equals(sampleModList.First()));

            Assert.AreEqual(newMod, sampleModList.First());
            Assert.AreEqual(sampleModList.First(), newMod);

            value.Add(new Tuple<int, Modification>(2, newMod));

            dictWithThisMod.Add("accession", value);
            var newModResEntries = ProteinDbWriter.WriteXmlDatabase(new List<Protein> { protein }, Path.Combine(TestContext.CurrentContext.TestDirectory, "test_modifications_with_proteins3.xml"));
            Assert.AreEqual(0, newModResEntries.Count);
            List<Protein> new_proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "test_modifications_with_proteins3.xml"),
                true, DecoyType.None, new List<Modification>(), false, new List<string>(), out Dictionary<string, Modification> um);
            Assert.AreEqual(1, new_proteins.Count);
            Assert.AreEqual(1, new_proteins[0].OneBasedPossibleLocalizedModifications.Count);
            Assert.AreEqual(1, new_proteins[0].OneBasedPossibleLocalizedModifications.SelectMany(kv => kv.Value).Count());
        }

        [Test]
        public void TestWritePtmWithNeutralLoss()
        {
            string filename = "test_neutral_loss_mod.xml";
            Dictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>>();

            ModificationMotif.TryGetMotif("T", out var motif);
            Modification m = new Modification(_originalId: "Phospho", _modificationType: "Test", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 80.0, _neutralLosses: new Dictionary<DissociationType, List<double>> { { DissociationType.HCD, new List<double> { 80.0, 0 } }, { DissociationType.ETD, new List<double> { 70.0, 0 } } });
            Assert.That(m.ValidModification);

            mods.Add(4, new List<Modification> { m });

            Protein protein = new Protein("PEPTIDE", "accession", oneBasedModifications: mods);
            Assert.That(protein.OneBasedPossibleLocalizedModifications.Count == 1);
            Assert.That(protein.OneBasedPossibleLocalizedModifications.First().Value.First().NeutralLosses.First().Value.Count == 2);

            ProteinDbWriter.WriteXmlDatabase(new List<Protein> { protein }, Path.Combine(TestContext.CurrentContext.TestDirectory, filename));

            // with passed-in mods
            List<Protein> new_proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, filename), true, DecoyType.None, new List<Modification> { m }, false, new List<string>(), out Dictionary<string, Modification> um);
            Assert.That(new_proteins.First().OneBasedPossibleLocalizedModifications.First().Value.First().NeutralLosses.First().Value.Count == 2);

            // should be able to read mod from top of database...
            new_proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, filename), true, DecoyType.None, new List<Modification>(), false, new List<string>(), out um);
            Assert.That(new_proteins.First().OneBasedPossibleLocalizedModifications.First().Value.First().NeutralLosses.First().Value.Count == 2);
        }

        [Test]
        public void TestWritePtmWithNeutralLoss_AsBioPolymer()
        {
            string filename = "test_neutral_loss_mod2.xml";
            Dictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>>();

            ModificationMotif.TryGetMotif("T", out var motif);
            Modification m = new Modification(_originalId: "Phospho", _modificationType: "Test", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 80.0, _neutralLosses: new Dictionary<DissociationType, List<double>> { { DissociationType.HCD, new List<double> { 80.0, 0 } }, { DissociationType.ETD, new List<double> { 70.0, 0 } } });
            Assert.That(m.ValidModification);

            mods.Add(4, new List<Modification> { m });

            IBioPolymer protein = new Protein("PEPTIDE", "accession", oneBasedModifications: mods);
            Assert.That(protein.OneBasedPossibleLocalizedModifications.Count == 1);
            Assert.That(protein.OneBasedPossibleLocalizedModifications.First().Value.First().NeutralLosses.First().Value.Count == 2);

            ProteinDbWriter.WriteXmlDatabase(new List<IBioPolymer> { protein }, Path.Combine(TestContext.CurrentContext.TestDirectory, filename));

            // with passed-in mods
            List<Protein> new_proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, filename), true, DecoyType.None, new List<Modification> { m }, false, new List<string>(), out Dictionary<string, Modification> um);
            Assert.That(new_proteins.First().OneBasedPossibleLocalizedModifications.First().Value.First().NeutralLosses.First().Value.Count == 2);

            // should be able to read mod from top of database...
            new_proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, filename), true, DecoyType.None, new List<Modification>(), false, new List<string>(), out um);
            Assert.That(new_proteins.First().OneBasedPossibleLocalizedModifications.First().Value.First().NeutralLosses.First().Value.Count == 2);
        }

        [Test]
        public void TestWritePtmWithDiagnosticIons()
        {
            string filename = "test_diagnostic_ion_mod.xml";
            Dictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>>();

            ModificationMotif.TryGetMotif("T", out var motif);
            Modification m = new Modification(_originalId: "Phospho", _modificationType: "Test", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 80.0, _diagnosticIons: new Dictionary<DissociationType, List<double>> { { DissociationType.HCD, new List<double> { 80.0, 0 } }, { DissociationType.ETD, new List<double> { 70.0, 0 } } });
            Assert.That(m.ValidModification);

            mods.Add(4, new List<Modification> { m });

            Protein protein = new Protein("PEPTIDE", "accession", oneBasedModifications: mods);
            Assert.That(protein.OneBasedPossibleLocalizedModifications.Count == 1);
            Assert.That(protein.OneBasedPossibleLocalizedModifications.First().Value.First().DiagnosticIons.First().Value.Count == 2);

            ProteinDbWriter.WriteXmlDatabase(new List<Protein> { protein }, Path.Combine(TestContext.CurrentContext.TestDirectory, filename));

            // with passed-in mods
            List<Protein> new_proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, filename), true, DecoyType.None, new List<Modification> { m }, false, new List<string>(), out Dictionary<string, Modification> um);
            Assert.That(new_proteins.First().OneBasedPossibleLocalizedModifications.First().Value.First().DiagnosticIons.First().Value.Count == 2);

            // should be able to read mod from top of database...
            new_proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, filename), true, DecoyType.None, new List<Modification>(), false, new List<string>(), out um);
            Assert.That(new_proteins.First().OneBasedPossibleLocalizedModifications.First().Value.First().DiagnosticIons.First().Value.Count == 2);
        }

        [Test]
        public void TestWritePtmWithNeutralLossAndDiagnosticIons()
        {
            string filename = "test_neutral_loss_diagnostic_ion_mod.xml";
            Dictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>>();

            ModificationMotif.TryGetMotif("T", out var motif);
            Modification m = new Modification(_originalId: "Phospho", _modificationType: "Test", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 80.0, _neutralLosses: new Dictionary<DissociationType, List<double>> { { DissociationType.HCD, new List<double> { 80.0, 0 } }, { DissociationType.ETD, new List<double> { 70.0, 0 } } }, _diagnosticIons: new Dictionary<DissociationType, List<double>> { { DissociationType.CID, new List<double> { 60.0, 0 } }, { DissociationType.EThcD, new List<double> { 40.0, 0 } } });
            Assert.That(m.ValidModification);

            mods.Add(4, new List<Modification> { m });

            Protein protein = new Protein("PEPTIDE", "accession", oneBasedModifications: mods);
            Assert.That(protein.OneBasedPossibleLocalizedModifications.Count == 1);
            Assert.That(protein.OneBasedPossibleLocalizedModifications.First().Value.First().NeutralLosses.First().Value.Count == 2);

            ProteinDbWriter.WriteXmlDatabase(new List<Protein> { protein }, Path.Combine(TestContext.CurrentContext.TestDirectory, filename));

            // with passed-in mods
            List<Protein> new_proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, filename), true, DecoyType.None, new List<Modification> { m }, false, new List<string>(), out Dictionary<string, Modification> um);
            Assert.That(new_proteins.First().OneBasedPossibleLocalizedModifications.First().Value.First().NeutralLosses.First().Value.Count == 2);
            Assert.That(new_proteins.First().OneBasedPossibleLocalizedModifications.First().Value.First().DiagnosticIons.First().Value.Count == 2);

            // should be able to read mod from top of database...
            new_proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, filename), true, DecoyType.None, new List<Modification>(), false, new List<string>(), out um);
            Assert.That(new_proteins.First().OneBasedPossibleLocalizedModifications.First().Value.First().NeutralLosses.First().Value.Count == 2);
            Assert.That(new_proteins.First().OneBasedPossibleLocalizedModifications.First().Value.First().DiagnosticIons.First().Value.Count == 2);
        }

        [Test]
        public static void IsoformReadTest()
        {
            string filepath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DatabaseTests\isoformTest.fasta");
            var proteinList = ProteinDbLoader.LoadProteinFasta(filepath, true, DecoyType.None, false, out var dbErrors, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex,
                ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex, ProteinDbLoader.UniprotOrganismRegex);

            Assert.AreNotEqual(proteinList[2].Accession, proteinList[4].Accession);
        }

        [Test]
        public static void TestRetrieveUniProtProteome()
        {
            string filepath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DatabaseTests");

            //UP000008595 is Uukuniemi virus (strain S23) (Uuk) which only has 4 proteins
            string returnedFilePath = ProteinDbRetriever.RetrieveProteome("UP000008595", filepath, ProteinDbRetriever.ProteomeFormat.fasta, ProteinDbRetriever.Reviewed.yes, ProteinDbRetriever.Compress.yes, ProteinDbRetriever.IncludeIsoforms.yes);

            filepath += "\\UP000008595_reviewed_isoform.fasta.gz";

            Assert.AreEqual(filepath, returnedFilePath);
            Assert.IsTrue(File.Exists(filepath));

            var proteinList = ProteinDbLoader.LoadProteinFasta(filepath, true, DecoyType.None, false, out var dbErrors, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex,
                ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex, ProteinDbLoader.UniprotOrganismRegex);

            Assert.AreEqual(4, proteinList.Count);
            Assert.IsTrue(proteinList.Select(p => p.Accession).ToList().Contains("P33453"));
            Assert.AreEqual("MLLAICSRTIRQQGLNCPPAVTFTSSHMRPPIPSFLLWTEGSDVLMDFDLDTIPAGSVTGSSIGPKFKIKTQAASSFVHDFTFAHWCDASDMPLRDHFPLVNDTFDHWTPDFISQRLDGSKVVVEFTTNRSDQEQSLISAFNTKVGKYEVALHNRSTTSSILFGVVV" +
                "VSETTVVTNLNLNQQEVDELCFRFLVARAVHLEMTTKMIIPEYDDEDEDKRSREVKAAFHSVQPDWNVTEANFAPFSRRMFSNFAQMEPDKEYLAHIILDSLKQAQADLDGNHYLNESLTEQARLDRNREESLNMVKDFERDFNNAAQRSAWSHKSTVPFPGVIPKVSGDTTSLSRLVEL" +
                "PVITGGSDATIRAWRSAYGSVSNGTVERCDEDVERERRAALCSLTVEELEESKALRMKYHRCKIDNGMMDKLDLAMQGVEAKEFKNHPSIIKKRSKSKKTFPLTADTRDIDLFLHHDDLMFNNEHSQTPPAAMIEAVKAGADAQSLHGLDKSANPWYASALWFLGLPIGLWLFMCTCIGVEL" +
                "SISLKQHCGRQKFIIKKLRFFDIFLLIKPTNSGSHVFYSIAFPESAILGKLHRSQCFKGLQFEDGWFWTEFSSFKMSKLTNVVKCLSTGFNLFWFWRDYYEVPFWAGNEKDFQTGKQRANKMFKFCLLMLLEDKARTEEIATLSRYVMMEGFVSPPCIPKPQKMIEKLPNLARTKFQVWLISR" +
                "MLQTIIRVSDYPFKITAGHKSANWTGMFNWVTGEPIESTQKLISLFYLGYLKNKEESPERNASIGMYKKILEYEDKHPGRYTYLGLGDPPSDDTRFHEYSISLLKHLCIHAEHDLRRNWGESFKAMISRDIVDAIASLDLERLATLKASSNFNEEWYQKRGDGKTYHRSKVLEKVSKYVKKSSSH" +
                "VHHIMEECLRKVESQGCMHVCLFKKPQHGGLREIYVLGFEERVVQLVIETIARQICKRFKSETLTNPKQKLAIPETHGLRAVKTCGIHHETVATSDDAAKWNQCHHVTKFALMLCHFTDPLFHGFIIRGCSMFMKKRIMIDQSLIDIIDSHTTLETSDAYLQKIHRGYHGSLDDQPRWISRGGAFVQ" +
                "TETGMMQGILHYTSSLLHTLLQEWLRTFSQRFIRTRVSVDQRPDVLVDVLQSSDDSGMMISFPSTDKGATGKYRYLSALIFKYKKVIGKYLGIYSSVKSTNNTLHLLEFNSEFFFHINHNRPLLRWITACDTISEQESLASRQEEMYNNLTSVLEGGGSFSLVSFCQFGQLLLHYTLLGMTVSPLFLEY" +
                "IKLVSEIKDPSLGYFLMDHPFGSGLSGFKYNVWVAVQNSILGSRYRSLLEAIQNSDSAAPKKTLDTTTSGTFVQSTIIRFGDRKKWQRLVDRLNLPEDWLDVIDKNPEIVYRRPRDGFEVSLRIAEKVHSPGVSNSLSKGNCIIRVISSSVYILSRSILSDGLAWLYDEEEEVKRPLLYKVMNQPELDLHSRLTPA" +
                "QLSTLFPMMAEFEKLQTHLRSYMKIEGEFISKKKVITQTRVNILETERFLRARPEDLIADKWFGFTRTRMTPRTFKEEWENLTSVFPWLTGNPSETLELSPFQHHVQLRNFFSRLDLKGRDIRIIGAPIKKSSGVSNVSTAIRDNFFPRFVLTHIPDEAAMERIEAAGILKHALFLTVTGPYTDQSKLDMCRDF" +
                "ITSSEPITLKPNHGKTRTNVLSLFQDYFSKRGPDIIFNRIQMANCGVIGGFTSPQKPKEVDGKIVYTGDGVWRGIVDGFQIQLVITYMPKQKSNELKSITVNSDRCISALSSFCQSWCKEMGVFNTEDFSKTQRFSKASFFMHKFKISGSKQTLGAPIFIVSEKIFRPICWDPSKLEFRVRGNTLNLTY" +
                "KEVNPGAGQRMFNILSYTVKDTDVSDENAFKLMSLSPRHKFHGREPSTSWICMRALPISTIDKLLERILNRERISGSIDNERLAECFKNVMESTLRRKGVFLSEFSRATQKMLDGLSRDMLDFFAEAGLNDDLLLEEEPWLSGLDTFMLDDEAYLEEYNLGPFGVFSVEQEMNTKYYHHLLLD" +
                "SLVEDVIQKLSLDGLRKLFQEEEAPLEYKKEVIRLLNILQRDASQIKWKSRDLLSENMGLDVDDDMFG", proteinList.Where(p => p.Accession == "P33453").FirstOrDefault().BaseSequence);

            File.Delete(filepath);

            //fasta; unreviewed; non-compressed; no isoforms
            filepath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DatabaseTests");
            returnedFilePath = ProteinDbRetriever.RetrieveProteome("UP000008595", filepath, ProteinDbRetriever.ProteomeFormat.fasta, ProteinDbRetriever.Reviewed.no, ProteinDbRetriever.Compress.no, ProteinDbRetriever.IncludeIsoforms.no);
            filepath += "\\UP000008595_unreviewed.fasta";
            Assert.AreEqual(filepath, returnedFilePath);
            Assert.IsTrue(File.Exists(filepath));
            File.Delete(filepath);

            //xml; reviewed; compresseded; no isoforms
            filepath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DatabaseTests");
            returnedFilePath = ProteinDbRetriever.RetrieveProteome("UP000008595", filepath, ProteinDbRetriever.ProteomeFormat.xml, ProteinDbRetriever.Reviewed.yes, ProteinDbRetriever.Compress.yes, ProteinDbRetriever.IncludeIsoforms.no);
            filepath += "\\UP000008595_reviewed.xml.gz";
            Assert.AreEqual(filepath, returnedFilePath);
            Assert.IsTrue(File.Exists(filepath));
            File.Delete(filepath);

            //xml; unreviewed; non-compresseded; no isoforms
            filepath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DatabaseTests");
            returnedFilePath = ProteinDbRetriever.RetrieveProteome("UP000008595", filepath, ProteinDbRetriever.ProteomeFormat.xml, ProteinDbRetriever.Reviewed.no, ProteinDbRetriever.Compress.no, ProteinDbRetriever.IncludeIsoforms.no);
            filepath += "\\UP000008595_unreviewed.xml";
            Assert.AreEqual(filepath, returnedFilePath);
            Assert.IsTrue(File.Exists(filepath));
            File.Delete(filepath);

            //junk null return
            filepath = "pathDoesNotExists";
            returnedFilePath = ProteinDbRetriever.RetrieveProteome("UP000008595", filepath, ProteinDbRetriever.ProteomeFormat.xml, ProteinDbRetriever.Reviewed.no, ProteinDbRetriever.Compress.no, ProteinDbRetriever.IncludeIsoforms.no);
            filepath += "\\UP000008595_unreviewed.xml";
            Assert.IsNull(returnedFilePath);

            //we don't support filetypes other than fasta or xml currently
            //requesting gff or other file formats will return null for now.
            filepath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DatabaseTests");
            returnedFilePath = ProteinDbRetriever.RetrieveProteome("UP000008595", filepath, ProteinDbRetriever.ProteomeFormat.gff, ProteinDbRetriever.Reviewed.no, ProteinDbRetriever.Compress.no, ProteinDbRetriever.IncludeIsoforms.no);
            filepath += "\\UP000008595_unreviewed.xml";
            Assert.IsNull(returnedFilePath);
        }

        [Test]
        public static void TestDownloadAvailableUniProtProteomes()
        {
            string filepath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DatabaseTests");
            string downloadedFilePath = ProteinDbRetriever.DownloadAvailableUniProtProteomes(filepath);
            Assert.AreEqual(Path.Combine(filepath, "availableUniProtProteomes.txt.gz"), downloadedFilePath);

            Dictionary<string, string> uniprotProteoms = ProteinDbRetriever.UniprotProteomesList(downloadedFilePath);

            Assert.IsTrue(uniprotProteoms.Keys.Contains("UP000005640"));
            Assert.AreEqual("Homo sapiens (Human)", uniprotProteoms["UP000005640"]);

            File.Delete(downloadedFilePath);
            uniprotProteoms.Clear();

            //test different types of compression
            uniprotProteoms = ProteinDbRetriever.UniprotProteomesList(Path.Combine(TestContext.CurrentContext.TestDirectory, @"DatabaseTests", @"b_availableUniProtProteomes.txt"));
            Assert.IsTrue(uniprotProteoms.Keys.Contains("UP000005640"));
            Assert.AreEqual("Homo sapiens (Human)", uniprotProteoms["UP000005640"]);
            uniprotProteoms.Clear();

            uniprotProteoms = ProteinDbRetriever.UniprotProteomesList(Path.Combine(TestContext.CurrentContext.TestDirectory, @"DatabaseTests", @"b_availableUniProtProteomes.txt.gz"));
            Assert.IsTrue(uniprotProteoms.Keys.Contains("UP000005640"));
            Assert.AreEqual("Homo sapiens (Human)", uniprotProteoms["UP000005640"]);
            uniprotProteoms.Clear();

            uniprotProteoms = ProteinDbRetriever.UniprotProteomesList(Path.Combine(TestContext.CurrentContext.TestDirectory, @"DatabaseTests", @"b_availableUniProtProteomes.zip"));
            Assert.IsTrue(uniprotProteoms.Keys.Contains("UP000005640"));
            Assert.AreEqual("Homo sapiens (Human)", uniprotProteoms["UP000005640"]);
            uniprotProteoms.Clear();

            //return null for bad filepath
            filepath = "bubba";
            downloadedFilePath = ProteinDbRetriever.DownloadAvailableUniProtProteomes(filepath);
            Assert.IsNull(downloadedFilePath);

            //bad file path returns null
            uniprotProteoms = ProteinDbRetriever.UniprotProteomesList("badFilePath");
            Assert.IsNull(uniprotProteoms);


            //wrong file extension returns null
            uniprotProteoms = ProteinDbRetriever.UniprotProteomesList(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"bad.fasta"));
            Assert.IsNull(uniprotProteoms);
        }

        [Test]
        public static void TestDownloadListOfColumnsAvailableAtUniProt()
        {
            var uniProtColumnDictionary = ProteinDbRetriever.UniprotColumnsList();
            Assert.IsTrue(uniProtColumnDictionary.Keys.Contains("Entry"));
            Assert.AreEqual("id", uniProtColumnDictionary["Entry"]);
        }

        [Test]
        public static void TestHyphenAccession()
        {
            string fastaFile = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "TestHypenAccession.fasta");

            string header = @">TTLL10-203.TTLL10-202|x|x|x|x|x|TTLL10|x";
            List<string> output = new List<string> { header };
            output.Add("PEPTIDE");
            File.WriteAllLines(fastaFile, output);

            var proteins = ProteinDbLoader.LoadProteinFasta(fastaFile, true, DecoyType.None, false, out var errors);

            Assert.That(proteins.First().Accession == "TTLL10-203.TTLL10-202");

            File.Delete(fastaFile);
        }

        [Test]
        public static void TestDifferentHeaderStyles()
        {
            // uniprot database
            string fastaFile = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "uniprot_aifm1.fasta");
            var proteins = ProteinDbLoader.LoadProteinFasta(fastaFile, true, DecoyType.Reverse, false, out var errors);
            Assert.That(proteins.Count == 2);

            var targetProtein = proteins.First(p => !p.IsDecoy);
            Assert.That(targetProtein.Accession == "Q9Z0X1");
            Assert.That(targetProtein.GeneNames.Count() == 1);
            Assert.That(targetProtein.GeneNames.First().Item2 == "Aifm1");
            Assert.That(targetProtein.FullName == "Apoptosis-inducing factor 1, mitochondrial");
            Assert.That(targetProtein.Name == "AIFM1_MOUSE");
            Assert.That(targetProtein.Organism == "Mus musculus");

            // gencode database
            fastaFile = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "gencode_mmp20.fa");
            proteins = ProteinDbLoader.LoadProteinFasta(fastaFile, true, DecoyType.Reverse, false, out errors);
            Assert.That(proteins.Count == 2);

            targetProtein = proteins.First(p => !p.IsDecoy);

            Assert.That(targetProtein.Accession == "ENSMUSP00000034487.2");
            Assert.That(targetProtein.GeneNames.Count() == 1);
            Assert.That(targetProtein.GeneNames.First().Item2 == "Mmp20");
            Assert.That(targetProtein.FullName == "Mmp20-201");

            // ensembl database
            fastaFile = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "ensembl_prrc2a.fa");
            proteins = ProteinDbLoader.LoadProteinFasta(fastaFile, true, DecoyType.Reverse, false, out errors);
            Assert.That(proteins.Count == 2);

            targetProtein = proteins.First(p => !p.IsDecoy);
            Assert.That(targetProtein.Accession == "ENSP00000372947.2");
            Assert.That(targetProtein.GeneNames.Count() == 1);
            Assert.That(targetProtein.GeneNames.First().Item2 == "ENSG00000206427.11");
        }
    }
}