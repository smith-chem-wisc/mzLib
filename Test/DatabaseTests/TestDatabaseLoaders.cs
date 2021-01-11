﻿// opyright 2016 Stefan Solntsev
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
using Easy.Common.Extensions;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Net;
using System.Text;
using UsefulProteomicsDatabases;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test
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
            Dictionary<string, HashSet<Tuple<int, Modification>>> mods = new Dictionary<string, HashSet<Tuple<int, Modification>>>();
            ProteinDbWriter.WriteXmlDatabase(mods, protein, Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "IsoformTest.xml"));
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
        public static void LoadOriginalMismatchedModifications()
        {
            var protein = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "oblm.xml"), true,
                DecoyType.Reverse, null, false, null, out var unknownModifications);
            Assert.AreEqual(0, protein[0].OneBasedPossibleLocalizedModifications.Count);
            var variant = protein[0].GetVariantProteins()[0];
            protein[0].NonVariantProtein.RestoreUnfilteredModifications();
            Assert.AreEqual(1, protein[0].NonVariantProtein.OneBasedPossibleLocalizedModifications.Count);
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
        public void FilesLoading() //delete mzLib\Test\bin\x64\Debug to update your local unimod list
        {
            Loaders.LoadElements();
            string uniModPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "unimod_tables2.xml");
            string psiModPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml");
            string uniProtPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt");

            // UniModPTMs
            var unimodMods = Loaders.LoadUnimod(uniModPath).ToList();
            Assert.AreEqual(2684, unimodMods.Count); // UniMod PTM list may be updated at some point, causing the unit test to fail

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
            Loaders.LoadElements();
            var sampleModList = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "z.txt"), out var errors).ToList();
            Assert.AreEqual(1, sampleModList.OfType<Modification>().Count());
            Protein protein = new Protein("MCSSSSSSSSSS", "accession", "organism", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>> { { 2, sampleModList.OfType<Modification>().ToList() } }, null, "name", "full_name", false, false, new List<DatabaseReference>(), new List<SequenceVariation>(), disulfideBonds: new List<DisulfideBond>());
            Assert.AreEqual(1, protein.OneBasedPossibleLocalizedModifications[2].OfType<Modification>().Count());
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { protein }, Path.Combine(TestContext.CurrentContext.TestDirectory, "test_modifications_with_proteins.xml"));
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
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinList, proteinDbFilePath);

            var lines = File.ReadAllLines(proteinDbFilePath);
            List<Protein> newProteinList = ProteinDbLoader.LoadProteinXML(proteinDbFilePath, true, DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out var um, -1);
        }

        [Test]
        public void DoNotWriteSameModTwiceAndDoNotWriteInHeaderSinceDifferent()
        {
            Loaders.LoadElements();
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
            var newModResEntries = ProteinDbWriter.WriteXmlDatabase(dictWithThisMod, new List<Protein> { protein }, Path.Combine(TestContext.CurrentContext.TestDirectory, "test_modifications_with_proteins3.xml"));
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

            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { protein }, Path.Combine(TestContext.CurrentContext.TestDirectory, filename));

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

            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { protein }, Path.Combine(TestContext.CurrentContext.TestDirectory, filename));

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

            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { protein }, Path.Combine(TestContext.CurrentContext.TestDirectory, filename));

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
            Assert.AreEqual(filepath + "\\availableUniProtProteomes.txt.gz", downloadedFilePath);

            Dictionary<string, string> uniprotProteoms = ProteinDbRetriever.UniprotProteomesList(downloadedFilePath);

            Assert.IsTrue(uniprotProteoms.Keys.Contains("UP000005640"));
            Assert.AreEqual("Homo sapiens (Human)", uniprotProteoms["UP000005640"]);

            File.Delete(downloadedFilePath);

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