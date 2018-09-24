// opyright 2016 Stefan Solntsev
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
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public class DatabaseLoaderTests
    {
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
        public void TestUpdateElements()
        {
            var elementLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "lal.dat");
            Loaders.UpdateElements(elementLocation);
            Loaders.UpdateElements(elementLocation);
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
            Loaders.UpdateElements(fake);
        }

        [Test]
        public void FilesLoading()
        {
            Loaders.LoadElements(Path.Combine(TestContext.CurrentContext.TestDirectory, "elements2.dat"));

            var unimodMods = Loaders.LoadUnimod(Path.Combine(TestContext.CurrentContext.TestDirectory, "unimod_tables2.xml")).ToList();
            Assert.AreEqual(2639, unimodMods.Count); // UniMod PTM list may be updated at some point, causing the unit test to fail

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
            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml"));
            
            // N6,N6,N6-trimethyllysine
            var trimethylLysine = psiModDeserialized.Items.OfType<UsefulProteomicsDatabases.Generated.oboTerm>().First(b => b.id.Equals("MOD:00083"));
            Assert.AreEqual("1+", trimethylLysine.xref_analog.First(b => b.dbname.Equals("FormalCharge")).name);

            // Phosphoserine
            Assert.IsFalse(psiModDeserialized.Items.OfType<UsefulProteomicsDatabases.Generated.oboTerm>().First(b => b.id.Equals("MOD:00046")).xref_analog.Any(b => b.dbname.Equals("FormalCharge")));

            Dictionary<string, int> formalChargesDictionary = Loaders.GetFormalChargesDictionary(psiModDeserialized);

            var uniprotPtms = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt"), formalChargesDictionary).ToList();
            Assert.AreEqual(334, uniprotPtms.Count()); // UniProt PTM list may be updated at some point, causing the unit test to fail

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

            var sampleModList = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "test.txt"), out var errors).ToList();

            Assert.AreEqual(2973, sampleModList.Count());
            string s = "";

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
            Assert.AreEqual(0,b.Count());
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
            Assert.AreEqual(0,b.Count());
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
            Loaders.LoadElements(Path.Combine(TestContext.CurrentContext.TestDirectory, "elements2.dat"));
            var sampleModList = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "z.txt"), out var errors).ToList();
            Assert.AreEqual(1, sampleModList.OfType<Modification>().Count());
            Protein protein = new Protein("MCSSSSSSSSSS", "accession", "organism", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>> { { 2, sampleModList.OfType<Modification>().ToList() } }, null, "name", "full_name", false, false, new List<DatabaseReference>(), new List<SequenceVariation>(), new List<DisulfideBond>());
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
            Loaders.LoadElements(Path.Combine(TestContext.CurrentContext.TestDirectory, "elements2.dat"));
            var sampleModList = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "z.txt"), out var errors).ToList();
            Protein protein = new Protein("MCSSSSSSSSSS", "accession", "organism", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>> { { 2, sampleModList.OfType<Modification>().ToList() } }, null, "name", "full_name", false, false, new List<DatabaseReference>(), new List<SequenceVariation>(), new List<DisulfideBond>());
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
    }
}
