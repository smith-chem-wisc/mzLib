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

        #region Public Methods

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

            Loaders.LoadPsiMod(Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml"));

            var uniprotPtms = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt")).ToList();

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

            var sampleModList = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "test.txt")).ToList();
            Console.WriteLine(sampleModList.First().ToString());
        }

        [Test]
        public void SampleModFileLoading()
        {
            var sampleModList = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "sampleModFile.txt")).ToList();
            Console.WriteLine(sampleModList.First().ToString());
        }

        [Test]
        public void SampleModFileLoadingFail1()
        {
            Assert.That(() => PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "sampleModFileFail1.txt")).ToList(),
                                            Throws.TypeOf<PtmListLoaderException>()
                                            .With.Property("Message")
                                            .EqualTo("Could not get motif from NxS"));
        }

        [Test]
        public void SampleModFileLoadingFail2()
        {
            Assert.That(() => PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "sampleModFileFail2.txt")).ToList(),
                                            Throws.TypeOf<PtmListLoaderException>()
                                            .With.Property("Message")
                                            .EqualTo("Could not get modification site from Anyplace."));
        }

        [Test]
        public void SampleModFileLoadingFail3()
        {
            Assert.That(() => PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "sampleModFileFail3.txt")).ToList(),
                                            Throws.TypeOf<FormatException>()
                                            .With.Property("Message")
                                            .EqualTo("Input string for chemical formula was in an incorrect format: $%#$%"));
        }

        [Test]
        public void SampleModFileLoadingFail4()
        {
            Assert.That(() => PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "m.txt")).ToList(),
                                            Throws.TypeOf<PtmListLoaderException>()
                                            .With.Property("Message")
                                            .EqualTo("0 or 238.229666 is not a valid monoisotopic mass"));
        }

        [Test]
        public void SampleModFileLoadingFail5()
        {
            Assert.That(() => PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "sampleModFileFail5.txt")).ToList(),
                                            Throws.TypeOf<PtmListLoaderException>()
                                            .With.Property("Message")
                                            .EqualTo("id is null"));
        }

        [Test]
        public void SampleModFileLoadingFail6()
        {
            Assert.That(() => PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "sampleModFileFail6.txt")).ToList(),
                                            Throws.TypeOf<PtmListLoaderException>()
                                            .With.Property("Message")
                                            .EqualTo("modificationType of lalaMod is null"));
        }

        [Test]
        public void CompactFormReading()
        {
            Assert.AreEqual(2, PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "sampleModFileDouble.txt")).Count());
        }

        [Test]
        public void CompactFormReading2()
        {
            Assert.AreEqual(2, PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "sampleModFileDouble2.txt")).Count());
        }

        [Test]
        public void Modification_read_write_into_proteinDb()
        {
            Loaders.LoadElements(Path.Combine(TestContext.CurrentContext.TestDirectory, "elements2.dat"));
            var sampleModList = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "z.txt")).ToList();
            Assert.AreEqual(1, sampleModList.OfType<ModificationWithMass>().Count());
            Protein protein = new Protein("MCSSSSSSSSSS", "accession", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>> { { 2, sampleModList.OfType<Modification>().ToList() } }, null, null, null, "name", "full_name", false, false, new List<DatabaseReference>());
            Assert.AreEqual(1, protein.OneBasedPossibleLocalizedModifications[2].OfType<ModificationWithMass>().Count());
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, ModificationWithMass>>>(), new List<Protein> { protein }, Path.Combine(TestContext.CurrentContext.TestDirectory, "test_modifications_with_proteins.xml"));
            List<Protein> new_proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "test_modifications_with_proteins.xml"), false, new List<Modification>(), false, null, null, out Dictionary<string, Modification> um);
            Assert.AreEqual(1, new_proteins.Count);
            Assert.AreEqual(1, new_proteins[0].OneBasedPossibleLocalizedModifications.Count);
            Assert.AreEqual(1, new_proteins[0].OneBasedPossibleLocalizedModifications.SelectMany(kv => kv.Value).Count());
            Assert.AreEqual("Type", new_proteins[0].OneBasedPossibleLocalizedModifications.SelectMany(kv => kv.Value).OfType<ModificationWithMass>().First().modificationType);
            Assert.AreEqual("Palmitoylation of C", new_proteins[0].OneBasedPossibleLocalizedModifications[2][0].id);
            Assert.AreEqual(1, new_proteins[0].OneBasedPossibleLocalizedModifications[2].OfType<ModificationWithMass>().Count());

            // Check that Modifications were saved after last load
            Assert.AreEqual(1, ProteinDbLoader.GetPtmListFromProteinXml(Path.Combine(TestContext.CurrentContext.TestDirectory, @"test_modifications_with_proteins.xml")).Count);
            Assert.True(ProteinDbLoader.GetPtmListFromProteinXml(Path.Combine(TestContext.CurrentContext.TestDirectory, @"test_modifications_with_proteins.xml"))[0] == new_proteins[0].OneBasedPossibleLocalizedModifications.SelectMany(kv => kv.Value).First());

            //But that we can still read modifications from other protein XMLs that exist
            Assert.AreEqual(0, ProteinDbLoader.GetPtmListFromProteinXml(Path.Combine(TestContext.CurrentContext.TestDirectory, "xml.xml")).Count);
        }

        #endregion Public Methods

    }
}