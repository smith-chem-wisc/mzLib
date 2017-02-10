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
using System;
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
            var sampleModList = PtmListLoader.ReadMods(Path.Combine(TestContext.CurrentContext.TestDirectory, "test.txt")).ToList();
            Console.WriteLine(sampleModList.First().ToString());
        }
        [Test]
        public void SampleModFileLoading()
        {
            var sampleModList = PtmListLoader.ReadMods(Path.Combine(TestContext.CurrentContext.TestDirectory, "sampleModFile.txt")).ToList();
            Console.WriteLine(sampleModList.First().ToString());
        }
        [Test]
        public void SampleModFileLoadingFail1()
        {
            Assert.That(() => PtmListLoader.ReadMods(Path.Combine(TestContext.CurrentContext.TestDirectory, "sampleModFileFail1.txt")).ToList(),
                                            Throws.TypeOf<PtmListLoaderException>()
                                            .With.Property("Message")
                                            .EqualTo("Could not get motif from NxS"));
        }
        [Test]
        public void SampleModFileLoadingFail2()
        {
            Assert.That(() => PtmListLoader.ReadMods(Path.Combine(TestContext.CurrentContext.TestDirectory, "sampleModFileFail2.txt")).ToList(),
                                            Throws.TypeOf<PtmListLoaderException>()
                                            .With.Property("Message")
                                            .EqualTo("Could not get modification site from Anyplace."));
        }
        [Test]
        public void SampleModFileLoadingFail3()
        {
            Assert.That(() => PtmListLoader.ReadMods(Path.Combine(TestContext.CurrentContext.TestDirectory, "sampleModFileFail3.txt")).ToList(),
                                            Throws.TypeOf<FormatException>()
                                            .With.Property("Message")
                                            .EqualTo("Input string for chemical formula was in an incorrect format"));
        }

        #endregion Public Methods
    }
}