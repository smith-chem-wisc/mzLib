using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using System.IO;
using System.Text;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;
using MzLibUtil;

namespace Test
{
    [TestFixture]
    public sealed class TestPtmListLoader
    {
        [Test]
        public static void Test_ReadAllModsFromFile()
        {
            string testModificationsFileLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", @"CommonArtifacts.txt");
            var a = PtmListLoaderGeneral.ReadModsFromFile(testModificationsFileLocation);
            Assert.AreEqual(33, a.Select(m=>m.Id).ToList().Count);


            testModificationsFileLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", @"CommonBiological.txt");
            var b = PtmListLoaderGeneral.ReadModsFromFile(testModificationsFileLocation);
            Assert.AreEqual(35, b.Select(m => m.Id).ToList().Count);
        }

        [Test]
        public static void Test_ModsFromFileAreSorted()
        {
            string testModificationsFileLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", @"CommonArtifacts.txt");
            var a = PtmListLoaderGeneral.ReadModsFromFile(testModificationsFileLocation);

            string id1 = a.First().Id.ToString();
            foreach (string modId in a.Select(m=>m.Id))
            {
                Assert.GreaterOrEqual(modId, id1);
            }


            testModificationsFileLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", @"CommonBiological.txt");
            var b = PtmListLoaderGeneral.ReadModsFromFile(testModificationsFileLocation);


            string id2 = b.First().Id.ToString();
            foreach (string modId in b.Select(m => m.Id))
            {
                Assert.GreaterOrEqual(modId, id2);
            }

        }
        [Test]
        public static void Test_ModsWithComments()
        {
            string testModificationsFileLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", @"ModsWithComments.txt");
            var a = PtmListLoaderGeneral.ReadModsFromFile(testModificationsFileLocation).ToList();
            Assert.AreEqual(4, a.Select(m => m.Id).ToList().Count);

            Assert.AreEqual("Deamidation", a[0].Id.ToString());
            Assert.AreEqual("Sodium", a[2].Id.ToString());//this has trailing whitespace that shouldn't be in the name

            //Make sure comments are okay on DR key and that key value pairs are still split correctly
            var someMod = a[2];
            var test = (Proteomics.ModificationGeneral)someMod;
            var residValueTest = test.DatabaseReference.First().Value.First();
            var residKeyTest = test.DatabaseReference.First().Key;
            Assert.AreEqual("RESID", residKeyTest);
            Assert.AreEqual("AA0441", residValueTest);
        }


        [Test]
        public static void Test_ReadAllModsFromFileGeneral()
        {
            string testModificationsFileLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", @"CommonArtifacts.txt");
            var a = PtmListLoaderGeneral.ReadModsFromFile(testModificationsFileLocation);
            Assert.AreEqual(33, a.Select(m => m.Id).ToList().Count);


            testModificationsFileLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", @"CommonBiological.txt");
            var b = PtmListLoaderGeneral.ReadModsFromFile(testModificationsFileLocation);
            Assert.AreEqual(35, b.Select(m => m.Id).ToList().Count);
        }

        [Test]
        public static void Test_ModsWithCommentsGeneral()
        {
            string testModificationsFileLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", @"ModsWithComments.txt");
            var a = PtmListLoaderGeneral.ReadModsFromFile(testModificationsFileLocation).ToList();
            Assert.AreEqual(4, a.Select(m => m.Id).ToList().Count);

            Assert.AreEqual("Deamidation", a[0].Id.ToString());
            Assert.AreEqual("Sodium", a[2].Id.ToString());//this has trailing whitespace that shouldn't be in the name

            //Make sure comments are okay on DR key and that key value pairs are still split correctly
            var someMod = a[2];
            var test = (Proteomics.ModificationGeneral)someMod;
            var residValueTest = test.DatabaseReference.First().Value.First();
            var residKeyTest = test.DatabaseReference.First().Key;
            Assert.AreEqual("RESID", residKeyTest);
            Assert.AreEqual("AA0441", residValueTest);
        }

        [Test]
        public void SampleModFileLoadingGeneral()
        {
            PtmListLoaderGeneral.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFile.txt"));
        }

        [Test]
        public void SampleModFileLoadingFail1General() //TG is not valide
        {
            var a = PtmListLoaderGeneral.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFileFail1.txt")).ToList();
            Assert.AreEqual(false, a.First().ValidModification);
        }


        [Test]
        public void SampleModFileLoadingFail2General()
        {
            var a = PtmListLoaderGeneral.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFileFail2.txt")).ToList();
            Assert.AreEqual(false, a.First().ValidModification); ;
        }

        [Test]
        public void PTMListLoader_ModWithComments_Equals_ModWithoutComments()
        {
            var a = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "SampleMod_Comments.txt")).ToList();
            var b = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "SampleMod_NoComments.txt")).ToList();
            Assert.IsTrue(a.First().Equals(b.First()));
        }

        [Test]
        public void PTMListLoaderGeneral_ModWithComments_Equals_ModWithoutComments()
        {
            var a = PtmListLoaderGeneral.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "SampleMod_Comments.txt")).ToList();
            var b = PtmListLoaderGeneral.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "SampleMod_NoComments.txt")).ToList();
            Assert.IsTrue(a.First().Equals(b.First()));
        }


        [Test]
        public void SampleModFileLoadingFail3General()
        {
            Assert.That(() => PtmListLoaderGeneral.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFileFail3.txt")).ToList(),
                                            Throws.TypeOf<MzLibException>()
                                            .With.Property("Message")
                                            .EqualTo("Input string for chemical formula was in an incorrect format: $%&$%"));
        }

        [Test]
        public void SampleModFileLoadingFail4General()
        {
            Assert.That(() => PtmListLoaderGeneral.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "m.txt")).ToList(),
                                            Throws.TypeOf<MzLibException>()
                                            .With.Property("Message")
                                            .EqualTo("0 or 238.229666 is not a valid monoisotopic mass"));
        }

        [Test]
        public void SampleModFileLoadingFail5General()
        {
            var a = PtmListLoaderGeneral.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFileFail5.txt")).ToList();
            Assert.AreEqual(false, a.First().ValidModification); // ID is missing

        }

        [Test]
        public void SampleModFileLoadingFail6General()
        {
            var a = PtmListLoaderGeneral.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFileFail6.txt")).ToList();
            Assert.AreEqual(false, a.First().ValidModification); // modification type is missing
        }
        [Test]
        public void SampleModFileLoadingFail5General_missingPosition()
        {
            var a = PtmListLoaderGeneral.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFileFail_missingPosition.txt")).ToList();
            Assert.AreEqual(false, a.First().ValidModification); // ID is missing

        }
        [Test]
        public void SampleModFileLoadingFail5General_missingMonoisotopicMassAndChemicalFormula()
        {
            var a = PtmListLoaderGeneral.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFileFail_missingChemicalFormulaAndMonoisotopicMass.txt")).ToList();
            Assert.AreEqual(false, a.First().ValidModification); // ID is missing

        }
        [Test]
        public void CompactFormReadingGeneral()
        {
            Assert.AreEqual(2, PtmListLoaderGeneral.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFileDouble.txt")).Count());
        }

        [Test]
        public void CompactFormReading2General()
        {
            Assert.AreEqual(2, PtmListLoaderGeneral.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFileDouble2.txt")).Count());
        }

    }
}
