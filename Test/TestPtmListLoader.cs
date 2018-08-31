using MzLibUtil;
using NUnit.Framework;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public sealed class TestPtmListLoader
    {
        [Test]
        public static void Test_ReadAllModsFromFile()
        {
            string testModificationsFileLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", @"CommonArtifacts.txt");
            var a = PtmListLoader.ReadModsFromFile(testModificationsFileLocation);
            Assert.AreEqual(33, a.Select(m => m.IdWithMotif).ToList().Count);

            testModificationsFileLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", @"CommonBiological.txt");
            var b = PtmListLoader.ReadModsFromFile(testModificationsFileLocation);
            Assert.AreEqual(35, b.Select(m => m.IdWithMotif).ToList().Count);
        }

        [Test]
        public static void Test_ModsFromFileAreSorted()
        {
            string testModificationsFileLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", @"CommonArtifacts.txt");
            var a = PtmListLoader.ReadModsFromFile(testModificationsFileLocation);

            string id1 = a.First().IdWithMotif.ToString();
            foreach (string modId in a.Select(m => m.IdWithMotif))
            {
                Assert.GreaterOrEqual(modId, id1);
            }

            testModificationsFileLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", @"CommonBiological.txt");
            var b = PtmListLoader.ReadModsFromFile(testModificationsFileLocation);

            string id2 = b.First().IdWithMotif.ToString();
            foreach (string modId in b.Select(m => m.IdWithMotif))
            {
                Assert.GreaterOrEqual(modId, id2);
            }
        }

        [Test]
        public static void Test_ModsWithComments()
        {
            string testModificationsFileLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", @"ModsWithComments.txt");
            var a = PtmListLoader.ReadModsFromFile(testModificationsFileLocation).ToList();
            Assert.AreEqual(4, a.Select(m => m.IdWithMotif).ToList().Count);

            Assert.AreEqual("Deamidation on N", a[0].IdWithMotif.ToString());
            Assert.AreEqual("Sodium on D", a[2].IdWithMotif.ToString());//this has trailing whitespace that shouldn't be in the name

            //Make sure comments are okay on DR key and that key value pairs are still split correctly
            var someMod = a[2];
            var test = (Proteomics.Modification)someMod;
            var residValueTest = test.DatabaseReference.First().Value.First();
            var residKeyTest = test.DatabaseReference.First().Key;
            Assert.AreEqual("RESID", residKeyTest);
            Assert.AreEqual("AA0441", residValueTest);
        }

        [Test]
        public static void Test_ReadAllModsFromFileGeneral()
        {
            string testModificationsFileLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", @"CommonArtifacts.txt");
            var a = PtmListLoader.ReadModsFromFile(testModificationsFileLocation);
            Assert.AreEqual(33, a.Select(m => m.IdWithMotif).ToList().Count);

            testModificationsFileLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", @"CommonBiological.txt");
            var b = PtmListLoader.ReadModsFromFile(testModificationsFileLocation);
            Assert.AreEqual(35, b.Select(m => m.IdWithMotif).ToList().Count);
        }

        [Test]
        public static void Test_ModsWithCommentsGeneral()
        {
            string testModificationsFileLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", @"ModsWithComments.txt");
            var a = PtmListLoader.ReadModsFromFile(testModificationsFileLocation).ToList();
            Assert.AreEqual(4, a.Select(m => m.IdWithMotif).ToList().Count);

            Assert.AreEqual("Deamidation on N", a[0].IdWithMotif.ToString());
            Assert.AreEqual("Sodium on D", a[2].IdWithMotif.ToString());//this has trailing whitespace that shouldn't be in the name

            //Make sure comments are okay on DR key and that key value pairs are still split correctly
            var someMod = a[2];
            var test = (Proteomics.Modification)someMod;
            var residValueTest = test.DatabaseReference.First().Value.First();
            var residKeyTest = test.DatabaseReference.First().Key;
            Assert.AreEqual("RESID", residKeyTest);
            Assert.AreEqual("AA0441", residValueTest);
        }

        [Test]
        public static void SampleModFileLoadingGeneral()
        {
            PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFile.txt"));
        }

        [Test]
        public static void SampleModFileLoadingFail1General() //TG is not valide
        {
            var a = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFileFail1.txt")).ToList();
            Assert.AreEqual(0, a.Count());
        }

        [Test]
        public static void SampleModFileLoadingFail2General()
        {
            var a = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFileFail2.txt")).ToList();
            Assert.AreEqual(0, a.Count()); ;
        }

        [Test]
        public static void PTMListLoader_ModWithComments_Equals_ModWithoutComments()
        {
            var a = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "SampleMod_Comments.txt")).ToList();
            var b = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "SampleMod_NoComments.txt")).ToList();
            Assert.IsTrue(a.First().Equals(b.First()));
        }

        [Test]
        public static void PTMListLoaderGeneral_ModWithComments_Equals_ModWithoutComments()
        {
            var a = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "SampleMod_Comments.txt")).ToList();
            var b = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "SampleMod_NoComments.txt")).ToList();
            Assert.IsTrue(a.First().Equals(b.First()));
        }

        [Test]
        public static void SampleModFileLoadingFail3General()
        {
            Assert.That(() => PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFileFail3.txt")).ToList(),
                                            Throws.TypeOf<MzLibException>()
                                            .With.Property("Message")
                                            .EqualTo("Input string for chemical formula was in an incorrect format: $%&$%"));
        }

        [Test]
        public static void SampleModFileLoadingFail4General()
        {
            Assert.That(() => PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "m.txt")).ToList(),
                                            Throws.TypeOf<MzLibException>()
                                            .With.Property("Message")
                                            .EqualTo("0 or 238.229666 is not a valid monoisotopic mass"));
        }

        [Test]
        public static void SampleModFileLoadingFail5General()
        {
            var a = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFileFail5.txt")).ToList();
            Assert.AreEqual(0, a.Count()); // ID is missing
        }

        [Test]
        public static void SampleModFileLoadingFail6General()
        {
            var a = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFileFail6.txt")).ToList();
            Assert.AreEqual(0, a.Count()); // modification type is missing
        }

        [Test]
        public static void SampleModFileLoadingFail5General_missingPosition()
        {
            var a = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFileFail_missingPosition.txt")).ToList();
            Assert.AreEqual(0, a.Count()); // ID is missing
        }

        [Test]
        public static void SampleModFileLoadingFail5General_missingMonoisotopicMassAndChemicalFormula()
        {
            var a = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFileFail_missingChemicalFormulaAndMonoisotopicMass.txt")).ToList();
            Assert.AreEqual(0, a.Count()); // ID is missing
        }

        [Test]
        public static void CompactFormReadingGeneral()
        {
            Assert.AreEqual(2, PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFileDouble.txt")).Count());
        }

        [Test]
        public static void CompactFormReading2General()
        {
            Assert.AreEqual(2, PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFileDouble2.txt")).Count());
        }

        [Test]
        public static void TestReadingIdWithMotif()
        {
            string modText = "ID   Detached EVK or XleDK\r\nPP   Peptide N-terminal.\r\nTG   evkX or vekX or ldkX or dlkX or idkX or dikX\r\nMT   Detached\r\nNL   C16H28N4O5\r\nCF   C16H28N4O5\r\n" + @"//";

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "detacher.txt");
            File.WriteAllLines(path, new string[] { modText });
            
            var mods = PtmListLoader.ReadModsFromFile(path).ToList();
            var motifs = mods.Select(p => p.Target.ToString()).Distinct().ToList();
            var ids = mods.Select(p => p.IdWithMotif).Distinct().ToList();

            Assert.That(mods.Count == 6);
            Assert.That(motifs.Count == 6);
            Assert.That(ids.Count == 6);
        }
    }
}