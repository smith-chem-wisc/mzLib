using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using System;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public sealed class TestPtmListLoader
    {
        [Test]
        public static void SampleModFileLoading()
        {
            PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFile.txt"), out var errors);
        }

        [Test]
        [TestCase("CommonArtifacts.txt", 33)]
        [TestCase("CommonBiological.txt", 35)]
        public static void Test_ReadAllModsFromFile(string filename, int modCount)
        {
            string testModificationsFileLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", filename);
            var a = PtmListLoader.ReadModsFromFile(testModificationsFileLocation, out var errors);
            Assert.AreEqual(modCount, a.Count());
        }

        [Test]
        [TestCase("CommonArtifacts.txt")]
        [TestCase("CommonBiological.txt")]
        public static void Test_ModsFromFileAreSorted(string filename)
        {
            string testModificationsFileLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", filename);
            var a = PtmListLoader.ReadModsFromFile(testModificationsFileLocation, out var errors);

            string id1 = a.First().IdWithMotif;
            foreach (string modId in a.Select(m => m.IdWithMotif))
            {
                Assert.GreaterOrEqual(modId, id1);
            }
        }

        [Test]
        public static void Test_ModsWithComments()
        {
            string testModificationsFileLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", @"ModsWithComments.txt");
            var a = PtmListLoader.ReadModsFromFile(testModificationsFileLocation, out var errors).ToList();
            Assert.AreEqual(4, a.Select(m => m.IdWithMotif).ToList().Count);

            Assert.AreEqual("Deamidation on N", a[0].IdWithMotif.ToString());
            Assert.AreEqual("Sodium on D", a[2].IdWithMotif.ToString());//this has trailing whitespace that shouldn't be in the name

            //Make sure comments are okay on DR key and that key value pairs are still split correctly
            Modification someMod = a[2];
            Modification test = someMod;
            var residValueTest = test.DatabaseReference.First().Value.First();
            var residKeyTest = test.DatabaseReference.First().Key;
            Assert.AreEqual("RESID", residKeyTest);
            Assert.AreEqual("AA0441", residValueTest);
        }

        [Test]
        [TestCase("sampleModFileFail1.txt")] // TG is not valid
        [TestCase("sampleModFileFail2.txt")]
        [TestCase("sampleModFileFail5.txt")] // ID is missing
        [TestCase("sampleModFileFail6.txt")] // modification type is missing
        [TestCase("sampleModFileFail_missingPosition.txt")] // missing position
        [TestCase("sampleModFileFail_missingChemicalFormulaAndMonoisotopicMass.txt")]
        public static void SampleModFileLoadingFail1General(string filename)
        {
            var a = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", filename), out var errors).ToList();
            Assert.AreEqual(0, a.Count);
        }

        [Test]
        public static void PTMListLoader_ModWithComments_Equals_ModWithoutComments()
        {
            var a = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "SampleMod_Comments.txt"), out var errors).ToList();
            var b = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "SampleMod_NoComments.txt"), out var errors2).ToList();
            Assert.IsTrue(a.First().Equals(b.First()));
        }

        [Test]
        [TestCase("sampleModFileFail3.txt", "Input string for chemical formula was in an incorrect format: $%&$%")]
        [TestCase("m.txt", "0 or 238.229666 is not a valid monoisotopic mass")]
        public static void SampleModFileLoadingFail3General(string filename, string errorMessage)
        {
            Assert.That(() =>
                PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", filename), out var errors).ToList(),
                Throws.TypeOf<MzLibException>().With.Property("Message").EqualTo(errorMessage));
        }

        [Test]
        [TestCase("sampleModFileDouble.txt")]
        [TestCase("sampleModFileDouble2.txt")]
        public static void CompactFormReadingGeneral(string filename)
        {
            Assert.AreEqual(2, PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", filename), out var errors).Count());
        }

        [Test]
        public static void TestReadingInvalidModifications()
        {
            string modText = "ID   mod\r\nMT   type\r\nPP   Anywhere.\r\nTG   X\r\nCF   H1\r\n" + @"//";
            var mods = PtmListLoader.ReadModsFromString(modText, out var errors1);
            Assert.AreEqual(1, mods.Count());
            Assert.AreEqual(0, errors1.Count);

            modText = "ID   mod\r\nMT   type\r\nPP   Anywhere.\r\nTG   X\r\n" + @"//"; // no monoisotopic mass, so invalid
            mods = PtmListLoader.ReadModsFromString(modText, out var errors2);
            Assert.AreEqual(0, mods.Count());
            Assert.AreEqual(1, errors2.Count);
            Assert.IsFalse(errors2.Single().Item1.ValidModification);
            Assert.IsTrue(errors2.Single().Item2.Split(new[] { "\r\n" }, StringSplitOptions.None).Any(x => x.StartsWith("#"))); // has an error comment
        }

        [Test]
        public static void TestReadingIdWithMotif()
        {
            string modText = "ID   Detached EVK or XleDK\r\nPP   Peptide N-terminal.\r\nTG   evkX or vekX or ldkX or dlkX or idkX or dikX\r\nMT   Detached\r\nNL   C16H28N4O5\r\nCF   C16H28N4O5\r\n" + @"//";

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "detacher.txt");
            File.WriteAllLines(path, new string[] { modText });

            var mods = PtmListLoader.ReadModsFromFile(path, out var errors).ToList();
            var motifs = mods.Select(p => p.Target.ToString()).Distinct().ToList();
            var ids = mods.Select(p => p.IdWithMotif).Distinct().ToList();

            Assert.That(mods.Count == 6);
            Assert.That(motifs.Count == 6);
            Assert.That(ids.Count == 6);
        }

        [Test]
        public static void TestInvalidModTypeError()
        {
            string mod = "ID   Deamidation\r\nTG   N or Q\r\nPP   Anywhere.\r\nMT   Mod:\r\nCF   H-1 N-1 O1\r\n//";
            string filepath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestInvalidModTypeError\ptmlist.txt");
            Directory.CreateDirectory(Directory.GetParent(filepath).FullName);
            File.WriteAllLines(filepath, new string[] { mod });

            var ptms = PtmListLoader.ReadModsFromFile(filepath, out var warnings).ToList();
            Assert.That(ptms.Count == 0);
            Assert.That(warnings.Count == 2);

            Assert.That(warnings.First().Item2.Contains("Modification type cannot contain ':'!"));
        }
    }
}