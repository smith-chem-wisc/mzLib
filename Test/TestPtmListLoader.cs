using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using System.IO;
using System.Text;
using System.Threading.Tasks;
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
            Assert.AreEqual(33, a.Select(m=>m.id).ToList().Count);


            testModificationsFileLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", @"CommonBiological.txt");
            var b = PtmListLoader.ReadModsFromFile(testModificationsFileLocation);
            Assert.AreEqual(35, b.Select(m => m.id).ToList().Count);
        }

        [Test]
        public static void Test_ModsFromFileAreSorted()
        {
            string testModificationsFileLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", @"CommonArtifacts.txt");
            var a = PtmListLoader.ReadModsFromFile(testModificationsFileLocation);

            string id1 = a.First().id.ToString();
            foreach (string modId in a.Select(m=>m.id))
            {
                Assert.GreaterOrEqual(modId, id1);
            }


            testModificationsFileLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", @"CommonBiological.txt");
            var b = PtmListLoader.ReadModsFromFile(testModificationsFileLocation);


            string id2 = b.First().id.ToString();
            foreach (string modId in b.Select(m => m.id))
            {
                Assert.GreaterOrEqual(modId, id2);
            }

        }
        [Test]
        public static void Test_ModsWithComments()
        {
            string testModificationsFileLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", @"ModsWithComments.txt");
            var a = PtmListLoader.ReadModsFromFile(testModificationsFileLocation).ToList();
            Assert.AreEqual(4, a.Select(m => m.id).ToList().Count);

            Assert.AreEqual("Deamidation on N", a[0].id.ToString());
            Assert.AreEqual("Sodium on D", a[2].id.ToString());//this has trailing whitespace that shouldn't be in the name

            //Make sure comments are okay on DR key and that key value pairs are still split correctly
            var someMod = a[2];
            var test = (Proteomics.ModificationWithMass)someMod;
            var residValueTest = test.linksToOtherDbs.First().Value.First();
            var residKeyTest = test.linksToOtherDbs.First().Key;
            Assert.AreEqual("RESID", residKeyTest);
            Assert.AreEqual("AA0441", residValueTest);
        }

    }
}
