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


    }
}
