using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Threading.Tasks;
using NUnit.Framework;
using Proteomics;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class TestProteinDbLoaderConcurrency
    {
        /// <summary>
        /// Several loads of one .gz at once used to fight over a single "temp.xml" beside the input,
        /// and the losers threw IOException from File.Create.
        /// </summary>
        [Test]
        public static void ConcurrentGzXmlLoadsDoNotShareATempFile()
        {
            string testDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "ConcurrentGzXmlLoads");
            if (Directory.Exists(testDirectory))
            {
                Directory.Delete(testDirectory, true);
            }
            Directory.CreateDirectory(testDirectory);
            string gzDatabase = Path.Combine(testDirectory, "concurrent.xml.gz");

            using (var uncompressed = new FileStream(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "xml.xml"), FileMode.Open, FileAccess.Read))
            using (var compressedFile = File.Create(gzDatabase))
            using (var compressor = new GZipStream(compressedFile, CompressionMode.Compress))
            {
                uncompressed.CopyTo(compressor);
            }

            int expected = ProteinDbLoader.LoadProteinXML(gzDatabase, true, DecoyType.None, null, false, null, out _).Count;
            Assert.That(expected, Is.GreaterThan(0));

            var counts = new int[4];
            var exceptions = new ConcurrentBag<Exception>();
            Parallel.For(0, counts.Length, i =>
            {
                try
                {
                    counts[i] = ProteinDbLoader.LoadProteinXML(gzDatabase, true, DecoyType.None, null, false, null, out _).Count;
                }
                catch (Exception e)
                {
                    exceptions.Add(e);
                }
            });

            Assert.That(exceptions.Select(e => e.ToString()).ToList(), Is.Empty);
            Assert.That(counts, Is.All.EqualTo(expected));

            // the decompressed file is cleaned up, so only the input is left
            Assert.That(Directory.GetFiles(testDirectory).Select(Path.GetFileName).ToList(), Is.EqualTo(new List<string> { "concurrent.xml.gz" }));

            Directory.Delete(testDirectory, true);
        }
    }
}
