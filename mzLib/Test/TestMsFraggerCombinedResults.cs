using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using FlashLFQ;
using MassSpectrometry;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;
using Proteomics.AminoAcidPolymer;
using System.IO;
using Test.FileReadingTests;
using UsefulProteomicsDatabases;
using ChromatographicPeak = FlashLFQ.ChromatographicPeak;
using Stopwatch = System.Diagnostics.Stopwatch;
using TopDownProteomics;
using FlashLFQ.ResultsReading;

namespace Test
{
    internal class TestMsFraggerCombinedResults
    {
        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\EditedMSFraggerResults")]
        public void TestLoadResultsCount(string path)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            MsFraggerCombinedResults ms = new MsFraggerCombinedResults(filePath);
            ms.LoadResults();

            Assert.That(ms.AllPsmFiles.Count.Equals(2));
            Assert.That(ms.Results.Count.Equals(8));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\EditedMSFraggerResults")]
        public void TestLoadResults(string path)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            MsFraggerCombinedResults ms = new MsFraggerCombinedResults(filePath);
            ms.LoadResults();

            List<string> results = ms.Results.Select(psm => psm.FileName).ToList();

            Assert.That((results.Count(s => s.Contains("A_1"))).Equals(4));
            Assert.That((results.Count(s => s.Contains("A_2"))).Equals(4));

        }

    }
}
