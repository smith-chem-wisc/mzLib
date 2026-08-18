using System.IO;
using System.Linq;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests.ProForma
{
    /// <summary>
    /// Reader side of the ProForma round-trip: the "ProForma" psmtsv column (written by
    /// MetaMorpheus's PsmTsvWriter) parses back into <see cref="SpectrumMatchFromTsv.ProForma"/>,
    /// and is optional so pre-ProForma result files still read.
    /// </summary>
    [TestFixture]
    internal class ProFormaReaderColumnTests
    {
        private const string SampleProForma = "EM[Oxidation]EVEES[Phospho]PEK";

        private static string SearchResult(string name) =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\SearchResults", name);

        [Test]
        public void ProForma_Column_IsParsed()
        {
            // Take a real result file and append a ProForma column (reader maps by header name, so position is irrelevant).
            var lines = File.ReadAllLines(SearchResult("BottomUpExample.psmtsv")).ToList();
            lines[0] += "\t" + SpectrumMatchFromTsvHeader.ProForma;
            for (int i = 1; i < lines.Count; i++)
                if (lines[i].Length > 0) lines[i] += "\t" + SampleProForma;

            string tmp = Path.Combine(TestContext.CurrentContext.TestDirectory,
                $"proforma_roundtrip_{TestContext.CurrentContext.Test.ID}.psmtsv");
            File.WriteAllLines(tmp, lines);
            try
            {
                var psms = SpectrumMatchTsvReader.ReadTsv<PsmFromTsv>(tmp, out _);
                Assert.That(psms, Is.Not.Empty);
                Assert.That(psms.All(p => p.ProForma == SampleProForma), Is.True);
            }
            finally
            {
                File.Delete(tmp);
            }
        }

        [Test]
        public void ProForma_Column_IsOptional_ForOlderFiles()
        {
            // BottomUpExample.psmtsv predates the column; reading must not throw and ProForma is null.
            var psms = SpectrumMatchTsvReader.ReadTsv<PsmFromTsv>(SearchResult("BottomUpExample.psmtsv"), out _);
            Assert.That(psms, Is.Not.Empty);
            Assert.That(psms.All(p => p.ProForma == null), Is.True);
        }
    }
}
