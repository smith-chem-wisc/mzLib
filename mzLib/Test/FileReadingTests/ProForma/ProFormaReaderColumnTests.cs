using System.IO;
using System.Linq;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests.ProForma
{
    /// <summary>
    /// Reader side of the ProForma round-trip: the "ProForma" psmtsv column (written by
    /// MetaMorpheus's PsmTsvWriter) parses back into <see cref="SpectrumMatchFromTsv.ProForma"/>.
    /// The column is optional: pre-ProForma result files (MetaMorpheus 1.1.11 and earlier) still read,
    /// and their ProForma is computed from the Full Sequence column instead.
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

                // A file's value is carried verbatim onto a disambiguated copy, never recomputed.
                var copy = new PsmFromTsv(psms[0], psms[0].FullSequence);
                Assert.That(copy.ProForma, Is.EqualTo(SampleProForma));
            }
            finally
            {
                File.Delete(tmp);
            }
        }

        [Test]
        public void ProForma_ColumnAbsent_IsComputedFromFullSequence()
        {
            // BottomUpExample.psmtsv predates the column; reading must not throw, and ProForma comes
            // from Full Sequence with UNIMOD accessions where the modification has one.
            var psms = SpectrumMatchTsvReader.ReadTsv<PsmFromTsv>(SearchResult("BottomUpExample.psmtsv"), out _);
            Assert.That(psms, Is.Not.Empty);
            Assert.That(psms.All(p => p.ProForma != null), Is.True);
            Assert.That(psms.Select(p => p.ProForma),
                Does.Contain("AYHEQLSVAEITNAC[UNIMOD:4]FEPANQMVK"));
        }

        // The MetaMorpheus notations the conversion has to survive, ported from the dataRepo project's
        // translator tests (datarepo/proforma.py), which was written because this field was never
        // populated for MetaMorpheus 1.1.11 output.
        [TestCase("PEPTIDEK", "PEPTIDEK", TestName = "FullSequence_Unmodified_IsUnchanged")]
        [TestCase("KLADQC[Common Fixed:Carbamidomethyl on C]TGLQ", "KLADQC[UNIMOD:4]TGLQ",
            TestName = "FullSequence_ResidueMod_BecomesUnimodAccession")]
        [TestCase("M[Common Variable:Oxidation on M]PEC[Common Fixed:Carbamidomethyl on C]K", "M[UNIMOD:35]PEC[UNIMOD:4]K",
            TestName = "FullSequence_TwoMods_KeepTheirPositions")]
        [TestCase("[Common Artifact:Ammonia loss on C]C[Common Fixed:Carbamidomethyl on C]AK", "[UNIMOD:385]-C[UNIMOD:4]AK",
            TestName = "FullSequence_LeadingBracket_IsNTerminal")]
        [TestCase("PEPD[Metal:Calcium on D]K", "PEPD[UNIMOD:951]K", TestName = "FullSequence_CalciumOnD_Converts")]
        // UniProt-named modifications resolve against the UniProt ptmlist and are written by their
        // UNIMOD accession (#1401); they used to keep the "UniProt:" name, which is not a ProForma vocabulary.
        [TestCase("[UniProt:N-acetylalanine on A]AAAGEAR", "[UNIMOD:1]-AAAGEAR",
            TestName = "FullSequence_UniProtNTerminalMod_BecomesUnimodAccession")]
        [TestCase("[UniProt:N-acetylserine on S]SEQK", "[UNIMOD:1]-SEQK",
            TestName = "FullSequence_UniProtNAcetylserine_BecomesUnimod1")]
        [TestCase("PEPK[UniProt:N6,N6-dimethyllysine on K]R", "PEPK[UNIMOD:36]R",
            TestName = "FullSequence_UniProtDimethyllysine_BecomesUnimod36")]
        [TestCase("PEPM[Common Variable:Oxidation on M]K", "PEPM[UNIMOD:35]K",
            TestName = "FullSequence_MetaMorpheusOxidation_StaysUnimod35")]
        // A UniProt modification with a PSI-MOD reference but no Unimod one takes the PSI-MOD accession,
        // written once: the loader stores the id already prefixed ("MOD:01694"), not as "MOD:MOD:01694".
        [TestCase("PEC[UniProt:S-(coelenterazin-3a-yl)cysteine on C]K", "PEC[MOD:01694]K",
            TestName = "FullSequence_UniProtModWithOnlyPsiMod_BecomesPsiModAccession")]
        // A UniProt modification with no ontology reference at all is written by name, as
        // ProFormaConverter.ToProFormaTerm writes the same modification.
        [TestCase("PED[UniProt:Cholesterol aspartate ester on D]K", "PED[Cholesterol aspartate ester]K",
            TestName = "FullSequence_UniProtModWithoutAccession_IsWrittenByName")]
        // MetaMorpheus writes a C-terminal modification after a '-', which is a terminus marker, not a residue.
        // The mod is MetaMorpheus's own (not in mzLib's set) and shares its name with the UniProt entry,
        // which carries Unimod 34 (methyl).
        [TestCase("KPVADYFL-[Common Artifact:Leucine methyl ester on L]", "KPVADYFL-[UNIMOD:34]",
            TestName = "FullSequence_CTerminalMod_StaysOnTheCTerminus")]
        [TestCase("PEPK[Custom:Nameless on K]R", "PEPK[Custom:Nameless on K]R", TestName = "FullSequence_UnknownMod_KeepsName")]
        public void ProFormaFromFullSequence_ConvertsMetaMorpheusNotation(string fullSequence, string expected)
        {
            Assert.That(SpectrumMatchFromTsv.ProFormaFromFullSequence(fullSequence), Is.EqualTo(expected));
        }

        [Test]
        public void ProForma_ColumnAbsent_UniProtModsAreWrittenAsAccessions()
        {
            // A real top-down G-PTM-D result without a ProForma column, full of UniProt-named modifications.
            var psms = SpectrumMatchTsvReader.ReadTsv<PsmFromTsv>(SearchResult("TDGPTMDSearchResults.psmtsv"), out _);
            var converted = psms.Where(p => p.FullSequence.Contains("[UniProt:") && p.ProForma != null).ToList();
            Assert.That(converted, Is.Not.Empty);

            // "UniProt:" is not a ProForma vocabulary, so it must never reach the output.
            Assert.That(converted.Where(p => p.ProForma!.Contains("UniProt:")).Select(p => p.ProForma), Is.Empty);
            var acetylSerine = converted.First(p => p.FullSequence.StartsWith("[UniProt:N-acetylserine on S]"));
            Assert.That(acetylSerine.ProForma, Does.StartWith("[UNIMOD:1]-S"));
            Assert.That(converted.Select(p => p.ProForma), Has.Some.Contains("K[UNIMOD:36]"),
                "N6,N6-dimethyllysine is Unimod 36");

            // Bottom-up file: Tele-methylhistidine is Unimod 34 (methyl).
            var bottomUp = SpectrumMatchTsvReader.ReadTsv<PsmFromTsv>(SearchResult("BottomUpExample.psmtsv"), out _);
            Assert.That(bottomUp.Select(p => p.ProForma), Does.Contain("YPIEH[UNIMOD:34]GIVTNWDDMEK"));
        }

        [TestCase(null)]
        [TestCase("")]
        [TestCase("PEPTIDE|PEPTLDE")]
        public void ProFormaFromFullSequence_EmptyOrAmbiguous_IsNull(string? fullSequence)
        {
            Assert.That(SpectrumMatchFromTsv.ProFormaFromFullSequence(fullSequence), Is.Null);
        }

        [Test]
        public void ProForma_ColumnAbsent_AmbiguousMatch_IsNull_ButEachCandidateConverts()
        {
            // One ProForma string cannot describe two candidates, so the ambiguous row has none, while a
            // disambiguated copy converts its own full sequence.
            var psms = SpectrumMatchTsvReader.ReadTsv<PsmFromTsv>(SearchResult("OneOverK0Example.psmtsv"), out _);
            var ambiguous = psms.Single(p => p.FullSequence.Contains('|'));
            Assert.That(ambiguous.ProForma, Is.Null);

            string firstCandidate = ambiguous.FullSequence.Split('|')[0];
            var disambiguated = new PsmFromTsv(ambiguous, firstCandidate, 0);
            Assert.That(disambiguated.ProForma, Is.Not.Null);
            Assert.That(disambiguated.ProForma, Is.EqualTo(SpectrumMatchFromTsv.ProFormaFromFullSequence(firstCandidate)));
        }

        [Test]
        public void ProForma_ColumnPresentButBlank_IsComputedLikeAnAbsentColumn()
        {
            // A blank cell carries no value from the file, so it must give the same answer as a file
            // without the column, including on a disambiguated candidate of an ambiguous row.
            var lines = File.ReadAllLines(SearchResult("OneOverK0Example.psmtsv")).ToList();
            lines[0] += "\t" + SpectrumMatchFromTsvHeader.ProForma;
            for (int i = 1; i < lines.Count; i++)
                if (lines[i].Length > 0) lines[i] += "\t";

            string tmp = Path.Combine(TestContext.CurrentContext.TestDirectory,
                $"proforma_blank_{TestContext.CurrentContext.Test.ID}.psmtsv");
            File.WriteAllLines(tmp, lines);
            try
            {
                var psms = SpectrumMatchTsvReader.ReadTsv<PsmFromTsv>(tmp, out _);
                Assert.That(psms, Is.Not.Empty);
                foreach (var psm in psms)
                    Assert.That(psm.ProForma, Is.EqualTo(SpectrumMatchFromTsv.ProFormaFromFullSequence(psm.FullSequence)));
                Assert.That(psms.Count(p => p.ProForma != null), Is.GreaterThan(0));

                var ambiguous = psms.Single(p => p.FullSequence.Contains('|'));
                string firstCandidate = ambiguous.FullSequence.Split('|')[0];
                var disambiguated = new PsmFromTsv(ambiguous, firstCandidate, 0);
                Assert.That(disambiguated.ProForma, Is.Not.Null);
                Assert.That(disambiguated.ProForma, Is.EqualTo(SpectrumMatchFromTsv.ProFormaFromFullSequence(firstCandidate)));
            }
            finally
            {
                File.Delete(tmp);
            }
        }

        [TestCase(@"FileReadingTests\SearchResults\TDGPTMDSearchResults.psmtsv")]
        [TestCase(@"FileReadingTests\SearchResults\XLink.psmtsv")]
        [TestCase(@"FileReadingTests\SearchResults\oglyco.psmtsv")]
        [TestCase(@"FileReadingTests\SearchResults\nglyco_f5.psmtsv")]
        [TestCase(@"Transcriptomics\TestData\OsmFileForTesting.osmtsv")]
        public void ProForma_ColumnAbsent_NeverThrows(string relativePath)
        {
            // Top-down, crosslink, glyco and RNA full sequences are not all expressible as flat ProForma;
            // the getter must answer null for those rather than throw from a property access.
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, relativePath);
            var matches = SpectrumMatchTsvReader.ReadTsv<SpectrumMatchFromTsv>(path, out _);
            Assert.That(matches, Is.Not.Empty);
            Assert.That(() => matches.Select(m => m.ProForma).ToList(), Throws.Nothing);
        }
    }
}
