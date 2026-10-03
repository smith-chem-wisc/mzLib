using System.Collections.Generic;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Omics.Modifications;
using Omics.SequenceConversion;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers;
using Readers.ProForma;

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
        // A UniProt modification carries its catalog entry's Unimod cross-reference.
        [TestCase("[UniProt:N-acetylalanine on A]AAAGEAR", "[UNIMOD:1]-AAAGEAR",
            TestName = "FullSequence_UniProtNTerminalMod_BecomesUnimodAccession")]
        // MetaMorpheus writes a C-terminal modification after a '-', which is a terminus marker, not a residue.
        // mzLib's catalogs define this name only as a UniProt modification, which cites UNIMOD:34.
        [TestCase("KPVADYFL-[Common Artifact:Leucine methyl ester on L]", "KPVADYFL-[UNIMOD:34]",
            TestName = "FullSequence_CTerminalMod_StaysOnTheCTerminus")]
        [TestCase("PEPK[Custom:Nameless on K]R", "PEPK[Custom:Nameless on K]R", TestName = "FullSequence_UnknownMod_KeepsName")]
        public void ProFormaFromFullSequence_ConvertsMetaMorpheusNotation(string fullSequence, string expected)
        {
            Assert.That(SpectrumMatchFromTsv.ProFormaFromFullSequence(fullSequence), Is.EqualTo(expected));
        }

        private static PeptideWithSetModifications DigestWithUniProtMod(string sequence, string idWithMotif, int position)
        {
            var mod = Mods.UniprotModifications.Single(m => m.IdWithMotif == idWithMotif);
            var localized = new Dictionary<int, List<Modification>> { [position] = new() { mod } };
            return new Protein(sequence, "P", oneBasedModifications: localized)
                .Digest(new DigestionParams(protease: "trypsin", maxMissedCleavages: 0, minPeptideLength: 1), new List<Modification>(), new List<Modification>())
                .Cast<PeptideWithSetModifications>()
                .First(p => p.AllModsOneIsNterminus.Count == 1);
        }

        private static double ReadBackMass(string proForma) =>
            ProFormaConverter.ToModificationDictionary(ProFormaReader.Read(proForma), Mods.AllKnownProteinModsDictionary)
                .Values.Single().MonoisotopicMass!.Value;

        [Test]
        public void ProForma_PsiModReference_IsWrittenWithOnePrefixAndReadsBack()
        {
            // UniProt's N6-crotonyllysine cites PSI-MOD (stored as "MOD:01892") and no UNIMOD record.
            var peptide = DigestWithUniProtMod("PEPKIDEK", "N6-crotonyllysine on K", 4);
            var mass = peptide.AllModsOneIsNterminus.Values.Single().MonoisotopicMass!.Value;
            Assert.That(peptide.FullSequence, Is.EqualTo("PEPK[UniProt:N6-crotonyllysine on K]"));

            foreach (var proForma in new[] { SpectrumMatchFromTsv.ProFormaFromFullSequence(peptide.FullSequence)!, peptide.ToProFormaString() })
            {
                Assert.That(proForma, Is.EqualTo("PEPK[MOD:01892]"));
                Assert.That(ReadBackMass(proForma), Is.EqualTo(mass).Within(1e-6));
            }

            // The accession with its prefix doubled reads as the same modification.
            Assert.That(ReadBackMass("PEPK[MOD:MOD:01892]"), Is.EqualTo(mass).Within(1e-6));
        }

        [Test]
        public void ProForma_UnimodReferenceOfAnotherMass_IsNotWrittenOrParsedAsUnimodId()
        {
            // UniProt's N,N-dimethylproline (+28.031) cites UNIMOD:529, whose record is +29.039.
            var peptide = DigestWithUniProtMod("PEPTIDEK", "N,N-dimethylproline on P", 1);
            var dimethylproline = peptide.AllModsOneIsNterminus.Values.Single();
            Assert.That(peptide.FullSequence, Is.EqualTo("[UniProt:N,N-dimethylproline on P]PEPTIDEK"));

            var parsed = MzLibSequenceParser.Instance.Parse(peptide.FullSequence)!.Value.Modifications.Single();
            Assert.That(parsed.MzLibModification, Is.SameAs(dimethylproline));
            Assert.That(parsed.UnimodId, Is.Null);

            foreach (var proForma in new[] { SpectrumMatchFromTsv.ProFormaFromFullSequence(peptide.FullSequence)!, peptide.ToProFormaString() })
            {
                Assert.That(proForma, Does.Not.Contain("UNIMOD:529"));
                Assert.That(ReadBackMass(proForma), Is.EqualTo(dimethylproline.MonoisotopicMass!.Value).Within(1e-6), proForma);
            }
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
