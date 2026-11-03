using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Readers.ProForma;
using Tdp = TopDownProteomics.ProForma;

namespace Test.FileReadingTests.ProForma
{
    /// <summary>
    /// Phase 1-3: Layer-1 round-trip (string &lt;-&gt; ProFormaTerm), delegated to the wrapped
    /// TopDownProteomics parser/writer. Verifies parsing succeeds and that the writer is
    /// canonically idempotent — a practical proxy for AST equivalence:
    /// write(parse(s)) == write(parse(write(parse(s)))).
    /// Covers every valid corpus example except the documented exclusion sets below.
    /// The SDK fully covers a single proteoform term; the remaining gaps are documented in
    /// deliverables/known-limitations.md.
    /// </summary>
    [TestFixture]
    internal class ProFormaLayer1RoundTripTests
    {
        /// <summary>
        /// v1.0-only notations the v2-only SDK does not parse. Intentionally unsupported
        /// (v2.0 removed them; user decision 2026-05-24). Asserted by
        /// <see cref="V1_DefaultSourcePrefix_IsUnsupported"/>.
        /// </summary>
        private static readonly HashSet<string> KnownUnsupportedV1 = new()
        {
            // v1.0 Rule 6 / best-practice ii: "default-source" prefix [SOURCE]+SEQ, removed in v2.0
            "v1-rule6-01", "v1-rule6-02", "v1-rule6-03", "v1-bp-ii",
        };

        /// <summary>
        /// Valid v2.0 features that live above a single ProFormaTerm and so are not handled by the
        /// SDK's single-term ParseString: charge/adducts (/z), chimeric (+), and multi-chain (//).
        /// Supporting these requires facade-level splitting/charge handling (planned, README
        /// Phases 7-8). Pinned by <see cref="MultiTerm_ChargeChimericMultichain_NotYetSupported"/>.
        /// </summary>
        private static readonly HashSet<string> RequiresMultiTermFacade = new()
        {
            "v2-4.2.3.2-01", "v2-4.2.3.2-02", "v2-4.2.3.3-03", // inter-chain // (crosslink)
            "v2-4.2.4-01", "v2-4.2.4-02",                       // branch //
            "v2-7.1-01", "v2-7.1-02", "v2-7.1-03", "v2-7.1-04", // charge /z [+adducts]
            "v2-7.1-05", "v2-7.1-06", "v2-7.1-07",
            "v2-7.2-01",                                         // chimeric +
        };

        private static IEnumerable<TestCaseData> RoundTripExamples()
        {
            foreach (var r in ProFormaTestCorpus.Load()
                         .Where(r => r.Valid
                                     && !KnownUnsupportedV1.Contains(r.Id)
                                     && !RequiresMultiTermFacade.Contains(r.Id)))
                yield return new TestCaseData(r.ProformaString)
                    .SetName($"{r.ComplianceLevel}_{ProFormaTestCorpus.ToTestName(r.Id)}");
        }

        [TestCaseSource(nameof(RoundTripExamples))]
        public void RoundTrip_IsCanonicallyIdempotent(string proForma)
        {
            string firstWrite = ProFormaWriter.Write(ProFormaReader.Read(proForma));
            string secondWrite = ProFormaWriter.Write(ProFormaReader.Read(firstWrite));
            Assert.That(secondWrite, Is.EqualTo(firstWrite),
                $"non-idempotent canonical form for '{proForma}' (first write = '{firstWrite}')");
        }

        /// <summary>
        /// Pins gap #2: the legacy v1.0 default-source prefix `[SOURCE]+SEQ` (dropped in v2.0) is
        /// not parsed by the v2-only SDK. Intentionally unsupported.
        /// </summary>
        [Test]
        public void V1_DefaultSourcePrefix_IsUnsupported()
        {
            var rows = ProFormaTestCorpus.Load().Where(r => KnownUnsupportedV1.Contains(r.Id)).ToList();
            Assert.That(rows, Has.Count.EqualTo(KnownUnsupportedV1.Count), "missing v1 default-source rows");
            foreach (var r in rows)
                Assert.Throws<Tdp.ProFormaParseException>(() => ProFormaReader.Read(r.ProformaString),
                    $"{r.Id} now parses — v1 default-source may be supported; revisit KnownUnsupportedV1.");
        }

        /// <summary>
        /// Pins gap #3: charge (/z), chimeric (+), and multi-chain (//) constructs are not handled
        /// by the SDK's single-term ParseString (it throws "/ is not an upper case letter"). When
        /// the multi-term facade is built (Phases 7-8), these move out of <see cref="RequiresMultiTermFacade"/>.
        /// </summary>
        [Test]
        public void MultiTerm_ChargeChimericMultichain_NotYetSupported()
        {
            var rows = ProFormaTestCorpus.Load().Where(r => RequiresMultiTermFacade.Contains(r.Id)).ToList();
            Assert.That(rows, Has.Count.EqualTo(RequiresMultiTermFacade.Count), "missing multi-term rows");
            foreach (var r in rows)
                Assert.Throws<Tdp.ProFormaParseException>(() => ProFormaReader.Read(r.ProformaString),
                    $"{r.Id} now parses via single-term Read — build the multi-term facade and move it out.");
        }

        /// <summary>
        /// ProForma 2.1 (section 6.3) allows several modifications on one terminus. Since TopDownProteomics 0.0.301 the term
        /// holds them all, in the order written (<see cref="Tdp.ProFormaTerm.NTerminalModifications"/>,
        /// <see cref="Tdp.ProFormaTerm.CTerminalModifications"/>), so the reader returns them. 0.0.299 kept only the last
        /// C-terminal one and read stacked N-terminal ones as unlocalized, with no error. Examples from the 2.1
        /// specification and the HUPO-PSI grammar tests.
        /// </summary>
        [TestCase("PEPTIDE-[UNIMOD:2][UNIMOD:35]", "C", "UNIMOD:2|UNIMOD:35")]
        [TestCase("PEPTIDEG-[Methyl][Amidated]", "C", "Methyl|Amidated")]
        [TestCase("PEPTIDEG-[Methyl][Amidated][INFO:A lot of C terminal mods]", "C", "Methyl|Amidated|A lot of C terminal mods")]
        [TestCase("[UNIMOD:1][UNIMOD:35]-PEPTIDE", "N", "UNIMOD:1|UNIMOD:35")]
        [TestCase("[Acetyl][Carbamyl]-QPEPTIDE", "N", "Acetyl|Carbamyl")]
        [TestCase("[Acetyl][Acetyl][Carbamyl]-QPEPTIDE", "N", "Acetyl|Acetyl|Carbamyl")]
        // After an unlocalized prefix: the unlocalized groups must stay unlocalized and the terminal ones terminal
        [TestCase("[Phospho]?[Acetyl]-PEPTIDE-[UNIMOD:2][UNIMOD:35]", "C", "UNIMOD:2|UNIMOD:35")]
        [TestCase("[Phospho]?[Acetyl][Carbamyl]-QPEPTIDE", "N", "Acetyl|Carbamyl")]
        [TestCase("[Phospho][Phospho]?[Acetyl][Carbamyl]-QPEPTIDE", "N", "Acetyl|Carbamyl")]
        [TestCase("[Phospho]^2?[Acetyl][Carbamyl]-QPEPTIDE", "N", "Acetyl|Carbamyl")]
        [TestCase("[Phospho]?[Methyl]^2?[Acetyl][Acetyl][Carbamyl]-QPEPTIDE", "N", "Acetyl|Acetyl|Carbamyl")]
        [TestCase("[Phospho]?{Glycan:Hex}[Acetyl][Carbamyl]-QPEPTIDE", "N", "Acetyl|Carbamyl")]
        [TestCase("{Glycan:Hex}[Phospho]?PEPTIDEG-[Methyl][Amidated]", "C", "Methyl|Amidated")]
        public void Read_StackedTerminalModifications_KeepEveryModification(string proForma, string terminus, string expected)
        {
            var term = ProFormaReader.Read(proForma);
            Assert.That(TerminalValues(term, terminus), Is.EqualTo(expected));

            var unlocalized = term.UnlocalizedTags?.SelectMany(t => t.Descriptors).Select(d => d.Value).ToList() ?? new List<string>();
            Assert.That(unlocalized, Has.None.AnyOf("Acetyl", "Carbamyl", "UNIMOD:1", "UNIMOD:35"),
                "a terminal modification was read as unlocalized");

            var reread = ProFormaReader.Read(ProFormaWriter.Write(term));
            Assert.That(TerminalValues(reread, terminus), Is.EqualTo(expected), "the written string reads back the same");
        }

        /// <summary>Several labile groups (ProForma 2.0, section 4.3.2) are each kept; 0.0.299 kept only the last.</summary>
        [Test]
        public void Read_SeveralLabileModifications_KeepEveryOne()
        {
            var term = ProFormaReader.Read("{Glycan:Hex}{Glycan:NeuAc}EMEVNESPEK");

            Assert.That(term.LabileModifications!.Select(m => m.Single().Value), Is.EqualTo(new[] { "Hex", "NeuAc" }));
            Assert.That(ProFormaWriter.Write(term), Is.EqualTo("{Glycan:Hex}{Glycan:NeuAc}EMEVNESPEK"));
        }

        private static string TerminalValues(Tdp.ProFormaTerm term, string terminus) => string.Join("|",
            (terminus == "N" ? term.NTerminalModifications : term.CTerminalModifications)!.Select(m => m.Single().Value));

        /// <summary>No residue may follow the C-terminal modification (ProForma 2.0 section 4.3.1). 0.0.299 read
        /// "PEPTIDE-[UNIMOD:2]K" as PEPTIDEK; 0.0.301 refuses it.</summary>
        [TestCase("PEPTIDE-[UNIMOD:2]K", 18)]
        [TestCase("[Phospho]?[Acetyl]-PEPTIDE-[UNIMOD:2]K", 37)]
        [TestCase("[Phospho]^2?PEPTIDE-[UNIMOD:2]K", 30)]
        public void Read_ResidueAfterCTerminalModification_IsNotValidProForma(string proForma, int position)
        {
            var ex = Assert.Throws<Tdp.ProFormaParseException>(() => ProFormaReader.Read(proForma));
            Assert.That(ex!.Message, Does.Contain($"at index {position} follows the C-terminal '-'"));
        }

        /// <summary>A string that is only modifications has no residue: a parse failure, never an index exception and
        /// never a term with an empty sequence.</summary>
        [TestCase("[UNIMOD:1]")]
        [TestCase("[UNIMOD:1][UNIMOD:35]")]
        [TestCase("[UNIMOD:1]-")]
        [TestCase("[Phospho]?")]
        [TestCase("[Phospho]?[Acetyl]")]
        [TestCase("[Phospho]?[Acetyl]-")]
        [TestCase("[Phospho]^2")]
        [TestCase("-[UNIMOD:2]")]
        public void Read_ModificationWithoutAResidue_IsAParseFailure(string proForma)
        {
            Assert.That(() => ProFormaReader.Read(proForma), Throws.TypeOf<Tdp.ProFormaParseException>());
        }

        /// <summary>The C-terminal '-' needs a residue before it. 0.0.299 read "-PEPTIDE" as PEPTIDE; 0.0.301 refuses it.</summary>
        [Test]
        public void Read_LeadingHyphen_IsNotValidProForma()
        {
            var ex = Assert.Throws<Tdp.ProFormaParseException>(() => ProFormaReader.Read("-PEPTIDE"));
            Assert.That(ex!.Message, Does.Contain("follows the C-terminal '-'"));
        }

        /// <summary>One modification per terminus, residue stacks, unlocalized and labile groups, nested brackets inside
        /// a modification, and a global fixed modification all still parse as before.</summary>
        [TestCase("[UNIMOD:1]-PEPTIDE-[UNIMOD:2]")]
        [TestCase("PEPC[UNIMOD:4][UNIMOD:35]TIDE")]
        [TestCase("[Phospho]?EMEVTSESPEK")]
        [TestCase("[Phospho][Phospho]?EMEVTSESPEK")]
        [TestCase("[Phospho]^2?EMEVTSESPEK")]
        [TestCase("[Phospho][Phospho]?[Acetyl]-EM[Oxidation]EVTSESPEK")]
        [TestCase("[Phospho]?[Acetyl]-PEPTIDE-[UNIMOD:2]")]
        [TestCase("<[MOD:01090]@C>[Phospho]?EM[Oxidation]EVTSECSPEK")]
        [TestCase("[Phospho]?{Glycan:Hex}[Acetyl]-PEPTIDE")]
        [TestCase("{Glycan:Hex}EM[Oxidation]EVNES[Phospho]PEK")]
        [TestCase("{Glycan:Hex}[iTRAQ4plex]-EM[Oxidation]EVNES[Phospho]PEK[iTRAQ4plex]-[Methyl]")]
        [TestCase("<[S-carboxamidomethyl-L-cysteine]@C>ATPEILTCNSIGCLK")]
        [TestCase("SEQUEN[Formula:[13C2][12C-2]H2N]CE")]
        [TestCase("PEPTIDE-[Formula:[13C2]C-2H2N]")]
        [TestCase("EM[-17.026549]EVEES[+79.966331]PEK")]
        public void Read_SupportedShapes_StillParse(string proForma)
        {
            Assert.That(() => ProFormaReader.Read(proForma), Throws.Nothing);
        }

        /// <summary>Empty brackets are still refused by the SDK, with its own message (pinned elsewhere).</summary>
        [TestCase("[]-PEPTIDE")]
        [TestCase("PEPTIDE-[]")]
        public void Read_EmptyTerminalBrackets_AreLeftToTheSdk(string proForma)
        {
            var ex = Assert.Throws<Tdp.ProFormaParseException>(() => ProFormaReader.Read(proForma));
            Assert.That(ex, Is.Not.InstanceOf<ProFormaUnsupportedException>());
        }
    }
}
