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
        /// ProForma 2.1 (section 6.3) allows several modifications on one terminus, but a ProFormaTerm holds one per
        /// terminus: the SDK kept only the last C-terminal one and misread stacked N-terminal ones as unlocalized. The
        /// reader refuses them rather than lose one. Examples from the 2.1 specification and the HUPO-PSI grammar tests.
        /// </summary>
        [TestCase("PEPTIDE-[UNIMOD:2][UNIMOD:35]", "C-terminus", 2)]
        [TestCase("PEPTIDEG-[Methyl][Amidated]", "C-terminus", 2)]
        [TestCase("PEPTIDEG-[Methyl][Amidated][INFO:A lot of C terminal mods]", "C-terminus", 3)]
        [TestCase("[UNIMOD:1][UNIMOD:35]-PEPTIDE", "N-terminus", 2)]
        [TestCase("[Acetyl][Carbamyl]-QPEPTIDE", "N-terminus", 2)]
        [TestCase("[Acetyl][Acetyl][Carbamyl]-QPEPTIDE", "N-terminus", 3)]
        public void Read_StackedTerminalModifications_AreRefusedNotLost(string proForma, string terminus, int count)
        {
            var ex = Assert.Throws<ProFormaUnsupportedException>(() => ProFormaReader.Read(proForma));
            Assert.That(ex!.Message, Does.Contain($"{count} modifications on the {terminus}").And.Contain("ProForma 2.1"));
            Assert.That(ex.IncompatibleItem,
                Is.EqualTo($"Multiple {(terminus == "N-terminus" ? "N" : "C")}-terminal modifications ({count} found)."));
            Assert.That(ex, Is.InstanceOf<Tdp.ProFormaParseException>(), "callers that catch parse failures still catch it");
        }

        /// <summary>No residue may follow the C-terminal modification (ProForma 2.0 section 4.3.1); the SDK read
        /// "PEPTIDE-[UNIMOD:2]K" as PEPTIDEK.</summary>
        [Test]
        public void Read_ResidueAfterCTerminalModification_IsNotValidProForma()
        {
            var ex = Assert.Throws<Tdp.ProFormaParseException>(() => ProFormaReader.Read("PEPTIDE-[UNIMOD:2]K"));
            Assert.That(ex!.Message, Is.EqualTo("Unexpected content at position 18 after the C-terminal modification."));
        }

        /// <summary>A string that is only a modification has no residue: a parse failure, never an index exception.</summary>
        [TestCase("[UNIMOD:1]")]
        [TestCase("[UNIMOD:1][UNIMOD:35]")]
        public void Read_ModificationWithoutAResidue_IsAParseFailure(string proForma)
        {
            Assert.That(() => ProFormaReader.Read(proForma), Throws.TypeOf<Tdp.ProFormaParseException>());
        }

        /// <summary>One modification per terminus, residue stacks, unlocalized and labile groups, nested brackets inside
        /// a modification, and a global fixed modification all still parse as before.</summary>
        [TestCase("[UNIMOD:1]-PEPTIDE-[UNIMOD:2]")]
        [TestCase("PEPC[UNIMOD:4][UNIMOD:35]TIDE")]
        [TestCase("[Phospho]?EMEVTSESPEK")]
        [TestCase("[Phospho][Phospho]?EMEVTSESPEK")]
        [TestCase("{Glycan:Hex}EM[Oxidation]EVNES[Phospho]PEK")]
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
