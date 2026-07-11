using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests.VariantCorpus
{
    /// <summary>
    /// Data-driven sequence-variant corpus. Each <see cref="CorpusCase"/> is a node in the possibility tree:
    /// a tiny protein spec (base sequence + mods + variants) paired with the proteoforms it MUST digest to,
    /// hand-derived from the domain "bible" (DatabaseTests/VariantCorpus/README.md) — never captured from the
    /// code. The cases live inline as records (not an external table) so the columns are compiler-checked and
    /// the input sits next to its expected answer, matching the house style in Transcriptomics/TestDigestion.cs.
    ///
    /// Two families run over the same oracle:
    ///   Corpus_Node           — build XML -> LoadProteinXML (variants expand at load) -> digest -> compare.
    ///   Corpus_Node_RoundTrip — same, but the consensus protein is serialized through the real ProteinDbWriter
    ///                           and re-read first, so encode+decode is exercised too (invariant 6 / #1083).
    /// L0 only (substitutions + PTMs, top-down, caps non-binding); indels/processing/digestion are later layers.
    /// </summary>
    [TestFixture]
    internal class VariantCorpusTests
    {
        /// <summary>
        /// One corpus node. <see cref="ExpectedForms"/> is the oracle (the exact proteoforms, in canonical
        /// least-modified-first order, in mzLib full-sequence notation). <see cref="ExpectedCount"/> is kept
        /// alongside as an independent over/under-generation guard (mirrors RnaDigestionTestCase, which keeps
        /// both DigestionProductCount and Sequences). Verdict/Reason document WHY, tied to a bible installment.
        /// </summary>
        internal record CorpusCase(
            string Id, string Layer, string Tests,
            string Base, string Mods, string Variants, string Protease,
            int MaxIsoforms, int MaxMods,
            int ExpectedCount, string Verdict, string Reason,
            string[] ExpectedForms,
            string Opts = "-");

        // Canonical test-mod registry (see README "Canonical test mods"). name -> (ModificationType, monoisotopicMass).
        private static readonly Dictionary<string, (string Type, double Mass)> ModRegistry = new()
        {
            ["Phosphorylation"] = ("Biological", 79.966331),
        };

        /// <summary>
        /// The corpus. L0 foundation: substitutions + PTMs only (encoding-unambiguous). Grow in layers per the
        /// README; simple nodes must keep passing as complexity is added.
        /// </summary>
        internal static IEnumerable<CorpusCase> GetCases()
        {
            // F00 — the floor: no variant, no mod -> one consensus proteoform. Caps non-binding.
            yield return new CorpusCase(
                Id: "F00", Layer: "L0a", Tests: "baseline-identity",
                Base: "PEPTIDE", Mods: "-", Variants: "-", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 1, Verdict: "applied",
                Reason: "No variant, no mod -> one consensus proteoform. Caps non-binding. The floor.",
                ExpectedForms: new[] { "PEPTIDE" });

            // F01 — variable PTM -> unmodified AND modified (0..MaxMods). 1 mod avail, cap 2 -> 2 forms.
            yield return new CorpusCase(
                Id: "F01", Layer: "L0b", Tests: "mod-only-variable",
                Base: "PEPTIDE", Mods: "Phosphorylation@4", Variants: "-", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 2, Verdict: "applied",
                Reason: "Variable PTM -> unmodified AND modified (0..MaxMods). 1 mod avail, cap 2 -> 2 forms.",
                ExpectedForms: new[] { "PEPTIDE", "PEPT[Biological:Phosphorylation on T]IDE" });

            // F02 — substitution -> consensus + variant sequence (installment 2).
            yield return new CorpusCase(
                Id: "F02", Layer: "L0a", Tests: "sub-only",
                Base: "PEPTIDE", Mods: "-", Variants: "OP=T VAR=V POS=4 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 2, Verdict: "applied",
                Reason: "Substitution -> consensus + variant sequence (installment 2).",
                ExpectedForms: new[] { "PEPTIDE", "PEPVIDE" });

            // F03 — sub x mod interaction: T->V drops the phospho at the edited position; no phosphorylated
            // variant produced. 3 forms NOT 4 (installment 2; invariants 1,2).
            yield return new CorpusCase(
                Id: "F03", Layer: "L0-int", Tests: "sub-x-mod-legal-consensus-only",
                Base: "PEPTIDE", Mods: "Phosphorylation@4", Variants: "OP=T VAR=V POS=4 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 3, Verdict: "pruned:illegal-mod",
                Reason: "T->V drops the phospho at the edited position; no phosphorylated variant produced. 3 forms NOT 4 (installment 2; invariants 1,2).",
                ExpectedForms: new[] { "PEPTIDE", "PEPT[Biological:Phosphorylation on T]IDE", "PEPVIDE" });

            // F04 — boundary: position 9 > length 7 -> variant pruned; only consensus (reader ~line 776;
            // outOfRangeVariant.xml).
            yield return new CorpusCase(
                Id: "F04", Layer: "L0a", Tests: "sub-out-of-range",
                Base: "PEPTIDE", Mods: "-", Variants: "OP=- VAR=A POS=9 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 1, Verdict: "pruned:out-of-range",
                Reason: "Position 9 > length 7 -> variant pruned; only consensus (reader ~line 776; outOfRangeVariant.xml).",
                ExpectedForms: new[] { "PEPTIDE" });
        }

        private static IEnumerable<TestCaseData> Cases()
        {
            foreach (var c in GetCases())
                yield return new TestCaseData(c).SetName($"Corpus_{c.Id}_{c.Tests}");
        }

        [TestCaseSource(nameof(Cases))]
        public void Corpus_Node(CorpusCase c)
        {
            AssertCaseSelfConsistent(c);

            List<string> produced = RunNode(c);
            List<string> producedAgain = RunNode(c);

            // Determinism (invariant 11): identical output, same order, across runs.
            Assert.That(producedAgain, Is.EqualTo(produced),
                $"{c.Id}: nondeterministic output across two runs.");

            // Count (fast over/under-generation check).
            Assert.That(produced.Count, Is.EqualTo(c.ExpectedCount),
                $"{c.Id} ({c.Tests}): expected {c.ExpectedCount} forms, got {produced.Count}.\n  got: {string.Join(" | ", produced)}");

            // Membership (order-independent for now; canonical-order assertion is a later tightening).
            Assert.That(produced.OrderBy(x => x, StringComparer.Ordinal),
                Is.EqualTo(c.ExpectedForms.OrderBy(x => x, StringComparer.Ordinal)),
                $"{c.Id} ({c.Tests}): form set mismatch.\n  expected: {string.Join(" | ", c.ExpectedForms)}\n  got:      {string.Join(" | ", produced)}");
        }

        /// <summary>
        /// Round-trip variant: same oracle, but the CONSENSUS protein is serialized through the REAL
        /// ProteinDbWriter and re-read before digestion. This exercises encode+decode (invariant 6 / #1083)
        /// WITHOUT the writer ever sitting upstream of the ground truth — the input XML is still the
        /// hand-authored, bible-grounded one; the writer only has to preserve what we already loaded.
        /// A writer that drops a variant or a variant-adjacent mod fails here while Corpus_Node stays green.
        /// </summary>
        [TestCaseSource(nameof(Cases))]
        public void Corpus_Node_RoundTrip(CorpusCase c)
        {
            AssertCaseSelfConsistent(c);

            List<string> produced = RunNodeRoundTrip(c);

            Assert.That(produced.Count, Is.EqualTo(c.ExpectedCount),
                $"{c.Id} ({c.Tests}) [round-trip]: expected {c.ExpectedCount} forms, got {produced.Count}.\n  got: {string.Join(" | ", produced)}");

            Assert.That(produced.OrderBy(x => x, StringComparer.Ordinal),
                Is.EqualTo(c.ExpectedForms.OrderBy(x => x, StringComparer.Ordinal)),
                $"{c.Id} ({c.Tests}) [round-trip]: form set mismatch after write/read cycle.\n  expected: {string.Join(" | ", c.ExpectedForms)}\n  got:      {string.Join(" | ", produced)}");
        }

        // The case's own two statements of size must agree — catches a typo where ExpectedForms and
        // ExpectedCount drift (mirrors the redundancy the TSV used to enforce across cases.tsv + expected/).
        private static void AssertCaseSelfConsistent(CorpusCase c) =>
            Assert.That(c.ExpectedForms.Length, Is.EqualTo(c.ExpectedCount),
                $"{c.Id}: case defines ExpectedCount={c.ExpectedCount} but lists {c.ExpectedForms.Length} ExpectedForms.");

        private static List<string> RunNodeRoundTrip(CorpusCase c)
        {
            List<Modification> knownMods = BuildMods(c);
            string xml = BuildXml(c);
            string tmp1 = Path.Combine(Path.GetTempPath(), $"corpus_{c.Id}_{Guid.NewGuid():N}.xml");
            string tmp2 = Path.Combine(Path.GetTempPath(), $"corpus_{c.Id}_rt_{Guid.NewGuid():N}.xml");
            try
            {
                // First load: the hand-authored (bible-grounded) XML -> consensus + expanded variants.
                File.WriteAllText(tmp1, xml);
                List<Protein> loaded = ProteinDbLoader.LoadProteinXML(tmp1,
                    generateTargets: true, decoyType: DecoyType.None,
                    allKnownModifications: knownMods, isContaminant: false,
                    modTypesToExclude: new List<string>(), out _);

                // Write ONLY the consensus protein(s): they retain SequenceVariations, so re-expansion
                // happens on reload exactly as it did on the first load. Writing the applied isoforms
                // instead would double-count and pin the expansion to write time.
                List<Protein> consensus = loaded.Where(p => !p.IsDecoy && !p.AppliedSequenceVariations.Any()).ToList();
                ProteinDbWriter.WriteXmlDatabase(
                    new Dictionary<string, HashSet<Tuple<int, Modification>>>(), consensus, tmp2);

                List<Protein> reread = ProteinDbLoader.LoadProteinXML(tmp2,
                    generateTargets: true, decoyType: DecoyType.None,
                    allKnownModifications: knownMods, isContaminant: false,
                    modTypesToExclude: new List<string>(), out _);

                return Digest(reread, c);
            }
            finally { File.Delete(tmp1); File.Delete(tmp2); }
        }

        private static List<string> RunNode(CorpusCase c)
        {
            List<Modification> knownMods = BuildMods(c);
            string xml = BuildXml(c);
            string tmp = Path.Combine(Path.GetTempPath(), $"corpus_{c.Id}_{Guid.NewGuid():N}.xml");
            try
            {
                File.WriteAllText(tmp, xml);
                List<Protein> proteins = ProteinDbLoader.LoadProteinXML(tmp,
                    generateTargets: true, decoyType: DecoyType.None,
                    allKnownModifications: knownMods, isContaminant: false,
                    modTypesToExclude: new List<string>(), out _);
                return Digest(proteins, c);
            }
            finally { File.Delete(tmp); }
        }

        private static List<string> Digest(List<Protein> proteins, CorpusCase c)
        {
            var dp = new DigestionParams(protease: c.Protease, maxMissedCleavages: 0, minPeptideLength: 1,
                maxModificationIsoforms: c.MaxIsoforms, maxModsForPeptides: c.MaxMods,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            return proteins
                .SelectMany(p => p.Digest(dp, new List<Modification>(), new List<Modification>()))
                .Select(pw => pw.FullSequence)
                .ToList();
        }

        // --- helpers: parse the case's mods/variants DSL and emit the matching mzLib XML ------------------

        private static List<Modification> BuildMods(CorpusCase c)
        {
            var mods = new List<Modification>();
            if (c.Mods == "-") return mods;
            foreach (var token in c.Mods.Split(';', StringSplitOptions.RemoveEmptyEntries))
            {
                var (name, pos) = ParseNameAtPos(token.Trim());
                char motifChar = c.Base[pos - 1];
                ModificationMotif.TryGetMotif(motifChar.ToString(), out ModificationMotif motif);
                var (type, mass) = ModRegistry[name];
                mods.Add(new Modification(_originalId: name, _modificationType: type, _target: motif,
                    _locationRestriction: "Anywhere.", _monoisotopicMass: mass));
            }
            return mods;
        }

        private static string BuildXml(CorpusCase c)
        {
            var features = new System.Text.StringBuilder();

            if (c.Mods != "-")
            {
                foreach (var token in c.Mods.Split(';', StringSplitOptions.RemoveEmptyEntries))
                {
                    var (name, pos) = ParseNameAtPos(token.Trim());
                    char motifChar = c.Base[pos - 1];
                    // The reader binds a "modified residue" feature by description == Modification.IdWithMotif.
                    string idWithMotif = $"{name} on {motifChar}";
                    features.Append(
                        $"    <feature type=\"modified residue\" description=\"{idWithMotif}\">\n" +
                        $"      <location><position position=\"{pos}\" /></location>\n" +
                        $"    </feature>\n");
                }
            }

            if (c.Variants != "-")
            {
                foreach (var v in c.Variants.Split(';', StringSplitOptions.RemoveEmptyEntries))
                {
                    var kv = ParseKeyVals(v.Trim());
                    string op = kv.TryGetValue("OP", out var o) && o != "-" ? o : "X";
                    string variation = kv["VAR"];
                    string pos = kv["POS"];
                    features.Append(
                        $"    <feature type=\"sequence variant\" description=\"{op}{pos}{variation} variant\">\n" +
                        $"      <original>{op}</original>\n      <variation>{variation}</variation>\n" +
                        $"      <location><position position=\"{pos}\" /></location>\n" +
                        $"    </feature>\n");
                }
            }

            return
                "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n" +
                "<mzLibProteinDb>\n" +
                "  <entry>\n" +
                $"    <accession>P_{c.Id}</accession>\n" +
                $"    <name>TEST_{c.Id}</name>\n" +
                "    <protein><recommendedName><fullName>Test protein</fullName></recommendedName></protein>\n" +
                features +
                $"    <sequence length=\"{c.Base.Length}\">{c.Base}</sequence>\n" +
                "  </entry>\n" +
                "</mzLibProteinDb>\n";
        }

        private static (string name, int pos) ParseNameAtPos(string token)
        {
            var parts = token.Split('@');
            return (parts[0], int.Parse(parts[1]));
        }

        private static Dictionary<string, string> ParseKeyVals(string s)
        {
            var d = new Dictionary<string, string>();
            foreach (var kv in s.Split(' ', StringSplitOptions.RemoveEmptyEntries))
            {
                var i = kv.IndexOf('=');
                if (i > 0) d[kv.Substring(0, i)] = kv.Substring(i + 1);
            }
            return d;
        }
    }
}
