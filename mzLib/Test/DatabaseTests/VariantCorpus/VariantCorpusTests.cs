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
    /// Data-driven sequence-variant corpus. Each row of cases.tsv is a node in the possibility tree:
    /// build a UniProt-style XML from the row, LoadProteinXML (which expands variants), digest each expanded
    /// protein, and compare the produced FullSequence set to the node's expected/&lt;id&gt;.txt.
    /// The expected files are the ORACLE (hand-derived from the domain "bible", never from the code).
    /// See DatabaseTests/VariantCorpus/README.md for the schema, notation, and the determinism contract.
    /// L0 only (substitutions + PTMs, top-down, caps non-binding); indels/processing/digestion are later layers.
    /// </summary>
    [TestFixture]
    internal class VariantCorpusTests
    {
        // Canonical test-mod registry (see README "Canonical test mods"). name -> (ModificationType, monoisotopicMass).
        private static readonly Dictionary<string, (string Type, double Mass)> ModRegistry = new()
        {
            ["Phosphorylation"] = ("Biological", 79.966331),
        };

        private static string CorpusDir => Path.Combine(TestContext.CurrentContext.TestDirectory,
            "DatabaseTests", "VariantCorpus");

        internal record Row(string Id, string Layer, string Tests, string Base, string Mods, string Variants,
            string Protease, int MaxIsoforms, int MaxMods, string Opts, int ExpectedCount, string Verdict, string Reason);

        private static IEnumerable<TestCaseData> Cases()
        {
            foreach (var r in LoadRows())
                yield return new TestCaseData(r).SetName($"Corpus_{r.Id}_{r.Tests}");
        }

        private static List<Row> LoadRows()
        {
            var rows = new List<Row>();
            foreach (var line in File.ReadAllLines(Path.Combine(CorpusDir, "cases.tsv")).Skip(1))
            {
                if (line.Trim().Length == 0) continue;
                var c = line.Split('\t');
                rows.Add(new Row(c[0], c[1], c[2], c[3], c[4], c[5], c[6],
                    int.Parse(c[7]), int.Parse(c[8]), c[9], int.Parse(c[10]), c[11], c[12]));
            }
            return rows;
        }

        [TestCaseSource(nameof(Cases))]
        public void Corpus_Node(Row row)
        {
            List<string> expected = File.ReadAllLines(Path.Combine(CorpusDir, "expected", row.Id + ".txt"))
                .Where(l => l.Trim().Length > 0).ToList();

            List<string> produced = RunNode(row);
            List<string> producedAgain = RunNode(row);

            // Determinism (invariant 11): identical output, same order, across runs.
            Assert.That(producedAgain, Is.EqualTo(produced),
                $"{row.Id}: nondeterministic output across two runs.");

            // Count (fast over/under-generation check).
            Assert.That(produced.Count, Is.EqualTo(row.ExpectedCount),
                $"{row.Id} ({row.Tests}): expected {row.ExpectedCount} forms, got {produced.Count}.\n  got: {string.Join(" | ", produced)}");

            // Membership (order-independent for now; canonical-order assertion is a later tightening).
            Assert.That(produced.OrderBy(x => x, StringComparer.Ordinal),
                Is.EqualTo(expected.OrderBy(x => x, StringComparer.Ordinal)),
                $"{row.Id} ({row.Tests}): form set mismatch.\n  expected: {string.Join(" | ", expected)}\n  got:      {string.Join(" | ", produced)}");
        }

        /// <summary>
        /// Round-trip variant: same oracle, but the CONSENSUS protein is serialized through the REAL
        /// ProteinDbWriter and re-read before digestion. This exercises encode+decode (invariant 6 / #1083)
        /// WITHOUT the writer ever sitting upstream of the ground truth — the input XML is still the
        /// hand-authored, bible-grounded one; the writer only has to preserve what we already loaded.
        /// A writer that drops a variant or a variant-adjacent mod fails here while Corpus_Node stays green.
        /// </summary>
        [TestCaseSource(nameof(Cases))]
        public void Corpus_Node_RoundTrip(Row row)
        {
            List<string> expected = File.ReadAllLines(Path.Combine(CorpusDir, "expected", row.Id + ".txt"))
                .Where(l => l.Trim().Length > 0).ToList();

            List<string> produced = RunNodeRoundTrip(row);

            Assert.That(produced.Count, Is.EqualTo(row.ExpectedCount),
                $"{row.Id} ({row.Tests}) [round-trip]: expected {row.ExpectedCount} forms, got {produced.Count}.\n  got: {string.Join(" | ", produced)}");

            Assert.That(produced.OrderBy(x => x, StringComparer.Ordinal),
                Is.EqualTo(expected.OrderBy(x => x, StringComparer.Ordinal)),
                $"{row.Id} ({row.Tests}) [round-trip]: form set mismatch after write/read cycle.\n  expected: {string.Join(" | ", expected)}\n  got:      {string.Join(" | ", produced)}");
        }

        private static List<string> RunNodeRoundTrip(Row row)
        {
            List<Modification> knownMods = BuildMods(row);
            string xml = BuildXml(row);
            string tmp1 = Path.Combine(Path.GetTempPath(), $"corpus_{row.Id}_{Guid.NewGuid():N}.xml");
            string tmp2 = Path.Combine(Path.GetTempPath(), $"corpus_{row.Id}_rt_{Guid.NewGuid():N}.xml");
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

                var dp = new DigestionParams(protease: row.Protease, maxMissedCleavages: 0, minPeptideLength: 1,
                    maxModificationIsoforms: row.MaxIsoforms, maxModsForPeptides: row.MaxMods,
                    initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

                return reread
                    .SelectMany(p => p.Digest(dp, new List<Modification>(), new List<Modification>()))
                    .Select(pw => pw.FullSequence)
                    .ToList();
            }
            finally { File.Delete(tmp1); File.Delete(tmp2); }
        }

        private static List<string> RunNode(Row row)
        {
            List<Modification> knownMods = BuildMods(row);
            string xml = BuildXml(row);
            string tmp = Path.Combine(Path.GetTempPath(), $"corpus_{row.Id}_{Guid.NewGuid():N}.xml");
            try
            {
                File.WriteAllText(tmp, xml);
                List<Protein> proteins = ProteinDbLoader.LoadProteinXML(tmp,
                    generateTargets: true, decoyType: DecoyType.None,
                    allKnownModifications: knownMods, isContaminant: false,
                    modTypesToExclude: new List<string>(), out _);

                var dp = new DigestionParams(protease: row.Protease, maxMissedCleavages: 0, minPeptideLength: 1,
                    maxModificationIsoforms: row.MaxIsoforms, maxModsForPeptides: row.MaxMods,
                    initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

                return proteins
                    .SelectMany(p => p.Digest(dp, new List<Modification>(), new List<Modification>()))
                    .Select(pw => pw.FullSequence)
                    .ToList();
            }
            finally { File.Delete(tmp); }
        }

        // --- helpers: parse the row's mods/variants DSL and emit the matching mzLib XML ------------------

        private static List<Modification> BuildMods(Row row)
        {
            var mods = new List<Modification>();
            if (row.Mods == "-") return mods;
            foreach (var token in row.Mods.Split(';', StringSplitOptions.RemoveEmptyEntries))
            {
                var (name, pos) = ParseNameAtPos(token.Trim());
                char motifChar = row.Base[pos - 1];
                ModificationMotif.TryGetMotif(motifChar.ToString(), out ModificationMotif motif);
                var (type, mass) = ModRegistry[name];
                mods.Add(new Modification(_originalId: name, _modificationType: type, _target: motif,
                    _locationRestriction: "Anywhere.", _monoisotopicMass: mass));
            }
            return mods;
        }

        private static string BuildXml(Row row)
        {
            var features = new System.Text.StringBuilder();

            if (row.Mods != "-")
            {
                foreach (var token in row.Mods.Split(';', StringSplitOptions.RemoveEmptyEntries))
                {
                    var (name, pos) = ParseNameAtPos(token.Trim());
                    char motifChar = row.Base[pos - 1];
                    // The reader binds a "modified residue" feature by description == Modification.IdWithMotif.
                    string idWithMotif = $"{name} on {motifChar}";
                    features.Append(
                        $"    <feature type=\"modified residue\" description=\"{idWithMotif}\">\n" +
                        $"      <location><position position=\"{pos}\" /></location>\n" +
                        $"    </feature>\n");
                }
            }

            if (row.Variants != "-")
            {
                foreach (var v in row.Variants.Split(';', StringSplitOptions.RemoveEmptyEntries))
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
                $"    <accession>P_{row.Id}</accession>\n" +
                $"    <name>TEST_{row.Id}</name>\n" +
                "    <protein><recommendedName><fullName>Test protein</fullName></recommendedName></protein>\n" +
                features +
                $"    <sequence length=\"{row.Base.Length}\">{row.Base}</sequence>\n" +
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
