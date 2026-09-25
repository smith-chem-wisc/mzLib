using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Security.Cryptography;
using Proteomics;

namespace UsefulProteomicsDatabases.GeneOntology
{
    /// <summary>
    /// go.obo read with its edges kept: the Gene Ontology as a DAG over is_a and part_of.
    ///
    /// Protein.GoTerms says which terms a UniProt entry was annotated with. This says what those terms
    /// imply: a protein annotated to "mitochondrial inner membrane" is in the mitochondrion only because
    /// the ontology says the one is part_of the other. ControlledVocabulary cannot answer that -- it keeps
    /// accession and name and drops every edge.
    ///
    /// go.obo is ~37 MB, so it is never embedded in the assembly (see ControlledVocabulary for the line).
    /// The caller supplies the file -- <see cref="Loaders.LoadGeneOntology"/> downloads and caches it -- and
    /// the graph records the file's name, sha256 and data-version, which travel with anything resolved
    /// against it.
    ///
    /// The file is read line by line rather than through TopDownProteomics' OboParser, because a
    /// malformed line has to be reported with its line number, and because the whole text never needs to
    /// be held in memory.
    /// </summary>
    public sealed class GeneOntologyGraph
    {
        private readonly Dictionary<string, GeneOntologyTerm> _terms;
        private readonly Dictionary<string, string> _altIdToPrimary;
        private readonly ConcurrentDictionary<string, IReadOnlySet<string>> _ancestors = new(StringComparer.Ordinal);

        private GeneOntologyGraph(Dictionary<string, GeneOntologyTerm> terms, Dictionary<string, string> altIdToPrimary,
            string sourceFileName, string sourceSha256, string release)
        {
            _terms = terms;
            _altIdToPrimary = altIdToPrimary;
            SourceFileName = sourceFileName;
            SourceSha256 = sourceSha256;
            Release = release;
            TermIds = terms.Keys.OrderBy(k => k, StringComparer.Ordinal).ToList();
        }

        /// <summary>The file name the graph was read from.</summary>
        public string SourceFileName { get; }

        /// <summary>Lower-case hex sha256 of the file's bytes as read.</summary>
        public string SourceSha256 { get; }

        /// <summary>
        /// The header's data-version, e.g. "releases/2026-07-26", or null when the file has none. Record it
        /// alongside anything resolved here: it is what makes the ontology a result was computed against
        /// auditable.
        /// </summary>
        public string Release { get; }

        /// <summary>Every primary term id, obsolete terms included, in ordinal order.</summary>
        public IReadOnlyList<string> TermIds { get; }

        /// <summary>Number of terms, obsolete terms included. Alternative ids are not counted.</summary>
        public int Count => _terms.Count;

        /// <summary>
        /// Looks a term up by its primary id or by one of its alternative ids; an alternative id returns the
        /// primary term. False for a null id or an id absent from this release.
        /// </summary>
        public bool TryGetTerm(string goId, out GeneOntologyTerm term)
        {
            term = null;
            if (goId == null)
            {
                return false;
            }
            if (_altIdToPrimary.TryGetValue(goId, out var primary))
            {
                goId = primary;
            }
            return _terms.TryGetValue(goId, out term);
        }

        /// <summary>
        /// Every term reachable from <paramref name="goId"/> over is_a and part_of, excluding the term
        /// itself, as primary ids. An alternative id is resolved first. An obsolete term has none.
        /// </summary>
        /// <exception cref="ArgumentException">
        /// The id is not in this release. A UniProt entry can cite a term newer than the pinned go.obo, and
        /// an empty set would read as "a root" rather than "not in this release".
        /// </exception>
        public IReadOnlySet<string> Ancestors(string goId)
        {
            if (!TryGetTerm(goId, out var term))
            {
                throw new ArgumentException($"{goId} is not a term in Gene Ontology release {Release ?? "(unversioned)"}.", nameof(goId));
            }
            return _ancestors.GetOrAdd(term.Id, ComputeAncestors);
        }

        private IReadOnlySet<string> ComputeAncestors(string id)
        {
            var found = new HashSet<string>(StringComparer.Ordinal);
            var pending = new Stack<string>(Parents(_terms[id]));
            while (pending.Count > 0)
            {
                string next = pending.Pop();
                if (!found.Add(next))
                {
                    continue;
                }
                foreach (string parent in Parents(_terms[next]))
                {
                    pending.Push(parent);
                }
            }
            // GO has no cycles, but a malformed file might; the visited set is what terminates the walk, and a
            // term that reaches itself is still not its own ancestor.
            found.Remove(id);
            return found;
        }

        private static IEnumerable<string> Parents(GeneOntologyTerm term) => term.IsAParents.Concat(term.PartOfParents);

        /// <summary>
        /// Reads a go.obo file. [Typedef] and [Instance] stanzas are skipped; every [Term] is kept, obsolete
        /// terms included. Tag values have their trailing "! comment" and "{qualifier}" removed.
        /// </summary>
        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="InvalidDataException">
        /// A line inside a stanza is not "tag: value"; a [Term] has no id; two terms share an id or an
        /// alternative id; or an is_a / part_of edge names a term absent from the file. The last is what a
        /// truncated download looks like, and it would otherwise make every ancestor query above it short.
        /// </exception>
        public static GeneOntologyGraph Load(string path)
        {
            if (!File.Exists(path))
            {
                throw new FileNotFoundException("Gene Ontology file not found.", path);
            }

            string fileName = Path.GetFileName(path);
            string sha256;
            using (var hashStream = File.OpenRead(path))
            {
                sha256 = Convert.ToHexString(SHA256.HashData(hashStream)).ToLowerInvariant();
            }

            var terms = new Dictionary<string, GeneOntologyTerm>(StringComparer.Ordinal);
            var altIds = new Dictionary<string, string>(StringComparer.Ordinal);
            string release = null;
            StanzaBuilder stanza = null;
            bool inHeader = true;
            // Set inside a stanza other than [Term]: its lines are still checked for well-formedness, then dropped.
            bool skipping = false;

            void Finish()
            {
                if (stanza == null)
                {
                    return;
                }
                var term = stanza.Build(fileName);
                if (!terms.TryAdd(term.Id, term))
                {
                    throw new InvalidDataException($"{fileName} line {stanza.StartLine}: term {term.Id} appears twice.");
                }
                stanza = null;
            }

            using (var reader = new StreamReader(path))
            {
                int lineNumber = 0;
                string line;
                while ((line = reader.ReadLine()) != null)
                {
                    lineNumber++;
                    string trimmed = line.Trim();
                    if (trimmed.Length == 0 || trimmed.StartsWith("!", StringComparison.Ordinal))
                    {
                        continue;
                    }

                    if (trimmed.StartsWith("[", StringComparison.Ordinal))
                    {
                        Finish();
                        inHeader = false;
                        stanza = trimmed == "[Term]" ? new StanzaBuilder(lineNumber) : null;
                        skipping = stanza == null;
                        continue;
                    }

                    int separator = trimmed.IndexOf(':');
                    if (separator <= 0)
                    {
                        throw new InvalidDataException($"{fileName} line {lineNumber}: expected 'tag: value', found '{trimmed}'.");
                    }
                    string tag = trimmed.Substring(0, separator);
                    string value = trimmed.Substring(separator + 1).Trim();

                    if (inHeader)
                    {
                        if (tag == "data-version")
                        {
                            release ??= value;
                        }
                        continue;
                    }
                    if (skipping)
                    {
                        continue;
                    }
                    stanza.Add(tag, StripTrailing(value));
                }
                Finish();
            }

            foreach (var term in terms.Values)
            {
                foreach (string alt in term.AltIds)
                {
                    if (terms.ContainsKey(alt) || !altIds.TryAdd(alt, term.Id))
                    {
                        throw new InvalidDataException($"{fileName}: alternative id {alt} of {term.Id} is already in use.");
                    }
                }
                foreach (string parent in Parents(term))
                {
                    if (!terms.ContainsKey(parent))
                    {
                        throw new InvalidDataException($"{fileName}: {term.Id} names parent {parent}, which is not in the file.");
                    }
                }
            }

            return new GeneOntologyGraph(terms, altIds, fileName, sha256, release);
        }

        /// <summary>Removes an OBO trailing comment ("! name") and trailing qualifier block ("{...}").</summary>
        private static string StripTrailing(string value)
        {
            int comment = value.IndexOf(" !", StringComparison.Ordinal);
            if (comment >= 0)
            {
                value = value.Substring(0, comment);
            }
            if (value.EndsWith("}", StringComparison.Ordinal))
            {
                int open = value.LastIndexOf(" {", StringComparison.Ordinal);
                if (open >= 0)
                {
                    value = value.Substring(0, open);
                }
            }
            return value.Trim();
        }

        private static GoAspect AspectOf(string ns) => ns switch
        {
            "biological_process" => GoAspect.BiologicalProcess,
            "cellular_component" => GoAspect.CellularComponent,
            "molecular_function" => GoAspect.MolecularFunction,
            _ => GoAspect.Unknown
        };

        private sealed class StanzaBuilder
        {
            private string _id;
            private string _name = "";
            private GoAspect _aspect = GoAspect.Unknown;
            private bool _obsolete;
            private readonly List<string> _altIds = new();
            private readonly List<string> _replacedBy = new();
            private readonly List<string> _isA = new();
            private readonly List<string> _partOf = new();

            public StanzaBuilder(int startLine) => StartLine = startLine;

            public int StartLine { get; }

            public void Add(string tag, string value)
            {
                switch (tag)
                {
                    case "id": _id = value; break;
                    case "name": _name = value; break;
                    case "namespace": _aspect = AspectOf(value); break;
                    case "is_obsolete": _obsolete = value == "true"; break;
                    case "alt_id": _altIds.Add(value); break;
                    case "replaced_by": _replacedBy.Add(value); break;
                    case "is_a": _isA.Add(value); break;
                    case "relationship":
                        // "part_of GO:0005740"; any other relationship type is not a containment edge.
                        string[] parts = value.Split(' ', 2, StringSplitOptions.RemoveEmptyEntries);
                        if (parts.Length == 2 && parts[0] == "part_of")
                        {
                            _partOf.Add(parts[1].Trim());
                        }
                        break;
                }
            }

            public GeneOntologyTerm Build(string fileName)
            {
                if (string.IsNullOrEmpty(_id))
                {
                    throw new InvalidDataException($"{fileName} line {StartLine}: [Term] has no id.");
                }
                // An obsolete term keeps its id but leaves the DAG: GO strips its parents, and a file that did not
                // would otherwise let an obsolete term be counted as somebody's ancestor.
                return new GeneOntologyTerm(_id, _name, _aspect, _obsolete, _altIds, _replacedBy,
                    _obsolete ? Array.Empty<string>() : _isA, _obsolete ? Array.Empty<string>() : _partOf);
            }
        }
    }
}
