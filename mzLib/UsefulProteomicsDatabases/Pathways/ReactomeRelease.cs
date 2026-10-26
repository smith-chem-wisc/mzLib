using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Security.Cryptography;
using System.Text.RegularExpressions;

namespace UsefulProteomicsDatabases.Pathways
{
    /// <summary>
    /// One Reactome release, read from Reactome's two flat files: ReactomePathways.txt (pathway id, name, species, for
    /// every species) and ReactomePathwaysRelation.txt (parent id, child id).
    ///
    /// UniProt cross-references only the lowest-level pathways a protein takes part in
    /// (Proteomics.Protein.ReactomePathways), so rolling a membership up to a broad area such as "Metabolism" needs this
    /// hierarchy. Like go.obo for GO, it is a small pinned input: the release records both files' sha256s, and the
    /// caller states the release number, because neither file carries one.
    ///
    /// The hierarchy is a DAG (388 pathways in release 96 have two or more parents), so <see cref="Ancestors"/> is a
    /// closure with each ancestor once. Loading refuses a cycle and an edge naming a pathway the pathways file lacks,
    /// rather than walking either.
    /// </summary>
    public sealed class ReactomeRelease
    {
        /// <summary>The species prefix of a human Reactome id.</summary>
        public const string HumanPrefix = "R-HSA-";

        /// <summary>A Reactome stable id: R-, a three-letter species code, -, a number.</summary>
        private static readonly Regex StableId = new(@"^R-[A-Z]{3}-(?<number>\d+)$", RegexOptions.Compiled);

        private readonly Dictionary<string, ReactomePathwayRow> _pathways;
        private readonly Dictionary<string, List<string>> _parents;
        private readonly Dictionary<string, IReadOnlySet<string>> _ancestors = new(StringComparer.Ordinal);

        private ReactomeRelease(Dictionary<string, (string Name, string Species)> pathways, List<ReactomeHierarchyRow> hierarchy,
            string release, string pathwaysSha256, string relationSha256)
        {
            Release = release;
            PathwaysSha256 = pathwaysSha256;
            RelationSha256 = relationSha256;

            _parents = new Dictionary<string, List<string>>(StringComparer.Ordinal);
            foreach (var edge in hierarchy)
            {
                if (!_parents.TryGetValue(edge.ChildId, out var parents))
                {
                    _parents[edge.ChildId] = parents = new List<string>();
                }
                parents.Add(edge.ParentId);
            }

            _pathways = new Dictionary<string, ReactomePathwayRow>(StringComparer.Ordinal);
            foreach (var (id, (name, species)) in pathways)
            {
                _pathways[id] = new ReactomePathwayRow(id, name, species, HumanEquivalentOf(id, name, pathways));
            }

            Pathways = _pathways.Keys.OrderBy(id => id, StringComparer.Ordinal).Select(id => _pathways[id]).ToList();
            Hierarchy = hierarchy.OrderBy(e => e.ParentId, StringComparer.Ordinal).ThenBy(e => e.ChildId, StringComparer.Ordinal).ToList();
            Roots = Pathways.Select(p => p.PathwayId).Where(id => !_parents.ContainsKey(id)).ToList();
        }

        /// <summary>The release number the caller stated, e.g. "96".</summary>
        public string Release { get; }

        /// <summary>Lower-case hex sha256 of the pathways file's bytes as read.</summary>
        public string PathwaysSha256 { get; }

        /// <summary>Lower-case hex sha256 of the relation file's bytes as read.</summary>
        public string RelationSha256 { get; }

        /// <summary>Every pathway of the release, every species, in ordinal order of id.</summary>
        public IReadOnlyList<ReactomePathwayRow> Pathways { get; }

        /// <summary>Every parent-child pair, ordered by parent then child (ordinal).</summary>
        public IReadOnlyList<ReactomeHierarchyRow> Hierarchy { get; }

        /// <summary>The pathways with no parent (the top-level areas, e.g. "Metabolism"), in ordinal order.</summary>
        public IReadOnlyList<string> Roots { get; }

        /// <summary>True when the release has this pathway.</summary>
        public bool Contains(string pathwayId) => pathwayId != null && _pathways.ContainsKey(pathwayId);

        /// <summary>
        /// Every pathway above <paramref name="pathwayId"/>, each once, the pathway itself excluded. Empty for a
        /// top-level pathway. Computed once per pathway and shared, so it is frozen.
        /// </summary>
        /// <exception cref="ArgumentException">The release has no such pathway.</exception>
        public IReadOnlySet<string> Ancestors(string pathwayId)
        {
            RequirePathway(pathwayId);
            if (_ancestors.TryGetValue(pathwayId, out var cached))
            {
                return cached;
            }
            var ancestors = new HashSet<string>(StringComparer.Ordinal);
            var stack = new Stack<string>(ParentsOf(pathwayId));
            while (stack.Count > 0)
            {
                string next = stack.Pop();
                if (ancestors.Add(next))
                {
                    foreach (string parent in ParentsOf(next))
                    {
                        stack.Push(parent);
                    }
                }
            }
            IReadOnlySet<string> frozen = ancestors;
            _ancestors[pathwayId] = frozen;
            return frozen;
        }

        /// <summary>
        /// The human pathway this one is, or null. A human pathway is its own equivalent. Another species' pathway
        /// (R-MMU-77289) is equivalent to R-HSA with the same number only when that pathway exists in this release with
        /// the identical name, which is how Reactome names the events it infers from human. A pathway curated in its
        /// own species has none.
        /// </summary>
        /// <exception cref="ArgumentException">The release has no such pathway.</exception>
        public string HumanEquivalentId(string pathwayId)
        {
            RequirePathway(pathwayId);
            return _pathways[pathwayId].HumanEquivalentId;
        }

        /// <summary>
        /// Reads a release from its two files. Blank lines are skipped; every other pathways line is
        /// id TAB name TAB species, and every relation line is parent TAB child. Names are kept verbatim.
        /// </summary>
        /// <param name="pathwaysPath">ReactomePathways.txt of the release.</param>
        /// <param name="relationPath">ReactomePathwaysRelation.txt of the same release.</param>
        /// <param name="release">The release number, e.g. "96". Reactome's files carry none.</param>
        /// <exception cref="ArgumentNullException">A path is null.</exception>
        /// <exception cref="ArgumentException">The release is null or empty.</exception>
        /// <exception cref="FileNotFoundException">A file does not exist.</exception>
        /// <exception cref="InvalidDataException">A line has the wrong number of cells or an empty id; a pathway id
        /// appears twice; a pair names a pathway the pathways file lacks; or the hierarchy has a cycle.</exception>
        public static ReactomeRelease Load(string pathwaysPath, string relationPath, string release)
        {
            ArgumentNullException.ThrowIfNull(pathwaysPath);
            ArgumentNullException.ThrowIfNull(relationPath);
            if (string.IsNullOrEmpty(release))
            {
                throw new ArgumentException("State the Reactome release number; its files carry none.", nameof(release));
            }

            var pathways = new Dictionary<string, (string Name, string Species)>(StringComparer.Ordinal);
            string pathwaysFile = Path.GetFileName(pathwaysPath);
            foreach (var (cells, lineNumber) in Lines(pathwaysPath, 3, "id, name, species"))
            {
                if (cells[0].Length == 0)
                {
                    throw new InvalidDataException($"{pathwaysFile} line {lineNumber}: the pathway id is empty.");
                }
                if (!pathways.TryAdd(cells[0], (cells[1], cells[2])))
                {
                    throw new InvalidDataException($"{pathwaysFile} line {lineNumber}: pathway {cells[0]} appears twice.");
                }
            }

            var hierarchy = new List<ReactomeHierarchyRow>();
            var seen = new HashSet<(string, string)>();
            string relationFile = Path.GetFileName(relationPath);
            foreach (var (cells, lineNumber) in Lines(relationPath, 2, "parent, child"))
            {
                foreach (string id in cells)
                {
                    if (!pathways.ContainsKey(id))
                    {
                        throw new InvalidDataException($"{relationFile} line {lineNumber}: {(id.Length == 0 ? "an empty id" : id)} is not a pathway in {pathwaysFile}.");
                    }
                }
                if (seen.Add((cells[0], cells[1])))
                {
                    hierarchy.Add(new ReactomeHierarchyRow(cells[0], cells[1]));
                }
            }

            var loaded = new ReactomeRelease(pathways, hierarchy, release, Sha256Of(pathwaysPath), Sha256Of(relationPath));
            loaded.RefuseCycles(relationFile);
            return loaded;
        }

        private IEnumerable<string> ParentsOf(string pathwayId) =>
            _parents.TryGetValue(pathwayId, out var parents) ? parents : Enumerable.Empty<string>();

        private void RequirePathway(string pathwayId)
        {
            if (!Contains(pathwayId))
            {
                throw new ArgumentException($"{pathwayId} is not a pathway in Reactome release {Release}.", nameof(pathwayId));
            }
        }

        /// <summary>
        /// Depth-first over parent links. A pathway is on the current path (false) or finished (true); meeting one that
        /// is still on the path closes a cycle.
        /// </summary>
        private void RefuseCycles(string relationFile)
        {
            var finished = new Dictionary<string, bool>(StringComparer.Ordinal);
            foreach (string start in _pathways.Keys)
            {
                if (finished.ContainsKey(start))
                {
                    continue;
                }
                var stack = new Stack<(string Id, IEnumerator<string> Parents)>();
                finished[start] = false;
                stack.Push((start, ParentsOf(start).GetEnumerator()));
                while (stack.Count > 0)
                {
                    var (id, parents) = stack.Peek();
                    if (!parents.MoveNext())
                    {
                        stack.Pop();
                        finished[id] = true;
                        continue;
                    }
                    string parent = parents.Current;
                    if (!finished.TryGetValue(parent, out bool parentFinished))
                    {
                        finished[parent] = false;
                        stack.Push((parent, ParentsOf(parent).GetEnumerator()));
                    }
                    else if (!parentFinished)
                    {
                        throw new InvalidDataException($"{relationFile}: the hierarchy has a cycle through {parent} and {id}.");
                    }
                }
            }
        }

        private static string HumanEquivalentOf(string id, string name, Dictionary<string, (string Name, string Species)> pathways)
        {
            if (id.StartsWith(HumanPrefix, StringComparison.Ordinal))
            {
                return id;
            }
            var match = StableId.Match(id);
            if (!match.Success)
            {
                return null;
            }
            string human = HumanPrefix + match.Groups["number"].Value;
            return pathways.TryGetValue(human, out var humanPathway) && string.Equals(humanPathway.Name, name, StringComparison.Ordinal)
                ? human
                : null;
        }

        private static IEnumerable<(string[] Cells, int LineNumber)> Lines(string path, int width, string expected)
        {
            if (!File.Exists(path))
            {
                throw new FileNotFoundException("Reactome release file not found.", path);
            }
            string fileName = Path.GetFileName(path);
            int lineNumber = 0;
            foreach (string line in File.ReadLines(path))
            {
                lineNumber++;
                if (line.Length == 0)
                {
                    continue;
                }
                string[] cells = line.Split('\t');
                if (cells.Length != width)
                {
                    throw new InvalidDataException($"{fileName} line {lineNumber}: expected {width} cells ({expected}), found {cells.Length}.");
                }
                yield return (cells, lineNumber);
            }
        }

        private static string Sha256Of(string path)
        {
            using var stream = File.OpenRead(path);
            return Convert.ToHexString(SHA256.HashData(stream)).ToLowerInvariant();
        }
    }
}
