using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace UsefulProteomicsDatabases.GeneOntology
{
    /// <summary>
    /// One category a term belongs to under a map.
    /// </summary>
    /// <param name="Category">The category label.</param>
    /// <param name="Subcategory">"category:subcategory" -- self-qualifying, so a subcategory can never be read
    /// against the wrong category when set-valued columns sit side by side -- or null when the term reaches
    /// the category only through a category-level anchor.</param>
    public sealed record GoCategory(string Category, string Subcategory);

    /// <summary>
    /// Applies a <see cref="GoCategoryMap"/> to terms of one <see cref="GeneOntologyGraph"/> release.
    ///
    /// The rule: a term belongs to a category when one of that category's anchors is the term itself or
    /// one of its ancestors (is_a + part_of). Within one category the most specific matching anchor wins --
    /// an anchor is dropped when another matching anchor of the same category lies below it -- so a term
    /// under "mitochondrion:inner_membrane" does not also report bare "mitochondrion". Across categories
    /// nothing is dropped: a term under two categories reports both.
    ///
    /// The map and the release are bound here, not at load, because whether an anchor exists depends on
    /// which release is loaded. A missing anchor is an error, never an empty category.
    /// </summary>
    public sealed class GoCategoryResolver
    {
        private readonly List<(GoCategoryAnchor Anchor, string PrimaryId)> _anchors;
        private readonly ConcurrentDictionary<string, IReadOnlyList<GoCategory>> _cache = new(StringComparer.Ordinal);

        /// <exception cref="ArgumentNullException">Either argument is null.</exception>
        /// <exception cref="InvalidDataException">
        /// An anchor is absent from the release, or is obsolete there. An obsolete term is no one's ancestor,
        /// so an anchor on it could only ever match itself.
        /// </exception>
        public GoCategoryResolver(GoCategoryMap map, GeneOntologyGraph ontology)
        {
            Map = map ?? throw new ArgumentNullException(nameof(map));
            Ontology = ontology ?? throw new ArgumentNullException(nameof(ontology));

            _anchors = new List<(GoCategoryAnchor, string)>();
            foreach (var anchor in map.Anchors)
            {
                if (!ontology.TryGetTerm(anchor.AnchorGoId, out var term))
                {
                    throw new InvalidDataException(
                        $"{map.SourceFileName}: anchor {anchor.AnchorGoId} is not a term in Gene Ontology release {ontology.Release ?? "(unversioned)"}.");
                }
                if (term.IsObsolete)
                {
                    throw new InvalidDataException(
                        $"{map.SourceFileName}: anchor {anchor.AnchorGoId} is obsolete in Gene Ontology release {ontology.Release ?? "(unversioned)"}.");
                }
                _anchors.Add((anchor, term.Id));
            }
        }

        /// <summary>The map being applied.</summary>
        public GoCategoryMap Map { get; }

        /// <summary>The release the map is applied against.</summary>
        public GeneOntologyGraph Ontology { get; }

        /// <summary>
        /// The categories <paramref name="goId"/> belongs to, ordered by category then subcategory (ordinal,
        /// category-level entries first). Empty when no anchor is the term or above it. An alternative id is
        /// resolved to its primary term first.
        /// </summary>
        /// <exception cref="ArgumentException">The id is not a term in this release.</exception>
        public IReadOnlyList<GoCategory> Categorize(string goId)
        {
            var ancestors = Ontology.Ancestors(goId);
            Ontology.TryGetTerm(goId, out var term);
            return _cache.GetOrAdd(term.Id, id => Compute(id, ancestors));
        }

        private IReadOnlyList<GoCategory> Compute(string termId, IReadOnlySet<string> ancestors)
        {
            var matched = _anchors
                .Where(a => a.PrimaryId == termId || ancestors.Contains(a.PrimaryId))
                .ToList();

            var kept = matched.Where(a => !matched.Any(other =>
                    other.Anchor.Category == a.Anchor.Category
                    && other.PrimaryId != a.PrimaryId
                    && Ontology.Ancestors(other.PrimaryId).Contains(a.PrimaryId)))
                .ToList();

            return kept
                .Select(a => new GoCategory(a.Anchor.Category,
                    a.Anchor.Subcategory == null ? null : a.Anchor.Category + ":" + a.Anchor.Subcategory))
                .Distinct()
                .OrderBy(c => c.Category, StringComparer.Ordinal)
                .ThenBy(c => c.Subcategory ?? "", StringComparer.Ordinal)
                .ToList();
        }
    }
}
