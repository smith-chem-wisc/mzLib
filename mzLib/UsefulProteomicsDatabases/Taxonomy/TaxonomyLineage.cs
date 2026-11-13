#nullable enable
namespace UsefulProteomicsDatabases.Taxonomy
{
    /// <summary>
    /// Where a taxonomy's taxa came from, carried in the header of every table written from them.
    /// </summary>
    /// <param name="Source">Where the taxa were read from: a URL, or a description of a hand-made file.</param>
    /// <param name="Release">The source's own release name, e.g. UniProt's "2026_04"; null when it has none.</param>
    /// <param name="Query">The query that selected the taxa, e.g. "reference:true"; null when there was none.</param>
    public sealed record TaxonomyProvenance(string Source, string? Release = null, string? Query = null);

    /// <summary>
    /// One taxon's place at each of the eight major ranks: the flat kingdom-to-species view of a
    /// <see cref="TaxonomyTree"/>, for a table column or a plot axis. Made by <see cref="TaxonomyTree.TryGetLineage"/>.
    /// </summary>
    /// <remarks>
    /// A rank is null when the taxon has no ancestor at it, which is an answer and not a gap: many viruses have no
    /// order, and a taxon above the rank of species has no species.
    /// </remarks>
    /// <param name="Taxon">The taxon the lineage is of.</param>
    public sealed record TaxonomyLineage(
        TaxonomyNode Taxon,
        TaxonomyNode? Domain,
        TaxonomyNode? Kingdom,
        TaxonomyNode? Phylum,
        TaxonomyNode? Class,
        TaxonomyNode? Order,
        TaxonomyNode? Family,
        TaxonomyNode? Genus,
        TaxonomyNode? Species);
}
