#nullable enable
namespace UsefulProteomicsDatabases.Taxonomy
{
    /// <summary>
    /// One taxon of the NCBI taxonomy as UniProt serves it: an organism, or a group of organisms at some rank.
    /// </summary>
    /// <param name="TaxonId">The NCBI taxonomy identifier, e.g. 9606 for Homo sapiens.</param>
    /// <param name="ParentTaxonId">The taxon directly above this one, or null for a root.</param>
    /// <param name="Rank">
    /// The rank exactly as the source states it, e.g. "family", "species", "no rank", "clade". Empty when the source
    /// states none, which is the case for the organism of a proteome: UniProt names it without a rank. Compare ranks
    /// through <see cref="TaxonomicRank"/>, which knows that "superkingdom" and "domain" are the same level.
    /// </param>
    /// <param name="ScientificName">The scientific name, e.g. "Enterobacteriaceae".</param>
    /// <param name="CommonName">The common name, e.g. "enterobacteria", or null when there is none.</param>
    /// <param name="IsHidden">
    /// True when the source marks the taxon as hidden from its displayed lineage ("cellular organisms", most unranked
    /// clades). A hidden taxon is still a real ancestor and still answers rank queries; the flag is for display only.
    /// </param>
    public sealed record TaxonomyNode(
        long TaxonId,
        long? ParentTaxonId,
        string Rank,
        string ScientificName,
        string? CommonName,
        bool IsHidden);
}
