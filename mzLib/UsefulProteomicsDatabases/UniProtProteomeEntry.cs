#nullable enable
namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// One row of UniProt's proteome catalogue, as <see cref="ProteinDbRetriever.UniprotProteomeCatalog"/> reads it.
    /// </summary>
    /// <param name="ProteomeId">The UniProt proteome identifier, e.g. "UP000005640".</param>
    /// <param name="Organism">The organism name as UniProt writes it, e.g. "Homo sapiens (Human)". Empty when the cell is blank.</param>
    /// <param name="TaxonId">
    /// The organism's NCBI taxonomy identifier, e.g. "9606", or null when the catalogue has no such column or the
    /// cell is blank. A string, like <c>Protein.NcbiTaxonomyId</c>, so the two compare directly.
    /// </param>
    /// <param name="ProteinCount">The number of proteins in the proteome, or null when absent or not a whole number.</param>
    /// <param name="Busco">UniProt's BUSCO completeness summary, verbatim, or null when absent.</param>
    /// <param name="Cpd">UniProt's Complete Proteome Detector verdict, e.g. "Standard", or null when absent.</param>
    public sealed record UniProtProteomeEntry(
        string ProteomeId,
        string Organism,
        string? TaxonId,
        long? ProteinCount,
        string? Busco,
        string? Cpd);
}
