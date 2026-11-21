namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// The regexes <see cref="ProteinDbLoader.LoadProteinFasta"/> uses for one FASTA header format.
    /// A null field is not parsed by that format; when the loader auto-detects, a null field leaves
    /// the caller's regex for that field in place.
    /// </summary>
    public sealed record FastaHeaderFieldRegexBundle(
        FastaHeaderFieldRegex Accession,
        FastaHeaderFieldRegex FullName,
        FastaHeaderFieldRegex Name,
        FastaHeaderFieldRegex GeneName,
        FastaHeaderFieldRegex Organism,
        FastaHeaderFieldRegex OrganismId);
}
