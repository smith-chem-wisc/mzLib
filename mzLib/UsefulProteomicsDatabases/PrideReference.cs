namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// A publication associated with a PRIDE Archive project, as returned by the PRIDE Archive
    /// REST API (v3 <c>projects/{accession}</c>). A plain data object populated by JSON
    /// deserialization.
    /// </summary>
    /// <remarks>
    /// PRIDE supplies the citation as a single pre-formatted line rather than as separate author,
    /// journal, volume and page fields, so it is preserved verbatim in <see cref="ReferenceLine"/>
    /// rather than parsed. A project may be published before its reference is registered, in which
    /// case <see cref="PubmedId"/> is 0 and <see cref="Doi"/> is empty.
    /// </remarks>
    public class PrideReference
    {
        /// <summary>The full citation exactly as PRIDE formats it.</summary>
        public string ReferenceLine { get; set; } = string.Empty;

        /// <summary>The PubMed identifier, or 0 if the publication has none.</summary>
        public int PubmedId { get; set; }

        /// <summary>The publication DOI, or empty if it has none.</summary>
        public string Doi { get; set; } = string.Empty;
    }
}
