using System;

namespace Proteomics
{
    /// <summary>
    /// Stores the UniProt XML entry-level attributes parsed from the &lt;entry&gt; element:
    /// dataset, created date, modified date, version, and XML namespace.
    /// </summary>
    public class UniProtEntryAttributes
    {
        public string Dataset { get; }
        public string Created { get; }
        public string Modified { get; }
        public string Version { get; }
        public string Xmlns { get; }

        public string CurrentDate => DateTime.Now.ToString("yyyy-MM-dd");

        public UniProtEntryAttributes(
            string dataset = "unknown",
            string? created = null,
            string? modified = null,
            string? version = null,
            string xmlns = "http://uniprot.org/uniprot")
        {
            Dataset = dataset;
            Created = created ?? CurrentDate;
            Modified = modified ?? CurrentDate;
            Version = version ?? "1";
            Xmlns = xmlns;
        }
    }
}
