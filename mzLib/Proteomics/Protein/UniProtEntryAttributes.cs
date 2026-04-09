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

        /// <summary>
        /// Helper property to get the current date formatted as yyyy-MM-dd for use in defaults.
        /// Note: This is computed on access and may return different values at different times.
        /// </summary>
        private static string CurrentDate => DateTime.Now.ToString("yyyy-MM-dd");

        /// <summary>
        /// Creates a new UniProtEntryAttributes instance.
        /// </summary>
        /// <param name="dataset">The dataset name. Defaults to "unknown" which is supported by ProSightPD and ProSight Annotator.</param>
        /// <param name="created">The created date in yyyy-MM-dd format. Defaults to current date.</param>
        /// <param name="modified">The modified date in yyyy-MM-dd format. Defaults to current date.</param>
        /// <param name="version">The version number. Defaults to "1".</param>
        /// <param name="xmlns">The XML namespace. Defaults to UniProt namespace.</param>
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
