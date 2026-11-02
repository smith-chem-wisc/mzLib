namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// How many times a project was downloaded in one calendar year, as reported by the PRIDE Archive
    /// search endpoint (v3 <c>search/projects</c>) in <see cref="PrideProjectSearchResult.YearlyDownloads"/>.
    /// </summary>
    /// <remarks>
    /// The year arrives as a string ("2026") rather than a number, and is kept as one: it is a label on
    /// a bucket, not a quantity to compute with, and parsing it here would invent a failure mode for a
    /// value no caller needs to do arithmetic on.
    /// </remarks>
    public class PrideYearlyDownloadCount
    {
        /// <summary>The calendar year the count covers, e.g. "2026".</summary>
        public string Year { get; set; } = string.Empty;

        /// <summary>Downloads recorded for the project in <see cref="Year"/>.</summary>
        public long Count { get; set; }
    }
}
