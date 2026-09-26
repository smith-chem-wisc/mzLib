using System;
using System.Collections.Generic;
using Newtonsoft.Json;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// One hit from a PRIDE Archive project search (v3 <c>search/projects</c>). A plain data object
    /// populated by JSON deserialization.
    /// </summary>
    /// <remarks>
    /// <para>
    /// This is deliberately NOT a <see cref="PrideProject"/>, and the two are not interchangeable.
    /// PRIDE serves search hits from a separate Elasticsearch projection (the API spec names it
    /// <c>ElasticPrideProject</c>) in which every controlled-vocabulary field has been FLATTENED TO A
    /// PLAIN STRING: the same project returns <c>instruments: [{ "accession": "MS:1001911", "name":
    /// "Q Exactive" }]</c> from <c>projects/{accession}</c> and <c>instruments: ["Q Exactive"]</c>
    /// here. Contacts collapse from ten-field objects to a display name, and references to a single
    /// pre-mangled citation string. Deserializing this payload into <see cref="PrideProject"/> would
    /// throw rather than silently degrade.
    /// </para>
    /// <para>
    /// So the collections below are <see cref="string"/>, not <see cref="MzLibUtil.CvParam"/>, and
    /// that is a property of the wire rather than a simplification chosen here. The accessions simply
    /// are not sent, and manufacturing one by looking a display name up in a vocabulary would put an
    /// identifier in a caller's hands that PRIDE never asserted. A caller who needs the controlled
    /// vocabulary should follow the hit's <see cref="Accession"/> to
    /// <see cref="PrideArchiveClient.GetProjectAsync"/>, which returns the full metadata object.
    /// </para>
    /// <para>
    /// Fields PRIDE omits for a given project arrive as empty collections, empty strings or zero,
    /// never null. Several are genuinely sparse: across 1 600 hits sampled from 16 different queries
    /// on 2026-08-21, <see cref="ProjectTags"/> was populated on 2.6%, <see cref="Sdrf"/> on 2.4%,
    /// <see cref="OtherOmicsLinks"/> on 18%, and <see cref="YearlyDownloads"/> with the
    /// <see cref="BotCount"/>/<see cref="HubCount"/>/<see cref="OrganicCount"/> trio on under half.
    /// Treat a zero or an empty list as "not reported", never as a measured zero. How sparse a given
    /// field looks tracks the RESULT SET rather than the archive, so it varies by query.
    /// </para>
    /// <para>
    /// One caveat the never-null guarantee does NOT cover: it applies to the fields, not to the
    /// elements inside them. PRIDE ships empty and whitespace-only strings inside
    /// <see cref="Keywords"/> on roughly 9% of hits, so filter before displaying or joining them.
    /// </para>
    /// </remarks>
    public class PrideProjectSearchResult
    {
        /// <summary>
        /// The project accession (e.g. "PXD012345") — the key to the full metadata object via
        /// <see cref="PrideArchiveClient.GetProjectAsync"/>.
        /// </summary>
        /// <remarks>
        /// Usually a ProteomeXchange "PXD" accession, but not always: search also returns legacy
        /// PRIDE-era "PRD" accessions and affinity "PAD" ones. Both resolve on the metadata endpoint,
        /// so treat this as an opaque key rather than assuming a prefix.
        /// </remarks>
        public string Accession { get; set; } = string.Empty;

        /// <summary>The project title.</summary>
        public string Title { get; set; } = string.Empty;

        /// <summary>The submitter's free-text description of the project.</summary>
        public string ProjectDescription { get; set; } = string.Empty;

        /// <summary>The submitter's description of how the sample was prepared.</summary>
        public string SampleProcessingProtocol { get; set; } = string.Empty;

        /// <summary>The submitter's description of how the data were searched and processed.</summary>
        public string DataProcessingProtocol { get; set; } = string.Empty;

        /// <summary>The dataset DOI, or empty if PRIDE has not minted one.</summary>
        public string Doi { get; set; } = string.Empty;

        /// <summary>
        /// How the project was submitted: "COMPLETE", "PARTIAL", "AFFINITY", or "PRIDE" for legacy
        /// pre-ProteomeXchange submissions. Not treated as a closed set — it is a plain string
        /// because PRIDE has added values before.
        /// </summary>
        public string SubmissionType { get; set; } = string.Empty;

        /// <summary>
        /// The project's SDRF sample metadata, flattened by PRIDE's search index into a single
        /// space-joined bag of the term VALUES — e.g. "Nt=label free sample;ac=ms:1002038
        /// Nt=trypsin;ac=ms:1001251 Mus musculus Liver Male 18 weeks". Empty when the project has no
        /// SDRF, which is the common case.
        /// </summary>
        /// <remarks>
        /// Not a file, a filename or a URL: nothing can be fetched with it, and the row/column
        /// structure of the real SDRF is gone. It is searchable text, useful for reading and matching
        /// rather than for parsing. Rare in practice — non-empty on roughly 2% of hits sampled
        /// 2026-08-21.
        /// </remarks>
        public string Sdrf { get; set; } = string.Empty;

        /// <summary>The date the project was submitted to PRIDE.</summary>
        /// <remarks>
        /// Typed as <see cref="DateTime"/>, not <see cref="DateTimeOffset"/> as on
        /// <see cref="PrideArchiveFile"/>: this endpoint sends a bare calendar date ("2026-08-15")
        /// with no time and no UTC offset, and parsing that into a <see cref="DateTimeOffset"/> would
        /// silently attach whatever offset the executing machine happens to be in. The split between
        /// the two types is per-endpoint and deliberate — it follows what each one actually sends.
        /// </remarks>
        public DateTime SubmissionDate { get; set; }

        /// <summary>The date the project became publicly visible. See <see cref="SubmissionDate"/> for why this is a <see cref="DateTime"/>.</summary>
        public DateTime PublicationDate { get; set; }

        /// <summary>
        /// The date the project was last modified. Has no counterpart on <see cref="PrideProject"/> —
        /// the search projection carries it and the metadata endpoint does not.
        /// </summary>
        public DateTime UpdatedDate { get; set; }

        /// <summary>PRIDE's coarse classification tags (e.g. "Human proteome project").</summary>
        public List<string> ProjectTags { get; set; } = new();

        /// <summary>The submitter's free-text keywords.</summary>
        public List<string> Keywords { get; set; } = new();

        /// <summary>The submitters' display names. Flattened from the ten-field contact objects the metadata endpoint returns.</summary>
        public List<string> Submitters { get; set; } = new();

        /// <summary>The principal investigators' display names. Flattened as <see cref="Submitters"/> is.</summary>
        public List<string> LabPIs { get; set; } = new();

        /// <summary>The submitting laboratories' affiliations. Has no counterpart on <see cref="PrideProject"/>.</summary>
        public List<string> Affiliations { get; set; } = new();

        /// <summary>The instruments used, by display name (e.g. "Orbitrap Exploris 480").</summary>
        public List<string> Instruments { get; set; } = new();

        /// <summary>The software used, by display name.</summary>
        public List<string> Softwares { get; set; } = new();

        /// <summary>The quantification methods used, by display name.</summary>
        public List<string> QuantificationMethods { get; set; } = new();

        /// <summary>
        /// The sample characteristics, by display value (e.g. "liver"). Flattened: the metadata
        /// endpoint returns these as key/value pairs of controlled-vocabulary terms, so which
        /// characteristic each value describes is not recoverable from a search hit.
        /// </summary>
        public List<string> SampleAttributes { get; set; } = new();

        /// <summary>The organisms studied, by display name (e.g. "Mus musculus (mouse)").</summary>
        public List<string> Organisms { get; set; } = new();

        /// <summary>The tissues or organism parts sampled.</summary>
        /// <remarks>
        /// PRIDE spells this "organismsPart" on the search endpoint and "organismParts" on
        /// <c>projects/{accession}</c>. The property keeps the metadata spelling so the two DTOs agree
        /// on the name of the same concept; the attribute carries the search endpoint's own. This is
        /// the one place a PRIDE DTO needs an explicit JSON name — every other member here binds by
        /// convention.
        /// </remarks>
        [JsonProperty("organismsPart")]
        public List<string> OrganismParts { get; set; } = new();

        /// <summary>The diseases studied, by display name.</summary>
        public List<string> Diseases { get; set; } = new();

        /// <summary>
        /// The associated publications, each as a single pre-formatted citation string. Flattened: the
        /// metadata endpoint returns structured <see cref="PrideReference"/> objects, so a PubMed ID or
        /// DOI cannot be read out of a search hit without parsing the string PRIDE assembled.
        /// </summary>
        public List<string> References { get; set; } = new();

        /// <summary>The experiment types, by display name (e.g. "Data-independent acquisition").</summary>
        public List<string> ExperimentTypes { get; set; } = new();

        /// <summary>
        /// The names of the project's files. A search-only convenience: it is not the manifest, and
        /// carries no sizes, categories or download locations. Use
        /// <see cref="PrideArchiveClient.GetProjectFilesAsync"/> — or
        /// <see cref="PrideArchiveClient.GetProjectFilesFromFtpAsync"/> for the complete list — to act
        /// on the files themselves.
        /// </summary>
        public List<string> ProjectFileNames { get; set; } = new();

        /// <summary>Links to related datasets in other omics repositories (e.g. "pride.project:PXD062751").</summary>
        public List<string> OtherOmicsLinks { get; set; } = new();

        /// <summary>
        /// The snippets that matched the query, keyed by the field each match was found in, with the
        /// matched terms wrapped in <c>&lt;em&gt;</c> markup.
        /// </summary>
        /// <remarks>
        /// A dynamic map, not a fixed shape: the keys are whichever fields the query happened to hit,
        /// so they vary per hit and per query. This is the one thing a search returns that
        /// <see cref="PrideArchiveClient.GetProjectAsync"/> cannot — it says WHY a project matched.
        /// The markup is PRIDE's; strip it if the value is going anywhere other than a highlighted view.
        /// </remarks>
        public Dictionary<string, List<string>> Highlights { get; set; } = new();

        /// <summary>Downloads recorded for the project, broken out by calendar year.</summary>
        public List<PrideYearlyDownloadCount> YearlyDownloads { get; set; } = new();

        /// <summary>Total downloads recorded for the project.</summary>
        public long DownloadCount { get; set; }

        /// <summary>Mean downloads per file in the project.</summary>
        public double AvgDownloadsPerFile { get; set; }

        /// <summary>The project's download-popularity percentile within PRIDE.</summary>
        public int Percentile { get; set; }

        /// <summary>Downloads attributed to automated crawlers.</summary>
        public long BotCount { get; set; }

        /// <summary>Downloads attributed to institutional or aggregating hubs.</summary>
        public long HubCount { get; set; }

        /// <summary>Downloads attributed to ordinary human traffic.</summary>
        public long OrganicCount { get; set; }
    }
}
