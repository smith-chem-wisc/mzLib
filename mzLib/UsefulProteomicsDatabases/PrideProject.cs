using System;
using System.Collections.Generic;
using MzLibUtil;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// The metadata describing one PRIDE Archive project, as returned by the PRIDE Archive REST API
    /// (v3 <c>projects/{accession}</c>). A plain data object populated by JSON deserialization;
    /// controlled-vocabulary fields are modeled with <see cref="CvParam"/>.
    /// </summary>
    /// <remarks>
    /// This is the "inspect before you download" companion to
    /// <see cref="PrideArchiveClient.GetProjectFilesAsync"/> — it lets a caller judge and cite a
    /// dataset without fetching gigabytes of spectra.
    /// <para>
    /// Match every controlled-vocabulary member on its <see cref="CvParam.Accession"/>, not on its
    /// <see cref="CvParam.Name"/>: the accession is stable across PRIDE releases and the display
    /// name is not. Members PRIDE omits for a given project arrive as empty collections or empty
    /// strings, never null.
    /// </para>
    /// </remarks>
    public class PrideProject
    {
        /// <summary>The ProteomeXchange accession (e.g. "PXD012345").</summary>
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

        /// <summary>Whether the submission is "COMPLETE" (identifications included) or "PARTIAL".</summary>
        public string SubmissionType { get; set; } = string.Empty;

        /// <summary>The license the data are released under.</summary>
        public string License { get; set; } = string.Empty;

        /// <summary>The date the project was submitted to PRIDE.</summary>
        /// <remarks>
        /// Typed as <see cref="DateTime"/>, not <see cref="DateTimeOffset"/> as on
        /// <see cref="PrideArchiveFile"/>: this endpoint sends a bare calendar date ("2019-01-15")
        /// with no time and no UTC offset, and parsing that into a <see cref="DateTimeOffset"/>
        /// would silently attach whatever offset the executing machine happens to be in.
        /// </remarks>
        public DateTime SubmissionDate { get; set; }

        /// <summary>The date the project became publicly visible. See <see cref="SubmissionDate"/> for why this is a <see cref="DateTime"/>.</summary>
        public DateTime PublicationDate { get; set; }

        /// <summary>PRIDE's coarse classification tags (e.g. "Technical", "Metaproteomics").</summary>
        public List<string> ProjectTags { get; set; } = new();

        /// <summary>The submitter's free-text keywords.</summary>
        public List<string> Keywords { get; set; } = new();

        /// <summary>The countries of the submitting laboratories.</summary>
        public List<string> Countries { get; set; } = new();

        /// <summary>The people who submitted the project.</summary>
        public List<PrideContact> Submitters { get; set; } = new();

        /// <summary>The principal investigators of the submitting laboratories.</summary>
        public List<PrideContact> LabPIs { get; set; } = new();

        /// <summary>The publications reporting this dataset.</summary>
        public List<PrideReference> References { get; set; } = new();

        /// <summary>The mass spectrometers used (e.g. accession "MS:1000449", name "LTQ Orbitrap").</summary>
        public List<CvParam> Instruments { get; set; } = new();

        /// <summary>The software used to process the data.</summary>
        public List<CvParam> Softwares { get; set; } = new();

        /// <summary>The kinds of experiment performed (e.g. "Shotgun proteomics").</summary>
        public List<CvParam> ExperimentTypes { get; set; } = new();

        /// <summary>The quantification methods used, if the project is quantitative.</summary>
        public List<CvParam> QuantificationMethods { get; set; } = new();

        /// <summary>The species the sample came from.</summary>
        public List<CvParam> Organisms { get; set; } = new();

        /// <summary>The tissues or anatomical parts the sample came from.</summary>
        public List<CvParam> OrganismParts { get; set; } = new();

        /// <summary>The diseases under study, if any.</summary>
        public List<CvParam> Diseases { get; set; } = new();

        /// <summary>
        /// The post-translational modifications reported in the dataset. A project with none carries
        /// the explicit term "PRIDE:0000398" ("No PTMs are included in the dataset") rather than an
        /// empty list.
        /// </summary>
        public List<CvParam> IdentifiedPTMStrings { get; set; } = new();

        /// <summary>Any further controlled-vocabulary annotations PRIDE attaches to the project.</summary>
        public List<CvParam> AdditionalAttributes { get; set; } = new();

        /// <summary>The sample characteristics, each a controlled-vocabulary key with one or more values.</summary>
        public List<PrideSampleAttribute> SampleAttributes { get; set; } = new();

        /// <summary>The number of times this project's files have been downloaded.</summary>
        public long TotalFileDownloads { get; set; }

        /// <summary>The portion of <see cref="TotalFileDownloads"/> attributed to automated crawlers.</summary>
        public long BotCount { get; set; }

        /// <summary>The portion of <see cref="TotalFileDownloads"/> attributed to aggregating hubs.</summary>
        public long HubCount { get; set; }

        /// <summary>The portion of <see cref="TotalFileDownloads"/> attributed to genuine human traffic.</summary>
        public long OrganicCount { get; set; }
    }
}
