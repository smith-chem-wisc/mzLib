using MzLibUtil;
using UsefulProteomicsDatabases;

namespace Readers
{
    /// <summary>The public repository that holds a deposit. It names the deposit's record in a drafted cell's provenance.</summary>
    internal enum SdrfRepository
    {
        /// <summary>PRIDE: the record is a <see cref="PrideProject"/>.</summary>
        Pride,

        /// <summary>MassIVE: the record is its PROXI dataset record, which has no protocols, organism parts or diseases.</summary>
        Massive
    }

    /// <summary>
    /// What <see cref="SdrfDrafter"/> reads from a deposit's repository record, whichever repository holds it (sdrf design
    /// MASSIVE-SDRF.md, MV1). PRIDE's <see cref="PrideProject"/> fills every field (<see cref="FromPride"/>). MassIVE's
    /// record has no sample or data processing protocol, organism part or disease, so those stay empty and the drafter
    /// writes those cells <c>not available</c>.
    ///
    /// <para>A value taken from the record is written with the record's own provenance word (<see cref="Source"/>:
    /// <c>pride project record</c>, <c>massive dataset record</c>), and every reason the drafter gives names the record
    /// (<see cref="Name"/>), so a MassIVE draft never claims PRIDE as its source.</para>
    /// </summary>
    internal sealed record SdrfDepositRecord(SdrfRepository Repository)
    {
        public string Title { get; init; } = string.Empty;
        public string Description { get; init; } = string.Empty;
        public string SampleProcessingProtocol { get; init; } = string.Empty;
        public string DataProcessingProtocol { get; init; } = string.Empty;
        public IReadOnlyList<string> Keywords { get; init; } = Array.Empty<string>();
        public IReadOnlyList<string> ProjectTags { get; init; } = Array.Empty<string>();
        public IReadOnlyList<CvParam> Organisms { get; init; } = Array.Empty<CvParam>();
        public IReadOnlyList<CvParam> OrganismParts { get; init; } = Array.Empty<CvParam>();
        public IReadOnlyList<CvParam> Diseases { get; init; } = Array.Empty<CvParam>();
        public IReadOnlyList<CvParam> Instruments { get; init; } = Array.Empty<CvParam>();

        /// <summary>PRIDE's project record, as the drafter reads it. Throws on null.</summary>
        internal static SdrfDepositRecord FromPride(PrideProject project)
        {
            if (project == null) throw new ArgumentNullException(nameof(project));
            return new SdrfDepositRecord(SdrfRepository.Pride)
            {
                Title = project.Title,
                Description = project.ProjectDescription,
                SampleProcessingProtocol = project.SampleProcessingProtocol,
                DataProcessingProtocol = project.DataProcessingProtocol,
                Keywords = project.Keywords ?? new List<string>(),
                ProjectTags = project.ProjectTags ?? new List<string>(),
                Organisms = project.Organisms ?? new List<CvParam>(),
                OrganismParts = project.OrganismParts ?? new List<CvParam>(),
                Diseases = project.Diseases ?? new List<CvParam>(),
                Instruments = project.Instruments ?? new List<CvParam>(),
            };
        }

        /// <summary>How a reason names the record: <c>PRIDE's project record</c>, <c>MassIVE's dataset record</c>.</summary>
        internal string Name => Repository switch
        {
            SdrfRepository.Pride => "PRIDE's project record",
            SdrfRepository.Massive => "MassIVE's dataset record",
            _ => throw new ArgumentOutOfRangeException(nameof(Repository), Repository, "no record name for this repository")
        };

        /// <summary>The provenance of a value taken from the record.</summary>
        internal SdrfDraftSource Source => Repository switch
        {
            SdrfRepository.Pride => SdrfDraftSource.PrideProjectRecord,
            SdrfRepository.Massive => SdrfDraftSource.MassiveDatasetRecord,
            _ => throw new ArgumentOutOfRangeException(nameof(Repository), Repository, "no provenance for this repository")
        };

        /// <summary>Whether a value came from a repository record: a deposit-wide summary, which fills a gap but never
        /// disputes a per-file statement, and which a row's provenance default prefers on a tie.</summary>
        internal static bool IsRecord(SdrfDraftSource source) =>
            source is SdrfDraftSource.PrideProjectRecord or SdrfDraftSource.MassiveDatasetRecord;
    }
}
