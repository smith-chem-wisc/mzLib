using MassSpectrometry;
using MzLibUtil;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;

namespace Readers
{
    /// <summary>
    /// The sample facts for one row — the half of an SDRF that no search can know.
    ///
    /// Nothing in a search knows what tissue was taken or from what condition, so these values come
    /// from a human: an SDRF sitting beside the spectra, or whatever the caller collected. That is
    /// why SDRF output is opt-in — opting in is the caller asserting it has these.
    ///
    /// <see cref="Characteristics"/> deliberately takes already-resolved <see cref="CvParam"/>s
    /// rather than free text. Resolving "liver" to an accession needs UBERON, which mzLib does not
    /// embed and should not; whoever supplies the sample metadata supplies the terms with it.
    /// </summary>
    public sealed record SdrfSample
    {
        /// <summary>The sample identifier. One source name may appear on several rows (one per label).</summary>
        public required string SourceName { get; init; }

        /// <summary>
        /// The organism, from the search database rather than from an ontology lookup: a UniProt
        /// database states it, and mzLib retains it as <c>Protein.NcbiTaxonomyId</c> plus
        /// <c>Protein.Organism</c>. Null only when the database supplied neither.
        /// </summary>
        public CvParam Organism { get; init; }

        /// <summary>
        /// Sample characteristics keyed by SDRF column name, e.g. "characteristics[organism part]".
        /// Values are terms, not free text (D13).
        /// </summary>
        public IReadOnlyDictionary<string, CvParam> Characteristics { get; init; }
            = new Dictionary<string, CvParam>();

        /// <summary>1-based, as SDRF writes them. Note mzLib's SpectraFileInfo stores these 0-based.</summary>
        public int BiologicalReplicate { get; init; } = 1;

        /// <summary>The label for this row: "label free sample", or a TMT channel.</summary>
        public CvParam Label { get; init; }

        /// <summary>Free-text factor value, or null. The experimental variable under study.</summary>
        public string FactorValue { get; init; }

        /// <summary>The SDRF column the factor value belongs under, e.g. "factor value[disease]".</summary>
        public string FactorValueColumn { get; init; }
    }

    /// <summary>
    /// The assay facts for one data file — the half a search knows exactly.
    ///
    /// PER FILE, not per run, because MetaMorpheus supports per-file parameter overrides and using
    /// the task-level values would silently flatten real differences.
    /// </summary>
    public sealed record SdrfAssay
    {
        /// <summary>
        /// The ORIGINAL raw file name. Not a calibrated or averaged derivative: those are
        /// intermediates this experiment produced, and an SDRF describes the data as deposited.
        /// </summary>
        public required string DataFileName { get; init; }

        /// <summary>The run identifier. Unique per file within the document.</summary>
        public required string AssayName { get; init; }

        /// <summary>
        /// Where the instrument comes from. mzML carries it already accessioned; a Thermo RAW gives
        /// a name with an empty accession, which the builder resolves against PSI-MS.
        /// </summary>
        public CvParam Instrument { get; init; }

        public Tolerance PrecursorMassTolerance { get; init; }
        public Tolerance ProductMassTolerance { get; init; }

        /// <summary>
        /// The cleavage agent. A <see cref="Protease"/> carries its own PSI-MS accession; any other
        /// <see cref="DigestionAgent"/> contributes a name with no accession, which the SDRF
        /// specification permits.
        /// </summary>
        public DigestionAgent CleavageAgent { get; init; }

        public IReadOnlyList<Modification> FixedModifications { get; init; } = new List<Modification>();
        public IReadOnlyList<Modification> VariableModifications { get; init; } = new List<Modification>();

        public DissociationType DissociationType { get; init; } = DissociationType.Unknown;

        /// <summary>DDA/DIA/PRM/SRM, as a PRIDE CV term. The corpus uses PRIDE here, not PSI-MS.</summary>
        public CvParam AcquisitionMethod { get; init; }

        /// <summary>1-based, as SDRF writes them.</summary>
        public int TechnicalReplicate { get; init; } = 1;

        /// <summary>1-based. 1 when the sample was not fractionated.</summary>
        public int Fraction { get; init; } = 1;
    }

    /// <summary>One row: a sample paired with the file it was acquired into.</summary>
    public sealed record SdrfRowInput(SdrfSample Sample, SdrfAssay Assay);

    /// <summary>
    /// Document-level facts and switches.
    /// </summary>
    public sealed record SdrfBuilderOptions
    {
        /// <summary>
        /// The accession of the dataset this experiment re-analysed, if any. Written to
        /// comment[proteomexchange accession number], and it is what ties our assay parameters back
        /// to somebody else's samples when the two halves come from different places.
        /// </summary>
        public string ProteomeXchangeAccession { get; init; }

        /// <summary>
        /// The software that produced the file, e.g. MetaMorpheus. MS:1002826 is MetaMorpheus's own
        /// PSI-MS accession. Null omits the column.
        /// </summary>
        public CvParam Software { get; init; }

        /// <summary>Version string of that software, recorded alongside it.</summary>
        public string SoftwareVersion { get; init; }

        /// <summary>
        /// Stamped into comment[sdrf version]. Defaults to the specification version this builder
        /// targets.
        /// </summary>
        public string SdrfVersion { get; init; } = "1.1.0";

        /// <summary>
        /// When true (the default), a sample fact that cannot be written as a controlled-vocabulary
        /// term throws rather than being silently replaced with "not available".
        ///
        /// This is the whole point of the opt-in design. A corpus quietly padded with reserved words
        /// is uniformly consistent, passes every validator, produces no drift findings, and cannot
        /// answer a single question — see <see cref="SdrfCoverage"/>. Failing loudly at build time,
        /// before a search starts, is the alternative.
        /// </summary>
        public bool RequireSampleMetadata { get; init; } = true;
    }
}
