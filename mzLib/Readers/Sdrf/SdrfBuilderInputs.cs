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

        /// <summary>
        /// The file the search actually read, when that is not <see cref="DataFileName"/> -- a
        /// calibrated or averaged derivative such as <c>X-calib.mzML</c>. Written to
        /// <c>comment[searched data file]</c> beside <c>comment[data file]</c>, so a row names both
        /// the acquisition it describes and the transformation its results came from: the results
        /// cannot be regenerated from the raw file without re-running that step, and naming only the
        /// raw file would hide it.
        ///
        /// Not a specification term. The spec's <c>comment[associated data file]</c> means a sidecar
        /// needed to READ the primary file (a <c>.wiff.scan</c>), and <c>comment[processed data
        /// file]</c> exists only in the metabolomics template, so borrowing either would mislead.
        /// Null when the search read the acquired file; the column is then omitted unless another row
        /// sets it.
        ///
        /// <para><b>THE JOIN RULE, for anyone reading an SDRF back into an experimental design.</b>
        /// An experimental design is keyed on the file the SEARCH read, and this column is the only
        /// place that name appears -- <c>comment[data file]</c> is the acquisition. So: key on
        /// <c>comment[searched data file]</c> when the document has it, fall back to
        /// <c>comment[data file]</c> when it does not, and match on the STEM, case-insensitively,
        /// because a converted <c>.mzML</c> matches no <c>.raw</c> by full name.</para>
        ///
        /// <para>Keying on <c>comment[data file]</c> instead fails on every CALIBRATED run, which is
        /// most of them -- and it fails by producing a design with zero matched files, which reads
        /// like a bad fixture rather than a mapping bug. That is why the rule is written here, beside
        /// the thing that creates the mismatch, rather than left for each consumer to rediscover.</para>
        ///
        /// <para><b>And a caveat for anyone MINING these documents.</b> The column is gated
        /// document-wide: one row that names a derivative gives every row the column, so a row where
        /// searched equals acquired asserts "no transformation" in one document and says nothing at
        /// all in another. Never read the column's ABSENCE as "nothing was transformed anywhere" --
        /// it means only that no row in that document set it.</para>
        /// </summary>
        public string SearchedDataFileName { get; init; }

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

        /// <summary>
        /// How <c>comment[label]</c> is written: <see cref="SdrfLabelForm.Accessioned"/> (the default)
        /// as <c>NT=TMT127N;AC=PRIDE:0000519</c>, or <see cref="SdrfLabelForm.Bare"/> as the term's
        /// name alone, <c>TMT127N</c>.
        ///
        /// Bare is a deliberate exception to "terms are resolved and written as terms", and the only
        /// one. The accessioned form is more precise, but about six in seven published SDRFs write the
        /// label bare, and sdrf-pipelines 0.1.6 mishandles the accessioned form: its OpenMS converter
        /// crashes on the second channel of every TMT file (bigbio/sdrf-pipelines#344). A caller whose
        /// output feeds quantms needs Bare until that is fixed.
        ///
        /// The label is still RESOLVED either way; only its spelling changes. A label with no name to
        /// write bare keeps its accessioned form rather than losing the accession.
        /// </summary>
        public SdrfLabelForm LabelForm { get; init; } = SdrfLabelForm.Accessioned;
    }
}
