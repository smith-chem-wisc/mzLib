using MassSpectrometry;
using MzLibUtil;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases;

namespace Readers
{
    /// <summary>
    /// Builds an SDRF-Proteomics document from mzLib-native inputs.
    ///
    /// Takes <see cref="Tolerance"/>, <see cref="Omics.Digestion.DigestionAgent"/>,
    /// <see cref="Modification"/>, <see cref="DissociationType"/> and <see cref="CvParam"/> —
    /// nothing application-specific — so a caller passes what it already holds rather than
    /// translating. That is what keeps the MetaMorpheus side of this a thin call rather than a
    /// second implementation.
    ///
    /// Two rules govern every value it writes:
    ///
    /// **Terms are resolved, never copied.** An accession comes from the pinned vocabulary
    /// (<see cref="ControlledVocabulary"/>), not from whatever an input document happened to say.
    /// The curated corpus is full of drift a copying builder would launder: bare "4" for UNIMOD:4
    /// in 39 documents, "Unimod:35" in 22, and two files that give one instrument a different
    /// accession on every row. Round-tripping someone else's file preserves their bytes; authoring
    /// a value is our decision.
    ///
    /// **Controlled vocabulary always, never free text.** Where a concept has a term, the term is
    /// written. Tolerating other people's free text is the pooling layer's job
    /// (<see cref="SdrfDriftLint"/>), not the writer's.
    ///
    /// What it will NOT do is invent a sample fact. If a caller cannot supply organism part or
    /// disease, the build fails rather than writing "not available" — see
    /// <see cref="SdrfBuilderOptions.RequireSampleMetadata"/>.
    /// </summary>
    public static class SdrfBuilder
    {
        // Column names, in the specification's block order: sample metadata, then data-file
        // metadata, then factor values.
        private const string SourceName = "source name";
        private const string Organism = "characteristics[organism]";
        private const string OrganismPart = "characteristics[organism part]";
        private const string BiologicalReplicate = "characteristics[biological replicate]";
        private const string AssayName = "assay name";
        private const string TechnologyType = "technology type";
        private const string AcquisitionMethod = "comment[proteomics data acquisition method]";
        private const string Label = "comment[label]";
        private const string Instrument = "comment[instrument]";
        private const string CleavageAgent = "comment[cleavage agent details]";
        private const string ModificationParameters = "comment[modification parameters]";
        private const string PrecursorTolerance = "comment[precursor mass tolerance]";
        private const string FragmentTolerance = "comment[fragment mass tolerance]";
        private const string DissociationMethod = "comment[dissociation method]";
        private const string FractionIdentifier = "comment[fraction identifier]";
        private const string TechnicalReplicate = "comment[technical replicate]";
        private const string DataFile = "comment[data file]";
        private const string PxAccession = "comment[proteomexchange accession number]";
        private const string SoftwareColumn = "comment[software]";
        private const string SdrfVersionColumn = "comment[sdrf version]";

        /// <summary>The one value SDRF defines for this column in an MS experiment.</summary>
        private const string TechnologyTypeValue = "proteomic profiling by mass spectrometry";

        /// <summary>
        /// Characteristics columns SDRF requires of every document, and which are therefore emitted
        /// whether or not the caller supplied a value.
        ///
        /// Every other column this builder writes is unconditional, so a caller cannot omit one by
        /// accident. These are the exception, because they arrive through
        /// <see cref="SdrfSample.Characteristics"/> — an open dictionary — and a caller that simply
        /// does not populate it produces a document with the column MISSING rather than empty.
        ///
        /// That distinction is the whole reason for this list. A missing column is invisible to
        /// every instrument built to detect thin metadata: <see cref="SdrfCoverage"/> measures the
        /// fill rate of columns that EXIST, so an absent column has no fill rate to report at all.
        /// An empty column is visible, honest, and permitted — the specification marks organism part
        /// required but allows the reserved words, and a fifth of the curated corpus uses one.
        ///
        /// So the rule is: always emit the column; let <see cref="Missing"/> decide what goes in it.
        /// </summary>
        private static readonly string[] RequiredCharacteristics = { OrganismPart };

        public static SdrfDocument Build(IEnumerable<SdrfRowInput> rows, SdrfBuilderOptions options = null)
        {
            if (rows is null) throw new ArgumentNullException(nameof(rows));
            options ??= new SdrfBuilderOptions();

            var inputs = rows.ToList();
            if (inputs.Count == 0)
                throw new ArgumentException("An SDRF document needs at least one row.", nameof(rows));

            // Build null-checks its argument, so the members it goes on to dereference are checked
            // here too. Without this, a null Assay or a null Characteristics dictionary surfaces as
            // a NullReferenceException thrown from inside a LINQ lambda, which tells the caller
            // neither which row was wrong nor what was missing -- and the entry point has already
            // established ArgumentException as the contract for a malformed input.
            for (int i = 0; i < inputs.Count; i++)
                RequireComplete(inputs[i], i);

            // Multi-cardinality columns are as wide as the widest row needs, and every row pads to
            // that width. Sizing per row would produce a ragged document.
            int modificationSlots = Math.Max(1, inputs.Max(r =>
                r.Assay.FixedModifications.Count + r.Assay.VariableModifications.Count));

            // Union with the required set, so a caller who supplied no characteristics at all still
            // gets a spec-conformant header. See RequiredCharacteristics for why absent is worse
            // than empty.
            var characteristicColumns = inputs
                .SelectMany(r => r.Sample.Characteristics.Keys)
                .Concat(RequiredCharacteristics)
                .Distinct(StringComparer.Ordinal)
                .OrderBy(c => c, StringComparer.Ordinal)
                .ToList();

            var factorColumns = inputs
                .Select(r => r.Sample.FactorValueColumn)
                .Where(c => !string.IsNullOrWhiteSpace(c))
                .Distinct(StringComparer.Ordinal)
                .OrderBy(c => c, StringComparer.Ordinal)
                .ToList();

            var header = new SdrfHeader(BuildHeader(
                characteristicColumns, factorColumns, modificationSlots, options));

            var built = inputs
                .Select(input => new SdrfRow(header,
                    BuildCells(input, characteristicColumns, factorColumns, modificationSlots, options)))
                .ToList();

            return new SdrfDocument(header, built);
        }

        /// <summary>
        /// Checks the one row for the nulls <see cref="Build"/> would otherwise dereference. The
        /// collection properties default to empty, so reaching this with a null means the caller set
        /// one to null explicitly; an empty list or dictionary is what it wants instead.
        /// </summary>
        private static void RequireComplete(SdrfRowInput input, int index)
        {
            if (input is null)
                throw new ArgumentException($"Row {index} is null.", nameof(input));
            if (input.Sample is null)
                throw new ArgumentException($"Row {index} has no sample.", nameof(input));
            if (input.Assay is null)
                throw new ArgumentException($"Row {index} has no assay.", nameof(input));
            if (input.Sample.Characteristics is null)
                throw new ArgumentException(
                    $"Row {index} has a null {nameof(SdrfSample.Characteristics)}; pass an empty " +
                    "dictionary for a sample with no characteristics.", nameof(input));
            if (input.Assay.FixedModifications is null)
                throw new ArgumentException(
                    $"Row {index} has a null {nameof(SdrfAssay.FixedModifications)}; pass an empty " +
                    "list for a search with none.", nameof(input));
            if (input.Assay.VariableModifications is null)
                throw new ArgumentException(
                    $"Row {index} has a null {nameof(SdrfAssay.VariableModifications)}; pass an " +
                    "empty list for a search with none.", nameof(input));
        }

        private static List<string> BuildHeader(
            IReadOnlyList<string> characteristics, IReadOnlyList<string> factors,
            int modificationSlots, SdrfBuilderOptions options)
        {
            var names = new List<string> { SourceName, Organism };
            names.AddRange(characteristics);
            names.Add(BiologicalReplicate);

            names.Add(AssayName);
            names.Add(TechnologyType);
            names.Add(AcquisitionMethod);
            names.Add(Label);
            names.Add(Instrument);
            names.Add(CleavageAgent);
            for (int i = 0; i < modificationSlots; i++) names.Add(ModificationParameters);
            names.Add(PrecursorTolerance);
            names.Add(FragmentTolerance);
            names.Add(DissociationMethod);
            names.Add(FractionIdentifier);
            names.Add(TechnicalReplicate);
            names.Add(DataFile);

            if (!string.IsNullOrWhiteSpace(options.ProteomeXchangeAccession)) names.Add(PxAccession);
            if (options.Software is not null) names.Add(SoftwareColumn);
            if (!string.IsNullOrWhiteSpace(options.SdrfVersion)) names.Add(SdrfVersionColumn);

            names.AddRange(factors);
            return names;
        }

        private static List<string> BuildCells(
            SdrfRowInput input, IReadOnlyList<string> characteristics, IReadOnlyList<string> factors,
            int modificationSlots, SdrfBuilderOptions options)
        {
            var sample = input.Sample;
            var assay = input.Assay;
            var cells = new List<string>
            {
                Required(sample.SourceName, SourceName, options),
                Term(sample.Organism, Organism, options)
            };

            foreach (var column in characteristics)
                cells.Add(sample.Characteristics.TryGetValue(column, out var value)
                    ? Term(value, column, options)
                    : Missing(column, options));

            cells.Add(Positive(sample.BiologicalReplicate, BiologicalReplicate));

            cells.Add(Required(assay.AssayName, AssayName, options));
            cells.Add(TechnologyTypeValue);
            cells.Add(Term(assay.AcquisitionMethod, AcquisitionMethod, options));
            cells.Add(Term(sample.Label, Label, options));
            cells.Add(InstrumentCell(assay.Instrument, options));
            cells.Add(CleavageAgentCell(assay.CleavageAgent, options));

            var modifications = ModificationCells(assay, options).ToList();
            for (int i = 0; i < modificationSlots; i++)
                cells.Add(i < modifications.Count ? modifications[i] : SdrfReserved.NotApplicable);

            cells.Add(ToleranceCell(assay.PrecursorMassTolerance, PrecursorTolerance, options));
            cells.Add(ToleranceCell(assay.ProductMassTolerance, FragmentTolerance, options));
            cells.Add(DissociationCell(assay.DissociationType, options));
            cells.Add(Positive(assay.Fraction, FractionIdentifier));
            cells.Add(Positive(assay.TechnicalReplicate, TechnicalReplicate));
            cells.Add(Required(assay.DataFileName, DataFile, options));

            if (!string.IsNullOrWhiteSpace(options.ProteomeXchangeAccession))
                cells.Add(options.ProteomeXchangeAccession);
            if (options.Software is not null)
                cells.Add(SoftwareCell(options));
            if (!string.IsNullOrWhiteSpace(options.SdrfVersion))
                // Either casing of the prefix is stripped before ours is added. Stripping only 'v'
                // turned a caller's "V1.1.0" into "vV1.1.0", and nothing downstream validates this
                // column, so the malformed value would have gone out silently. The specification
                // asks for vMAJOR.MINOR.PATCH.
                cells.Add("v" + options.SdrfVersion.TrimStart('v', 'V'));

            foreach (var column in factors)
                cells.Add(string.Equals(sample.FactorValueColumn, column, StringComparison.Ordinal)
                          && !string.IsNullOrWhiteSpace(sample.FactorValue)
                    ? sample.FactorValue
                    : SdrfReserved.NotApplicable);

            return cells;
        }

        /// <summary>
        /// The instrument. mzML supplies it already accessioned; a Thermo RAW supplies a name with
        /// an empty accession, which is resolved here against PSI-MS — the one place a name-to-
        /// accession lookup genuinely belongs, since the raw reader must not carry an ontology.
        /// </summary>
        private static string InstrumentCell(CvParam instrument, SdrfBuilderOptions options)
        {
            if (instrument is null) return Missing(Instrument, options);

            if (string.IsNullOrEmpty(instrument.Accession)
                && !string.IsNullOrEmpty(instrument.Name)
                && ControlledVocabulary.PsiMs.TryGetByName(instrument.Name, out var resolved))
                instrument = resolved;

            // Neither a name nor an accession is not a term, and SdrfCell.ToCell rightly refuses
            // it. That refusal must not be how a lenient build ends: RequireSampleMetadata = false
            // exists precisely so a finished search is never thrown away, and an ArgumentException
            // raised from a formatting helper would defeat it. Route it through Missing, which is
            // the one place that decides what an absent value costs.
            if (string.IsNullOrEmpty(instrument.Accession) && string.IsNullOrEmpty(instrument.Name))
                return Missing(Instrument, options);

            return string.IsNullOrEmpty(instrument.Accession)
                // Named but unresolvable. The specification allows NT= with no AC=, and that is
                // honest: it says which instrument without claiming a term we could not find.
                ? SdrfCell.ToCell(new CvParam("", "", instrument.Name, ""))
                : SdrfCell.ToCell(instrument);
        }

        private static string CleavageAgentCell(Omics.Digestion.DigestionAgent agent, SdrfBuilderOptions options)
        {
            if (agent is null) return Missing(CleavageAgent, options);

            // Protease carries its own PSI-MS accession, from the embedded proteases.tsv. Nine
            // MetaMorpheus proteases have none; the specification permits NT= alone, and inventing
            // an accession for them would be worse than omitting one.
            string accession = agent is Protease protease ? protease.PsiMsAccessionNumber : "";
            string name = agent is Protease named && !string.IsNullOrEmpty(named.PsiMsName)
                ? named.PsiMsName
                : agent.Name;

            // Same reason as InstrumentCell: an agent that names nothing cannot be written as a
            // term, and a lenient build must get a reserved word rather than an exception.
            if (string.IsNullOrEmpty(accession) && string.IsNullOrEmpty(name))
                return Missing(CleavageAgent, options);

            return SdrfCell.ToCell(new CvParam(
                string.IsNullOrEmpty(accession) ? "" : "MS", accession ?? "", name, ""));
        }

        /// <summary>
        /// One cell per modification, fixed first. The accession comes from the modification's own
        /// UNIMOD database reference — mzLib already resolves it when loading — and the target
        /// residues from its motif.
        /// </summary>
        private static IEnumerable<string> ModificationCells(SdrfAssay assay, SdrfBuilderOptions options)
        {
            foreach (var mod in assay.FixedModifications) yield return ModificationCell(mod, "Fixed", options);
            foreach (var mod in assay.VariableModifications) yield return ModificationCell(mod, "Variable", options);
        }

        /// <summary>
        /// NT= carries the BARE modification name and TA= carries the residue, which is what the
        /// specification's own examples show — "NT=Oxidation;AC=UNIMOD:35;TA=M;MT=Variable".
        ///
        /// So the name comes from <see cref="Modification.OriginalId"/>, not
        /// <see cref="Modification.IdWithMotif"/>. mzLib's Modification constructor sets
        /// IdWithMotif to "Oxidation on M" whenever a target is present (Modification.cs:78), so
        /// taking the name from there wrote the motif twice, and once in a field that is not for it.
        /// OriginalId is the same choice ProFormaConverter.cs:159 makes when it emits a modification
        /// name into an external format.
        /// </summary>
        private static string ModificationCell(Modification mod, string modificationType, SdrfBuilderOptions options)
        {
            string accession = "";
            if (mod.DatabaseReference is not null
                && mod.DatabaseReference.TryGetValue("Unimod", out var unimod)
                && unimod.Count > 0)
                accession = "UNIMOD:" + unimod.First();

            var extras = new List<(string, string)>();
            if (mod.Target is not null) extras.Add(("TA", mod.Target.ToString()));
            extras.Add(("MT", modificationType));

            // IdWithMotif is only ever populated when OriginalId is, so this fallback order is the
            // safe one; the reverse could never reach OriginalId.
            string name = mod.OriginalId ?? mod.IdWithMotif ?? "";

            // As in InstrumentCell: a modification naming nothing and accessioned to nothing cannot
            // be written as a term, and must not throw out of a lenient build.
            if (string.IsNullOrEmpty(name) && string.IsNullOrEmpty(accession))
                return Missing(ModificationParameters, options);

            return SdrfCell.ToCell(
                new CvParam(string.IsNullOrEmpty(accession) ? "" : "UNIMOD", accession, name, ""),
                extras.ToArray());
        }

        /// <summary>
        /// Tolerances are written as a value and a unit ("10 ppm", "0.02 Da"), which is what the
        /// corpus does and what the specification's examples show — not as the MS:1001412/1001413
        /// cvParam pair the mzIdentML writer uses.
        ///
        /// INVARIANT CULTURE, because this is a number in a file format, not a number shown to a
        /// person: on a de-DE machine the default would write "0,02 Da", which no consumer parses.
        /// PpmTolerance and AbsoluteTolerance already format their own values this way
        /// (PpmTolerance.cs:38, AbsoluteTolerance.cs:38), as does <see cref="Positive"/> below.
        /// Their ToString is not reused here only because it renders "±10.0000 PPM", which is not
        /// the grammar SDRF asks for.
        /// </summary>
        private static string ToleranceCell(Tolerance tolerance, string column, SdrfBuilderOptions options)
        {
            if (tolerance is null) return Missing(column, options);
            string value = tolerance.Value.ToString(System.Globalization.CultureInfo.InvariantCulture);
            return tolerance is PpmTolerance ? value + " ppm" : value + " Da";
        }

        private static string DissociationCell(DissociationType dissociation, SdrfBuilderOptions options)
        {
            if (DissociationTypeCvTerms.TryGetTerm(dissociation, out var term))
                return SdrfCell.ToCell(term);

            // PSI-MS has no term for aEPD, and Custom/Autodetect/AnyActivationType do not name a
            // method at all. Writing a near neighbour would silently join this run to the wrong
            // population, so say nothing instead.
            return Missing(DissociationMethod, options);
        }

        private static string SoftwareCell(SdrfBuilderOptions options)
        {
            var software = options.Software;
            return string.IsNullOrWhiteSpace(options.SoftwareVersion)
                ? SdrfCell.ToCell(software)
                : SdrfCell.ToCell(software, ("VV", options.SoftwareVersion));
        }

        private static string Term(CvParam term, string column, SdrfBuilderOptions options) =>
            term is null ? Missing(column, options) : SdrfCell.ToCell(term);

        private static string Required(string value, string column, SdrfBuilderOptions options) =>
            string.IsNullOrWhiteSpace(value) ? Missing(column, options) : value;

        /// <summary>
        /// 1-based, and validated. SDRF replicate and fraction identifiers start at 1; mzLib's
        /// SpectraFileInfo stores them 0-based, so a caller that forwards those directly would write
        /// a 0 that every consumer reads as an error.
        /// </summary>
        private static string Positive(int value, string column) =>
            value >= 1
                ? value.ToString(System.Globalization.CultureInfo.InvariantCulture)
                : throw new MzLibException(
                    $"'{column}' must be 1-based and at least 1, but was {value}. mzLib's " +
                    "SpectraFileInfo stores replicates and fractions 0-based; add 1 before passing them.");

        /// <summary>
        /// What to do when a value is not available.
        ///
        /// Under the default (<see cref="SdrfBuilderOptions.RequireSampleMetadata"/>) this THROWS.
        /// Writing a reserved word would produce a document that passes every validator, generates
        /// no drift findings, and says nothing — the failure mode this whole design exists to
        /// prevent. A caller that genuinely wants a partial document opts out explicitly.
        ///
        /// NOTE that the switch is named for the columns it exists for, but it guards EVERY column
        /// routed through here — comment[label] and comment[dissociation method] among them. So the
        /// message names the column that has no value rather than asserting that sample metadata is
        /// what is missing, which for those two would be untrue. Splitting the switch in two would
        /// be a design change rather than a wording fix, and no caller needs the halves separately
        /// yet.
        /// </summary>
        private static string Missing(string column, SdrfBuilderOptions options) =>
            options.RequireSampleMetadata
                ? throw new MzLibException(
                    $"No value for '{column}'. An SDRF that writes a reserved word wherever it " +
                    "cannot state a fact is uniformly consistent, passes validation, and cannot be " +
                    "mined. Supply the value, or set RequireSampleMetadata = false to accept a " +
                    "document that says nothing here.")
                : SdrfReserved.NotAvailable;
    }

    /// <summary>
    /// The SDRF reserved words, which the specification requires lowercase.
    ///
    /// Internal, per D19: nothing outside this assembly calls it. The builder is the only thing that
    /// writes a reserved word, and a consumer reading one back compares against the literal text the
    /// specification fixes, not against a constant this library owns.
    /// </summary>
    internal static class SdrfReserved
    {
        public const string NotAvailable = "not available";
        public const string NotApplicable = "not applicable";
        public const string Anonymized = "anonymized";
        public const string Pooled = "pooled";
    }
}
