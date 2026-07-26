using System.Globalization;
using Chemistry;
using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;
using Easy.Common.Extensions;

namespace Readers
{
    /// <summary>
    /// A single row of DIA-NN's main report (the long-format report.tsv written by DIA-NN 1.8/1.9).
    /// Each row is one precursor quantified in one run, so the same peptide appears once per run
    /// it was detected in.
    /// </summary>
    public class DiaNnPrecursor : IQuantifiableRecord
    {
        public static CsvConfiguration CsvConfiguration = new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
            HasHeaderRecord = true,
            IgnoreBlankLines = true,
            TrimOptions = TrimOptions.Trim,
            BadDataFound = null,
        };

        #region DIA-NN Fields

        /// <summary>
        /// Full path of the spectra file as DIA-NN saw it, which is frequently a path on a
        /// different machine. <see cref="Run"/> is the reliable identifier; use that instead.
        /// </summary>
        [Name("File.Name")]
        public string SpectraFilePath { get; set; }

        /// <summary>
        /// The spectra file name with no directory and no extension, e.g. "20240316_PBMCs_24min_3Th-1_centroid".
        /// </summary>
        [Name("Run")]
        public string Run { get; set; }

        [Name("Protein.Group")]
        public string ProteinGroup { get; set; }

        /// <summary>
        /// Semicolon-delimited UniProt accessions of every protein in the group.
        /// Runs parallel to <see cref="Genes"/> and <see cref="ProteinNames"/>.
        /// </summary>
        [Name("Protein.Ids")]
        public string ProteinIds { get; set; }

        /// <summary>
        /// Semicolon-delimited UniProt entry names, e.g. "TPM1_HUMAN;TPM2_HUMAN".
        /// </summary>
        [Name("Protein.Names")]
        public string ProteinNames { get; set; }

        /// <summary>
        /// Semicolon-delimited gene names. Blank when DIA-NN's database has no gene for the group.
        /// </summary>
        [Name("Genes")]
        public string Genes { get; set; }

        [Name("PG.Quantity")]
        public double ProteinGroupQuantity { get; set; }

        [Name("PG.Normalised")]
        public double ProteinGroupNormalized { get; set; }

        [Name("PG.MaxLFQ")]
        public double ProteinGroupMaxLfq { get; set; }

        // The Genes.* aggregates are blank for every row whose Genes column is blank,
        // hence nullable rather than defaulted to zero -- a blank is "not reported", not "zero".
        [Name("Genes.Quantity")]
        public double? GenesQuantity { get; set; }

        [Name("Genes.Normalised")]
        public double? GenesNormalized { get; set; }

        [Name("Genes.MaxLFQ")]
        public double? GenesMaxLfq { get; set; }

        [Name("Genes.MaxLFQ.Unique")]
        public double? GenesMaxLfqUnique { get; set; }

        /// <summary>
        /// Sequence with modifications in DIA-NN's UniMod accession notation, e.g. "AAAAC(UniMod:4)LDK".
        /// </summary>
        [Name("Modified.Sequence")]
        public string ModifiedSequence { get; set; }

        [Name("Stripped.Sequence")]
        public string BaseSequence { get; set; }

        /// <summary>
        /// The modified sequence with the charge appended, e.g. "AAAAC(UniMod:4)LDK2".
        /// This is DIA-NN's unique key for a precursor.
        /// </summary>
        [Name("Precursor.Id")]
        public string PrecursorId { get; set; }

        [Name("Precursor.Charge")]
        public int PrecursorCharge { get; set; }

        [Name("Q.Value")]
        public double QValue { get; set; }

        [Name("PEP")]
        public double PosteriorErrorProbability { get; set; }

        [Name("Global.Q.Value")]
        public double GlobalQValue { get; set; }

        [Name("Protein.Q.Value")]
        public double ProteinQValue { get; set; }

        [Name("PG.Q.Value")]
        public double ProteinGroupQValue { get; set; }

        [Name("Global.PG.Q.Value")]
        public double GlobalProteinGroupQValue { get; set; }

        [Name("GG.Q.Value")]
        public double GeneGroupQValue { get; set; }

        [Optional]
        [Name("Translated.Q.Value")]
        public double? TranslatedQValue { get; set; }

        [Name("Proteotypic")]
        [TypeConverter(typeof(IntegerBooleanConverter))]
        public bool Proteotypic { get; set; }

        [Name("Precursor.Quantity")]
        public double PrecursorQuantity { get; set; }

        [Name("Precursor.Normalised")]
        public double PrecursorNormalized { get; set; }

        [Optional]
        [Name("Precursor.Translated")]
        public double? PrecursorTranslated { get; set; }

        [Optional]
        [Name("Translated.Quality")]
        public double? TranslatedQuality { get; set; }

        [Optional]
        [Name("Ms1.Translated")]
        public double? Ms1Translated { get; set; }

        [Optional]
        [Name("Quantity.Quality")]
        public double? QuantityQuality { get; set; }

        /// <summary>
        /// Apex retention time in minutes. DIA-NN writes retention times in minutes already,
        /// so this needs no conversion to satisfy <see cref="IQuantifiableRecord.RetentionTime"/>.
        /// </summary>
        [Name("RT")]
        public double RetentionTime { get; set; }

        [Name("RT.Start")]
        public double RetentionTimeStart { get; set; }

        [Name("RT.Stop")]
        public double RetentionTimeStop { get; set; }

        [Name("iRT")]
        public double IndexedRetentionTime { get; set; }

        [Name("Predicted.RT")]
        public double PredictedRetentionTime { get; set; }

        [Name("Predicted.iRT")]
        public double PredictedIndexedRetentionTime { get; set; }

        [Name("Lib.Q.Value")]
        public double LibraryQValue { get; set; }

        [Name("Lib.PG.Q.Value")]
        public double LibraryProteinGroupQValue { get; set; }

        [Optional]
        [Name("Ms1.Profile.Corr")]
        public double? Ms1ProfileCorrelation { get; set; }

        [Optional]
        [Name("Ms1.Area")]
        public double? Ms1Area { get; set; }

        [Optional]
        [Name("Evidence")]
        public double? Evidence { get; set; }

        [Optional]
        [Name("Spectrum.Similarity")]
        public double? SpectrumSimilarity { get; set; }

        [Optional]
        [Name("Averagine")]
        public double? Averagine { get; set; }

        [Optional]
        [Name("Mass.Evidence")]
        public double? MassEvidence { get; set; }

        [Optional]
        [Name("CScore")]
        public double? CScore { get; set; }

        [Optional]
        [Name("Decoy.Evidence")]
        public double? DecoyEvidence { get; set; }

        [Optional]
        [Name("Decoy.CScore")]
        public double? DecoyCScore { get; set; }

        /// <summary>
        /// Semicolon-delimited per-fragment intensities, in the same order as <see cref="FragmentInfo"/>.
        /// Kept as the raw string so the file round-trips on write.
        /// </summary>
        [Optional]
        [Name("Fragment.Quant.Raw")]
        public string FragmentQuantRaw { get; set; }

        [Optional]
        [Name("Fragment.Quant.Corrected")]
        public string FragmentQuantCorrected { get; set; }

        [Optional]
        [Name("Fragment.Correlations")]
        public string FragmentCorrelations { get; set; }

        [Name("MS2.Scan")]
        public int Ms2ScanNumber { get; set; }

        [Name("Precursor.Mz")]
        public double PrecursorMz { get; set; }

        /// <summary>
        /// Semicolon-delimited fragment annotations and m/z values, e.g. "y7^1/601.3427124;b3^1/214.1191711;".
        /// </summary>
        [Optional]
        [Name("Fragment.Info")]
        public string FragmentInfo { get; set; }

        [Optional]
        [Name("Lib.Index")]
        public int LibraryIndex { get; set; }

        // The ion mobility columns are written by every DIA-NN version but are zero for
        // instruments without ion mobility, and are absent from some older reports.
        [Optional]
        [Name("IM")]
        public double? IonMobility { get; set; }

        [Optional]
        [Name("iIM")]
        public double? IndexedIonMobility { get; set; }

        [Optional]
        [Name("Predicted.IM")]
        public double? PredictedIonMobility { get; set; }

        [Optional]
        [Name("Predicted.iIM")]
        public double? PredictedIndexedIonMobility { get; set; }

        #endregion

        #region IQuantifiableRecord Implementation

        /// <summary>
        /// The Run column, which is the spectra file name with no directory and no extension.
        /// The File.Name column is not used because it records the path on whichever machine
        /// DIA-NN ran on, which is usually not where the spectra file sits now.
        /// </summary>
        [Ignore]
        public string FileName => Run.IsNullOrEmpty() ? StripPathAndExtension(SpectraFilePath) : Run;

        [Ignore]
        public string FullSequence => ModifiedSequence.IsNullOrEmpty() ? BaseSequence : ModifiedSequence;

        [Ignore]
        public int ChargeState => PrecursorCharge;

        /// <summary>
        /// DIA-NN's main report contains only target identifications; decoy scores are reported
        /// in the Decoy.CScore column rather than as separate rows.
        /// </summary>
        [Ignore]
        public bool IsDecoy => false;

        /// <summary>
        /// DIA-NN reports no peptide mass, so it is calculated from the precursor m/z and charge.
        /// </summary>
        [Ignore]
        public double MonoisotopicMass => PrecursorMz.ToMass(PrecursorCharge);

        [Ignore] private List<(string, string, string)> _proteinGroupInfos;

        [Ignore]
        public List<(string proteinAccessions, string geneName, string organism)> ProteinGroupInfos =>
            _proteinGroupInfos ??= AddProteinGroupInfos();

        /// <summary>
        /// Builds one tuple per protein in the group.
        /// <para>
        /// Protein.Ids, Genes, and Protein.Names look like parallel lists but are not. DIA-NN writes
        /// each as an independently sorted, de-duplicated set over the group, so position i of one
        /// does not correspond to position i of another, and zipping them by index silently invents
        /// wrong pairings. In a 9-run PBMC report that mispaired 1.5% of rows: for the group
        /// "P06753;P07951;P09493;P67936" with genes "TPM1;TPM2;TPM3;TPM4" it claims P06753 is TPM1,
        /// when P06753 is TPM3 and P09493 is TPM1. Three quarters of the mispaired groups have
        /// equal-length lists, so no length check can detect them.
        /// </para>
        /// <para>
        /// The report alone therefore cannot say which gene belongs to which accession, and a wrong
        /// gene is worse than none: FlashLFQ keys its protein groups on the accession and keeps
        /// whichever gene it saw first, so one mispaired row poisons that accession everywhere. A
        /// gene is assigned only when the group names exactly one, which is always true of a
        /// single-protein group. Organism is treated separately because it survives the ambiguity
        /// more often -- when every entry name in the group ends in the same suffix, that suffix is
        /// right no matter how the accessions pair up, so it stays populated for single-organism
        /// searches. Recovering per-accession genes for multi-protein groups would take the FASTA or
        /// spectral library that DIA-NN searched, which this reader does not have.
        /// </para>
        /// </summary>
        private List<(string, string, string)> AddProteinGroupInfos()
        {
            string[] accessions = SplitList(ProteinIds);

            // Fall back to the group identifier when Protein.Ids is empty so that a record always
            // yields at least one protein group -- FlashLFQ keys its ProteinGroup dictionary on this.
            if (accessions.Length == 0)
                accessions = SplitList(ProteinGroup);

            // Only a group naming exactly one accession and exactly one gene ties them together
            // beyond doubt. A group of two accessions reporting a single gene does not: that gene
            // belongs to one of them, and the report does not say which.
            string[] genes = SplitList(Genes);
            string gene = accessions.Length == 1 && genes.Length == 1 ? genes[0] : "";

            // Organism is coarse enough to survive the ambiguity: if every entry name the group
            // reports ends in the same suffix, that suffix holds however the accessions pair up,
            // which keeps it populated for the single-organism searches that dominate.
            string[] organisms = SplitList(ProteinNames)
                .Select(entryName =>
                {
                    int underscore = entryName.LastIndexOf('_');
                    return underscore >= 0 ? entryName.Substring(underscore + 1) : "";
                })
                .Distinct()
                .ToArray();
            string organism = organisms.Length == 1 ? organisms[0] : "";

            // A row naming no protein at all still needs one group for FlashLFQ to key on
            if (accessions.Length == 0)
                return new List<(string, string, string)> { ("", gene, organism) };

            return accessions.Select(accession => (accession, gene, organism)).ToList();
        }

        private static string[] SplitList(string value) =>
            value.IsNullOrEmpty() ? Array.Empty<string>() : value.Split(';', StringSplitOptions.RemoveEmptyEntries);

        /// <summary>
        /// Strips directory and extension from a path written by DIA-NN. Path.GetFileNameWithoutExtension
        /// is not used because DIA-NN commonly runs on Linux and writes '/'-delimited paths that
        /// .NET on Linux would not split.
        /// </summary>
        internal static string StripPathAndExtension(string path)
        {
            if (path.IsNullOrEmpty()) return "";

            int lastSeparator = path.LastIndexOfAny(new[] { '/', '\\' });
            string fileName = lastSeparator >= 0 ? path.Substring(lastSeparator + 1) : path;

            int lastDot = fileName.LastIndexOf('.');
            return lastDot > 0 ? fileName.Substring(0, lastDot) : fileName;
        }

        #endregion
    }
}
