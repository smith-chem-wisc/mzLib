using System.Globalization;
using Omics;
using Omics.SpectrumMatch;

namespace Readers
{
    /// <summary>
    /// A lightweight PSM record that implements only <see cref="Omics.SpectralMatch.ISpectralMatch"/> and <see cref="IQuantifiableRecord"/>.
    /// Parses only the ~15 columns required by those interfaces plus QValue and PEP_QValue,
    /// making it significantly faster to construct than <see cref="SpectrumMatchFromTsv"/>.
    /// </summary>
    public class LightWeightSpectralMatch : Omics.SpectralMatch.ISpectralMatch, IQuantifiableRecord
    {
        private static readonly string[] AcceptedSpectraFormats = SpectrumMatchFromTsvHeader.AcceptedSpectraFormats;

        #region Properties

        // Shared by both interfaces
        public string FullSequence { get; }
        public string BaseSequence { get; }
        public string Accession { get; }
        public bool IsDecoy { get; }

        // ISpectralMatch
        public string FullFilePath { get; }
        public int OneBasedScanNumber { get; }
        public double Score { get; }
        public double[]? Intensities { get; }

        // IQuantifiableRecord
        public string FileName { get; }
        public double RetentionTime { get; }
        public int ChargeState { get; }
        public double MonoisotopicMass { get; }

        // Additional fields for filtering / downstream use
        public double QValue { get; }
        public double PepQValue { get; }

        // Backing fields for lazy ProteinGroupInfos
        private readonly string _geneName;
        private readonly string _organismName;
        private List<(string proteinAccessions, string geneName, string organism)>? _proteinGroupInfos;

        public List<(string proteinAccessions, string geneName, string organism)> ProteinGroupInfos
        {
            get
            {
                _proteinGroupInfos ??= ConstructProteinGroupInfo();
                return _proteinGroupInfos;
            }
        }

        #endregion

        /// <summary>
        /// Constructs a lightweight spectral match from a single TSV line.
        /// Only parses the columns required by <see cref="ISpectralMatch"/> and <see cref="IQuantifiableRecord"/>.
        /// </summary>
        public LightWeightSpectralMatch(string[] spl, Dictionary<string, int> parsedHeader)
        {
            // File name — strip known spectra file extensions
            string rawFileName = spl[parsedHeader[SpectrumMatchFromTsvHeader.FileName]].Trim();
            if (rawFileName.Contains('.'))
            {
                foreach (var ext in AcceptedSpectraFormats)
                {
                    rawFileName = Path.GetFileName(rawFileName.Replace(ext, string.Empty, StringComparison.InvariantCultureIgnoreCase));
                }
            }
            FileName = rawFileName;
            FullFilePath = rawFileName;

            // Scan number
            OneBasedScanNumber = int.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.Ms2ScanNumber]]);

            // Sequences
            FullSequence = spl[parsedHeader[SpectrumMatchFromTsvHeader.FullSequence]];
            BaseSequence = SpectrumMatchFromTsv.RemoveParentheses(spl[parsedHeader[SpectrumMatchFromTsvHeader.BaseSequence]].Trim());

            // Score
            Score = double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.Score]].Trim(), CultureInfo.InvariantCulture);

            // Charge
            ChargeState = (int)double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.PrecursorCharge]].Trim(), CultureInfo.InvariantCulture);

            // Decoy/Contaminant/Target
            string dct = spl[parsedHeader[SpectrumMatchFromTsvHeader.DecoyContaminantTarget]].Trim();
            IsDecoy = dct.Contains('D');

            // Accession (already resolved to canonical key by ParseHeader)
            Accession = parsedHeader[SpectrumMatchFromTsvHeader.Accession] >= 0
                ? spl[parsedHeader[SpectrumMatchFromTsvHeader.Accession]].Trim()
                : string.Empty;

            // Gene name and organism (for ProteinGroupInfos)
            _geneName = parsedHeader[SpectrumMatchFromTsvHeader.GeneName] >= 0
                ? spl[parsedHeader[SpectrumMatchFromTsvHeader.GeneName]].Trim()
                : string.Empty;
            _organismName = parsedHeader[SpectrumMatchFromTsvHeader.OrganismName] >= 0
                ? spl[parsedHeader[SpectrumMatchFromTsvHeader.OrganismName]].Trim()
                : string.Empty;

            // Monoisotopic mass
            string monoMassStr = parsedHeader[SpectrumMatchFromTsvHeader.MonoisotopicMass] >= 0
                ? spl[parsedHeader[SpectrumMatchFromTsvHeader.MonoisotopicMass]].Trim()
                : string.Empty;
            MonoisotopicMass = double.TryParse(monoMassStr.Split('|')[0], CultureInfo.InvariantCulture, out double monoMass)
                ? monoMass
                : -1;

            // Retention time
            RetentionTime = parsedHeader[SpectrumMatchFromTsvHeader.Ms2ScanRetentionTime] >= 0
                && double.TryParse(spl[parsedHeader[SpectrumMatchFromTsvHeader.Ms2ScanRetentionTime]].Trim(), CultureInfo.InvariantCulture, out double rt)
                ? rt
                : -1;

            // QValue
            QValue = parsedHeader[SpectrumMatchFromTsvHeader.QValue] >= 0
                ? double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.QValue]].Trim(), CultureInfo.InvariantCulture)
                : double.NaN;

            // PEP_QValue
            PepQValue = parsedHeader[SpectrumMatchFromTsvHeader.PEP_QValue] >= 0
                ? double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.PEP_QValue]].Trim(), CultureInfo.InvariantCulture)
                : double.NaN;

            // TMT/isobaric reporter ion columns
            Intensities = ParseReporterIonColumns(spl, parsedHeader);
        }

        #region Interface Methods

        public IEnumerable<IBioPolymerWithSetMods> GetIdentifiedBioPolymersWithSetMods()
        {
            throw new NotSupportedException(
                $"{nameof(LightWeightSpectralMatch)} does not support biopolymer resolution. " +
                $"Use {nameof(SpectrumMatchFromTsv)} for full parsing.");
        }

        public int CompareTo(Omics.SpectralMatch.ISpectralMatch? other)
        {
            if (other is null) return 1;
            // Higher score is better — sort descending
            return other.Score.CompareTo(Score);
        }

        public bool Equals(Omics.SpectralMatch.ISpectralMatch? other)
        {
            if (other is null) return false;
            if (ReferenceEquals(this, other)) return true;

            return string.Equals(FullFilePath, other.FullFilePath, StringComparison.Ordinal)
                && OneBasedScanNumber == other.OneBasedScanNumber
                && string.Equals(FullSequence, other.FullSequence, StringComparison.Ordinal);
        }

        public override bool Equals(object? obj)
        {
            if (obj is Omics.SpectralMatch.ISpectralMatch sm) return Equals(sm);
            return false;
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(FullFilePath, OneBasedScanNumber, FullSequence);
        }

        #endregion

        #region Private Helpers

        private List<(string proteinAccessions, string geneName, string organism)> ConstructProteinGroupInfo()
        {
            string[] accessions = Accession.Split('|');
            string[] genes = _geneName.Split('|');
            string[] organisms = _organismName.Split('|');
            var list = new List<(string, string, string)>(accessions.Length);
            for (int i = 0; i < accessions.Length; i++)
            {
                var gene = genes.Length > i ? genes[i] : genes[0];
                var organism = organisms.Length > i ? organisms[i] : organisms[0];
                list.Add((accessions[i], gene, organism));
            }
            return list;
        }

        private static double[]? ParseReporterIonColumns(string[] spl, Dictionary<string, int> parsedHeader)
        {
            var presentChannels = new List<int>();
            foreach (var channelName in SpectrumMatchFromTsvHeader.TmtChannelNames)
            {
                if (parsedHeader.TryGetValue(channelName, out int colIndex) && colIndex >= 0 && colIndex < spl.Length)
                {
                    presentChannels.Add(colIndex);
                }
            }

            if (presentChannels.Count == 0)
                return null;

            double[] values = new double[presentChannels.Count];
            for (int i = 0; i < presentChannels.Count; i++)
            {
                if (double.TryParse(spl[presentChannels[i]].Trim(),
                    NumberStyles.Any, CultureInfo.InvariantCulture, out double val))
                {
                    values[i] = val;
                }
            }
            return values;
        }

        #endregion
    }
}
