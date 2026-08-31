using Easy.Common.Extensions;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using System.Globalization;

namespace Readers
{
    public class PsmFromTsv : SpectrumMatchFromTsv
    {

        public string ProteinAccession => Accession;
        public string ProteinName => Name;
        public string PeptideMonoMass => MonoisotopicMassString;
        public string PeptideDescription => Description;
        public string PreviousAminoAcid => PreviousResidue;
        public string NextAminoAcid => NextResidue;
        public int PsmCount => SpectrumMatchCount;
        public string StartAndEndResiduesInProtein => StartAndEndResiduesInParentSequence;

        //For crosslink
        public string CrossType { get; }
        public string LinkResidues { get; }
        public int? ProteinLinkSite { get; }
        public int? Rank { get; }
        public string BetaPeptideProteinAccession { get; }
        public int? BetaPeptideProteinLinkSite { get; }
        public string BetaPeptideBaseSequence { get; }
        public string BetaPeptideFullSequence { get; }
        public string BetaPeptideTheoreticalMass { get; }
        public double? BetaPeptideScore { get; }
        public int? BetaPeptideRank { get; }
        public List<MatchedFragmentIon> BetaPeptideMatchedIons { get; }
        public Dictionary<int, List<MatchedFragmentIon>> BetaPeptideChildScanMatchedIons { get; }
        /// <summary>
        /// If Crosslink, this contains the alpha and beta sequences. Otherwise, it contains the full sequence
        /// </summary>
        public string UniqueSequence { get; }
        public double? XLTotalScore { get; }
        public string ParentIons { get; }

        public PsmFromTsv(string line, char[] split, Dictionary<string, int> parsedHeader, SpectrumMatchParsingParameters? parsingParams = null)
            : base (line, split, parsedHeader, parsingParams ??= new())
        {
            var spl = line.Split(split).Select(p => p.Trim('\"')).ToArray();

            //For crosslinks
            CrossType = GetOptionalValue(SpectrumMatchFromTsvHeader.CrossTypeLabel, parsedHeader, spl);
            LinkResidues = GetOptionalValue(SpectrumMatchFromTsvHeader.LinkResiduesLabel, parsedHeader, spl);
            ProteinLinkSite = GetOptionalValue<int>(SpectrumMatchFromTsvHeader.ProteinLinkSiteLabel, parsedHeader, spl);
            Rank = GetOptionalValue<int>(SpectrumMatchFromTsvHeader.RankLabel, parsedHeader, spl);
            BetaPeptideProteinAccession = GetOptionalValue(SpectrumMatchFromTsvHeader.BetaPeptideProteinAccessionLabel, parsedHeader, spl);
            BetaPeptideProteinLinkSite = GetOptionalValue<int>(SpectrumMatchFromTsvHeader.BetaPeptideProteinLinkSiteLabel, parsedHeader, spl);
            BetaPeptideBaseSequence = GetOptionalValue(SpectrumMatchFromTsvHeader.BetaPeptideBaseSequenceLabel, parsedHeader, spl);
            BetaPeptideFullSequence = GetOptionalValue(SpectrumMatchFromTsvHeader.BetaPeptideFullSequenceLabel, parsedHeader, spl);
            BetaPeptideTheoreticalMass = GetOptionalValue(SpectrumMatchFromTsvHeader.BetaPeptideTheoreticalMassLabel, parsedHeader, spl);
            BetaPeptideScore = GetOptionalValue<double>(SpectrumMatchFromTsvHeader.BetaPeptideScoreLabel, parsedHeader, spl);
            BetaPeptideRank = GetOptionalValue<int>(SpectrumMatchFromTsvHeader.BetaPeptideRankLabel, parsedHeader, spl);
            BetaPeptideMatchedIons = parsingParams.ParseMatchedFragmentIons
                ? (parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonsLabel] < 0) 
                    ? null 
                    : ((spl[parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonsLabel]].StartsWith("{")) 
                        ? ReadChildScanMatchedIons(spl[parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonsLabel]].Trim(), spl[parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonIntensitiesLabel]].Trim(), BetaPeptideBaseSequence, parsingParams).First().Value 
                        : ReadFragmentIonsFromString(spl[parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonsLabel]].Trim(), spl[parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonIntensitiesLabel]].Trim(), BetaPeptideBaseSequence, parsingParams))
                : [];
            XLTotalScore = GetOptionalValue<double>(SpectrumMatchFromTsvHeader.XLTotalScoreLabel, parsedHeader, spl);
            ParentIons = GetOptionalValue(SpectrumMatchFromTsvHeader.ParentIonsLabel, parsedHeader, spl);
            // This ensures backwards compatibility with old Crosslink Search Results
            // This works because the alpha and beta peptide full sequences are written to tsv with their crosslink site included (e.g., PEPTIDEK(4))
            if (UniqueSequence == null && BetaPeptideFullSequence != null)
            {
                UniqueSequence = FullSequence + BetaPeptideFullSequence;
            }

            // child scan matched ions for xlink and glyco. we are getting them all above and then deleting primary scan ions here.
            ChildScanMatchedIons = parsingParams.ParseMatchedFragmentIons 
                ? (!spl[parsedHeader[SpectrumMatchFromTsvHeader.MatchedIonMzRatios]].StartsWith("{")) 
                    ? null 
                    : ReadChildScanMatchedIons(spl[parsedHeader[SpectrumMatchFromTsvHeader.MatchedIonMzRatios]].Trim(), spl[parsedHeader[SpectrumMatchFromTsvHeader.MatchedIonIntensities]].Trim(), BaseSeq, parsingParams)
                : [];
            if (ChildScanMatchedIons != null && ChildScanMatchedIons.ContainsKey(Ms2ScanNumber))
            {
                ChildScanMatchedIons.Remove(Ms2ScanNumber);
            }

            // beta peptide child scan matched ions (for crosslinks)
            BetaPeptideChildScanMatchedIons = parsingParams.ParseMatchedFragmentIons
                ? (parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonsLabel] < 0) 
                    ? null 
                    : ((!spl[parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonsLabel]].StartsWith("{")) 
                        ? null 
                        : ReadChildScanMatchedIons(spl[parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonsLabel]].Trim(), spl[parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonIntensitiesLabel]].Trim(), BetaPeptideBaseSequence, parsingParams))
                : [];
            if (BetaPeptideChildScanMatchedIons != null && BetaPeptideChildScanMatchedIons.ContainsKey(Ms2ScanNumber))
            {
                BetaPeptideChildScanMatchedIons.Remove(Ms2ScanNumber);
            }
        }

        /// <summary>
        /// Constructor used to disambiguate PsmFromTsv to a single psm object
        /// </summary>
        /// <param name="psm">psm to disambiguate</param>
        /// <param name="fullSequence">sequence of ambiguous psm to use</param>
        public PsmFromTsv(PsmFromTsv psm, string fullSequence, int index = 0, string baseSequence = "")
            : base(psm, fullSequence, index, baseSequence)
        {
            // XL and/or Glyco Specific Fields
            VariantCrossingIons = psm.VariantCrossingIons?.ToList();
            CrossType = psm.CrossType;
            LinkResidues = psm.LinkResidues;
            ProteinLinkSite = psm.ProteinLinkSite;
            Rank = psm.Rank;
            BetaPeptideProteinAccession = psm.BetaPeptideProteinAccession;
            BetaPeptideProteinLinkSite = psm.BetaPeptideProteinLinkSite;
            BetaPeptideBaseSequence = psm.BetaPeptideBaseSequence;
            BetaPeptideFullSequence = psm.BetaPeptideFullSequence;
            BetaPeptideTheoreticalMass = psm.BetaPeptideTheoreticalMass;
            BetaPeptideScore = psm.BetaPeptideScore;
            BetaPeptideRank = psm.BetaPeptideRank;
            BetaPeptideMatchedIons = psm.BetaPeptideMatchedIons?.ToList();
            BetaPeptideChildScanMatchedIons = psm.BetaPeptideChildScanMatchedIons;
            XLTotalScore = psm.XLTotalScore;
            ParentIons = psm.ParentIons;
            RetentionTime = psm.RetentionTime;
        }

        /// <summary>
        /// Override library spectrum for cross link library spectrum implict conversion
        /// </summary>
        /// <returns></returns>
        public override LibrarySpectrum ToLibrarySpectrum()
        {
            bool isDecoy = this.DecoyContamTarget == "D";

            List<MatchedFragmentIon> fragments = new List<MatchedFragmentIon>();

            double matchedIonIntensitySum = Math.Max(1.0, this.MatchedIons.Select(i => i.Intensity).Sum());

            foreach (MatchedFragmentIon ion in this.MatchedIons)
            {
                fragments.Add(new MatchedFragmentIon(ion.NeutralTheoreticalProduct, ion.Mz, ion.Intensity / matchedIonIntensitySum, ion.Charge));
            }

            if (BetaPeptideMatchedIons.IsNotNullOrEmpty())
            {
                List<MatchedFragmentIon> betaFragments = new();
                foreach (var ion in BetaPeptideMatchedIons)
                {
                    betaFragments.Add(new MatchedFragmentIon(ion.NeutralTheoreticalProduct, ion.Mz, ion.Intensity / matchedIonIntensitySum, ion.Charge));
                }
                string uniqueSequence = UniqueSequence ?? FullSequence + BetaPeptideFullSequence;
                return new CrosslinkLibrarySpectrum(uniqueSequence, PrecursorMz, PrecursorCharge, fragments, RetentionTime, betaFragments);
            }

            return (new(this.FullSequence, this.PrecursorMz, this.PrecursorCharge, fragments, RetentionTime, isDecoy));
        }
    }
}
