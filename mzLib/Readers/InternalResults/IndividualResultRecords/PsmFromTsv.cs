using System.Globalization;
using Easy.Common.Extensions;
using Omics.Fragmentation;
using Omics.SpectrumMatch;

namespace Readers
{
    public class PsmFromTsv : SpectrumMatchFromTsv
    {

        public string ProteinAccession => Accession;
        public string ProteinName => Name;
        public string PeptideMonoMass => MonoisotopicMass;
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

        //For Glyco
        public string GlycanStructure { get; set; }
        public double? GlycanMass { get; set; }
        public string GlycanComposition { get; set; }
        public LocalizationLevel? GlycanLocalizationLevel { get; set; }
        public string LocalizedGlycan { get; set; }

        public PsmFromTsv(string line, char[] split, Dictionary<string, int> parsedHeader)
            : base (line, split, parsedHeader)
        {
            var spl = line.Split(split).Select(p => p.Trim('\"')).ToArray();

            //For crosslinks
            CrossType = (parsedHeader[SpectrumMatchFromTsvHeader.CrossTypeLabel] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.CrossTypeLabel]].Trim();
            LinkResidues = (parsedHeader[SpectrumMatchFromTsvHeader.LinkResiduesLabel] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.LinkResiduesLabel]].Trim();
            ProteinLinkSite = (parsedHeader[SpectrumMatchFromTsvHeader.ProteinLinkSiteLabel] < 0) ? null : (spl[parsedHeader[SpectrumMatchFromTsvHeader.ProteinLinkSiteLabel]] == "" ? null : (int?)int.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.ProteinLinkSiteLabel]].Trim()));
            Rank = (parsedHeader[SpectrumMatchFromTsvHeader.RankLabel] < 0) ? null : (int?)int.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.RankLabel]].Trim());
            BetaPeptideProteinAccession = (parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideProteinAccessionLabel] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideProteinAccessionLabel]].Trim();
            BetaPeptideProteinLinkSite = (parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideProteinLinkSiteLabel] < 0) ? null : (spl[parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideProteinLinkSiteLabel]] == "" ? null : (int?)int.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideProteinLinkSiteLabel]].Trim()));
            BetaPeptideBaseSequence = (parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideBaseSequenceLabel] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideBaseSequenceLabel]].Trim();
            BetaPeptideFullSequence = (parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideFullSequenceLabel] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideFullSequenceLabel]].Trim();
            BetaPeptideTheoreticalMass = (parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideTheoreticalMassLabel] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideTheoreticalMassLabel]].Trim();
            BetaPeptideScore = (parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideScoreLabel] < 0) ? null : (double?)double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideScoreLabel]].Trim(), CultureInfo.InvariantCulture);
            BetaPeptideRank = (parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideRankLabel] < 0) ? null : (int?)int.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideRankLabel]].Trim());
            BetaPeptideMatchedIons = (parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonsLabel] < 0) ? null :
                ((spl[parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonsLabel]].StartsWith("{")) ? ReadChildScanMatchedIons(spl[parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonsLabel]].Trim(), spl[parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonIntensitiesLabel]].Trim(), BetaPeptideBaseSequence).First().Value : ReadFragmentIonsFromString(spl[parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonsLabel]].Trim(), spl[parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonIntensitiesLabel]].Trim(), BetaPeptideBaseSequence));
            XLTotalScore = (parsedHeader[SpectrumMatchFromTsvHeader.XLTotalScoreLabel] < 0) ? null : (double?)double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.XLTotalScoreLabel]].Trim(), CultureInfo.InvariantCulture);
            ParentIons = (parsedHeader[SpectrumMatchFromTsvHeader.ParentIonsLabel] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.ParentIonsLabel]].Trim();
            // This ensures backwards compatibility with old Crosslink Search Results
            // This works because the alpha and beta peptide full sequences are written to tsv with their crosslink site included (e.g., PEPTIDEK(4))
            if (UniqueSequence == null && BetaPeptideFullSequence != null)
            {
                UniqueSequence = FullSequence + BetaPeptideFullSequence;
            }

            // child scan matched ions for xlink and glyco. we are getting them all above and then deleting primary scan ions here.
            ChildScanMatchedIons = (!spl[parsedHeader[SpectrumMatchFromTsvHeader.MatchedIonMzRatios]].StartsWith("{")) ? null : ReadChildScanMatchedIons(spl[parsedHeader[SpectrumMatchFromTsvHeader.MatchedIonMzRatios]].Trim(), spl[parsedHeader[SpectrumMatchFromTsvHeader.MatchedIonIntensities]].Trim(), BaseSeq);
            if (ChildScanMatchedIons != null && ChildScanMatchedIons.ContainsKey(Ms2ScanNumber))
            {
                ChildScanMatchedIons.Remove(Ms2ScanNumber);
            }

            // beta peptide child scan matched ions (for crosslinks)
            BetaPeptideChildScanMatchedIons = (parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonsLabel] < 0) ? null :
                ((!spl[parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonsLabel]].StartsWith("{")) ? null : ReadChildScanMatchedIons(spl[parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonsLabel]].Trim(), spl[parsedHeader[SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonIntensitiesLabel]].Trim(), BetaPeptideBaseSequence));
            if (BetaPeptideChildScanMatchedIons != null && BetaPeptideChildScanMatchedIons.ContainsKey(Ms2ScanNumber))
            {
                BetaPeptideChildScanMatchedIons.Remove(Ms2ScanNumber);
            }

            //For Glyco            
            try // Try is so that glyco and non-glyco psms can be read from the same file
            {
                GlycanMass = (parsedHeader[PsmTsvHeader_Glyco.GlycanMass] < 0) ? null : (double?)double.Parse(spl[parsedHeader[PsmTsvHeader_Glyco.GlycanMass]], CultureInfo.InvariantCulture);
                GlycanComposition = (parsedHeader[PsmTsvHeader_Glyco.GlycanComposition] < 0) ? null : spl[parsedHeader[PsmTsvHeader_Glyco.GlycanComposition]];
                GlycanStructure = (parsedHeader[PsmTsvHeader_Glyco.GlycanStructure] < 0) ? null : spl[parsedHeader[PsmTsvHeader_Glyco.GlycanStructure]];
                var localizationLevel = (parsedHeader[PsmTsvHeader_Glyco.GlycanLocalizationLevel] < 0) ? null : spl[parsedHeader[PsmTsvHeader_Glyco.GlycanLocalizationLevel]];
                if (localizationLevel != null)
                {
                    if (localizationLevel.Equals("NA"))
                        GlycanLocalizationLevel = null;
                    else
                        GlycanLocalizationLevel = (LocalizationLevel)Enum.Parse(typeof(LocalizationLevel), localizationLevel);
                }
                LocalizedGlycan = (parsedHeader[PsmTsvHeader_Glyco.LocalizedGlycan] < 0) ? null : spl[parsedHeader[PsmTsvHeader_Glyco.LocalizedGlycan]];

            }
            catch
            {
                GlycanMass = null;
                GlycanComposition = null;
                GlycanStructure = null;
                GlycanLocalizationLevel = null;
                LocalizedGlycan = null;
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
            GlycanStructure = psm.GlycanStructure;
            GlycanMass = psm.GlycanMass;
            GlycanComposition = psm.GlycanComposition;
            GlycanLocalizationLevel = psm.GlycanLocalizationLevel;
            LocalizedGlycan = psm.LocalizedGlycan;
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
                Product product = new Product(ion.NeutralTheoreticalProduct.ProductType, ion.NeutralTheoreticalProduct.Terminus, ion.NeutralTheoreticalProduct.NeutralMass, ion.NeutralTheoreticalProduct.FragmentNumber, ion.NeutralTheoreticalProduct.AminoAcidPosition, ion.NeutralTheoreticalProduct.NeutralLoss);
                fragments.Add(new MatchedFragmentIon(product, ion.Mz, ion.Intensity / matchedIonIntensitySum, ion.Charge));
            }
            double retentionTime = RetentionTime ?? -1;

            if (BetaPeptideMatchedIons.IsNotNullOrEmpty())
            {
                List<MatchedFragmentIon> betaFragments = new();
                foreach (var ion in BetaPeptideMatchedIons)
                {
                    Product product = new Product(ion.NeutralTheoreticalProduct.ProductType, ion.NeutralTheoreticalProduct.Terminus, ion.NeutralTheoreticalProduct.NeutralMass, ion.NeutralTheoreticalProduct.FragmentNumber, ion.NeutralTheoreticalProduct.AminoAcidPosition, ion.NeutralTheoreticalProduct.NeutralLoss);
                    betaFragments.Add(new MatchedFragmentIon(product, ion.Mz, ion.Intensity / matchedIonIntensitySum, ion.Charge));
                }
                string uniqueSequence = UniqueSequence ?? FullSequence + BetaPeptideFullSequence;
                return new CrosslinkLibrarySpectrum(uniqueSequence, PrecursorMz, PrecursorCharge, fragments, retentionTime, betaFragments);
            }

            return (new(this.FullSequence, this.PrecursorMz, this.PrecursorCharge, fragments, retentionTime, isDecoy));
        }
    }
}
