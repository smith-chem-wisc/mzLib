using Chemistry;
using Omics.SpectrumMatch;
using Transcriptomics;
using Transcriptomics.Digestion;

namespace Readers
{
    public class OsmFromTsv : SpectrumMatchFromTsv
    {
        public IHasChemicalFormula FivePrimeTerminus { get; set; }
        public IHasChemicalFormula ThreePrimeTerminus { get; set; }

        public OsmFromTsv(string line, char[] split, Dictionary<string, int> parsedHeader)
            : base(line, split, parsedHeader)
        {
            var spl = line.Split(split).Select(p => p.Trim('\"')).ToArray();

            if (parsedHeader[SpectrumMatchFromTsvHeader.FivePrimeTerminus] >= 0)
                FivePrimeTerminus = ChemicalFormula.ParseFormula(spl[parsedHeader[SpectrumMatchFromTsvHeader.FivePrimeTerminus]]);
            else if (PreviousResidue == "-")
                FivePrimeTerminus = NucleicAcid.DefaultFivePrimeTerminus;
            else
                FivePrimeTerminus = Rnase.DefaultFivePrimeTerminus;

            if (parsedHeader[SpectrumMatchFromTsvHeader.ThreePrimeTerminus] >= 0)
                ThreePrimeTerminus = ChemicalFormula.ParseFormula(spl[parsedHeader[SpectrumMatchFromTsvHeader.ThreePrimeTerminus]]);
            else if (NextResidue == "-")
                ThreePrimeTerminus = NucleicAcid.DefaultThreePrimeTerminus;
            else
                ThreePrimeTerminus = Rnase.DefaultThreePrimeTerminus;
        }

        /// <summary>
        /// Constructor used to disambiguate PsmFromTsv to a single osm object
        /// </summary>
        /// <param name="osm">osm to disambiguate</param>
        /// <param name="fullSequence">sequence of ambiguous osm to use</param>
        public OsmFromTsv(OsmFromTsv osm, string fullSequence, int index = 0, string baseSequence = "",
            IHasChemicalFormula? fivePrimeTerm = null, IHasChemicalFormula? threePrimeTerm = null)
            : base(osm, fullSequence, index, baseSequence)
        {
            FivePrimeTerminus = fivePrimeTerm ?? osm.FivePrimeTerminus;
            ThreePrimeTerminus = threePrimeTerm ?? osm.ThreePrimeTerminus;
        }
    }
}
