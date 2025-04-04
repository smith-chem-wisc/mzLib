using Omics.SpectrumMatch;

namespace Readers
{
    public class OsmFromTsv : SpectrumMatchFromTsv
    {
        public OsmFromTsv(string line, char[] split, Dictionary<string, int> parsedHeader)
            : base(line, split, parsedHeader)
        {
            // TODO: Parse Oligo specific columns
        }

        /// <summary>
        /// Constructor used to disambiguate PsmFromTsv to a single psm object
        /// </summary>
        /// <param name="psm">psm to disambiguate</param>
        /// <param name="fullSequence">sequence of ambiguous psm to use</param>
        public OsmFromTsv(OsmFromTsv psm, string fullSequence, int index = 0, string baseSequence = "")
            : base(psm, fullSequence, index, baseSequence)
        {
            // TODO: Parse Oligo specific columns
        }
    }
}
