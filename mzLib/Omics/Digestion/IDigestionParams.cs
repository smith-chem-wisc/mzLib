using Omics.Fragmentation;

namespace Omics.Digestion
{
    public interface IDigestionParams
    {
        int MaxMissedCleavages { get; set; }
        int MinLength { get; set; }
        int MaxLength { get; set; }
        int MaxModificationIsoforms { get; set; }
        int MaxMods { get; set; }
        DigestionAgent DigestionAgent { get; }
        FragmentationTerminus FragmentationTerminus { get; }
        /// <summary>
        /// Search mode type refers to the CleavageSpecificity enum and is used for MetaMorpheus to determine if it should perform a non-specific, semi-specific, or fully specific search.
        /// For the initial implementation, RNA will have it hardcoded to a fully specific search.
        /// </summary>
        CleavageSpecificity SearchModeType { get; }

        /// <summary>
        /// new terminus parameter is for non and semi specific searches
        /// The Fragmentation terminus will never be null, if a null value is inputted, the fragmentation terminus will be set to the fragmentation terminus of the object being cloned.
        /// </summary>
        /// <param name="newTerminus"></param>
        /// <returns></returns>
        IDigestionParams Clone(FragmentationTerminus? newTerminus = null);
    }
}
