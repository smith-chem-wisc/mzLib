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
        CleavageSpecificity SearchModeType { get; }

        /// <summary>
        /// new terminus parameter is for non and semi specific searches
        /// </summary>
        /// <param name="newTerminus"></param>
        /// <returns></returns>
        IDigestionParams Clone(FragmentationTerminus? newTerminus = null);
    }
}
