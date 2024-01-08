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
    }
}
