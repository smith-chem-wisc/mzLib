using Omics.Fragmentation;

namespace Transcriptomics;

public class RnaFragmentationParams : FragmentationParams
{
    public static readonly RnaFragmentationParams Default = new();
    public RnaFragmentationParams()
    {
        GenerateMIon = true;
        MIonLosses = new List<MIonLoss>();
    }
}
