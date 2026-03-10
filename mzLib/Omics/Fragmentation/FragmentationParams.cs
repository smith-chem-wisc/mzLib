namespace Omics.Fragmentation;

public class FragmentationParams : IFragmentationParams
{
    public bool GenerateMIon { get; set; } = false;
    public List<MIonLoss> MIonLosses { get; set; } = new();
}