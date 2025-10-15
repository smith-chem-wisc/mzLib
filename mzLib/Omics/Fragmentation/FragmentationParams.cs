namespace Omics.Fragmentation;

public class FragmentationParams
{
    public bool GenerateMIon { get; set; } = false;
    public List<MIonLoss> MIonLosses { get; set; } = new();
}