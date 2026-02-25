namespace Omics.Fragmentation;

public class FragmentationParams
{
    public bool GenerateMIon { get; set; } = false;
    public List<MIonLoss> MIonLosses { get; set; } = new();
    public int MaxModsForCumulativeNeutralLosses { get; set; } = 1; // Default to 1, which means only consider individual modifications without combinations. Keep the current searching not crashing.
}