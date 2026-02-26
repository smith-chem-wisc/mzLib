namespace Omics.Fragmentation;

public interface IFragmentationParams
{
    /// <summary>
    /// Whether to generate M ions (the intact molecule with a charge state of 1)
    /// </summary>
    bool GenerateMIon { get; set; }

    /// <summary>
    /// The types of M ion losses to generate, if GenerateMIon is true
    /// </summary>
    List<MIonLoss> MIonLosses { get; set; }
}