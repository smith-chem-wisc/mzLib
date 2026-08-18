namespace Omics.Fragmentation;

public class FragmentationParams : IFragmentationParams, IEquatable<FragmentationParams>
{
    public bool GenerateMIon { get; set; } = false;
    public List<MIonLoss> MIonLosses { get; set; } = new();

    #region Equality

    public override bool Equals(object? obj)
        => obj is FragmentationParams fp && Equals(fp);

    bool IEquatable<IFragmentationParams>.Equals(IFragmentationParams? other)
        => other is FragmentationParams fp && Equals(fp);

    public bool Equals(FragmentationParams? other)
    {
        if (other is null) return false;
        return GenerateMIon == other.GenerateMIon
               && MIonListComparer.Instance.Equals(MIonLosses, other.MIonLosses);
    }

    public override int GetHashCode() => HashCode.Combine(GenerateMIon, MIonListComparer.Instance.GetHashCode(MIonLosses));

    #endregion
}