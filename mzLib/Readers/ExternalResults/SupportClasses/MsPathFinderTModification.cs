namespace Readers;

public class MsPathFinderTModification(string modName, int modLocation, char modifiedResidue, double monoisotopicMass)
    : ILocalizedModification
{
    public string Name { get; } = modName;
    public int OneBasedLocalization { get; } = modLocation;
    public char ModifiedResidue { get; } = modifiedResidue;
    public double MonoisotopicMass { get; } = monoisotopicMass;

    public bool Equals(ILocalizedModification? other)
    {
        if (ReferenceEquals(null, other)) return false;
        if (ReferenceEquals(this, other)) return true;
        return OneBasedLocalization == other.OneBasedLocalization
               && Name == other.Name
               && ModifiedResidue == other.ModifiedResidue;
    }

    public override bool Equals(object? obj)
    {
        if (ReferenceEquals(null, obj)) return false;
        if (ReferenceEquals(this, obj)) return true;
        if (obj.GetType() != this.GetType()) return false;
        return Equals((ILocalizedModification)obj);
    }

    public override int GetHashCode()
    {
        return HashCode.Combine(OneBasedLocalization, Name, ModifiedResidue, MonoisotopicMass);
    }

    public override string ToString()
    {
        return $"{OneBasedLocalization}{ModifiedResidue}-{Name}";
    }
}