using Chemistry;

namespace Omics.Fragmentation;

public class MIonLoss(string name, string annotation, ChemicalFormula chemicalFormula) : IHasChemicalFormula, IEquatable<MIonLoss>
{
    public static MIonLoss PhosphoLoss;
    public static MIonLoss WaterLoss;
    public static MIonLoss PhosphoWaterLoss;
    public static MIonLoss AmmoniaLoss;
    public static readonly Dictionary<string, MIonLoss> AllMIonLosses;

    static MIonLoss()
    {
        PhosphoLoss = new MIonLoss("Phosphate Loss", "-P", ChemicalFormula.ParseFormula("H1P1O3"));
        WaterLoss = new MIonLoss("Water Loss", "-H2O", ChemicalFormula.ParseFormula("H2O1"));
        PhosphoWaterLoss = new MIonLoss("Phosphate and Water Loss", "-P-H2O", ChemicalFormula.ParseFormula("H3P1O4"));
        AmmoniaLoss = new MIonLoss("Ammonia Loss", "-NH3", ChemicalFormula.ParseFormula("H3N1"));
        AllMIonLosses  = new Dictionary<string, MIonLoss>
        {
            { "-P", PhosphoLoss },
            { "-H2O", WaterLoss },
            { "-P-H2O", PhosphoWaterLoss },
            { "-NH3", AmmoniaLoss }
        };
    }

    public string Name { get; init; } = name;
    public string Annotation { get; init; } = annotation;
    public double MonoisotopicMass => ThisChemicalFormula.MonoisotopicMass;
    public ChemicalFormula ThisChemicalFormula { get; init; } = chemicalFormula;

    public override string ToString() => Annotation;

    public bool Equals(MIonLoss? other)
    {
        if (other is null) return false;
        return Name == other.Name
            && Annotation == other.Annotation
            && ThisChemicalFormula.Equals(other.ThisChemicalFormula);
    }

    public override bool Equals(object? obj) => obj is MIonLoss m && Equals(m);

    public override int GetHashCode() => HashCode.Combine(Name, Annotation, ThisChemicalFormula.GetHashCode());
}

/// <summary>
/// Compares MIonLoss lists for equality, ignoring order. Two lists are considered equal if they contain the same MIonLosses, regardless of their order in the list.
/// </summary>
public class MIonListComparer : IEqualityComparer<List<MIonLoss>>
{
    public static MIonListComparer Instance { get; } = new();

    public bool Equals(List<MIonLoss>? x, List<MIonLoss>? y)
    {
        if (x == null || y == null) return false;
        if (x.Count != y.Count) return false;
        var temp = new List<MIonLoss>(y);
        foreach (var item in x)
        {
            int idx = temp.FindIndex(item.Equals);
            if (idx < 0) return false;
            temp.RemoveAt(idx);
        }
        return true;
    }

    public int GetHashCode(List<MIonLoss> obj)
    {
        int hash = 17;
        foreach (var item in obj.OrderBy(i => i.Annotation))
        {
            hash = hash * 31 + item.GetHashCode();
        }
        return hash;
    }
}