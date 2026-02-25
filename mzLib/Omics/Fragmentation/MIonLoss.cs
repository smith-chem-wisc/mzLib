using Chemistry;

namespace Omics.Fragmentation;

public class MIonLoss(string name, string annotation, ChemicalFormula chemicalFormula) : IHasChemicalFormula
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
    public ChemicalFormula ThisChemicalFormula { get; set; } = chemicalFormula;

    public override string ToString() => Annotation;
}