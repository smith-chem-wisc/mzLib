using Chemistry;

namespace Omics.Fragmentation;

public class MIonLoss(string name, string annotation, ChemicalFormula chemicalFormula) : IHasChemicalFormula
{
    public static MIonLoss PhosphoLoss;
    public static MIonLoss WaterLoss;
    public static MIonLoss PhosphoWaterLoss;
    public static readonly Dictionary<string, MIonLoss> AllMIonLosses;

    static MIonLoss()
    {
        PhosphoLoss = new MIonLoss("Phosphate Loss", "-P", ChemicalFormula.ParseFormula("H-1P-1O-3"));
        WaterLoss = new MIonLoss("Water Loss", "-H2O", ChemicalFormula.ParseFormula("H-2O-1"));
        PhosphoWaterLoss = new MIonLoss("Phosphate and Water Loss", "-P-H2O", ChemicalFormula.ParseFormula("H-3P-1O-4"));
        AllMIonLosses  = new Dictionary<string, MIonLoss>
        {
            { "-P", PhosphoLoss },
            { "-H2O", WaterLoss },
            { "-P-H2O", PhosphoWaterLoss }
        };
    }

    public string Name { get; init; } = name;
    public string Annotation { get; init; } = annotation;
    public double MonoisotopicMass => ThisChemicalFormula.MonoisotopicMass;
    public ChemicalFormula ThisChemicalFormula { get; set; } = chemicalFormula;
}