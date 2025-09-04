using Chemistry;
using Omics.Fragmentation;

namespace Transcriptomics;

public class RnaFragmentationParams : FragmentationParams
{
    static RnaFragmentationParams()
    {
        Default = new();

        // Initialize common RNA M-Ion Losses
        var phosphoLossFomula = MIonLoss.PhosphoLoss.ThisChemicalFormula;
        var waterLossFormula = MIonLoss.WaterLoss.ThisChemicalFormula;

        // Adenine Base Loss
        var aLossFormula = Nucleotide.AdenineBase.BaseChemicalFormula - 2 * Nucleotide.AdenineBase.BaseChemicalFormula;
        var aLoss = new MIonLoss("Adenine Base Loss", "-A", aLossFormula);
        MIonLoss.AllMIonLosses.Add(aLoss.Annotation, aLoss);

        // Cytosine Base Loss
        var cLossFormula = Nucleotide.CytosineBase.BaseChemicalFormula - 2 * Nucleotide.CytosineBase.BaseChemicalFormula;
        var cLoss = new MIonLoss("Cytosine Base Loss", "-C", cLossFormula);
        MIonLoss.AllMIonLosses.Add(cLoss.Annotation, cLoss);

        // Guanine Base Loss
        var gLossFormula = Nucleotide.GuanineBase.BaseChemicalFormula - 2 * Nucleotide.GuanineBase.BaseChemicalFormula;
        var gLoss = new MIonLoss("Guanine Base Loss", "-G", gLossFormula);
        MIonLoss.AllMIonLosses.Add(gLoss.Annotation, gLoss);

        // Uracil Base Loss
        var uLossFormula = Nucleotide.UracilBase.BaseChemicalFormula - 2 * Nucleotide.UracilBase.BaseChemicalFormula;
        var uLoss = new MIonLoss("Uracil Base Loss", "-U", uLossFormula);
        MIonLoss.AllMIonLosses.Add(uLoss.Annotation, uLoss);

        // Water Losses
        var aWaterLossFormula = aLossFormula + waterLossFormula;
        var aWaterLoss = new MIonLoss("Adenine Base Water Loss", "-A-H2O", aWaterLossFormula);
        MIonLoss.AllMIonLosses.Add(aWaterLoss.Annotation, aWaterLoss);

        var cWaterLossFormula = cLossFormula + waterLossFormula;
        var cWaterLoss = new MIonLoss("Cytosine Base Water Loss", "-C-H2O", cWaterLossFormula);
        MIonLoss.AllMIonLosses.Add(cWaterLoss.Annotation, cWaterLoss);

        var gWaterLossFormula = gLossFormula + waterLossFormula;
        var gWaterLoss = new MIonLoss("Guanine Base Water Loss", "-G-H2O", gWaterLossFormula);
        MIonLoss.AllMIonLosses.Add(gWaterLoss.Annotation, gWaterLoss);

        var uWaterLossFormula = uLossFormula + waterLossFormula;
        var uWaterLoss = new MIonLoss("Uracil Base Water Loss", "-U-H2O", uWaterLossFormula);
        MIonLoss.AllMIonLosses.Add(uWaterLoss.Annotation, uWaterLoss);

        // Phospho Losses
        var aPhosphoLossFormula = aLossFormula + phosphoLossFomula;
        var aPhosphoLoss = new MIonLoss("Adenine Base Phospho Loss", "-A-P", aPhosphoLossFormula);
        MIonLoss.AllMIonLosses.Add(aPhosphoLoss.Annotation, aPhosphoLoss);

        var cPhosphoLossFormula = cLossFormula + phosphoLossFomula;
        var cPhosphoLoss = new MIonLoss("Cytosine Base Phospho Loss", "-C-P", cPhosphoLossFormula);
        MIonLoss.AllMIonLosses.Add(cPhosphoLoss.Annotation, cPhosphoLoss);

        var gPhosphoLossFormula = gLossFormula + phosphoLossFomula;
        var gPhosphoLoss = new MIonLoss("Guanine Base Phospho Loss", "-G-P", gPhosphoLossFormula);
        MIonLoss.AllMIonLosses.Add(gPhosphoLoss.Annotation, gPhosphoLoss);

        var uPhosphoLossFormula = uLossFormula + phosphoLossFomula;
        var uPhosphoLoss = new MIonLoss("Uracil Base Phospho Loss", "-U-P", uPhosphoLossFormula);
        MIonLoss.AllMIonLosses.Add(uPhosphoLoss.Annotation, uPhosphoLoss);

        // Phospho Water Losses
        var aPhosphoWaterLossFormula = aLossFormula + phosphoLossFomula + waterLossFormula;
        var aPhosphoWaterLoss = new MIonLoss("Adenine Base Phospho Water Loss", "-A-P-H2O", aPhosphoWaterLossFormula);
        MIonLoss.AllMIonLosses.Add(aPhosphoWaterLoss.Annotation, aPhosphoWaterLoss);

        var cPhosphoWaterLossFormula = cLossFormula + phosphoLossFomula + waterLossFormula;
        var cPhosphoWaterLoss = new MIonLoss("Cytosine Base Phospho Water Loss", "-C-P-H2O", cPhosphoWaterLossFormula);
        MIonLoss.AllMIonLosses.Add(cPhosphoWaterLoss.Annotation, cPhosphoWaterLoss);

        var gPhosphoWaterLossFormula = gLossFormula + phosphoLossFomula + waterLossFormula;
        var gPhosphoWaterLoss = new MIonLoss("Guanine Base Phospho Water Loss", "-G-P-H2O", gPhosphoWaterLossFormula);
        MIonLoss.AllMIonLosses.Add(gPhosphoWaterLoss.Annotation, gPhosphoWaterLoss);

        var uPhosphoWaterLossFormula = uLossFormula + phosphoLossFomula + waterLossFormula;
        var uPhosphoWaterLoss = new MIonLoss("Uracil Base Phospho Water Loss", "-U-P-H2O", uPhosphoWaterLossFormula);
        MIonLoss.AllMIonLosses.Add(uPhosphoWaterLoss.Annotation, uPhosphoWaterLoss);
    }

    public static readonly RnaFragmentationParams Default;
    public RnaFragmentationParams()
    {
        GenerateMIon = true;
        MIonLosses = new List<MIonLoss>();
    }
}
