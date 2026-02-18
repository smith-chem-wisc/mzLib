using Chemistry;
using MassSpectrometry;
using System.Text;

namespace Omics.Modifications;
public enum BaseLossBehavior
{
    Default,      // Normal base loss
    Suppressed,   // No base loss (e.g., 2'-O-methyl)
    Modified      // Base loss includes modification (e.g., m6A)
}

/// <summary>
/// Represents a modification that modified the base in some form, which can affect the base loss behavior during fragmentation. 
/// For example, 2-o methylation adds its mass to the sugar and would be present in the residual after base loss whereas m6 methyl adds its mass to the base and would be lost during base loss. 
/// In addition, 2-o stabilizes the bond and suppresses base loss, whereas m6A does not suppress base loss. 
/// 
/// Cases:
/// 2-O methylation: BaseLossType = Suppressed, BaseLossModification = null
/// m6A methylation: BaseLossType = Modified, BaseLossModification = ChemicalFormula.ParseFormula("CH3")
/// 2-Om6 methylation: BaseLossType = Suppressed, BaseLossModification = ChemicalFormula.ParseFormula("CH3") (the base loss includes the m6A modification, but the 2-O methylation is not included in the base loss because it is on the sugar and not the base)
/// </summary>
public class BaseModification : Modification
{
    public BaseLossBehavior BaseLossType { get; }
    public ChemicalFormula? BaseLossModification { get; }

    public BaseModification(string _originalId = null, string _accession = null, string _modificationType = null, string _featureType = null,
        ModificationMotif _target = null, string _locationRestriction = "Unassigned.", ChemicalFormula _chemicalFormula = null,
        double? _monoisotopicMass = null, Dictionary<string, IList<string>> _databaseReference = null,
        Dictionary<string, IList<string>> _taxonomicRange = null, List<string> _keywords = null,
        Dictionary<DissociationType, List<double>> _neutralLosses = null, Dictionary<DissociationType, List<double>> _diagnosticIons = null,
        string _fileOrigin = null, BaseLossBehavior baseLossType = BaseLossBehavior.Default, ChemicalFormula baseLossModification = null)
        : base(_originalId, _accession, _modificationType, _featureType, _target, _locationRestriction, _chemicalFormula, _monoisotopicMass, _databaseReference, _taxonomicRange, _keywords, _neutralLosses, _diagnosticIons, _fileOrigin)
    {
        BaseLossType = baseLossType;
        BaseLossModification = baseLossModification;
    }

    public override bool ValidModification
    {
        get
        {
            if (!base.ValidModification) return false;

            // If Modified type, BaseLossModification should be provided
            if (BaseLossType == BaseLossBehavior.Modified && BaseLossModification == null)
                return false;

            return true;
        }
    }

    public override string ToString()
    {
        var sb = new StringBuilder(base.ToString());

        if (BaseLossType != BaseLossBehavior.Default)
        {
            sb.Append("BL   ");
            sb.Append(BaseLossType.ToString());
            if (BaseLossModification != null)
            {
                sb.Append(":");
                sb.Append(BaseLossModification);
            }
        }
        return sb.ToString();
    }
}