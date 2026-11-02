using Chemistry;
using MassSpectrometry;
using Omics.Fragmentation;

namespace Omics.Modifications;

/// <summary>
/// Represents a modification that occurs at the cleavage site during enzymatic digestion.
/// </summary>
public class CleavageModification : Modification
{
    public bool IsFixedMod { get; }
    public bool IsVariableMod { get; }

    public CleavageModification(bool isFixed, bool isVariable, Modification mod)
        : this(isFixed, isVariable,mod.OriginalId, mod.Accession, mod.ModificationType, mod.FeatureType, mod.Target, mod.LocationRestriction, mod.ChemicalFormula, mod.MonoisotopicMass, mod.DatabaseReference, mod.TaxonomicRange, mod.Keywords, mod.NeutralLosses, mod.DiagnosticIons, mod.FileOrigin)
    {
    }

    public CleavageModification(bool isFixed, bool isVariable, string _originalId = null, string _accession = null, string _modificationType = null, string _featureType = null,
    ModificationMotif _target = null, string _locationRestriction = "Unassigned.", ChemicalFormula _chemicalFormula = null,
    double? _monoisotopicMass = null, Dictionary<string, IList<string>> _databaseReference = null,
    Dictionary<string, IList<string>> _taxonomicRange = null, List<string> _keywords = null,
    Dictionary<DissociationType, List<double>> _neutralLosses = null, Dictionary<DissociationType, List<double>> _diagnosticIons = null,
    string _fileOrigin = null)
    : base(_originalId, _accession, _modificationType, _featureType, _target, _locationRestriction, _chemicalFormula, _monoisotopicMass, _databaseReference, _taxonomicRange, _keywords, _neutralLosses, _diagnosticIons, _fileOrigin)
    {
        if (isFixed && isVariable)
            throw new ArgumentException("A modification cannot be both fixed and variable.");
        if (!isFixed && !isVariable)
            throw new ArgumentException("A modification must be either fixed or variable.");

        IsFixedMod = isFixed;
        IsVariableMod = isVariable;
    }
}
