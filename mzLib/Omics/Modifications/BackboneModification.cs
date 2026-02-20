using System.Text;
using Chemistry;
using MassSpectrometry;
using Omics.Fragmentation;

namespace Omics.Modifications;    

/// <summary>
/// For backbone modifications: fragment types affected by the modification and their corresponding mass shifts. For example, phosphothiolate replaces the O with an S in the RNA backbone, so a, b, y, and z ions will not have a mass shift imparted by the modification, but c, d, w, x ions will have a mass shift O-1S1.
/// </summary>
public class BackboneModification : Modification
{
    /// <summary>
    /// The sorted product types in which the modification mass shift applies. 
    /// </summary>
    public ProductType[] ProductsContainingModMass { get; }

    public BackboneModification(string _originalId = null, string _accession = null, string _modificationType = null, string _featureType = null,
        ModificationMotif _target = null, string _locationRestriction = "Unassigned.", ChemicalFormula _chemicalFormula = null,
        double? _monoisotopicMass = null, Dictionary<string, IList<string>> _databaseReference = null,
        Dictionary<string, IList<string>> _taxonomicRange = null, List<string> _keywords = null,
        Dictionary<DissociationType, List<double>> _neutralLosses = null, Dictionary<DissociationType, List<double>> _diagnosticIons = null,
        string _fileOrigin = null, ProductType[] productsContainingShift = null) : base(_originalId, _accession, _modificationType, _featureType, _target, _locationRestriction, _chemicalFormula, _monoisotopicMass, _databaseReference, _taxonomicRange, _keywords, _neutralLosses, _diagnosticIons, _fileOrigin)
    {
        ProductsContainingModMass = productsContainingShift;
        Array.Sort(ProductsContainingModMass);
    }

    public override string ToString()
    {
        var sb =  new StringBuilder(base.ToString());

        if (this.ProductsContainingModMass is { Length: > 0 })
        {
            sb.Append("BM   ");
            sb.Append(string.Join(",", this.ProductsContainingModMass.Select(k => k.ToString())));
            sb.Append(":");
        }

        return sb.ToString();
    }

    private static readonly (ProductType, ProductType)[] _forbiddenPairs = 
        [
            (ProductType.a, ProductType.w), 
            (ProductType.b, ProductType.x), 
            (ProductType.c, ProductType.y), 
            (ProductType.d, ProductType.z)
        ];

    public override bool ValidModification
    {
        get
        {
            if (ProductsContainingModMass == null || ProductsContainingModMass.Length == 0)
            {
                // If ProductsContainingModMass is null or empty, return false to indicate this is not a valid backbone modification because we need to know which product types are affected by the modification to determine if it is a valid backbone modification.
                return false;
            }

            if (!base.ValidModification)
            {
                return false;
            }


            // No backbone mod should have the same mass shift on both types of ions in a pair (e.g., a and w, b and x, c and y, d and z) because that would indicate the modification is not actually on the backbone. If both types of ions in a pair have the same mass shift, return false to indicate this is not a valid backbone modification.
            foreach (var pair in _forbiddenPairs)
            {
                if (ProductsContainingModMass.Contains(pair.Item1) && ProductsContainingModMass.Contains(pair.Item2))
                {
                    return false;
                }
            }
            return true;
        }
    }
}
