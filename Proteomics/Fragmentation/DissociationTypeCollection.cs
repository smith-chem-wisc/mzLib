using Chemistry;
using MassSpectrometry;
using System.Collections.Generic;

namespace Proteomics.Fragmentation
{
    public class DissociationTypeCollection
    {
        public static Dictionary<DissociationType, List<ProductType>> ProductsFromDissociationType = new Dictionary<DissociationType, List<ProductType>>
        {
            { DissociationType.Unknown, new List<ProductType>() },
            { DissociationType.CID, new List<ProductType>{ ProductType.b, ProductType.y } },
            { DissociationType.IRMPD, new List<ProductType>{ ProductType.b, ProductType.y } },
            { DissociationType.ECD, new List<ProductType>{ ProductType.c, ProductType.y, ProductType.zPlusOne } },
            { DissociationType.PQD, new List<ProductType>() },
            { DissociationType.ETD, new List<ProductType>{ ProductType.c, ProductType.y, ProductType.zPlusOne } },
            { DissociationType.HCD, new List<ProductType>{ ProductType.b, ProductType.y } },//HCD often creates a-, aStar, and aDegree-ions and we should examine what other prominent algoroithms do to see if that would benefit our search results
            { DissociationType.AnyActivationType, new List<ProductType>{ ProductType.b, ProductType.y } },
            { DissociationType.EThcD, new List<ProductType>{ ProductType.b, ProductType.y, ProductType.c, ProductType.zPlusOne } },
            { DissociationType.Custom, new List<ProductType>() },
            { DissociationType.ISCID, new List<ProductType>() }
        };

        private static Dictionary<ProductType, double?> NeutralMassShiftFromProductType = new Dictionary<ProductType, double?>
        {
            { ProductType.a, null},//-C -O
            { ProductType.aStar, null},//-C -O -N -H3
            { ProductType.aDegree, null},//-C -O2 -H2
            { ProductType.b, null},//no change
            { ProductType.bStar, null},//-N -H3
            { ProductType.bDegree, null},//-H2 -O1
            { ProductType.c, null},//+N1 +H3
            { ProductType.x, null},//+C1 +O2
            { ProductType.y, null},//+O +H2
            { ProductType.yStar, null},//+O -H -N
            { ProductType.yDegree, null},//no change
            { ProductType.zPlusOne, null},//+O +H -N: A Zdot ion is also known as z+1. It is not a z-ion in the Biemann nomenclature. It differs from a y-ion by N-1 H-1;
            { ProductType.M, null},// neutral Molecular product can be used with neutral loss as fragment
            { ProductType.D, null},// diagnostic ions are not shifted but added sumarily
        };

        public static double GetMassShiftFromProductType(ProductType productType)
        {
            if (NeutralMassShiftFromProductType.TryGetValue(productType, out double? shift))
            {
                if (!shift.HasValue)
                {
                    // compute formula
                    switch (productType)
                    {
                        case ProductType.a: NeutralMassShiftFromProductType[productType] = ChemicalFormula.ParseFormula("C-1O-1").MonoisotopicMass; break;
                        case ProductType.aStar: NeutralMassShiftFromProductType[productType] = ChemicalFormula.ParseFormula("C-1O-1N-1H-3").MonoisotopicMass; break;
                        case ProductType.aDegree: NeutralMassShiftFromProductType[productType] = ChemicalFormula.ParseFormula("C-1O-2H-2").MonoisotopicMass; break; // -46.0054793036,-C -O2 -H2
                        case ProductType.b: NeutralMassShiftFromProductType[productType] = 0; break;// 0, no change
                        case ProductType.bStar: NeutralMassShiftFromProductType[productType] = ChemicalFormula.ParseFormula("N-1H-3").MonoisotopicMass; break;// -17.02654910112, -N -H3
                        case ProductType.bDegree: NeutralMassShiftFromProductType[productType] = ChemicalFormula.ParseFormula("H-2O-1").MonoisotopicMass; break;// -18.01056468403, -H2 -O1
                        case ProductType.c: NeutralMassShiftFromProductType[productType] = ChemicalFormula.ParseFormula("N1H3").MonoisotopicMass; break;// 17.02654910112, +N1 +H3
                        case ProductType.x: NeutralMassShiftFromProductType[productType] = ChemicalFormula.ParseFormula("C1O2").MonoisotopicMass; break;// 43.98982923914, +C1 +O2
                        case ProductType.y: NeutralMassShiftFromProductType[productType] = ChemicalFormula.ParseFormula("H2O1").MonoisotopicMass; break;// 18.01056468403, +O +H2
                        case ProductType.yStar: NeutralMassShiftFromProductType[productType] = ChemicalFormula.ParseFormula("O1H-1N-1").MonoisotopicMass; break;// 0.98401558291000057, +O -H -N
                        case ProductType.yDegree: NeutralMassShiftFromProductType[productType] = 0; break;// 0, no change
                        case ProductType.zPlusOne: NeutralMassShiftFromProductType[productType] = ChemicalFormula.ParseFormula("O1H1N-1").MonoisotopicMass; break;//; 2.9996656473699996, +O +H -N:
                        case ProductType.M: NeutralMassShiftFromProductType[productType] = 0; break;// no change
                        case ProductType.D: NeutralMassShiftFromProductType[productType] = 0; break;// no change
                    }
                }

                return NeutralMassShiftFromProductType[productType].Value;
            }
            else
            {
                throw new MzLibUtil.MzLibException("Unknown product type!");
            }
        }

        public static double ProductTypeSpecificFragmentNeutralMass(double mass, ProductType p)
        {
            return (double)ClassExtensions.RoundedDouble(mass + GetMassShiftFromProductType(p), 9);
        }
    }
}