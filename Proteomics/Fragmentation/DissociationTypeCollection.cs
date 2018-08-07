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
            { DissociationType.CID, new List<ProductType>{ ProductType.B, ProductType.Y } },
            { DissociationType.MPD, new List<ProductType>() },
            { DissociationType.ECD, new List<ProductType>{ ProductType.C, ProductType.Y, ProductType.Zdot } },
            { DissociationType.PQD, new List<ProductType>() },
            { DissociationType.ETD, new List<ProductType>{ ProductType.C, ProductType.Y, ProductType.Zdot } },
            { DissociationType.HCD, new List<ProductType>{ ProductType.B, ProductType.Y } },
            { DissociationType.AnyActivationType, new List<ProductType>{ ProductType.B, ProductType.Y } },
            { DissociationType.EThCD, new List<ProductType>{ ProductType.B, ProductType.Y, ProductType.C, ProductType.Zdot } },
            { DissociationType.Custom, new List<ProductType>() },
            { DissociationType.ISCID, new List<ProductType>() }
        };

        private static Dictionary<ProductType, double?> NeutralMassShiftFromProductType = new Dictionary<ProductType, double?>
        {
            { ProductType.A, null},//-C -0
            { ProductType.Astar, null},//-C -O -N -H3
            { ProductType.Adot, null},//-C -O2 -H2
            { ProductType.B, null},//no change
            { ProductType.Bstar, null},//-N -H3
            { ProductType.Bdot, null},//-H2 -O1
            { ProductType.C, null},//+N1 +H3
            { ProductType.X, null},//+C1 +O2
            { ProductType.Y, null},//+O +H2
            { ProductType.Ystar, null},//+O -H -N
            { ProductType.Ydot, null},//no change
            { ProductType.Zdot, null},//+O +H -N: A Zdot ion is also known as z+1. It is not a z-ion in the Biemann nomenclature. It differs from a y-ion by N-1 H-1;
            { ProductType.M, null},// neutral Molecular product can be used with neutral loss as fragment
            { ProductType.D, null},// diagnostic ions are not shifted but added sumarily
        };

        private static double GetMassShiftFromProductType(ProductType productType)
        {
            if (NeutralMassShiftFromProductType.TryGetValue(productType, out double? shift))
            {
                if (!shift.HasValue)
                {
                    // compute formula
                    switch (productType)
                    {
                        case ProductType.A: NeutralMassShiftFromProductType[productType] = ChemicalFormula.ParseFormula("C-1O-1").MonoisotopicMass; break;
                        case ProductType.Astar: NeutralMassShiftFromProductType[productType] = ChemicalFormula.ParseFormula("C-1O-1N-1H-3").MonoisotopicMass; break;
                        case ProductType.Adot: NeutralMassShiftFromProductType[productType] = ChemicalFormula.ParseFormula("C-1O-2H-2").MonoisotopicMass; break; // -46.0054793036,-C -O2 -H2
                        case ProductType.B: NeutralMassShiftFromProductType[productType] = 0; break;// 0, no change
                        case ProductType.Bstar: NeutralMassShiftFromProductType[productType] = ChemicalFormula.ParseFormula("N-1H-3").MonoisotopicMass; break;// -17.02654910112, -N -H3
                        case ProductType.Bdot: NeutralMassShiftFromProductType[productType] = ChemicalFormula.ParseFormula("H-2O-1").MonoisotopicMass; break;// -18.01056468403, -H2 -O1
                        case ProductType.C: NeutralMassShiftFromProductType[productType] = ChemicalFormula.ParseFormula("N1H3").MonoisotopicMass; break;// 17.02654910112, +N1 +H3
                        case ProductType.X: NeutralMassShiftFromProductType[productType] = ChemicalFormula.ParseFormula("C1O2").MonoisotopicMass; break;// 43.98982923914, +C1 +O2
                        case ProductType.Y: NeutralMassShiftFromProductType[productType] = ChemicalFormula.ParseFormula("H2O1").MonoisotopicMass; break;// 18.01056468403, +O +H2
                        case ProductType.Ystar: NeutralMassShiftFromProductType[productType] = ChemicalFormula.ParseFormula("O1H-1N-1").MonoisotopicMass; break;// 0.98401558291000057, +O -H -N
                        case ProductType.Ydot: NeutralMassShiftFromProductType[productType] = 0; break;// 0, no change
                        case ProductType.Zdot: NeutralMassShiftFromProductType[productType] = ChemicalFormula.ParseFormula("O1H1N-1").MonoisotopicMass; break;//; 2.9996656473699996, +O +H -N: A Zdot ion is also known as z+1. It is not a z-ion in the Biemann nomenclature. It differs from a y-ion by N-1 H-1;

                        case ProductType.M: NeutralMassShiftFromProductType[productType] = 0; break;// no change
                        case ProductType.D: NeutralMassShiftFromProductType[productType] = 0; break;// no change

                        default: throw new MzLibUtil.MzLibException("Unknown product type!");
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

        //Ion Type      Neutral Mr
        //a             [N]+[M]-CHO
        //a*	        a-NH3
        //a°	        a-H2O
        //b             [N]+[M]-H
        //b*	        b-NH3
        //b°	        b-H2O
        //c             [N]+[M]+NH2
        //d             a – partial side chain
        //v             y – complete side chain
        //w             z – partial side chain
        //x             [C]+[M]+CO-H
        //y             [C]+[M]+H
        //y*	        y-NH3
        //y°	        y-H2O
        //z             [C]+[M]-NH2
    }
}