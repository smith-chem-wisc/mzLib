using Chemistry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public enum ProductType
    {
        a,
        aStar,
        aDegree,
        adot,
        aBase,
        b,
        bAmmoniaLoss,
        bWaterLoss,
        //BnoB1ions,
        bdot,
        bBase,
        c,
        cdot,
        cBase,
        d,
        ddot,
        dBase,
        dH2O, // d-H20
        w,
        wdot,
        wBase,
        x,
        xdot,
        xBase,
        y,
        yAmmoniaLoss,
        yWaterLoss,
        ydot,
        yBase,
        z,
        zPlusOne,//This is zDot plus H
        zDot,
        zBase,
        M, //this is the molecular ion // [M]
        D, //this is a diagnostic ion // Modification loss mass
        Ycore, //Glyco core Y ions // [pep] + Neutral core Glycan mass (such as: [pep] + [N]) //Which already consider the loss of H2O and H-transfer
        Y //Glyco Y ions // [pep] + other Glycan mass 
    }

    public static class RnaProductTypeExtensions
    {
        private static readonly Dictionary<ProductType, ChemicalFormula> FragmentIonCaps =
            new Dictionary<ProductType, ChemicalFormula>
            {
                { ProductType.a, ChemicalFormula.ParseFormula("H") },
                { ProductType.adot, ChemicalFormula.ParseFormula("H2") },
                { ProductType.b, ChemicalFormula.ParseFormula("OH") },
                { ProductType.bdot, ChemicalFormula.ParseFormula("OH2") },
                { ProductType.c, ChemicalFormula.ParseFormula("O3H2P") },
                { ProductType.cdot, ChemicalFormula.ParseFormula("O3H3P") },
                { ProductType.d, ChemicalFormula.ParseFormula("O4H2P") },
                { ProductType.ddot, ChemicalFormula.ParseFormula("O4H3P") },

                { ProductType.w, ChemicalFormula.ParseFormula("H") },
                { ProductType.wdot, ChemicalFormula.ParseFormula("H2") },
                { ProductType.x, ChemicalFormula.ParseFormula("O-1H") },
                { ProductType.xdot, ChemicalFormula.ParseFormula("O-1H2") },
                { ProductType.y, ChemicalFormula.ParseFormula("O-3P-1") },
                { ProductType.ydot, ChemicalFormula.ParseFormula("O-3HP-1") },
                { ProductType.z, ChemicalFormula.ParseFormula("O-4P-1") },
                { ProductType.zDot, ChemicalFormula.ParseFormula("O-4HP-1") },
                //fragment - Base chemical formula is the corresponding fragment chemical formula subtracing 1 H as H is lost when base is removed
                { ProductType.aBase, ChemicalFormula.ParseFormula("H-2") }, // "H-1" -H 
                { ProductType.bBase, ChemicalFormula.ParseFormula("O1H-2") }, //"OH1" -H
                { ProductType.cBase, ChemicalFormula.ParseFormula("O3H-1P") }, //"O3P" -H
                { ProductType.dBase, ChemicalFormula.ParseFormula("O4H-1P") }, //"O4H2P" -H

                { ProductType.wBase, ChemicalFormula.ParseFormula("H-2") }, //"H"-H
                { ProductType.xBase, ChemicalFormula.ParseFormula("O-1H-2") }, //"O-1H" -H
                { ProductType.yBase, ChemicalFormula.ParseFormula("O-3H-2P-1") }, //"O-3P-1" -H
                { ProductType.zBase, ChemicalFormula.ParseFormula("O-4H-3P-1") }, //"O-4H-1P-1" -1
                //d-H2O
                { ProductType.dH2O, ChemicalFormula.ParseFormula("O3P") },
            };

        private static readonly Dictionary<DissociationType, List<ProductType>> ProductsFromDissociationType =
            new Dictionary<DissociationType, List<ProductType>>()
            {
                { DissociationType.Unknown, new List<ProductType>() },
                {
                    DissociationType.CID,
                    new List<ProductType> { ProductType.w, ProductType.y, ProductType.aBase, ProductType.dH2O }
                },
                { DissociationType.LowCID, new List<ProductType>() { } },
                { DissociationType.IRMPD, new List<ProductType>() { } },
                { DissociationType.ECD, new List<ProductType> { } },
                { DissociationType.PQD, new List<ProductType> { } },
                { DissociationType.ETD, new List<ProductType> { } },
                {
                    DissociationType.HCD,
                    new List<ProductType> { ProductType.w, ProductType.y, ProductType.aBase, ProductType.dH2O }
                },
                { DissociationType.AnyActivationType, new List<ProductType> { } },
                { DissociationType.EThcD, new List<ProductType> { } },
                { DissociationType.Custom, new List<ProductType> { } },
                { DissociationType.ISCID, new List<ProductType> { } }
            };

        public static List<ProductType> GetProductTypesFromDissociationType(this DissociationType dissociationType) =>
            ProductsFromDissociationType[dissociationType];

        public static FragmentationTerminus GetRnaTerminusType(this ProductType fragmentType)
        {
            // Super handy: http://stackoverflow.com/questions/4624248/c-logical-riddle-with-bit-operations-only-one-bit-is-set

            switch (fragmentType)
            {
                case ProductType.a:
                case ProductType.adot:
                case ProductType.aBase:
                case ProductType.b:
                case ProductType.bdot:
                case ProductType.bBase:
                case ProductType.c:
                case ProductType.cdot:
                case ProductType.cBase:
                case ProductType.d:
                case ProductType.ddot:
                case ProductType.dBase:
                case ProductType.dH2O:
                    return FragmentationTerminus.FivePrime;

                case ProductType.w:
                case ProductType.wdot:
                case ProductType.wBase:
                case ProductType.x:
                case ProductType.xdot:
                case ProductType.xBase:
                case ProductType.y:
                case ProductType.ydot:
                case ProductType.yBase:
                case ProductType.z:
                case ProductType.zDot:
                case ProductType.zBase:
                    return FragmentationTerminus.ThreePrime;

                case ProductType.aStar:
                case ProductType.aDegree:
                case ProductType.bAmmoniaLoss:
                case ProductType.bWaterLoss:
                case ProductType.yAmmoniaLoss:
                case ProductType.yWaterLoss:
                case ProductType.zPlusOne:
                case ProductType.M:
                case ProductType.D:
                case ProductType.Ycore:
                case ProductType.Y:
                default:
                    throw new ArgumentOutOfRangeException(nameof(fragmentType), fragmentType, null);
            }
        }

        public static ChemicalFormula GetIonCap(this ProductType fragmentType)
        {
            return FragmentIonCaps[fragmentType];
        }
    }
}
