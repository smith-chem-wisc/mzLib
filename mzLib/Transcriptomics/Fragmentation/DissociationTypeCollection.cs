using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using MassSpectrometry;

namespace Transcriptomics
{
    /// <summary>
    /// Methods dealing with specific product type for RNA molecules
    /// </summary>
    public static class DissociationTypeCollection
    {
        /// <summary>
        /// Product Ion types by dissociation method
        /// </summary>
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

        /// <summary>
        /// Returns list of products types based upon the dissociation type
        /// </summary>
        /// <param name="dissociationType"></param>
        /// <returns></returns>
        public static List<ProductType> GetRnaProductTypesFromDissociationType(this DissociationType dissociationType) =>
            ProductsFromDissociationType[dissociationType];


        /// <summary>
        /// Mass to be added or subtracted
        /// </summary>
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

        /// <summary>
        /// Returns mass shift by product type
        /// </summary>
        /// <param name="type"></param>
        /// <returns></returns>
        public static double GetRnaMassShiftFromProductType(this ProductType type) => FragmentIonCaps[type].MonoisotopicMass;

        public static FragmentationTerminus GetRnaTerminusType(this ProductType fragmentType)
        {
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

        /// <summary>
        /// Product ion types by Fragmentation Terminus
        /// </summary>
        private static readonly Dictionary<FragmentationTerminus, List<ProductType>>
            ProductIonTypesFromSpecifiedTerminus = new Dictionary<FragmentationTerminus, List<ProductType>>
            {
                {
                    FragmentationTerminus.FivePrime, new List<ProductType>
                    {
                        ProductType.a, ProductType.adot, ProductType.aBase,
                        ProductType.b, ProductType.bdot, ProductType.bBase,
                        ProductType.c, ProductType.cdot, ProductType.cBase,
                        ProductType.d, ProductType.ddot, ProductType.dBase, ProductType.dH2O,
                    }
                },
                {
                    FragmentationTerminus.ThreePrime, new List<ProductType>
                    {
                        ProductType.w, ProductType.wdot, ProductType.wBase,
                        ProductType.x, ProductType.xdot, ProductType.xBase,
                        ProductType.y, ProductType.ydot, ProductType.yBase,
                        ProductType.z, ProductType.zDot, ProductType.zBase,
                    }
                },
                {
                    FragmentationTerminus.Both, new List<ProductType>
                    {

                        ProductType.a, ProductType.adot, ProductType.aBase,
                        ProductType.b, ProductType.bdot, ProductType.bBase,
                        ProductType.c, ProductType.cdot, ProductType.cBase,
                        ProductType.d, ProductType.ddot, ProductType.dBase, ProductType.dH2O,
                        ProductType.w, ProductType.wdot, ProductType.wBase,
                        ProductType.x, ProductType.xdot, ProductType.xBase,
                        ProductType.y, ProductType.ydot, ProductType.yBase,
                        ProductType.z, ProductType.zDot, ProductType.zBase,
                    }
                }
            };


        public static List<ProductType> GetRnaTerminusSpecificProductTypes(
            this FragmentationTerminus fragmentationTerminus)
        {
            return ProductIonTypesFromSpecifiedTerminus[fragmentationTerminus];
        }

        /// <summary>
        /// Returns all product ion types based upon specified terminus
        /// </summary>
        /// <param name="dissociationType"></param>
        /// <param name="fragmentationTerminus"></param>
        /// <returns></returns>
        public static List<ProductType> GetRnaTerminusSpecificProductTypesFromDissociation(
            this DissociationType dissociationType, FragmentationTerminus fragmentationTerminus)
        {
            var terminusSpecific = fragmentationTerminus.GetRnaTerminusSpecificProductTypes();
            var dissociationSpecific = dissociationType.GetRnaProductTypesFromDissociationType();
            return terminusSpecific.Intersect(dissociationSpecific).ToList();
        }
    }
}
