using Chemistry;
using MassSpectrometry;

namespace Omics.Fragmentation.Oligo
{
    /// <summary>
    /// Methods dealing with specific product type for RNA molecules
    /// </summary>
    public static class DissociationTypeCollection
    {
        /// <summary>
        /// Product Ion types by dissociation method
        /// </summary>
        /// <remarks>
        /// HCD ions were taken from the following paper: https://www.nature.com/articles/s41598-023-36193-2
        /// Ion types below here should be validated with experimental results.
        /// Base and water losses occur very frequently and may also be present in these activation types.
        /// CID, UVPD, and aEPD ions were taken from the following paper: https://pubs.acs.org/doi/10.1021/acs.analchem.3c05428?ref=PDF
        /// NETD ions were taken from the following paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7161943/
        /// lowCID ions were taken from this Thermo Poster: https://assets.thermofisher.com/TFS-Assets/CMD/Flyers/fl-489263-asms23-optimized-fragmentation-oligonucleotides-suppresses-undesired-fragmentation-fl489263-en.pdf
        /// </remarks>
        public static Dictionary<DissociationType, List<ProductType>> ProductsFromDissociationType =
            new Dictionary<DissociationType, List<ProductType>>()
            {
                { DissociationType.Unknown, new List<ProductType>() },
                { DissociationType.Custom, new List<ProductType>() },
                {
                    DissociationType.AnyActivationType, new List<ProductType>
                    {
                        ProductType.a, ProductType.aBaseLoss, ProductType.aWaterLoss,
                        ProductType.b, ProductType.bBaseLoss, ProductType.bWaterLoss,
                        ProductType.c, ProductType.cBaseLoss, ProductType.cWaterLoss,
                        ProductType.d, ProductType.dBaseLoss, ProductType.dWaterLoss,
                        ProductType.w, ProductType.wBaseLoss, ProductType.wWaterLoss,
                        ProductType.x, ProductType.xBaseLoss, ProductType.xWaterLoss,
                        ProductType.y, ProductType.yBaseLoss, ProductType.yWaterLoss,
                        ProductType.z, ProductType.zBaseLoss, ProductType.zWaterLoss,
                        ProductType.M
                    }
                },
                {
                    DissociationType.CID, new List<ProductType>
                    {
                        ProductType.a, ProductType.aBaseLoss, ProductType.c, ProductType.dWaterLoss, ProductType.w,
                        ProductType.y, ProductType.yWaterLoss, ProductType.M
                    }
                },
                {
                    DissociationType.HCD, new List<ProductType>
                    {
                        ProductType.a, ProductType.aBaseLoss, ProductType.b, ProductType.c, ProductType.d,
                        ProductType.dWaterLoss, ProductType.w, ProductType.x, ProductType.y, ProductType.z,
                        ProductType.M
                    }
                },
                {
                    DissociationType.UVPD, new List<ProductType>
                    {
                        ProductType.a, ProductType.c, ProductType.d, ProductType.w, ProductType.M
                    }
                },
                {
                    DissociationType.aEPD, new List<ProductType>
                    {
                        ProductType.a, ProductType.c, ProductType.d, ProductType.w, ProductType.x, ProductType.z, ProductType.M
                    }
                },
                {
                    DissociationType.NETD, new List<ProductType>
                    {
                        ProductType.w, ProductType.d, ProductType.M
                    }
                },
                {
                    DissociationType.LowCID, new List<ProductType>()
                    {
                        ProductType.aBaseLoss, ProductType.c, ProductType.dWaterLoss, ProductType.w,
                        ProductType.y, ProductType.yWaterLoss, ProductType.M
                    }
                },
                { DissociationType.IRMPD, new List<ProductType>() { } },
                { DissociationType.ECD, new List<ProductType> { } },
                { DissociationType.PQD, new List<ProductType> { } },
                { DissociationType.ETD, new List<ProductType> { } },
                { DissociationType.EThcD, new List<ProductType> { } },
            };

        /// <summary>
        /// Returns all dissociation types with implemented product type collections
        /// </summary>
        public static IEnumerable<DissociationType> AllImplementedDissociationTypes =>
            ProductsFromDissociationType.Where(p => p.Value.Any())
                .Select(p => p.Key);

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
                { ProductType.aWaterLoss, ChemicalFormula.ParseFormula("H-1O-1") },
                { ProductType.b, ChemicalFormula.ParseFormula("OH") },
                { ProductType.bWaterLoss, ChemicalFormula.ParseFormula("H-1") },
                { ProductType.c, ChemicalFormula.ParseFormula("O3H2P") },
                { ProductType.cWaterLoss, ChemicalFormula.ParseFormula("O2P") },
                { ProductType.d, ChemicalFormula.ParseFormula("O4H2P") },
                { ProductType.dWaterLoss, ChemicalFormula.ParseFormula("O3P") },

                { ProductType.w, ChemicalFormula.ParseFormula("H") },
                { ProductType.wWaterLoss, ChemicalFormula.ParseFormula("H-1O-1") },
                { ProductType.x, ChemicalFormula.ParseFormula("O-1H") },
                { ProductType.xWaterLoss, ChemicalFormula.ParseFormula("O-2H-1") },
                { ProductType.y, ChemicalFormula.ParseFormula("O-3P-1") },
                { ProductType.yWaterLoss, ChemicalFormula.ParseFormula("O-4H-2P-1") },
                { ProductType.z, ChemicalFormula.ParseFormula("O-4P-1") },
                { ProductType.zWaterLoss, ChemicalFormula.ParseFormula("O-5H-2P-1") },
                //fragment - Base chemical formula is the corresponding fragment chemical formula subtracing 1 H as H is lost when base is removed
                { ProductType.aBaseLoss, ChemicalFormula.ParseFormula("H-2") }, // "H-1" -H 
                { ProductType.bBaseLoss, ChemicalFormula.ParseFormula("O1H-2") }, //"OH1" -H
                { ProductType.cBaseLoss, ChemicalFormula.ParseFormula("O3H-1P") }, //"O3P" -H
                { ProductType.dBaseLoss, ChemicalFormula.ParseFormula("O4H-1P") }, //"O4H2P" -H

                { ProductType.wBaseLoss, ChemicalFormula.ParseFormula("H-2") }, //"H"-H
                { ProductType.xBaseLoss, ChemicalFormula.ParseFormula("O-1H-2") }, //"O-1H" -H
                { ProductType.yBaseLoss, ChemicalFormula.ParseFormula("O-3H-2P-1") }, //"O-3P-1" -H
                { ProductType.zBaseLoss, ChemicalFormula.ParseFormula("O-4H-3P-1") }, //"O-4H-1P-1" -1

                { ProductType.M, new ChemicalFormula() }
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
                case ProductType.aWaterLoss:
                case ProductType.aBaseLoss:
                case ProductType.b:
                case ProductType.bWaterLoss:
                case ProductType.bBaseLoss:
                case ProductType.c:
                case ProductType.cWaterLoss:
                case ProductType.cBaseLoss:
                case ProductType.d:
                case ProductType.dWaterLoss:
                case ProductType.dBaseLoss:
                    return FragmentationTerminus.FivePrime;

                case ProductType.w:
                case ProductType.wWaterLoss:
                case ProductType.wBaseLoss:
                case ProductType.x:
                case ProductType.xWaterLoss:
                case ProductType.xBaseLoss:
                case ProductType.y:
                case ProductType.yWaterLoss:
                case ProductType.yBaseLoss:
                case ProductType.z:
                case ProductType.zWaterLoss:
                case ProductType.zBaseLoss:
                    return FragmentationTerminus.ThreePrime;

                case ProductType.M:
                    return FragmentationTerminus.None;

                case ProductType.aStar:
                case ProductType.aDegree:
                case ProductType.bAmmoniaLoss:
                case ProductType.yAmmoniaLoss:
                case ProductType.zPlusOne:
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
                        ProductType.a, ProductType.aWaterLoss, ProductType.aBaseLoss,
                        ProductType.b, ProductType.bWaterLoss, ProductType.bBaseLoss,
                        ProductType.c, ProductType.cWaterLoss, ProductType.cBaseLoss,
                        ProductType.d, ProductType.dWaterLoss, ProductType.dBaseLoss, 
                    }
                },
                {
                    FragmentationTerminus.ThreePrime, new List<ProductType>
                    {
                        ProductType.w, ProductType.wWaterLoss, ProductType.wBaseLoss,
                        ProductType.x, ProductType.xWaterLoss, ProductType.xBaseLoss,
                        ProductType.y, ProductType.yWaterLoss, ProductType.yBaseLoss,
                        ProductType.z, ProductType.zWaterLoss, ProductType.zBaseLoss,
                    }
                },
                {
                    FragmentationTerminus.Both, new List<ProductType>
                    {

                        ProductType.a, ProductType.aWaterLoss, ProductType.aBaseLoss,
                        ProductType.b, ProductType.bWaterLoss, ProductType.bBaseLoss,
                        ProductType.c, ProductType.cWaterLoss, ProductType.cBaseLoss,
                        ProductType.d, ProductType.dWaterLoss, ProductType.dBaseLoss, 
                        ProductType.w, ProductType.wWaterLoss, ProductType.wBaseLoss,
                        ProductType.x, ProductType.xWaterLoss, ProductType.xBaseLoss,
                        ProductType.y, ProductType.yWaterLoss, ProductType.yBaseLoss,
                        ProductType.z, ProductType.zWaterLoss, ProductType.zBaseLoss,
                        ProductType.M
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
