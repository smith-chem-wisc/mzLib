namespace Omics.Fragmentation.Peptide
{
    public class TerminusSpecificProductTypes
    {
        /// <summary>
        /// The types of ions that can be generated from a peptide fragment, based on the terminus of the fragment
        /// </summary>
        public static Dictionary<FragmentationTerminus, List<ProductType>> ProductIonTypesFromSpecifiedTerminus = new Dictionary<FragmentationTerminus, List<ProductType>>
        {
            {FragmentationTerminus.N, new List<ProductType>{ ProductType.a, ProductType.aDegree, ProductType.aStar, ProductType.b, ProductType.bWaterLoss, ProductType.bAmmoniaLoss, ProductType.c } }, //all ion types that include the N-terminus
            {FragmentationTerminus.C, new List<ProductType>{ ProductType.x, ProductType.y, ProductType.yWaterLoss, ProductType.yAmmoniaLoss, ProductType.zDot, ProductType.zPlusOne } }, //all ion types that include the C-terminus
            {FragmentationTerminus.Both, new List<ProductType>{ ProductType.a, ProductType.aDegree, ProductType.aStar, ProductType.b, ProductType.bWaterLoss, ProductType.bAmmoniaLoss, ProductType.c, ProductType.x, ProductType.y, ProductType.yWaterLoss, ProductType.yAmmoniaLoss, ProductType.zDot, ProductType.zPlusOne} },
            {FragmentationTerminus.None, new List<ProductType>() }
        };

        /// <summary>
        /// The terminus of the peptide fragment that the product ion is generated from
        /// </summary>
        public static Dictionary<ProductType, FragmentationTerminus> ProductTypeToFragmentationTerminus = new Dictionary<ProductType, FragmentationTerminus>
        {
            { ProductType.a, FragmentationTerminus.N },
            { ProductType.aDegree, FragmentationTerminus.N },
            { ProductType.aStar, FragmentationTerminus.N },
            { ProductType.b, FragmentationTerminus.N },
            { ProductType.bWaterLoss, FragmentationTerminus.N },
            { ProductType.bAmmoniaLoss, FragmentationTerminus.N },
            { ProductType.c, FragmentationTerminus.N },
            { ProductType.x, FragmentationTerminus.C },
            { ProductType.y, FragmentationTerminus.C },
            { ProductType.yWaterLoss, FragmentationTerminus.C },
            { ProductType.yAmmoniaLoss, FragmentationTerminus.C },
            { ProductType.zDot, FragmentationTerminus.C },
            { ProductType.zPlusOne, FragmentationTerminus.C },
        };

    }
}
