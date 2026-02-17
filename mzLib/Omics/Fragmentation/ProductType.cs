namespace Omics.Fragmentation
{
    public enum ProductType
    {
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

        // Base ions are for Nucleic acids in which the base is cleaved as a neutral loss during fragmentation
        // schematic for RNA fragmentation modes can be found below
        // https://www.researchgate.net/figure/The-standard-nomenclature-for-oligonucleotide-fragmentation-during-collisioninduced_fig6_271536997
        // Base loss ions are for Nucleic acids in which the base is cleaved as a neutral loss during fragmentation
        //     These base losses have only been explicetly confirmed for 3' fragments (a,b,c,d)
        //     The base loss ions for 5' fragments (w,x,y,z) are theoretical and have not been confirmed

        a,
        aStar,
        aDegree,
        aWaterLoss,
        aBaseLoss,
        b,
        bAmmoniaLoss,
        bWaterLoss,
        //BnoB1ions,
        bBaseLoss,
        c,
        cWaterLoss,
        cBaseLoss,
        d,
        dWaterLoss,
        dBaseLoss,
        w,
        wWaterLoss,
        wBaseLoss,
        x,
        xWaterLoss,
        xBaseLoss,
        y,
        yAmmoniaLoss,
        yWaterLoss,
        yBaseLoss,
        z,
        zPlusOne,       //This is zDot plus H
        zDot,
        zWaterLoss,
        zBaseLoss,
        M,              //this is the molecular ion // [M]
        D,              //this is a diagnostic ion // Modification loss mass
        Ycore,          //Glyco core Y ions // [pep] + Neutral core Glycan mass (such as: [pep] + [N]) //Which already consider the loss of H2O and H-transfer
        Y               //Glyco Y ions // [pep] + other Glycan mass 
    }

    public static class ProductTypeExtensions 
    {
        public static bool IsBaseLoss(this ProductType productType)
        {
            return productType is ProductType.aBaseLoss or ProductType.bBaseLoss or ProductType.cBaseLoss or ProductType.dBaseLoss or ProductType.wBaseLoss or ProductType.xBaseLoss or ProductType.yBaseLoss or ProductType.zBaseLoss;
        }

        public static bool IsWaterLoss(this ProductType productType)
        {
            return productType is ProductType.aWaterLoss or ProductType.bWaterLoss or ProductType.cWaterLoss or ProductType.dWaterLoss or ProductType.wWaterLoss or ProductType.xWaterLoss or ProductType.yWaterLoss or ProductType.zWaterLoss;
        }

        /// <summary>
        /// Returns all product types that are in the same "family" as the input product type. This means they share a cleavage site and direction. For example, a, a*, a°, a-H2O, and a-BaseLoss are all in the same family because they are all a-type ions. 
        /// </summary>
        /// <param name="productType"></param>
        /// <returns></returns>
        public static ProductType[] GetFragmentFamilyMembers(this ProductType productType)
        {
            return productType switch
            {
                ProductType.a => [ProductType.aStar, ProductType.aDegree, ProductType.aWaterLoss, ProductType.aBaseLoss],
                ProductType.aStar => [ProductType.a, ProductType.aDegree, ProductType.aWaterLoss, ProductType.aBaseLoss],
                ProductType.aDegree => [ProductType.a, ProductType.aStar, ProductType.aWaterLoss, ProductType.aBaseLoss],
                ProductType.aWaterLoss => [ProductType.a, ProductType.aStar, ProductType.aDegree, ProductType.aBaseLoss],
                ProductType.aBaseLoss => [ProductType.a, ProductType.aStar, ProductType.aDegree, ProductType.aWaterLoss],

                ProductType.b => [ProductType.bAmmoniaLoss, ProductType.bWaterLoss, ProductType.bBaseLoss],
                ProductType.bAmmoniaLoss => [ProductType.b, ProductType.bWaterLoss, ProductType.bBaseLoss],
                ProductType.bWaterLoss => [ProductType.b, ProductType.bAmmoniaLoss, ProductType.bBaseLoss],
                ProductType.bBaseLoss => [ProductType.b, ProductType.bAmmoniaLoss, ProductType.bWaterLoss],

                ProductType.c => [ProductType.cWaterLoss, ProductType.cBaseLoss],
                ProductType.cWaterLoss => [ProductType.c, ProductType.cBaseLoss],
                ProductType.cBaseLoss => [ProductType.c, ProductType.cWaterLoss],

                ProductType.d => [ProductType.dWaterLoss, ProductType.dBaseLoss],
                ProductType.dWaterLoss => [ProductType.d, ProductType.dBaseLoss],
                ProductType.dBaseLoss => [ProductType.d, ProductType.dWaterLoss],

                ProductType.w => [ProductType.wWaterLoss, ProductType.wBaseLoss],
                ProductType.wWaterLoss => [ProductType.w, ProductType.wBaseLoss],
                ProductType.wBaseLoss => [ProductType.w, ProductType.wWaterLoss],

                ProductType.x => [ProductType.xWaterLoss, ProductType.xBaseLoss],
                ProductType.xWaterLoss => [ProductType.x, ProductType.xBaseLoss],
                ProductType.xBaseLoss => [ProductType.x, ProductType.xWaterLoss],

                ProductType.y => [ProductType.yAmmoniaLoss, ProductType.yWaterLoss, ProductType.yBaseLoss],
                ProductType.yAmmoniaLoss => [ProductType.y, ProductType.yWaterLoss, ProductType.yBaseLoss],
                ProductType.yWaterLoss => [ProductType.y, ProductType.yAmmoniaLoss, ProductType.yBaseLoss],
                ProductType.yBaseLoss => [ProductType.y, ProductType.yAmmoniaLoss, ProductType.yWaterLoss],

                ProductType.z => [ProductType.zPlusOne, ProductType.zDot, ProductType.zWaterLoss, ProductType.zBaseLoss],
                ProductType.zPlusOne => [ProductType.z, ProductType.zDot, ProductType.zWaterLoss, ProductType.zBaseLoss],
                ProductType.zDot => [ProductType.z, ProductType.zPlusOne, ProductType.zWaterLoss, ProductType.zBaseLoss],
                ProductType.zWaterLoss => [ProductType.z, ProductType.zPlusOne, ProductType.zDot, ProductType.zBaseLoss],
                ProductType.zBaseLoss => [ProductType.z, ProductType.zPlusOne, ProductType.zDot, ProductType.zWaterLoss],

                ProductType.Ycore => [ProductType.Y],
                ProductType.Y => [ProductType.Ycore],

                ProductType.M => [],
                ProductType.D => [],
                _ => []
            };
        }
    }
}
