using System;
using System.Collections.Generic;
using System.Text;

namespace Proteomics.Fragmentation
{
    public class TerminusSpecificProductTypes
    {
        public static Dictionary<FragmentationTerminus, List<ProductType>> ProductIonTypesFromSpecifiedTerminus = new Dictionary<FragmentationTerminus, List<ProductType>>
        {
            {FragmentationTerminus.N, new List<ProductType>{ ProductType.a, ProductType.aDegree, ProductType.aStar, ProductType.b, ProductType.bDegree, ProductType.bStar, ProductType.c } }, //all ion types that include the N-terminus
            {FragmentationTerminus.C, new List<ProductType>{ ProductType.x, ProductType.y, ProductType.yDegree, ProductType.yStar, ProductType.zPlusOne } }, //all ion types that include the C-terminus
            {FragmentationTerminus.Both, new List<ProductType>{ ProductType.a, ProductType.aDegree, ProductType.aStar, ProductType.b, ProductType.bDegree, ProductType.bStar, ProductType.c, ProductType.x, ProductType.y, ProductType.yDegree, ProductType.yStar, ProductType.zPlusOne} },
            {FragmentationTerminus.None, new List<ProductType>() }
        };

        public static Dictionary<ProductType, FragmentationTerminus> ProductTypeToFragmentationTerminus = new Dictionary<ProductType, FragmentationTerminus>
        {
            { ProductType.a, FragmentationTerminus.N },
            { ProductType.aDegree, FragmentationTerminus.N },
            { ProductType.aStar, FragmentationTerminus.N },
            { ProductType.b, FragmentationTerminus.N },
            { ProductType.bDegree, FragmentationTerminus.N },
            { ProductType.bStar, FragmentationTerminus.N },
            { ProductType.c, FragmentationTerminus.N },
            { ProductType.x, FragmentationTerminus.C },
            { ProductType.y, FragmentationTerminus.C },
            { ProductType.yDegree, FragmentationTerminus.C },
            { ProductType.yStar, FragmentationTerminus.C },
            { ProductType.zPlusOne, FragmentationTerminus.C },
        };

    }
}
